//! The actual simulation

use crate::{
    float, vec_utils, Branch, BranchId, BranchKind, BranchTree, Float, Tube, ATMOSPHERIC_PRESSURE,
};
use sparse21::Matrix;

/// The environment of the simulation, at a given tick
///
/// All pressures are in units of Pascals.
pub struct SimulationEnvironment {
    /// The pressure exerted on the lungs by the pleural cavity -- the portion of the pressure that
    /// the diaphragm is responsible for. Units of pascals (kg/ms²)
    ///
    /// TODO: This could have an attached gradient to take gravity into account.
    pub pleural_pressure: Float,

    /// The pressure at the tracheal boundary. In Foy's paper, this is always equal to atmospheric
    /// pressure.
    pub tracheal_pressure: Float,

    /// The viscosity of the air inside the lungs. Units of kg/m*s
    pub air_viscosity: Float,

    /// The density of the air (assumed incompressible) inside the lungs. Units of kg/m³
    pub air_density: Float,
}

impl SimulationEnvironment {
    /// Resets all of the airflow and pressure in the tree
    ///
    /// The airflow is all set to zero and the volume of each acinar region is calculated to sum to
    /// the target volume, distributed into each acinar region according to its compliance relative
    /// to the others. All of the pressure values are set to match the volume of the acinar regions
    /// and the current pleural pressure.
    ///
    /// If the target volume is not given, the tree will be set with atmospheric pressure.
    pub fn reset(&self, tree: &mut BranchTree, target_volume: Option<Float>) {
        // Returns the total compliance (i.e. summed) of the subtree rooted at the given branch
        //
        // We use this instead of e.g., looping through all branches, because adding together terms
        // with similar magnitude will reduce floating-point inaccuracy. Because the outer function
        // (`SimulationEnvironment::reset`) is only run once, at startup, we're happy to pay the
        // penalty of running in O(n log n) time instead of O(n) to get the benefit of more
        // accurate compliance measurements.
        fn subtree_compliance(tree: &BranchTree, id: BranchId) -> Float {
            match &tree[id] {
                Branch::Bifurcation(b) => {
                    subtree_compliance(tree, b.left_child) + subtree_compliance(tree, b.right_child)
                }
                Branch::Acinar(b) => b.compliance,
            }
        }

        let internal_pressure: Float;
        let acinar_pressure_differential: Float;

        if let Some(target) = target_volume {
            let total_compliance = subtree_compliance(tree, tree.root_id());

            // Using the total compliance, determine the ratio we use to set the volume of an acinar
            // region from its compliance.
            //
            // Essentially, because `volume = compliance * pressure`, we're determining the pressure
            // we'd need to have at all acinar regions in order to get the desired volume.
            acinar_pressure_differential = target / total_compliance;
            // From this and the pleural pressure, we can then derive the actual pressure inside the
            // lungs.
            //
            // This doesn't matter so much, because only volume is carried over from one simulation
            // tick to the next.
            internal_pressure = acinar_pressure_differential - self.pleural_pressure;
        } else {
            internal_pressure = ATMOSPHERIC_PRESSURE;
            acinar_pressure_differential = ATMOSPHERIC_PRESSURE - self.pleural_pressure;
        };

        // The acinar pressure differential is passed separately to reduce floating point
        // inaccuracy
        fn set_all(
            tree: &mut BranchTree,
            pressure: Float,
            acinar_pressure_differential: Float,
            id: BranchId,
        ) {
            let branch = &mut tree[id];

            let tube = branch.tube_mut();
            tube.flow_in = 0.0;
            tube.end_pressure = pressure;

            match branch {
                &mut Branch::Bifurcation(b) => {
                    set_all(tree, pressure, acinar_pressure_differential, b.left_child);
                    set_all(tree, pressure, acinar_pressure_differential, b.right_child);
                }
                Branch::Acinar(b) => {
                    b.volume = b.compliance * acinar_pressure_differential;
                }
            }
        }

        let root_id = tree.root_id();
        set_all(
            tree,
            internal_pressure,
            acinar_pressure_differential,
            root_id,
        );
    }

    /// Updates the state of the tree, given the external pressures
    ///
    /// The update corresponds to the change over the duration in seconds specified by `timestep`.
    ///
    /// This method may return an error if some kind of unexpected numerical hiccup occurs. If this
    /// is the case, the tree will not have been updated but it is unlikely to be possible to
    /// continue.
    pub fn do_tick(&self, tree: &mut BranchTree, timestep: Float) -> Result<(), &'static str> {
        // The relationship between individual state values in the tree is nonlinear. So in order
        // to perform the actual update from a given timestep, we need to solve the system of
        // nonlinear equations.
        //
        // Because analytical solutions are difficult, we'll use Newton's method to determine the
        // state of the tree at the current timestep.
        //
        // For more information about the way we represent the current state of the tree, refer to
        // the `state_vector` function towards the bottom of this file.

        const MAX_ITERS: usize = 100;
        const TOL: Float = 1e-6;

        // The current state is probably pretty close to where we want to be.
        let mut state = state_vector(tree);
        // Pre-compute the value of the optimization function, because we reuse this across
        // iterations.
        let mut opt_func = self.optimization_func_at(tree, &state, timestep);

        for iter in 0.. {
            if iter >= MAX_ITERS {
                return Err("Reached maximum iteration count without solution");
            }

            // Standard vector newton's method is:
            //
            //   dx = -[J(x)]^-1 * f(x)
            //
            // which requires finding the inverse of the jacobian (possibly very large!). The
            // normal way to resolve this is by instead solving for dx with this system of
            // equations:
            //
            //   J(x)(-dx) = f(x)
            //
            // ... which is what we do below, preserving the negative sign around 'dx'.

            let mut jacobian = self.jacobian_at(tree, &state, timestep);

            // precompute this, so we can pass it in by value.
            let opt_func_norm_squared = vec_utils::norm_squared(&opt_func);

            // solve for dx:
            let minus_dx = jacobian
                .solve(opt_func)
                .ok()
                .ok_or("failed to solve system of equations for dx in newton's method")?;

            vec_utils::sub_assign(&mut state, &minus_dx);
            opt_func = self.optimization_func_at(tree, &state, timestep);

            // And then we're done if the state is "close enough" to the optimal solution.
            if vec_utils::norm_squared(&minus_dx) <= TOL && opt_func_norm_squared <= TOL {
                break;
            }
        }

        set_state_from_vector(tree, &state);
        Ok(())
    }

    /// Returns the value of the optimization function at the point given by the state vector
    ///
    /// This function is set up so that the system of simultaneous equations has a solution when
    /// the optimization function equals zero this (this is pretty standard). It's used above in
    /// [`do_tick`], inside the loop.
    ///
    /// ## Layout of returned values
    ///
    /// Each returned value represents a single equation. These equations are, in order:
    ///
    /// * ΔP - flow_resistance(Q)Q = 0; for each branch bifurcations..., acinar regions...
    /// * Qᵢ - Qₗ - Qᵣ = 0; for each bifurcation
    /// * V - dt * Q - V₀ = 0; for each acinar region (V₀ is the volume at the previous timestep)
    /// * P - plural_pressure - V/C = 0; for each acinar region (C is compliance)
    ///
    /// where `dt` is the given timestep. Counting carefully gives a total of `2 *
    /// tree.count_bifurcations() + 3 * tree.count_acinar_regions()` equations; or that many values
    /// in the returned vector.
    ///
    /// [`do_tick`]: Self::do_tick
    fn optimization_func_at(
        &self,
        tree: &BranchTree,
        state: &[Float],
        timestep: Float,
    ) -> Vec<Float> {
        let n_outputs = tree.optimization_func_output_size();
        let mut output = vec![0.0; n_outputs];

        let root_id = tree.root_id();
        assert_eq!(root_id.deconstruct().0, BranchKind::Bifurcation);

        for idx in 0..tree.count_bifurcations() {
            let id = BranchId::new(BranchKind::Bifurcation, idx);
            let b = match &tree[id] {
                Branch::Bifurcation(b) => b,
                Branch::Acinar(_) => unreachable!(),
            };

            let this_pressure = state[tree.state_pressure_index(id)];
            let this_flow = state[tree.state_flow_index(id)];

            // (possible) root pressure differential:
            if id == root_id {
                let this_resistance = self.flow_resistance(tree[id].tube(), this_flow);
                // Convert `self.tracheal_pressure` from P to P̂:
                let tracheal_pressure = self.tracheal_pressure - ATMOSPHERIC_PRESSURE;
                let this_pressure_eqn =
                    tracheal_pressure - this_pressure - this_resistance * this_flow;

                output[tree.pressure_delta_eqn_index(id)] = this_pressure_eqn;
            }

            // Pressure differential of the children:
            //
            // Left first:
            let left_pressure = state[tree.state_pressure_index(b.left_child)];
            let left_flow = state[tree.state_flow_index(b.left_child)];
            let left_resistance = self.flow_resistance(tree[b.left_child].tube(), left_flow);
            let left_pressure_eqn = this_pressure - left_pressure - left_resistance * left_flow;

            output[tree.pressure_delta_eqn_index(b.left_child)] = left_pressure_eqn;

            // And then right:
            let right_pressure = state[tree.state_pressure_index(b.right_child)];
            let right_flow = state[tree.state_flow_index(b.right_child)];
            let right_resistance = self.flow_resistance(tree[b.right_child].tube(), right_flow);
            let right_pressure_eqn = this_pressure - right_pressure - right_resistance * right_flow;

            output[tree.pressure_delta_eqn_index(b.right_child)] = right_pressure_eqn;

            // Sum of flow from children:
            output[tree.flow_conservation_eqn_index(id)] = this_flow - (left_flow + right_flow);
        }

        for idx in 0..tree.count_acinar_regions() {
            let id = BranchId::new(BranchKind::Acinar, idx);
            let b = match &tree[id] {
                Branch::Acinar(b) => b,
                Branch::Bifurcation(_) => unreachable!(),
            };

            let this_pressure = state[tree.state_pressure_index(id)];
            let this_flow = state[tree.state_flow_index(id)];

            // Air conservation:
            let old_volume = b.volume - v_atmospheric(b.compliance);
            let new_volume = state[tree.state_volume_index(id)];
            output[tree.air_conservation_eqn_index(id)] =
                new_volume - timestep * this_flow - old_volume;

            // Acinar pressure:
            output[tree.acinar_pressure_eqn_index(id)] =
                this_pressure - self.pleural_pressure - new_volume / b.compliance;
        }

        output
    }

    /// Returns the jacobian of the optimization function at the point given by the state vector
    ///
    /// This is used with the standard vector Newton's method to find the state vector that sets
    /// the optimization function to zero.
    ///
    /// For general information about the optimization function, refer to [`optimization_func_at`].
    ///
    /// [`optimization_func_at`]: Self::optimization_func_at
    fn jacobian_at(&self, tree: &BranchTree, state: &[Float], timestep: Float) -> Matrix {
        // Make it a matrix right off the bat, so we can use index pairs to set values.
        let mut jacobian = Matrix::new();

        // Because the equations are *mostly* linear, we're mostly filling slots with 1 or -1.
        // There are a few that aren't, the pressure delta equation has the nasty resistance(Q) * Q
        // term in it. The only derivative we care about there is w.r.t. the flow. That's given by
        // the `flow_resistance_term_derivative` method.
        //
        // To actually add all of the equations, we'll mostly follow the way that the output in
        // `optimization_func_at` is set. That way, it should be easier to track (and it gives us a
        // consistent way to go through all of the values).

        let root_id = tree.root_id();
        assert_eq!(root_id.deconstruct().0, BranchKind::Bifurcation);

        for idx in 0..tree.count_bifurcations() {
            let id = BranchId::new(BranchKind::Bifurcation, idx);
            let b = match &tree[id] {
                Branch::Bifurcation(b) => b,
                Branch::Acinar(_) => unreachable!(),
            };

            let p_idx = tree.state_pressure_index(id);
            let q_idx = tree.state_flow_index(id);

            // (possible) root pressure differential:
            if id == root_id {
                let eqn_idx = tree.pressure_delta_eqn_index(id);
                let flow_rate = state[q_idx];

                jacobian.add_element(eqn_idx, p_idx, -1.0);
                jacobian.add_element(
                    eqn_idx,
                    q_idx,
                    -self.flow_resistance_term_derivative(&b.tube, flow_rate),
                );
            }

            // Child pressure differentials:
            let left_p_idx = tree.state_pressure_index(b.left_child);
            let left_q_idx = tree.state_flow_index(b.left_child);
            let left_flow = state[left_q_idx];
            let left_eqn_idx = tree.pressure_delta_eqn_index(b.left_child);

            jacobian.add_element(left_eqn_idx, p_idx, 1.0);
            jacobian.add_element(left_eqn_idx, left_p_idx, -1.0);
            jacobian.add_element(
                left_eqn_idx,
                left_q_idx,
                -self.flow_resistance_term_derivative(tree[b.left_child].tube(), left_flow),
            );

            let right_p_idx = tree.state_pressure_index(b.right_child);
            let right_q_idx = tree.state_flow_index(b.right_child);
            let right_flow = state[right_q_idx];
            let right_eqn_idx = tree.pressure_delta_eqn_index(b.right_child);

            jacobian.add_element(right_eqn_idx, p_idx, 1.0);
            jacobian.add_element(right_eqn_idx, right_p_idx, -1.0);
            jacobian.add_element(
                right_eqn_idx,
                right_q_idx,
                -self.flow_resistance_term_derivative(tree[b.right_child].tube(), right_flow),
            );

            // Sum of flow:
            let sum_eqn_idx = tree.flow_conservation_eqn_index(id);
            jacobian.add_element(sum_eqn_idx, q_idx, 1.0);
            jacobian.add_element(sum_eqn_idx, right_q_idx, -1.0);
            jacobian.add_element(sum_eqn_idx, left_q_idx, -1.0);
        }

        for idx in 0..tree.count_acinar_regions() {
            let id = BranchId::new(BranchKind::Acinar, idx);
            let b = match &tree[id] {
                Branch::Acinar(b) => b,
                Branch::Bifurcation(_) => unreachable!(),
            };

            // Air conservation:
            let v_idx = tree.state_volume_index(id);
            let q_idx = tree.state_flow_index(id);
            let air_cons_eqn_idx = tree.air_conservation_eqn_index(id);

            jacobian.add_element(air_cons_eqn_idx, v_idx, 1.0);
            jacobian.add_element(air_cons_eqn_idx, q_idx, -timestep);

            // Acinar pressure:
            let p_idx = tree.state_pressure_index(id);
            let acinar_p_eqn_idx = tree.acinar_pressure_eqn_index(id);

            jacobian.add_element(acinar_p_eqn_idx, p_idx, 1.0);
            jacobian.add_element(acinar_p_eqn_idx, v_idx, -1.0 / b.compliance);
        }

        jacobian
    }

    /// Returns the the flow resistance of the tube, in units of kg/m⁴*s
    ///
    /// This method is only allowed to use the properties of the tube that *aren't* tracked
    /// separately in the system state. Hence why `flow_rate` is passed in separately.
    fn flow_resistance(&self, tube: &Tube, flow_rate: Float) -> Float {
        // Normally would be 4.0 / ... diameter; we're using radius, so 2.0 at the top
        let reynolds = (2.0 * self.air_density * flow_rate.abs())
            / (self.air_viscosity * float::PI * tube.radius);

        // Similarly to above. normally it's diameter ^ 4, with 32 in the numerator, but since
        // we're using the radius, we just use a factor of 2.
        (2.0 * self.air_viscosity * tube.length * RESISTANCE_CORRECTION)
            / (float::PI * tube.radius.powi(4))
            * (reynolds * 2.0 * tube.radius / tube.length).sqrt()
        // 1.0
    }

    /// Returns the derivative of the flow resistance TERM w.r.t. the flow rate
    ///
    /// N.B. This is more than just what the `flow_resistance` method covers; we're producing the
    /// derivative of the full `flow_resistance(Q) * Q` term.
    fn flow_resistance_term_derivative(&self, tube: &Tube, flow_rate: Float) -> Float {
        // 1.0
        1.5 * self.flow_resistance(tube, flow_rate)
    }
}

// Correction factor that Foy mentions, taken from Pedley et al. 1970, to compensate for energy
// dissipation in a branching network.
const RESISTANCE_CORRECTION: Float = 1.85;

/// Returns the volume than an acinar region *would have* under atmospheric pressure with no
/// pleural pressure
fn v_atmospheric(compliance: Float) -> Float {
    // In general, V = compliance * (pressure - pleural pressure), so we get this simple formula
    // because we're taking pleural pressure to be zero.
    compliance * ATMOSPHERIC_PRESSURE
}

/// Produces the "state vector" for the tree -- i.e. a vector representation of all of the state
/// values for the tree
///
/// This is used in combination with [`set_state_from_vector`] to set the state of the tree from a
/// similar vector once the new values have been determined.
///
/// The state vector is formatted so that it is filled in order of `(P..., Q..., V...)`, with
/// tube pressures `P`, tube flow rates `Q`, and acinar volumes `V`. Additionally, all of the
/// bifurcations are placed in order of `BranchId`, before the acinar regions within `P...` and
/// `Q...`.
///
/// The expected size of the vector is then exactly `2 * tree.count_branches() +
/// tree.count_acinar_regions()`.
///
/// [`set_state_from_vector`]: Self::set_state_from_vector
fn state_vector(tree: &BranchTree) -> Vec<Float> {
    let var_count = tree.state_size();
    let mut state_vec: Vec<Float> = vec![0.0; var_count];

    // Crawl the entire tree, setting each value as appropriate
    fn crawl(tree: &BranchTree, vec: &mut [Float], id: BranchId) {
        let branch = &tree[id];
        let (kind, index) = id.deconstruct();
        debug_assert_eq!(kind, branch.kind());

        // For determining the index of the values to set, we might need to offset the index within
        // the 'P' or 'Q' region, depending on whether the branch is a bifurcation is not.
        //
        // Because `BanchKind::Bifurcation` has a discriminant of 0 (and `Acinar` is 1), we can
        // multiply the discriminant for conditional addition.
        let p_base = kind as usize * tree.count_bifurcations();
        let q_base = tree.count_branches() + kind as usize * tree.count_bifurcations();

        let tube = branch.tube();
        vec[p_base + index] = tube.end_pressure - ATMOSPHERIC_PRESSURE;
        vec[q_base + index] = tube.flow_in;

        match branch {
            Branch::Bifurcation(b) => {
                crawl(tree, vec, b.left_child);
                crawl(tree, vec, b.right_child);
            }
            Branch::Acinar(b) => {
                let v_base = 2 * tree.count_branches();
                vec[v_base + index] = b.volume - v_atmospheric(b.compliance);
            }
        }
    }

    crawl(tree, &mut state_vec, tree.root_id());

    state_vec
}

/// Sets the entire state of the tree from a state vector
///
/// For more information, see: [`state_vector`](Self::state_vector).
///
/// ## Panics
///
/// Panics if the length of the vector doesn't match what's expected.
fn set_state_from_vector(tree: &mut BranchTree, vec: &[Float]) {
    // We're expecting it to be a column vector with the length specified in `state_vector`.
    let expected_len = 2 * tree.count_branches() + tree.count_acinar_regions();
    assert_eq!(vec.len(), expected_len);

    fn crawl(tree: &mut BranchTree, vec: &[Float], id: BranchId) {
        // prefetch these so we don't have borrow conflicts with `branch`
        let n_bifurcations = tree.count_bifurcations();
        let n_branches = tree.count_branches();

        let branch = &mut tree[id];
        let (kind, index) = id.deconstruct();
        debug_assert_eq!(kind, branch.kind());

        // Taken from `crawl` in `state_vector` above. Refer there for more info.
        let p_base = kind as usize * n_bifurcations;
        let q_base = n_branches + kind as usize * n_bifurcations;

        let tube = branch.tube_mut();
        tube.end_pressure = vec[p_base + index] + ATMOSPHERIC_PRESSURE;
        tube.flow_in = vec[q_base + index];
        drop(tube);

        match branch {
            &mut Branch::Bifurcation(b) => {
                crawl(tree, vec, b.left_child);
                crawl(tree, vec, b.right_child);
            }
            Branch::Acinar(b) => {
                let v_base = 2 * n_branches;
                b.volume = vec[v_base + index] + v_atmospheric(b.compliance);
            }
        }
    }

    let root = tree.root_id();
    crawl(tree, vec, root);
}

// Helper methods for indexing state vectors and the optimization function
impl BranchTree {
    /// Returns the size of the state vector for this tree
    fn state_size(&self) -> usize {
        2 * self.count_branches() + self.count_acinar_regions()
    }

    /// Returns the number of outputs (i.e. equations represented) in the optimization function for
    /// this tree
    fn optimization_func_output_size(&self) -> usize {
        2 * self.count_bifurcations() + 3 * self.count_acinar_regions()
    }

    /// Returns the index of the branch's pressure in the state vector for this tree
    ///
    /// N.B.: The pressure is from the `end_pressure` field of the branch.
    fn state_pressure_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        // offset by bifurcations if it's an acinar region.
        idx + kind as usize * self.count_bifurcations()
    }

    /// Returns the index of the flow through this branch in the state vector for this tree
    fn state_flow_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        self.count_branches() + idx + kind as usize * self.count_bifurcations()
    }

    /// Returns the index of the volume of this acinar region in the state vector for this tree
    ///
    /// ## Panics
    ///
    /// This method panics if the provided `BranchId` doesn't correspond to an acinar region.
    fn state_volume_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        assert_eq!(kind, BranchKind::Acinar);

        2 * self.count_branches() + idx
    }

    /// Returns the index in the output vector of [`optimization_func_at`] that corresponds to the
    /// pressure delta equation for the given branch
    ///
    /// The equation is of the form: ΔP - flow_resistance(Q)Q = 0.
    fn pressure_delta_eqn_index(&self, id: BranchId) -> usize {
        // This equation is first, with all bifurcations and then acinar regions.
        let (kind, idx) = id.deconstruct();

        // offset by bifurcations if it's an acinar region.
        idx + kind as usize * self.count_bifurcations()
    }

    /// Returns the index in the output vector of [`optimization_func_at`] that corresponds to the
    /// flow conservation equation for the given bifurcation
    ///
    /// The equation is of the form: Qᵢ - Qₗ - Qᵣ = 0, where Qᵢ is the provided branch, with
    /// children Qₗ and Qᵣ.
    ///
    /// ## Panics
    ///
    /// This method panics if the given `BranchId` doesn't correspond to a bifurcation
    fn flow_conservation_eqn_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        assert_eq!(kind, BranchKind::Bifurcation);

        self.count_branches() + idx
    }

    /// Returns the index in the output vector of [`optimization_func_at`] that corresponds to the
    /// air conservation (flow into volume) equation for the given acinar region
    ///
    /// The equation is of the form: Q - dt * V = 0.
    ///
    /// ## Panics
    ///
    /// This method panics if the given `BranchId` doesn't correspond to an acinar region
    fn air_conservation_eqn_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        assert_eq!(kind, BranchKind::Acinar);

        self.count_branches() + self.count_bifurcations() + idx
    }

    /// Returns the index in the output vector of [`optimization_func_at`] that corresponds to the
    /// acinar pressure equation for the given acinar region
    ///
    /// The equation balances the internal pressure with the pleural pressure on the region (from
    /// the diaphragm). It is of the form: P - plural_pressure - V/C = 0.
    ///
    /// ## Panics
    ///
    /// This method panics if the given `BranchId` doesn't correspond to an acinar region
    fn acinar_pressure_eqn_index(&self, id: BranchId) -> usize {
        let (kind, idx) = id.deconstruct();
        assert_eq!(kind, BranchKind::Acinar);

        2 * self.count_branches() + idx
    }
}

// Just some simple tests that setting state vectors works correctly.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::gen::equal::EqualChildGenerator;
    use crate::gen::ParentInfo;
    use crate::Point;
    use rand::Rng;

    // check that (tree -> state vector -> tree) doesn't change the values
    #[test]
    fn test_state_vector_idempotent() {
        let generator = EqualChildGenerator {
            max_depth: 3,
            compliance: 2e-12,
        };
        let start = ParentInfo {
            id: 0,
            pos: Point { x: 0.01, y: 0.02 },
            total_angle: -float::FRAC_PI_2,
            length: 0.01,
            tube_radius: 0.001,
        };

        let mut tree = BranchTree::from_generator(start, &generator);
        let mut rng = rand::thread_rng();

        // Randomly set values in the tree:
        for branch in tree.branches_mut() {
            if let Branch::Acinar(b) = branch {
                b.volume = rng.gen();
            }

            let tube = branch.tube_mut();
            tube.flow_in = rng.gen();
            tube.end_pressure = rng.gen();
        }

        // Generate the state vector
        let state_vector = state_vector(&tree);

        // Clone this copy for later and un-set all of the values:
        let original = tree.clone();
        for branch in tree.branches_mut() {
            if let Branch::Acinar(b) = branch {
                b.volume = rng.gen();
            }

            let tube = branch.tube_mut();
            tube.flow_in = 0.0;
            tube.end_pressure = 0.0;
        }

        // Re-set the tree from the state vector and compare equality.
        set_state_from_vector(&mut tree, &state_vector);
        assert_eq!(original, tree);
    }
}
