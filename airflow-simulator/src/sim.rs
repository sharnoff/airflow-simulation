//! The actual simulation

use crate::{Branch, BranchId, BranchTree, Float};

/// The environment of the simulation, at a given tick
///
/// All pressures are in units of MegaPascals.
pub struct SimulationEnvironment {
    /// The pressure on the lungs -- the portion of the pressure that the diaphragm is responsible
    /// for
    ///
    /// TODO: This should actually be distributed, because we need to take the effects of gravity
    /// into account.
    pub pleural_pressure: Float,

    /// The pressure at the trachel boundary. In Foy's paper, this is always equal to atmospheric
    /// pressure.
    pub tracheal_pressure: Float,
}

impl SimulationEnvironment {
    /// Resets all of the airflow and pressure in the tree, ignoring the pleural pressure
    ///
    /// The airflow is all set to zero and the pressure is set to equal the tracheal pressure. This
    /// will only be stable if the pleural pressure is not zero -- pleural pressure *is* taken into
    /// account for determining acinar region volumes, but does not affect flow rate.
    pub fn reset(&self, tree: &mut BranchTree) {
        // This is generally pretty simple; we're just applying a pre-set value across the entire
        // tree.
        //
        // The only *mildly* complicated bit is that we need to calculate the volumes of acinar
        // regions, but there's a pretty simple formula for that.
        fn set_all(env: &SimulationEnvironment, tree: &mut BranchTree, id: BranchId) {
            let branch = &mut tree[id];

            let tube = branch.tube_mut();
            tube.flow_rate = 0.0;
            tube.end_pressure = env.tracheal_pressure;

            match branch {
                &mut Branch::Bifurcation(b) => {
                    set_all(env, tree, b.left_child);
                    set_all(env, tree, b.right_child);
                }
                Branch::Acinar(b) => {
                    // The volume of the acinar region is essentially just given by this formula,
                    // no fancy work required:
                    b.volume = b.compliance * (b.tube.end_pressure - env.pleural_pressure);
                }
            }
        }

        let root_id = tree.root_id;
        set_all(self, tree, root_id);
    }

    /// Updates the state of the tree, given the external pressures
    ///
    /// The update corresponds to the change over the duration in seconds specified by `timestep`.
    pub fn do_tick(&self, tree: &mut BranchTree, timestep: Float) {
        todo!()
    }
}
