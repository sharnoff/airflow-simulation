#!/bin/python
#
# Helper script for generating the TikZ figures for sparse matrices. Called by the paper's makefile
# for each JSON "build" description in turn.
#
# The JSON build format is pretty simple for the most part. It consists of two fields:
#  * 'model' defines the model to use. Currently, only { "kind": "symmetric", "depth": <N> }
#  * 'layout' defines the layout of the jacobean, with lists 'vars' and 'eqns' to define the order
#     of the variables and equations, respectively.
#
# Variables are named "<branch>.<name>" where <branch> is the index of the branch, starting from
# zero, and <name> is one of "P", "Q", or "V" (corresponding to pressure, flow, or volume).
#
# Equations are named "<branch>.<eq-id>" where <eq-id> gives the number of the equation, as
# indicated in the paper. For reference, they are:
# 
#   1:  P_parent - P_i - R(i)Q_i     = 0
#   2:  Q_i - Q_left - Q_right       = 0
#   3:  V_i^t - V_i^{t-1} - dt Q_i^t = 0
#   4:  P_i - P_pl(t) - 1/C_i V_i    = 0
#
# Existing examples should hopefully make this clear.

import json
import argparse
from typing import Callable

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("buildfile", nargs=1,
            help="Sets the build file to use for generating the output file. Must end in .build.json")

    args = parser.parse_args()

    buildfile: str = args.buildfile[0]
    if not buildfile.endswith('.build.json'):
        raise Exception(f"Build file must end in .build.json, got: {buildfile}")

    # outfile = buildfile : s/.build.json$/.tikz.tex/
    outfile = buildfile[:len(buildfile) - len('.build.json')] + '.tikz.tex'

    # Read the build file
    with open(buildfile, 'r') as f:
        build_desc = json.load(f)

    extra_fields = get_extra_fields(build_desc, ok = ['model', 'layout'])
    if extra_fields:
        raise Exception(f"unexpected JSON fields {list(extra_fields)}")

    if "model" not in build_desc or not isinstance(build_desc["model"], dict):
        raise Exception("expected 'model' field as an object")
    if "layout" not in build_desc or not isinstance(build_desc["layout"], dict):
        raise Exception("expected 'layout' field as an object")

    model = Model(build_desc["model"])
    layout = Layout(build_desc["layout"])

    tikz_str = tikz_to_string(generate_values(model, layout))
    with open(outfile, 'w+') as f:
        f.write(tikz_str)

# Helper function for returning the unused fields in a JSON dict, so that we can appropriately
# return an error if there are any.
def get_extra_fields(json_dict: dict, ok: list[str] = []) -> set[str]:
    keys = set(json_dict.keys())
    for k in ok:
        if k in keys:
            keys.remove(k)
    return keys


# Information about a particular branch
class Branch:
    def __init__(self, index: int, parent: int | None, children: tuple[int, int] | None):
        self._index = index
        self._parent = parent
        self._children = children

    @staticmethod
    def acinar(index: int, parent_index: int | None) -> 'Branch':
        return Branch(index, parent_index, None)

    @staticmethod
    def bifurcation(index: int, parent_index: int | None, left_child: int, right_child: int) -> 'Branch':
        return Branch(index, parent_index, (left_child, right_child))

    def index(self) -> int:
        return self._index

    def parent(self) -> int | None:
        return self._parent

    def children(self) -> tuple[int, int] | None:
        return self._children

# Abstraction over the types of models used for 
class Model:
    def __init__(self, json: dict):
        if "kind" not in json or not isinstance(json["kind"], str):
            raise Exception("Model kind must be supplied")

        kind = json["kind"]

        assert isinstance(kind, str)

        class Inner:
            def __init__(self, count_branches: Callable[[], int], branch: Callable[[int], Branch]):
                self.count_branches = count_branches
                self.branch = branch

            @staticmethod
            def symmetric(depth: int) -> 'Inner':
                assert depth >= 1

                # For depth 1, n = 1; depth 2, n = 3; depth 3, n = 7; ...
                #
                # Essentially, the *height* of the tree.
                n_branches = (2 ** depth) - 1

                # the index of the first acinar node.
                fst_acinar = (2 ** (depth - 1)) - 1

                def branch(index: int) -> Branch:
                    assert 0 <= index < n_branches

                    parent = (index - 1) // 2 if index != 0 else None

                    # acinar
                    if index >= fst_acinar:
                        return Branch.acinar(index, parent)

                    # bifurcation
                    left_child = index * 2 + 1
                    right_child = left_child + 1

                    return Branch.bifurcation(index, parent, left_child, right_child)

                return Inner(lambda: n_branches, branch)

        used_fields = ['kind']

        match kind:
            case 'symmetric':
                if "depth" not in json or not isinstance(json["depth"], int):
                    raise Exception("symmetric model requires integer 'depth' parameter")

                used_fields += ['depth']
                self.inner = Inner.symmetric(json["depth"])
            case _:
                raise Exception(f"Unknown model kind '{kind}'")

        # Check for any extra fields that were provided to the model but unused
        extra_fields = get_extra_fields(json, ok=used_fields)
        if extra_fields:
            raise Exception(f"unexpected model fields {list(extra_fields)}")

    def count_branches(self) -> int:
        return self.inner.count_branches()

    def branch(self, index: int) -> Branch:
        return self.inner.branch(index)

Variable = tuple[int, str]

# Class to represent equation types
class Equation:
    def __init__(self, branch: int, variables: Callable[[Branch], list[Variable]]):
        self._branch = branch
        self._variables = variables

    def variables(self, model: Model) -> list[Variable]:
        return self._variables(model.branch(self._branch))

# The layout of the matrix image we generate
class Layout:
    def __init__(self, json: dict):
        extra_fields = get_extra_fields(json, ok = ['eqns', 'vars'])
        if extra_fields:
            raise Exception(f"unexpected layout fields {list(extra_fields)}")

        if "eqns" not in json or not isinstance(json["eqns"], list):
            raise Exception("Expected a list in 'layout.eqns'")
        if "vars" not in json or not isinstance(json["vars"], list):
            raise Exception("Expected a list in 'layout.vars'")


        self.eqns = [Layout.parse_eqn(eq) for eq in json["eqns"]]

        variables = [Layout.parse_var(v) for v in json["vars"]]
        self._var_idxs = dict([(v, index) for index, v in enumerate(variables)])

    @staticmethod
    def parse_eqn(eq: str) -> Equation:
        def eq_1_vars(branch: Branch) -> list[Variable]:
            i = branch.index()

            vs = [(i, 'P'), (i, 'Q')]
            p = branch.parent()
            if p is not None:
                vs += [(p, 'P')]

            return vs

        def eq_2_vars(branch: Branch) -> list[Variable]:
            children = branch.children()
            if children is None:
                raise Exception(f"Equation 2 used for non-bifurcation branch #{branch.index()}")

            return [(branch.index(), 'Q'), (children[0], 'Q'), (children[1], 'Q')]

        def eq_3_vars(branch: Branch) -> list[Variable]:
            assert branch.children() is None

            i = branch.index()
            return [(i, 'Q'), (i, 'V')]

        def eq_4_vars(branch: Branch) -> list[Variable]:
            assert branch.children() is None

            i = branch.index()
            return [(i, 'P'), (i, 'V')]

        eq_fns = [None, eq_1_vars, eq_2_vars, eq_3_vars, eq_4_vars]

        try:
            branch_number, eq_number = [int(s) for s in eq.split('.')]

            assert 1 <= eq_number <= 4
            return Equation(branch_number, eq_fns[eq_number])
        except Exception as e:
            raise Exception(f'could not parse equation {eq}') from e

    @staticmethod
    def parse_var(var: str) -> Variable:
        branch_number, varname = var.split('.')
        assert varname in ['P', 'Q', 'V']
        assert int(branch_number) >= 0

        return (int(branch_number), varname)

    # Returns an iterator of the equations
    def equations(self) -> list[Equation]:
        return self.eqns

    # Returns a mapping from variable names to their index
    def variable_indexes(self) -> dict[Variable, int]:
        return self._var_idxs

def generate_values(model: Model, layout: Layout) -> list[list[int]]:
    equations = layout.equations()
    var_idxs = layout.variable_indexes()

    width = len(var_idxs)

    def row(eqn: Equation) -> list[int]:
        r = [0 for _ in range(width)]
        for v in eqn.variables(model):
            r[var_idxs[v]] = 1
        return r

    return [row(eqn) for eqn in equations]

def tikz_to_string(values: list[list[int]]) -> str:
    def pad(v: int) -> str:
        return f'|[{v}]|'

    def row(vs: list[int]) -> str:
        return '    ' + ' & '.join([pad(v) for v in vs]) + r' \\'

    comment = "% Adapted from https://tex.stackexchange.com/a/493402"
    begin_tikz = r"\begin{tikzpicture}[0/.style={draw,ultra thin},1/.style={0,fill=black}]"
    begin_matrix = r"\matrix[matrix of nodes,cells={minimum size=1.5em,anchor=center}]{"
    content = '\n'.join([row(vs) for vs in values])
    end="};\n\\end{tikzpicture}"

    return f"{comment}\n{begin_tikz}\n{begin_matrix}\n{content}\n{end}"

if __name__ == '__main__':
    main()
