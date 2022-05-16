#!/bin/python
#
# Runs a healthy model, then constricts one side for a duration & lets it recover.

import base
from runner import bin_hash, Runner, run_parallel
import csv
from hashlib import sha1
from pathlib import Path

DATADIR = 'data/asymmetric'

def main():
    total_time = 28
    depth = 10
    onsets = [8]
    method = 'tanh'
    transition_time = 4 # number of seconds to transition on either end
    stay_degraded_for = 4 # stays degraded for 4 seconds
    constrictions = [i for i in range(91)] # % to constrict by

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path(DATADIR).mkdir(parents=False, exist_ok=True)
    Path(f'{DATADIR}/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    runners = []
    for start in onsets:
        for constrict_by in constrictions:
            r_tuple = make_runner(template_hash, depth, total_time, method, start, transition_time,
                    stay_degraded_for, constrict_by)
            if r_tuple is not None:
                runners.append(r_tuple)

    if len(runners) != 0:
        results = run_parallel(runners)
        for returncode, filename in results:
            if returncode != 0:
                print(f'warning: trial {filename} exited with code {returncode}')

    print('done generating initial data')
    print('producing summary...')
    for onset in onsets:
        make_summary(template_hash, depth, total_time, method, onset, transition_time,
                stay_degraded_for, constrictions)

    print('all jobs complete.', end='\n\n')

    print(f'check \'{DATADIR}/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
	"schedule": {{
		"interpolate": "{interpolate_method}",
		"keyframes": [
			{{ "time": 0, "is_degraded": false }},
			{{ "time": {start_degrade}, "is_degraded": false }},
			{{ "time": {start_degraded}, "is_degraded": true }},
			{{ "time": {end_degraded}, "is_degraded": true }},
			{{ "time": {end_degrade}, "is_degraded": false }}
		]
	}},
	"config": {{
		"branch_length_decrease": """+str(base.BRANCH_LENGTH_DECREASE)+""",
		"branch_radius_decrease": """+str(base.BRANCH_RADIUS_DECREASE)+""",
        "split_angle": 0.8,
		"max_depth": {depth}
    }},
    "root": {{
        "type": "manual",
        "left_override": {{ "relative_radius_abnormal": {constricted_radius_decrease} }}
	}}
}}
"""

def trial_filename(depth: int, total_time: int, method: str, onset: float, transition_time: float,
        degraded_for: float, constrict_by: int) -> str:
    return f'{depth}d@{constrict_by}%-{total_time}s-{method}:{onset:.2f}s-{transition_time}s-{degraded_for}s.csv'

def make_runner(template_hash: str, depth: int, total_time: int, method: str, onset: float,
        transition_time: float, degraded_for: float, constrict_by: int) -> tuple[Runner, str] | None:
    filename = trial_filename(depth, total_time, method, onset, transition_time, degraded_for,
            constrict_by)
    path = f'{DATADIR}/{template_hash}/{filename}'

    if Path(path).exists():
        return None

    fmt_args = {
        'interpolate_method': method,
        'start_degrade': onset,
        'start_degraded': onset + transition_time,
        'end_degraded': onset + transition_time + degraded_for,
        'end_degrade': onset + transition_time + degraded_for + transition_time,
        'depth': depth,
        'constricted_radius_decrease': base.BRANCH_RADIUS_DECREASE * (1 - 0.01*constrict_by),
    }
    config = TEMPLATE.format_map(fmt_args)
    return (Runner(config, total_time, path, additional_subtrees=['.left','.right']), filename)

def make_summary(template_hash: str, depth: int, total_time: int, method: str, onset: float,
        transition_time: float, degraded_for: float, constrictions: list[int]) -> None:
    if constrictions[0] != 0:
        raise Exception('expected constrictions[0] == 0')

    def get_flows(constriction: int) -> list[float]:
        filename = trial_filename(depth, total_time, method, onset, transition_time,
                degraded_for, constriction)
        per_tick_flow = []
        with open(f'{DATADIR}/{template_hash}/{filename}', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                flow = float(row['.right flow in'])
                per_tick_flow.append(flow)
        return per_tick_flow

    flows = [(c, get_flows(c)) for c in constrictions]
    _, zero_flows = flows[0] # because constrictions[0] == 0
    zero_max = max(zero_flows)

    per_constriction_results = []
    for c, fs in flows:
        max_flow__diff = max(fs) - zero_max
        max__flow_diff = max([abs(t - z) for t, z in zip(fs, zero_flows)])
        per_constriction_results.append((c, max_flow__diff, max__flow_diff))

    per_constriction_results.sort(key = lambda t: t[0])

    summary_filename = f'summary-{depth}d-{total_time}s-{method}:{onset:.2f}s-{transition_time}s-{degraded_for}s.csv'
    with open(f'{DATADIR}/{template_hash}/{summary_filename}', 'w+') as f:
        f.write('constriction,(max flow) difference,max (flow difference)\n')

        for c, m1, m2 in per_constriction_results:
            f.write(f'{c},{m1:.6f},{m2:.6f}\n')

if __name__ == '__main__':
    main()
