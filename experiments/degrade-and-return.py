#!/bin/python
#
# Runs the following experiment:
#
#  1. Start with healthy lungs
#  2. After some time, slowly constrict the radius of the two bronchi (and all their children)
#  3. Let the unhealthy lungs sit for a bit
#  4. After some time, slowly bring the lungs back to healthy
#
# The experiment is configured by a set of values below, and the resulting data is placed into the
# 'data/degrade-and-return/' directory in the repository root, using file names with the following
# template:
#
#   {TEMPLATE_HASH}/{DEPTH}d@{UNHEALTHY}%-{TOTAL_TIME}s-{METHOD}:{ONSET}s-{CHANGE}s-{UNHEALTHY}s

import base
from runner import bin_hash, Runner, run_parallel
import csv
from hashlib import sha1
from pathlib import Path

DATADIR = 'data/degrade-and-return'

def main():
    total_time = 10
    depth = 10
    onsets = [6.1]
    method = 'tanh'
    transition_time = 0.1
    stay_degraded_for = 1.5
    constrict_by = 80

    # total_time = 28
    # depth = 10
    # onsets = [8 + 0.05*i for i in range(81)] # 8, 8.05, 8.1, ... 12
    # method = 'tanh' # possible: ['linear','tanh','cubic','trig']
    # transition_time = 4 # number of seconds to transition on either end
    # stay_degraded_for = 4 # stays degraded for 4 seconds
    # constrict_by = 50 # 50%

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
    make_summary(template_hash, depth, total_time, method, onsets, transition_time,
            stay_degraded_for, constrict_by)
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
		"max_depth": {depth},
        "trachea_radius": """+str(base.TRACHEA_RADIUS)+""",
        "trachea_radius_abnormal": {abnormal_trachea_radius}
    }},
    "root": {{
        "type": "auto"
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
        'abnormal_trachea_radius': base.TRACHEA_RADIUS * (1 - 0.01*constrict_by),
    }
    config = TEMPLATE.format_map(fmt_args)
    return (Runner(config, total_time, path), filename)

def make_summary(template_hash: str, depth: int, total_time: int, method: str, onsets: list[float],
        transition_time: float, degraded_for: float, constrict_by: int) -> None:
    results = []

    for onset in onsets:
        start_recovery = onset + transition_time + degraded_for
        end_recovery = start_recovery + transition_time

        filename = trial_filename(depth, total_time, method, onset, transition_time, degraded_for,
                constrict_by)
        path = f'{DATADIR}/{template_hash}/{filename}'

        max_flow = None # maximum MAGNITUDE flow -- may be negative or positive
        max_time = None # time at which max flow occured
        max_constriction = None # constriction level at which max flow occured

        with open(path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                time = float(row['time'])
                flow = float(row['flow in'])

                if time < start_recovery:
                    continue
                elif time > end_recovery:
                    break

                if max_flow is None or abs(flow) > abs(max_flow):
                    max_flow = flow
                    max_time = time
                    max_constriction = 0.01 * (constrict_by * float(row['degradation ratio']))

        if max_flow is None:
            raise Exception(f"trial {filename}: simulation didn't include recovery time range")

        results.append((onset, max_flow, max_time, max_constriction))

    summary_filename = f'summary-{depth}d@{constrict_by}%-{total_time}s-{method}:-{transition_time}s-{degraded_for}s.csv'
    summary_path = f'{DATADIR}/{template_hash}/{summary_filename}'

    with open(summary_path, 'w+') as f:
        f.write('onset,max flow,abs max flow,max flow time,constriction at max flow\n')

        for onset, max_flow, max_time, max_constriction in results:
            # scale onset and max_time so that they're relative to the breath cycles:
            onset %= 4.0
            max_time %= 4.0

            f.write(f'{onset:.2f},{max_flow},{abs(max_flow)},{max_time:.2f},{max_constriction}\n')

if __name__ == '__main__':
    main()
