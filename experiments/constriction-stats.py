#!/bin/python
#
# Runs many models at various levels of constriction for long enough for it to even out

import base
from runner import bin_hash, Runner, run_parallel
from hashlib import sha1
from pathlib import Path
import csv

TIMESTEPS = [0.01, 0.0025]
DATADIR = 'data/constriction-stats'

def main():
    total_time = 28 # 20/24s to even out, +4s for a full breath cycle
    depth = 10
    constrictions = [i for i in range(91)] + [-i for i in range(1,51)]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path(DATADIR).mkdir(parents=False, exist_ok=True)
    Path(f'{DATADIR}/{template_hash}').mkdir(exist_ok=True)

    runners = []
    for timestep in TIMESTEPS:
        for c in constrictions:
            r_tuple = make_runner(template_hash, timestep, total_time, depth, c)
            if r_tuple is not None:
                runners.append(r_tuple)

    print(f'running {len(runners)} trials...')

    if len(runners) != 0:
        results = run_parallel(runners)
        for returncode, filename in results:
            if returncode != 0:
                print(f'warning: trial {filename} exited with code {returncode}')

    print('done generating initial data')
    if len(TIMESTEPS) == 1:
        print('producing summary...')
    else:
        print('producing summaries...')
    for timestep in TIMESTEPS:
        make_summary(template_hash, timestep, total_time, depth, constrictions, (20, 24))
    print('all jobs complete.', end='\n\n')

    print(f'check \'{DATADIR}/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": """+str(base.BRANCH_LENGTH_DECREASE)+""",
        "branch_radius_decrease": """+str(base.BRANCH_RADIUS_DECREASE)+""",
        "split_angle": 0.8,
        "max_depth": {depth},
        "trachea_radius": {trachea_radius}
    }},
    "root": {{
        "type": "auto"
    }}
}}
"""

def trial_filename(timestep: float, total_time: float, depth: int, constriction: int) -> str:
    return f'{timestep}x{total_time}s-{depth}d@{constriction}%'

def make_runner(template_hash: str, timestep: float, total_time: float, depth: int,
        constriction: int) -> tuple[Runner, str] | None:
    filename = trial_filename(timestep, total_time, depth, constriction);
    path = f'{DATADIR}/{template_hash}/{filename}.csv'

    # Check if this trial has already been run
    if Path(path).exists():
        return

    cons_frac = 0.01 * constriction
    fmt_args = {
        'depth': depth,
        'trachea_radius': (1 - cons_frac) * base.TRACHEA_RADIUS, # trachea is 0.01 meters by default
    }
    config = TEMPLATE.format_map(fmt_args)
    return (Runner(config, total_time, path, timestep), filename)

def make_summary(template_hash: str, timestep: float, total_time: float, depth: int,
        constrictions: list[int], valid_time_range: tuple[float, float]) -> None:
    results = []

    assert constrictions[0] == 0

    for c in constrictions:
        filename = trial_filename(timestep, total_time, depth, c)
        path = f'{DATADIR}/{template_hash}/{filename}.csv'

        # Read the data from this trial
        max_flow = 0
        max_time = None # time at which max_flow occured
        max_acceleration = 0
        last_pair = None

        with open(path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                time = float(row['time'])
                flow = float(row['flow in'])

                if valid_time_range[0] <= time < valid_time_range[1]:
                    if last_pair is not None:
                        acc = (flow - last_pair['flow']) / (time - last_pair['time'])
                        if acc > max_acceleration:
                            max_acceleration = acc

                    if flow > max_flow:
                        max_flow = flow
                        max_time = time

                    pass

                last_pair = { 'time': time, 'flow': flow }

                pass
        
        assert 0 not in [max_flow, max_acceleration]
        assert max_time is not None

        # shift max_time back by the start of valid_time_range so that it starts at zero, not
        # whatever temporal offset we're using
        max_time -= valid_time_range[0]

        results.append((c, max_flow, max_time, max_acceleration))

    base_time = results[0][2] # remember: results[0] is still 0% constriction
    results.sort(key=lambda r: r[0]) # ... but now it might not be

    # Write the results to the summary file
    summary_path = f'{DATADIR}/{template_hash}/summary-{timestep}x{total_time}-{depth}d.csv'
    with open(summary_path, 'w+') as f:
        f.write('constriction,max flow,max flow time,max acceleration,scaled max acceleration\n')

        for constriction, max_flow, max_time, max_acceleration in results:
            # center max_time so that it's within a 2s window around base_time
            if max_time > base_time + 2:
                max_time -= 4
            elif max_time < base_time - 2:
                max_time += 4

            max_time -= base_time

            scaled = max_acceleration / max_flow

            f.write(f'{constriction},{max_flow:.4},{max_time:.4},{max_acceleration:.4},{scaled:.4}\n')

if __name__ == '__main__':
    main()
