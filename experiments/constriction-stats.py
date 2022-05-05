#!/bin/python
#
# Runs many models at various levels of constriction for long enough for it to even out

from runner import bin_hash, Runner
from hashlib import sha1
from pathlib import Path
import csv

def main():
    total_time = 24 # 20s to even out, +4s for a full breath cycle
    depth = 8
    constrictions = [i for i in range(96)]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/constriction-stats').mkdir(parents=False, exist_ok=True)
    Path(f'data/constriction-stats/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    count = len(constrictions)

    i = 0
    for c in constrictions:
        # Status indicator
        print(f'{i}/{count}...\r', end='', flush=True)

        run(template_hash, total_time, depth, c)
        i += 1

    print('done generating initial data')
    print('producing summary...')
    make_summary(template_hash, total_time, depth, constrictions, (20, 24))
    print('all jobs complete.', end='\n\n')

    print(f'check \'data/constriction-stats/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": 0.75,
        "branch_radius_decrease": 0.75,
        "split_angle": 0.8,
        "max_depth": {depth},
        "trachea_radius": {trachea_radius}
    }},
    "root": {{
        "type": "auto"
    }}
}}
"""

def run(template_hash: str, total_time: float, depth: int, constriction: int) -> None:
    filename = f'{total_time}-{depth}d@{int(constriction)}%'
    path = f'data/constriction-stats/{template_hash}/{filename}.csv'

    # Check if this trial has already been run
    if Path(path).exists():
        return

    cons_frac = 0.01 * constriction
    fmt_args = {
        'depth': depth,
        'trachea_radius': (1 - cons_frac) * 0.01, # trachea is 0.01 meters by default
    }
    config = TEMPLATE.format_map(fmt_args)

    output, exit_code = Runner(config, total_time).run()
    with open(path, 'w+') as f:
        f.write(output)

    if exit_code != 0:
        print(f'warning: trial {filename} exited with code {exit_code}')

def make_summary(template_hash: str, total_time: float, depth: int, constrictions: list[int],
        valid_time_range: tuple[float, float]) -> None:
    results = []

    constrictions.sort()
    assert constrictions[0] == 0

    for c in constrictions:
        filename = f'{total_time}-{depth}d@{c}%'
        path = f'data/constriction-stats/{template_hash}/{filename}.csv'

        # Read the data from this trial
        max_flow = 0
        max_time = None # time at which max_flow occured
        max_acceleration = 0
        last_pair = None

        with open(path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                time = float(row['time'])
                flow = float(row['flow out'])

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

    # Write the results to the summary file
    summary_path = f'data/constriction-stats/{template_hash}/summary.csv'
    with open(summary_path, 'w+') as f:
        base_time = results[0][2] # shift all the base times to be relative to constriction=0%

        f.write('constriction,max flow,max flow time,max acceleration,scaled max acceleration\n')

        for constriction, max_flow, max_time, max_acceleration in results:
            max_time -= base_time
            scaled = max_acceleration / max_flow

            f.write(f'{constriction},{max_flow:.4},{max_time:.4},{max_acceleration:.4},{scaled:.4}\n')

if __name__ == '__main__':
    main()
