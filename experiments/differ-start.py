#!/bin/python
#
# Runs a healthy model, but with a (possibly outlandish) starting volume

from runner import bin_hash, Runner
from hashlib import sha1
from pathlib import Path

def main():
    total_time = 8
    depth = 6
    constrictions = [0, 0.75, 0.8, 0.85, 0.9] # percent constricted, at the trachea and lower
    # Volumes ranging from 0.5L to 4L, in increments of 0.25L
    initial_volumes = [x * 0.25 for x in range(2, 17)]
    assert initial_volumes[0] == 0.5 and initial_volumes[-1] == 4.0

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/differ-start').mkdir(parents=False, exist_ok=True)
    Path(f'data/differ-start/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    count = len(constrictions) * len(initial_volumes)

    i = 0
    for c in constrictions:
        for v in initial_volumes:
            # Status indicator
            print(f'{i}/{count}...\r', end='', flush=True)

            run(template_hash, total_time, depth, c, v)
            i += 1

    print('done.', end='\n\n')
    print(f'check \'data/differ-start/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": 0.75,
        "branch_radius_decrease": 0.75,
        "split_angle": 0.8,
        "max_depth": {depth},
        "trachea_radius": {trachea_radius}
    }},
    "env": {{
        "initial_volume": {initial_volume}
    }},
    "root": {{
        "type": "auto"
    }}
}}
"""

def run(template_hash: str, total_time: float, depth: int, constriction: float, initial_volume: float) -> None:
    fmt_args = {
        'depth': depth,
        'trachea_radius': (1 - constriction) * 0.01, # trachea is 0.01 meters by default
        'initial_volume': initial_volume * 0.001, # convert L -> m³
    }
    config = TEMPLATE.format_map(fmt_args)

    output, exit_code = Runner(config, total_time).run()
    name = f'{initial_volume}L-{int(constriction*100)}%'
    path = f'data/differ-start/{template_hash}/{name}.csv'
    with open(path, 'w+') as f:
        f.write(output)

    if exit_code != 0:
        print(f'warning: trial {name} exited with code {exit_code}')

if __name__ == '__main__':
    main()
