#!/bin/python
#
# Runs a healthy model, but with a (possibly outlandish) starting volume

import base
from runner import bin_hash, Runner, run_parallel
from hashlib import sha1
from pathlib import Path

def main():
    total_time = 8
    depth = 10
    constrictions = [0, 0.4] # percent constricted, at the trachea and lower
    # Volumes ranging from 0.5L to 6L, in increments of 0.25L
    initial_volumes = [x * 0.25 for x in range(2, 24)]
    assert initial_volumes[0] == 0.5 and initial_volumes[-1] == 5.75

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/differ-start').mkdir(parents=False, exist_ok=True)
    Path(f'data/differ-start/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    runners = []
    for c in constrictions:
        for v in initial_volumes:
            r_tuple = make_runner(template_hash, total_time, depth, c, v)
            if r_tuple is not None:
                runners.append(r_tuple)

    if len(runners) != 0:
        results = run_parallel(runners)
        for returncode, filename in results:
            if returncode != 0:
                print(f'warning: trial {filename} exited with code {returncode}')

    print('done.', end='\n\n')
    print(f'check \'data/differ-start/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": """+str(base.BRANCH_LENGTH_DECREASE)+""",
        "branch_radius_decrease": """+str(base.BRANCH_RADIUS_DECREASE)+""",
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

def make_runner(template_hash: str, total_time: float, depth: int, constriction: float, initial_volume: float) -> tuple[Runner, str] | None:
    name = f'{initial_volume}L-{int(constriction*100)}%'
    path = f'data/differ-start/{template_hash}/{name}.csv'

    if Path(path).exists():
        return

    fmt_args = {
        'depth': depth,
        'trachea_radius': (1 - constriction) * base.TRACHEA_RADIUS, # trachea is 0.01 meters by default
        'initial_volume': initial_volume * 0.001, # convert L -> m³
    }
    config = TEMPLATE.format_map(fmt_args)
    return (Runner(config, total_time, path), name)

if __name__ == '__main__':
    main()
