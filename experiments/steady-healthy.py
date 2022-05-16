#!/bin/python
#
# Runs a model without any degradation for a duration.

import base
from runner import bin_hash, Runner, run_parallel
from hashlib import sha1
from pathlib import Path

def main():
    total_times = [1004] # +4 so that we include something _starting_ at 1k
    depths = [10]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/steady-healthy').mkdir(parents=False, exist_ok=True)
    Path(f'data/steady-healthy/{template_hash}').mkdir(exist_ok=True)

    runners = []
    for t in total_times:
        for d in depths:
            r_tuple = make_runner(template_hash, t, d)
            if r_tuple is not None:
                runners.append(r_tuple)

    print(f'running {len(runners)} trials...')

    if len(runners) != 0:
        results = run_parallel(runners)
        for returncode, filename in results:
            if returncode != 0:
                print(f'warning: trial {filename} exited with code {returncode}')

    print('done.', end='\n\n')
    print(f'check \'data/steady-healthy/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": """+str(base.BRANCH_LENGTH_DECREASE)+""",
        "branch_radius_decrease": """+str(base.BRANCH_RADIUS_DECREASE)+""",
        "split_angle": 0.8,
        "max_depth": {depth}
    }},
    "root": {{
        "type": "auto"
    }}
}}
"""

def make_runner(template_hash: str, total_time: float, depth: int) -> tuple[Runner, str] | None:
    name = f'{total_time}@{depth}d'
    path = f'data/steady-healthy/{template_hash}/{name}.csv'

    if Path(path).exists():
        return

    fmt_args = {
        'depth': depth,
    }
    config = TEMPLATE.format_map(fmt_args)
    return (Runner(config, total_time, path), name)

if __name__ == '__main__':
    main()
