#!/bin/python
#
# Runs many models at various levels of constriction for long enough for it to even out

from runner import bin_hash, Runner
from hashlib import sha1
from pathlib import Path

def main():
    total_time = 28 # 20s to even out, +8s for a couple full breaths
    depth = 6
    constrictions = [0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/constricted-flow').mkdir(parents=False, exist_ok=True)
    Path(f'data/constricted-flow/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    count = len(constrictions)

    i = 0
    for c in constrictions:
        # Status indicator
        print(f'{i}/{count}...\r', end='', flush=True)

        run(template_hash, total_time, depth, c)
        i += 1

    print('done.', end='\n\n')
    print(f'check \'data/constricted-flow/{template_hash}/*.csv\' for results')

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

def run(template_hash: str, total_time: float, depth: int, constriction: float) -> None:
    fmt_args = {
        'depth': depth,
        'trachea_radius': (1 - constriction) * 0.01, # trachea is 0.01 meters by default
    }
    config = TEMPLATE.format_map(fmt_args)

    output, exit_code = Runner(config, total_time).run()
    name = f'{total_time}-{depth}d@{int(constriction*100)}%'
    path = f'data/constricted-flow/{template_hash}/{name}.csv'
    with open(path, 'w+') as f:
        f.write(output)

    if exit_code != 0:
        print(f'warning: trial {name} exited with code {exit_code}')

if __name__ == '__main__':
    main()
