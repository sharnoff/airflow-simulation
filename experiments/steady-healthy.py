#!/bin/python
#
# Runs a model without any degradation for a duration.

from runner import bin_hash, Runner
from hashlib import sha1
from pathlib import Path

def main():
    total_times = [1004] # +4 so that we include something _starting_ at 1k
    depths = [6]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    # Ensure that our directory is there:
    Path('data/steady-healthy').mkdir(parents=False, exist_ok=True)
    Path(f'data/steady-healthy/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    count = len(total_times) * len(depths)

    i = 0
    for t in total_times:
        for d in depths:
            # Status indicator
            print(f'{i}/{count}...\r', end='', flush=True)

            run(template_hash, t, d)
            i += 1

    print('done.', end='\n\n')
    print(f'check \'data/steady-healthy/{template_hash}/*.csv\' for results')

TEMPLATE = """
{{
    "config": {{
        "branch_length_decrease": 0.75,
        "branch_radius_decrease": 0.75,
        "split_angle": 0.8,
        "max_depth": {depth}
    }},
    "root": {{
        "type": "auto"
    }}
}}
"""

def run(template_hash: str, total_time: float, depth: int) -> None:
    fmt_args = {
        'depth': depth,
    }
    config = TEMPLATE.format_map(fmt_args)

    output, exit_code = Runner(config, total_time).run()
    name = f'{total_time}@{depth}d'
    path = f'data/steady-healthy/{template_hash}/{name}.csv'
    with open(path, 'w+') as f:
        f.write(output)

    if exit_code != 0:
        print(f'warning: trial {name} exited with code {exit_code}')

if __name__ == '__main__':
    main()
