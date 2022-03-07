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
#   $TEMPLATE_HASH/$TOTAL_TIME-$ONSET-$METHOD:$CHANGE-$UNHEALTHY@$UNHEALTHY_RADIUS

from runner import bin_hash, Runner
from hashlib import sha1
from pathlib import Path

def main():
    total_times = [16]
    onsets = [4 + 0.5*i for i in range(9)] # 4, 4.5, ... 8
    methods = ['linear','tanh','cubic','trig']
    changes = [4]
    unhealthies = [4]
    unhealthy_radii = [0.075]

    hasher = sha1()
    hasher.update(bin_hash())
    hasher.update(TEMPLATE.encode('utf-8'))
    template_hash = hasher.hexdigest()

    count = len(total_times) * len(onsets) * len(methods) * len(changes) * len(unhealthies) * len(unhealthy_radii)

    # Ensure that our directory is there:
    Path('data/degrade-and-return').mkdir(parents=False, exist_ok=True)
    Path(f'data/degrade-and-return/{template_hash}').mkdir(exist_ok=True)

    print('running...')

    # Most of these lists will always only have a single item, but it's still worth doing the
    # iteration like this.
    i = 0
    for t in total_times:
        for o in onsets:
            for m in methods:
                for c in changes:
                    for u in unhealthies:
                        for r in unhealthy_radii:
                            # Status indicator
                            print(f'{i}/{count}...\r', end='', flush=True)

                            run(template_hash, t, o, m, c, u, r)
                            i += 1

    print('done.', end='\n\n')
    print(f'check \'data/degrade-and-return/{template_hash}/*.csv\' for results')

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
		"branch_length_decrease": 0.75,
		"branch_radius_decrease": 0.75,
		"split_angle": 0.8,
		"max_depth": 6
	}},
	"root": {{
		"type": "manual",
		"left_override": {{ "relative_angle": 0.7, "relative_radius_abnormal": {radius} }},
		"left": {{
			"type": "auto",
			"child_angles": [1.2, 0.3, {{ "reset_after": 2 }}],
			"branch_length_decrease": [0.65, 0.9, {{ "reset_after": 1 }}]
		}},
		"right_override": {{ "relative_angle": 0.7, "relative_radius_abnormal": {radius} }},
		"right": {{
			"type": "auto",
			"child_angles": [0.3, 1.2, {{ "reset_after": 2 }}],
			"branch_length_decrease": [0.9, 0.65, {{ "reset_after": 1 }}]
		}}
	}}
}}
"""

def run(template_hash: str, total_time: float, onset: float, method: str, change: float,
        unhealthy: float, unhealthy_radius: float) -> None:
    fmt_args = {
        'interpolate_method': method,
        'start_degrade': onset,
        'start_degraded': onset + change,
        'end_degraded': onset + change + unhealthy,
        'end_degrade': onset + change + unhealthy + change,
        'radius': unhealthy_radius,
    }
    config = TEMPLATE.format_map(fmt_args)

    output, exit_code = Runner(config, total_time).run()
    name = f'{total_time}-{onset}-{method}:{change}-{unhealthy}@{unhealthy_radius}'
    path = f'data/degrade-and-return/{template_hash}/{name}.csv'
    with open(path, 'w+') as f:
        f.write(output)

    if exit_code != 0:
        print(f'warning: trial {name} exited with code {exit_code}')

if __name__ == '__main__':
    main()
