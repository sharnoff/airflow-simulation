#!/bin/python
#
# Matches the timestep-by-timestep constriction in a degrade-and-return trial with the steady data
# on varied constriction levels. Linearly interpolates data from neighboring constriction levels.
#
# Places the result in the same directory, at 'matching-{original_filename}'

RECOVERY_DATADIR = 'data/degrade-and-return/3ed8c5cf2279df60ea4b7b353a8c91b6aa202a51'
CONSTRICTION_DATDIR = 'data/constriction-stats/d2b951ba10c28e9cfafdb3a67574bd882a649ad3'

ORIGINAL_FILENAMES = ['10d@50%-28s-tanh:8.00s-4s-4s.csv']

# parameters for the constriction runs:
TIMESTEP = '0.01'
TOTAL_TIME = 28

import csv
import re
import math

def main():
    cache = Cache(CONSTRICTION_DATDIR, TIMESTEP, TOTAL_TIME)

    for filename in ORIGINAL_FILENAMES:
        run(cache, filename)

def run(cache: 'Cache', filename: str):
    # extract some bits of information from the filename
    pat = re.compile(r'^(?P<depth>\d+)d@(?P<constriction>\d+)%.*$')
    if (m := pat.match(filename)) is None:
        raise Exception(f'malformed filename {filename}')

    depth = int(m.group('depth'))
    constriction = int(m.group('constriction'))

    # loop through all the lines in the file & generate the "expected" state at each time
    results = []
    with open(f'{RECOVERY_DATADIR}/{filename}', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            time = row['time']
            c_frac = float(row['degradation ratio'])
            # c_frac = 0 is at constriction = 0, so we remove the "(1 - c_frac) * 0" term.
            actual_constriction = c_frac * constriction

            lower_c = math.floor(actual_constriction)
            upper_c = math.ceil(actual_constriction)
            frac = actual_constriction - lower_c # fractional piece of the constriction

            lower_state = cache.get(depth, lower_c, time)
            upper_state = cache.get(depth, upper_c, time)
            lerped = lower_state.lerp(upper_state, frac)

            results.append((time, lerped))

    results_path = f'{RECOVERY_DATADIR}/matched-{filename}'
    with open(results_path, 'w+') as f:
        f.write('time,flow,volume\n')

        for time, state in results:
            f.write(f'{time},{state.flow},{state.volume}\n')
    
    print(f'output calculation to {results_path}')

class InstantState:
    def __init__(self, volume: float, flow: float):
        self.volume = volume
        self.flow = flow

    def lerp(self, other: 'InstantState', self_bias: float) -> 'InstantState':
        print(self_bias)
        v = self.volume * self_bias + other.volume * (1 - self_bias)
        f = self.flow * self_bias + other.flow * (1 - self_bias)
        return InstantState(v, f)

class Cache:
    def __init__(self, constriction_datadir: str, timestep: str, total_time: int):
        self.constriction_datadir = constriction_datadir
        self.timestep = timestep
        self.total_time = total_time
        # filename -> "time" -> state
        self.cached: dict[str, dict[str, InstantState]] = {}

    def get(self, depth: int, constriction: int, time: str) -> InstantState:
        fn = self.filename(depth, constriction)

        if (file := self.cached.get(fn)) is None:
            file = self.cache(fn)

        if (state := file.get(time)) is None:
            raise Exception(f'time {time} not found in file {fn}')

        return state

    # returns the filename associated with the given depth & constriction level
    def filename(self, depth: int, constriction: int) -> str:
        return f'{self.timestep}x{self.total_time}s-{depth}d@{constriction}%.csv'

    # read 'filename' and cache all of its values
    def cache(self, filename: str) -> dict[str, InstantState]:
        path = f'{self.constriction_datadir}/{filename}'
        
        mapping = {}

        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                time = row['time']
                volume = float(row['total volume'])
                flow = float(row['flow in'])
                mapping[time] = InstantState(volume, flow)

        self.cached[filename] = mapping
        return mapping

if __name__ == '__main__':
    main()
