# Helper module for running the airflow simulator
#
# To make sure that things are up to date, it calls `cargo build --release` on import. Subsequent
# calls to the program use `target/release/airflow-simulator`.
#
# The primary export is the `Runner` class.

import asyncio
import os
import sys
from hashlib import sha1
from typing import TypeVar, Any

print('\x1b[34m::\x1b[0m building airflow-simulator binary...')

status = os.system('RUSTFLAGS="-C target-cpu=native" cargo build --release')
# if building failed, the script calling this should fail as well.
exit_code = os.waitstatus_to_exitcode(status)
if exit_code != 0:
    print("\x1b[31m::\x1b[0m cannot import runner.py; build failed")
    sys.exit(exit_code)

print('\x1b[34m::\x1b[0m done building')

os.nice(1) # make this process a little bit lower-priority

_build_hash = None

# Returns the SHA-1 hash of the build binary
def bin_hash() -> bytes:
    global _build_hash

    if _build_hash is not None:
        return _build_hash

    hasher = sha1()
    with open('target/release/airflow-simulator', 'br') as f:
        hasher.update(f.read())
    _build_hash = hasher.digest()
    return _build_hash

class Runner:
    def __init__(self, config: str, total_time: float, out_file: str, timestep: float | None = None,
            additional_subtrees: list[str] | None = None):
        self.config = config
        self.total_time = total_time
        self.out_file = out_file
        self.timestep = timestep
        self.subtrees = additional_subtrees

T = TypeVar('T')

# Runs all "trials" in parallel, using at most numproc concurrent subprocesses
#
# Results are equivalent to [t.run() for t in trials].
def run_parallel(trials: list[tuple[Runner, T]], numproc: int | None = None) -> list[tuple[int, T]]:
    if numproc is None:
        if (numproc_env := os.getenv('MAX_PROCS')) is not None:
            if numproc_env.isnumeric():
                numproc = int(numproc_env)
            else:
                raise Exception(f'environment variable MAX_PROCS must be integer (got "{numproc_env}")')
        else:
            numproc = os.cpu_count()
            if numproc is None:
                print('could not determine number of CPUs, defaulting to running 4x parallel')
                numproc = 4
    elif numproc < 1:
        raise Exception('numproc must be >= 1')
    
    results: list[tuple[int, T] | None] = [None for _ in range(len(trials))]

    sema = asyncio.Semaphore(numproc)
    completed = [0] # use a list so that we pass a reference

    async def run(completed: list[int], idx: int):
        async with sema:
            this, passthru = trials[idx]

            if os.path.exists(this.out_file):
                raise Exception(f'cannot run: file {this.out_file} already exists')

            model_file = f'.tmp-model.{idx}.json'

            with open(model_file, 'w+') as f:
                f.write(this.config)

            cmd = ['target/release/airflow-simulator',
                    '-m', 'from-json',
                    '--input-file', model_file,
                    '--output', 'csv',
                    '-t', str(this.total_time)]
            if this.timestep is not None:
                cmd += ['--timestep', str(this.timestep)]
            if this.subtrees is not None:
                for subtree_path in this.subtrees:
                    cmd += ['-a', subtree_path]

            with open(this.out_file, 'w+') as f:
                child = await asyncio.create_subprocess_exec(*cmd, stdout=f)
                return_code = await child.wait()

            # regardless of result, cleanup the model file
            os.remove(model_file)

            results[idx] = (return_code, passthru)
            completed[0] += 1
            print(f'\r{completed[0]}/{len(trials)} done', end='')

    print(f'0/{len(trials)} done', end='', flush=True)
    tasks = asyncio.gather(*[run(completed, i) for i in range(len(trials))])

    loop = asyncio.get_event_loop()
    try:
        loop.run_until_complete(tasks)
    finally:
        loop.run_until_complete(loop.shutdown_asyncgens())
        loop.close()

    print('') # add a newline

    actual_results: Any = results # escape hatch for type checking.
    return actual_results
