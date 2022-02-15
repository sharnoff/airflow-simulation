# Helper module for running the airflow simulator
#
# To make sure that things are up to date, it calls `cargo build --release` on import. Subsequent
# calls to the program use `target/release/airflow-simulator`.
#
# The primary export is the `Runner` class.

import os
import sys
import subprocess
from hashlib import sha1

print('\x1b[34m::\x1b[0m building airflow-simulator binary...')

status = os.system('cargo build --release')
# if building failed, the script calling this should fail as well.
exit_code = os.waitstatus_to_exitcode(status)
if exit_code != 0:
    print("\x1b[31m::\x1b[0m cannot import runner.py; build failed")
    sys.exit(exit_code)

print('\x1b[34m::\x1b[0m done building')

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
    def __init__(self, config: str, total_time: float, timestep: float | None = None):
        self.config = config
        self.total_time = total_time
        self.timestep = timestep
        pass
    
    # Runs the program with all of the configuration provided
    def run(self) -> tuple[str, int]:
        MODEL_FILE = '.tmp-model.json'

        # Write the configuration to a temporary file
        with open(MODEL_FILE, 'w+') as f:
            f.write(self.config)

        cmd = ['target/release/airflow-simulator',
                '-m', 'from-json',
                '--input-file', MODEL_FILE,
                '--output', 'csv',
                '-t', str(self.total_time)]

        if self.timestep is not None:
            cmd += ['--timestep', str(self.timestep)]

        child = subprocess.run(cmd, stdout=subprocess.PIPE)

        # cleanup the temp file
        os.remove(MODEL_FILE)

        return (child.stdout.decode('utf-8'), child.returncode)
