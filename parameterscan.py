import numpy as np
import subprocess
import os

# -----------------------------------------------------------------------
# Parameter scan script for the electrostatics exercise.
#
# Default use: convergence study of phi(0) vs N (= N1 = N2) in the
# trivial (uniform) case.  Change 'paramstr' and 'variable_array' to
# scan any other parameter.
# -----------------------------------------------------------------------

# Path to compiled executable (adjust if needed)
repertoire     = ''
executable     = 'engine.exe'
input_filename = 'trivial.in'   # base configuration file

# Base parameters (values here are overwritten by the scan below)
input_parameters = {
    'b'      : 0.05,   # Inner radius [m]
    'R'      : 0.1,    # Outer radius [m]
    'V0'     : 0,      # Boundary potential [V]
    'a0'     : 1,      # Charge density scale [V/m^2]  (unused when trivial=true)
    'trivial': 'true', # true: uniform test case
    'N1'     : 5,      # Intervals in [0, b]
    'N2'     : 5,      # Intervals in [b, R]
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
paramstr       = 'N1'                        # parameter name in engine
variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

# Build a label for output directories / filenames
outstr = (f"electrostatics_b_{input_parameters['b']:.2g}"
          f"_R_{input_parameters['R']:.2g}"
          f"_trivial_{input_parameters['trivial']}")

# -----------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

# -----------------------------------------------------------------------
# Run the scan
# -----------------------------------------------------------------------
for val in variable_array:

    params = input_parameters.copy()
    params[paramstr] = val
    # For a convergence study keep N1 = N2
    if paramstr == 'N1':
        params['N2'] = val

    output_file = f"{outstr}_{paramstr}_{val}"
    output_path = os.path.join(outdir, output_file)

    # Build the command-line parameter string
    param_string = " ".join(f"{k}={v}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")
