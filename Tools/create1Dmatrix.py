import numpy as np
from Tools import run_cmd

def create_1D_matrix(filename, matrix,diary_file):
    if matrix.shape != (3, 4):
        nl = 'ERROR: Matrix must be 3x4 for a valid AFNI transformation.'
        raise Exception(run_cmd.error(nl, diary_file))

    # Flatten the matrix to 1D and save it
    flat_matrix = matrix.flatten()
    np.savetxt(filename, [flat_matrix], fmt='%.6f')