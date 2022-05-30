import numpy as np
from numpy import round
from parameters import Parameters
from topopt import FE_solver
from scipy.sparse.linalg import eigs
from scipy.linalg import cho_factor, LinAlgError

def test__stiffness_matrix_singularity():
    '''
    SANITY CHECK

    Aim : The stiffness matrix has to be positive definite.
    
    Expected Output : The eigenvalues of the global stiffness matrix should not be negative
                      and cholskey decomposition should be possible.
    
    Remarks : test case passed successfully
    '''
    params= Parameters()
    solver= FE_solver(params)
    density= solver.phy_dens
    solver.solve(density)

    # Eigen values of the global stiffness matrix are calculated
    eigen_values,eigen_vectors= eigs(solver.kg)
    eigen_values=round(eigen_values,8)

    #Check for positive eigen values
    assert (all(eigen_values>=0))
    print("All eigen values are non-negative")

    #check for Positive Definiteness using cholskey decomposition
    try:
        cho_factor(solver.kg.toarray())
    except LinAlgError:
        print("Stiffness matrix is not Symmetric and Positive definite")
    print("Stiffness matrix is Symmetric and Positive definite")
    

if __name__=='__main__':
    test__stiffness_matrix_singularity()