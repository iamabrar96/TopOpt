import numpy as np
from geometry import Rectangle_beam
from parameters import Parameters
from topopt import FE_solver

def test__comparision_with_analytical_soln():
    '''
    Aim :   The total error of Numerical solution (wrt. analytical solution) 
            should be within the tolerance provided.
    
    Expected Output :   Test passes successfully
    
    Remarks :   The case considered here is the famous Timoschenko beam problem 
                who's analytical solution is known.
                
                Test case passed successfully
    '''
    params= Parameters()


    #########################    PreProcesing    ##########################
    # select the geometry
    params.geometry_type= 'Rectangle_beam'

    # select the discretisation
    # length, width and height are chosen same as the no. of elements in the 
    # respective dimension. This is just a design choice.
    params.nelx=30
    params.nely=10
    params.nelz=1

    # select force b.c.
    params.force= np.array([[0, -1, 0]])
    
    # select the integration scheme
    params.integration_type="CompositeGauss4"   # 4th order gauss quadrature 

    params.volfrac= 1.0
    ######################    Perform Simulation    ########################
    # setup the solver with the parameters given above
    # Upon initialization this does three things outof the box. They are
    #   - Creates the geometry corresponding using gmsh
    #   - Initializes gauss quadrature (i.e nunts and weights)
    #   - Prepares System of Equations (i.e creates Constitutive(C), Strain-displacement(B), 
    #     Stiffness(K) and Force(f) matrices)
    solver= FE_solver(params)
    
    # solve the system of equations to get the numerical displacements
    U, _, _ = solver.solve()

    #########################    PostProcesing    #########################
    # since for this cantilever beam case there is no variation in displacement along 
    # z-dimension, lets just consider the nodes along the central axis of beam
    center_node_coords= np.zeros((params.n_dim, params.nelx+1), dtype= 'float')
    center_node_coords[0]= np.arange(0, params.nelx+1, dtype= 'float')  # x-coords
    U_exact= timoschenko__exact_diplacement(center_node_coords, params)

    # calculate the corresponding dofs
    center_node_tags= solver.geometry.centerNodeTags
    dof_center= solver.helper.getDofsForNodeTags(center_node_tags)[0].T
    U_num= U[dof_center,0]

    # indicate some absolute tolerance
    params.tol= 1e-3

    error= np.abs(U_exact[1]*8 -U_num[1]).mean()    #compare y-components of displacements
    assert error<params.tol


def timoschenko__exact_diplacement(x, params):
        x, y, z = x[0], x[1], x[2]
        
        length, width, height= params.nelx, params.nely, params.nelz
        P= params.force[0,1]
        c= P/(6*params.I* params.E0)
        
        u = np.zeros((params.n_dim, len(x)))
        u[0]= -c*y*(6*length-3*x)*x + (2+params.nu)*(y**2 -(width/2)**2)
        
        temp = 3*params.nu*(length-x)*y**2
        temp+= (4 + 5*params.nu)*x*(width/2)**2 + (3*length - x)*x**2
        u[1]= c*temp*8

        return u   

if __name__=='__main__':
    test__comparision_with_analytical_soln()