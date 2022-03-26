from dataclasses import dataclass
import numpy as np

@dataclass
class Parameters:
    """
    Parameters necessary for running the simulation.
    Attributes
    ----------

    n_dim : int
        number of dimension of the problem
        Default: 3

    nu : float
        poison ratio
        Default : 0.3
    
    nelx: int
        number of element in x direction
    nely : int
        number of element in y direction
    nelz : int
        number of element in z direction
    
    volfrac : float                     
        Volume Fraction Limit   
        Default : 0.3 

    p : int 
        Penalisation power for SIMP interpolaion     
         
    rmin : float
        Minimum size of the selected domain 
        Default : 1.5

    E0 : np.array   
        Youngs modulus of solid material                              
        
    Emin : float                        
        Youngs modulus of void material 
        Default : 1e-9 

  
    max_loop : int
        maximum nuber of iterations 

    tol =0.001
        termination criteria 

    force : np.array
        Force Boundary condition  
        Default : np.array([0, -1, 0])
    
    num_load : int                  
        Specify which load_case. 1 if single load, n if n multiple loads   
    
    integration_type : String
        Select the Gauss quadrature used for integration
        Default : 'CompositeGauss4'

    geometry_type : String
        Select the Geometry for simulation. Currently supported geometries are  
            - Rectangle_beam (default)
            - Michell_beam 
            - Mid_cantilever 
            - multiple_load_case  
    
    density_cutoff : float
        Minimum density value for elements to be considered in the final topology
        Default : 0.5 

    filter : int 
        Select filter types. Currently supported filter types are
            - 1 : density filter (default)
            - 2 : sensitivity filter
            - 3 : grey scale filter

    node_per_ele : int 
        number of nodes per element
        Default : 8
    
    dof_per_node : int
        degrees of freedom  per node
        Default: 3

    dof_per_ele : int
        degrees of freedom per element

    num_elems : int  
        Total number of elements 
    tnodes : int
        Total number of nodes      
    tdof : int
        Total number of degrees of freedom
    """

    n_dim=3
    nu=0.3

    nelx=30
    nely=10
    nelz=4

    volfrac=0.3                      
    p=3           
    E0=1   
    Emin=1e-9                         
    rmin=1.5

    max_loop=60
    tol =0.001

    geometry_type= "multiple_load_case"
    # force = np.array([[0, -1, 0]])
    force = np.array([[0, -1, 0],[0, 1, 0]]) #for multiple_load_case
    num_load= len(force)            
    
    integration_type="CompositeGauss4"
    density_cutoff= 0.4
    filter=2 
    
    node_per_ele=8
    dof_per_node=3
    dof_per_ele= node_per_ele*dof_per_node  
    
    num_elems=nelx*nely*nelz  
    tnodes=(nelx+1)*(nely+1)*(nelz+1)
    tdof=dof_per_node*tnodes      