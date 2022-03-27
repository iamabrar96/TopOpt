from dataclasses import dataclass
import numpy as np

@dataclass
class Parameters:
    """
    Parameters necessary for running the simulation.
    Attributes
    ----------
    

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

    n_dim : int
        number of dimension of the problem
        Default: 3
    n_out : int
        number of dimension of the output: density function
        Default: 1
    epoch : int
        number of epochs for the model to run
        Default: 50
    size_hidden : int
        size of the hidden layer
        Default: 500
    num_hidden_layer : int
        number of hidden layer in the neural network
        Default: 3
    seed_item: int
        initialisation for manual seeds for reproducable results
    """

    nu=0.3

    nelx=30
    nely=10
    nelz=2

    volfrac=0.3                      
    p=3           
    E0=1   
    Emin=1e-9                         
    rmin=1.5

    max_loop=60
    tol =0.001

    geometry_type= "Rectangle_beam"
    force = np.array([[0, -1, 0]])              #for Rectangle_beam
    # force = np.array([[0, -1, 0],[0, 1, 0]])  #for multiple_load_case
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

    n_dim = 3
    out_dim =1
    epoch = 50
    size_hidden = 500
    num_hidden_layer = 3
    seed_item = 9867985
    optimizer_name = 'Adam'
    dropout = 0.20
    I= (nelz* nely**3)/12    

    c2= np.array(  [1.0, 18.45074218, 12.28496332, 10.18983032,  9.20225366,  8.67839379,
                    8.37293717,  8.18218208,  8.05635729,  7.96962493,  7.90764666,  7.86201934,
                    7.82757955,  7.80102058,  7.78014569 , 7.76344685 , 7.74985785,  7.73860515,
                    7.7291165,  7.720966,    7.71384449,  7.70755057,  7.7020022 ,  7.69727332,
                    7.69367585,  7.691901,    7.693454,    7.70106127,  7.72390748,  7.77382679,
                    8.01467146]) 