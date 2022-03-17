from dataclasses import dataclass


@dataclass
class Parameters:
    ''' number of nodes per element'''
    node_per_ele=8
    '''  degrees of freedom  per node'''
    dof_per_node=3
    ''' degrees of freedom per element '''
    dof_per_ele= node_per_ele*dof_per_node
    ''' number of dimension of the problem'''
    n_dim=3
    ''' poison ratio '''
    nu=0.3
    ''' number of element in x direction'''
    nelx=30
    ''' number of element in y direction'''
    nely=10
    ''' number of element in z direction '''
    nelz=2
    ''' volume fraction limit '''
    volfrac=0.3                      
    ''' penalisation power for SIMP interpolaion ''' 
    p=3           
    ''' minimum size of the selected domain '''                    
    rmin=1.5
    '''  youngs modulus of solid material'''                           
    E0=1   
    '''  youngs modulus of void material '''                          
    Emin=1e-9                         
    '''Total number of elements'''
    num_elems=nelx*nely*nelz  
    ''' Total number of nodes'''   
    tnodes=(nelx+1)*(nely+1)*(nelz+1)
    ''' Total number of degrees of freedom'''
    tdof=dof_per_node*tnodes
    ''' maximum nuber of iterations'''
    max_loop=200
    ''' termination criteria'''
    tol =0.01
    '''Force magnitude'''
    force = 1
    '''Integration type'''
    integration_type="CompositeGauss4"

    geometry_type= "Rectangle_beam"

    density_cutoff= 0.5
    '''multiple load_case'''
    mul_load=0                   # 0 default single load, 2 multiload case
    ''' filter types'''
    filter=1            # 1 default density filter , 2 for sensitivity filter
    grey_filter=1   # 1 defualt topology optimization case , 2 for grey scale filter
    density_cutoff= 0.4
