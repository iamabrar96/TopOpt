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
    nelx=4
    ''' number of element in y direction'''
    nely=1
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
    force = -1
    '''Integration type'''
    integration_type="CompositeGauss4"
