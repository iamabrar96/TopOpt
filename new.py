from functools import cache
import gmsh
import numpy as np
from numpy.matlib import repmat
import sys
import math
import itertools
from numpy import exp
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
#######################################################################
''' Input Parameters'''
#######################################################################

''' number of nodes per element'''
node_per_ele=8
'''  degrees of freedom  per node'''
dof_per_node=3
''' degrees of freedom per element '''
dof_per_ele= node_per_ele*dof_per_node
''' number of dimension of the problem'''
nd_per_element=3
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
no_of_ele=nelx*nely*nelz  
''' Total number of nodes'''   
tnodes=(nelx+1)*(nely+1)*(nelz+1)
''' Total number of degrees of freedom'''
tdof=dof_per_node*tnodes
''' maximum nuber of iterations'''
max_loop=200
''' termination criteria'''
tol =0.01

##############################################################################


class Rectangle_beam:
   

    def __init__(self,point1,point2):
        self.point1=point1
        self.point2=point2
        gmsh.initialize()
        gmsh.model.add("rectangular beam")
        gmsh.option.setNumber("Mesh.RecombineAll",1)
    '''
    This function takes dimensions(x,y,z) and the number of divisions on each axis as input 
    '''
    def create_geometry(self):
        self.box =   gmsh.model.occ.addBox(*self.point1,*self.point2)
        gmsh.model.occ.synchronize()
        for i in [9,10,11,12] :
            gmsh.model.mesh.setTransfiniteCurve(i,nelx+1)
        for i in [1,3,5,7] : 
            gmsh.model.mesh.setTransfiniteCurve(i,nelz+1)
        for i in [2,4,6,8]:
            gmsh.model.mesh.setTransfiniteCurve(i,nely+1)
        for i in range(1,7):
            gmsh.model.mesh.setTransfiniteSurface(i)

    
    def create_mesh(self):
        gmsh.model.mesh.setTransfiniteVolume(self.box)
        gmsh.model.mesh.generate(3)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.recombine() 

    '''
    This function gives the nodes and coordinates of the entire geometry
    '''
    def get_node_coord(self, dimTag=(-1,-1)): 
        return gmsh.model.mesh.getNodes(*dimTag, includeBoundary=True,returnParametricCoord=True)
    
    
    '''
    It defines the point load  in y direction on specified nodal positions.
    It returns the tag number which in  this case  is 2
    '''

    def add_forcebc(self):
        self.force_bc=gmsh.model.addPhysicalGroup(1,[5],2)
        #gmsh.model.setPhysicalName(1, self.force_bc, "Forced boundary condition")

    '''
    It defines the nodal points where there is no displacement .i.e. the boundary condition are fixed .
    similarly it also returns the tag number which in this case is 3
    '''
    def add_fixedbc(self):
        self.fixed_bc= gmsh.model.addPhysicalGroup(1,[4,3,1,2],3)
        #gmsh.model.setPhysicalName(1, self.fixed_bc, "Fixed boundary condition")
            
      
    def getNodesForPhysicalGroup(self,dimTag=(1,2)):
        self.p=gmsh.model.mesh.getNodesForPhysicalGroup(*dimTag)
        self.forcedofy=((self.p[0]-1 )*nd_per_element)                                          # node numbers of forcedof
        self.forcedofyc= self.forcedofy+1                                                       # node coordinate of forcedof
        self.fixeddof=((self.p[0]-1)*nd_per_element)                                             # node numbers of fixddof
        self.fixeddofc=np.concatenate((self.fixeddof,self.fixeddof+1,self.fixeddof+2))           # node coordinates of fixeddof 

    def free_dof(self):
        self.setdiff=np.arange(0,tdof)
        self.freedof=[]                                                                # it contains all the dofs except the fixeddof
        for i in self.setdiff:
            if i not in self.fixeddofc:
                self.freedof.append(i)

    def elementtype(self):
        self.my_element= 5 #gmsh.model.mesh.getElementType(familyName="Hexahedron", order=1, serendip = False)
       
    def elemprop(self):
        self.prop=gmsh.model.mesh.getElementProperties(elementType=self.my_element)
    '''
    this function gives the nodes and coordinates of each hexahydron element
    '''
    def element_nodes(self):
        self.nodetags, self.coord, self.parametricCoord=gmsh.model.mesh.getNodesByElementType(elementType=self.my_element,tag = -1, returnParametricCoord = True)
        self.my_nodes=np.split(self.nodetags,no_of_ele)
        self.my_coord=np.split(self.coord,no_of_ele)
        self.my_coord_edof=[]      # degrees of freedom of each element along x,y,z direction
        
        # self.my_coord=np.split(self.coord,no_of_ele)
        # self.my_coord_edof=[]      # degrees of freedom of each element along x,y,z direction 
        
        for i in range(no_of_ele):

            self.my_coord_edof.append(np.vstack((self.my_nodes[i]*3-3,self.my_nodes[i]*3-2,self.my_nodes[i]*3-1)).flatten('F'))
        self.my_coord_edof=np.vstack(self.my_coord_edof)
        # obtaining center point(x,y,z) of each element 
        self.my_coord1=np.array(self.my_coord)
        self.coord2=[]
        for i in (self.my_coord1):
            self.coord2.append(i.reshape(-1,3))
        self.center_point=[]
        for j in range(no_of_ele):
            self.center_point.append((np.sum(self.coord2[j],axis=0)/node_per_ele))
        self.centroid=np.vstack(self.center_point)
    '''
    this function takes the output of elementtype function as input  and fourth order gauss quadrature rule is applied
    '''
    def gauss_points(self):
        self.Integration_points, self.weights=gmsh.model.mesh.getIntegrationPoints(elementType=self.my_element,integrationType="CompositeGauss4")

    def len_quad_pts(self):
        self.total_quad_pt=len(self.weights)
   
    def basis_function(self):
        _, self.shapefunc,_ = gmsh.model.mesh.getBasisFunctions(elementType=self.my_element,localCoord=self.Integration_points,functionSpaceType="Lagrange")
        _,self.shapefunc_dertv,_= gmsh.model.mesh.getBasisFunctions(elementType=self.my_element,localCoord=self.Integration_points,functionSpaceType="GradLagrange")
    
    def split_shapefunc(self):
        self.spl=np.split(self.shapefunc_dertv,self.total_quad_pt)
        self.spl_ele=np.array(self.spl)                               # conversion of list to array of shapefunction derivative 

    def Jacob(self):
        self.jacobians, self.determinants, self.coord=gmsh.model.mesh.getJacobians(elementType=self.my_element,localCoord=self.Integration_points,tag = -1, task = 0, numTasks = 1)
        self.determinants=np.array(self.determinants)
        self.determinants=self.determinants.reshape(-1,self.total_quad_pt)
    
    '''
    input : shape function derivatives with respect to x, y ,z axis
    output: all stacked strain displacmenet matrix (B) together i.e. No. of B matrix is realted to  number of gauss points 
    '''
    def Bmat(self):
        self.a = []                        # it consists of all the B matrices(27) stacked together 
        for row in self.shapefunc_dertv.reshape(-1, 24):
            B=[]
            row = row.reshape(-1,3)
            for row2 in row:
                B_i= np.zeros((6,3))
                B_i[0,0]=row2[0]
                B_i[1,1]=row2[1]
                B_i[2,2]=row2[2]
                B_i[3,0]=row2[1]
                B_i[3,1]=row2[0]
                B_i[4,1]=row2[2]
                B_i[4,2]=row2[1]
                B_i[5,0]=row2[2]
                B_i[5,2]=row2[0]
                B.append(B_i)
            self.a.append(np.hstack(B))
        self.a=np.array(self.a)  

    '''
    it defines the threee dimensional constitutive matrix for an isotropic element
    '''

    def constitutive_matrix(self):
        
        C=np.zeros((6,6))
        C[0,0]=1-nu
        C[0,1]=nu
        C[0,2]=C[0,1]
        C[1,0]=C[0,1]
        C[1,1]=C[0,0]
        C[1,2]=C[0,1]
        C[3,3]=(1-2*nu)/2
        C[2,2]=C[0,0]
        C[4,4]=C[3,3]
        C[5,5]=C[3,3]
        C[2,0]=C[0,1]
        C[2,1]=C[0,1] 
        self.C= 1/((1+nu)*(1-2*nu))*C
    
    # def ref_element(self):
    #     k= np.zeros((24,24))
    #     for i,weight in enumerate(self.weights):
    #         k+= weight*np.matmul(self.a[i].T, np.matmul(self.C, self.a[i]))
    #     print(k[10,9])
    #     exit()
    '''
    obtaining individual element stiffeness matrices from Bmat and consitutive matrix functions
    '''

    def element_stiffness_matrix(self):
        self.ke = []
        for j in range(self.determinants.shape[0]):
            self.k=[]
            for i in range (self.a.shape[0]):
                self.k.append(self.determinants[j][i]*self.weights[i]*(np.matmul(np.transpose(self.a[i]),np.matmul(self.C, self.a[i]))))
            self.k= np.sum(self.k,axis = 0) 
            #self.k=np.array(self.k)
            self.ke.append(self.k)
        self.ke = np.array(self.ke)

    '''
    forming of global stiffness matrix by assembling all the element matrices which is obtained from the element_stiffness_matrix function
    '''
    def globalstiffness_matrix(self):
        self.kg=np.zeros((tdof,tdof))
        for i in range(self.ke.shape[0]):
            self.Nodes=np.vstack((self.my_nodes[i]*3-3,self.my_nodes[i]*3-2,self.my_nodes[i]*3-1)).flatten('F')
            self.Nodes=np.array(self.Nodes)
            x,y=np.meshgrid(self.Nodes,self.Nodes)
            self.kg[y,x]+= self.msimp[:,i]*self.ke[i]
    
    '''
    preparing the filter function
    '''
    def density_filter(self):
        self.ih=[1]*no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2
        self.jh=[1]*no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2
        self.sh=[0]*len(self.ih)
        self.cn=0                  # counter 
        for k1 in range(1,nelz+1):
            for i1 in range(1,nelx+1):
                for j1 in range(1,nely+1):
                    self.e1=(k1-1)*nelx*nely + (i1-1)*nely + (j1 -1)
                    for k2 in range  (np.maximum(k1-math.floor(rmin),1),np.minimum(k1+math.floor(rmin),nelz)+1):
                        for i2 in range(np.maximum(i1-math.floor(rmin),1),np.minimum(i1+math.floor(rmin),nelx)+1):
                            for j2 in range(np.maximum(j1-math.floor(rmin),1),np.minimum(j1+math.floor(rmin),nely)+1):
                               self.e2=((k2-1)*nelx*nely) +(i2-1)*nely+ (j2 -1)
                               if self.cn<no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2:  
                                    self.ih[self.cn]=self.e1                        # row indexing
                                    self.jh[self.cn]=self.e2                        # column indexing
                                    self.sh[self.cn]=np.maximum(0,rmin-math.sqrt((i1-i2)**2 + (j1-j2)**2 +(k1-k2)**2))
                               else:
                            
                                    self.ih.append(self.e1)
                                    self.jh.append(self.e2)
                                    self.sh.append(np.maximum(0,rmin-math.sqrt((i1-i2)**2 + (j1-j2)**2 +(k1-k2)**2)))

                               self.cn=self.cn+1

        self.row=np.array(self.ih)
        self.column=np.array(self.jh)
        self.val=np.array(self.sh)                       
        self.H=coo_matrix((self.val,(self.row,self.column)),shape=(nelx*nely*nelz,nelx*nely*nelz)).tocsc()
        self.HS=np.sum(self.H,axis=1)
    '''
    initialisation of design and physical variables (densities)
    '''
    def densities(self):
        self.des_dens=np.ones((no_of_ele,1))*volfrac
        self.pyh_dens=self.des_dens                  # physical densities are assigned a constant and unifrom values(initially)
    '''
    defining the formula for the modified SIMP method(relation betwen density and youngs modulus)
    '''
    def simp_formula(self):
        self.msimp=Emin+(np.transpose(self.pyh_dens.reshape(nelx*nely*nelz,1))**p)*(E0-Emin)

    '''
    defining the nodal forces and nodal displacements
    '''
    def nodal_forces(self):
        self.F=np.zeros((tdof,1))
        self.F[self.forcedofyc[0],0]=-1
        self.F[self.forcedofyc[1],0]=-1
        self.F[self.forcedofyc[2],0]=-1

    def nodal_displacements(self):
        self.U=np.zeros((tdof,1))
        x,y=np.meshgrid(self.freedof,self.freedof)
        self.U[self.freedof]=np.linalg.solve(self.kg[y,x],self.F[self.freedof])
    
##############################################################################
        ''' 
            iteration  process

        '''
############################################################################

    def densityestimation(self):
 
        loop=0
        difference=1
        #while difference>tol and loop < max_loop:
        if difference>tol:
            if loop<max_loop:

                # update global stiffness using x.globalstiffness()
                # solve for nodal disp using x.nodaldisp()
                # complicance
                # flter
                # Update densities
                ''''minimum compliance objective function and sensitivity analysis'''
                comp = []

                for i in range(no_of_ele):
                    temp= np.dot(self.U[self.my_coord_edof[i]].T,np.dot(self.ke[i],self.U[self.my_coord_edof[i]]))
                    comp.append(temp)
                comp= np.array(comp).reshape(no_of_ele,1)
                #print(comp)
                der_comp=(comp*( -p*(E0-Emin)*self.pyh_dens**(p-1) ))
                der_vol=np.ones((no_of_ele,1))
                ''' use of filtering function to improve the sensitivity analysis  '''
                der_comp=self.H*( der_comp/self.HS)
                der_vol= self.H*( der_vol/self.HS)
                ''' Optimality criteria update scheme'''
                ##### implementing the bisection algorithm to predict the lambda value ######
                l1=0 
                l2=1e9
                forward=0.2
                while (l2-l1)/(l1+l2)>1e-3:
                    lmid=0.5*(l2+l1)
                    new_density= np.maximum(0.0,np.maximum(self.des_dens-forward,np.minimum(1.0,np.minimum(self.des_dens+forward,np.multiply(self.des_dens,np.sqrt(-(der_comp)/ der_vol/lmid))))))
                    self.pyh_dens= (self.H*new_density)/self.HS
                    if np.sum(self.pyh_dens)>volfrac*no_of_ele:
                        l1=lmid
                    else :
                        l2=lmid
                difference=abs(new_density-self.des_dens)
                self.des_dens=new_density
                loop=loop+1

    def visualize(self):
        gmsh.fltk.run()
    
    def finalize(self):
        gmsh.finalize()


x=Rectangle_beam(point1=[0,0,0],point2=[nelx,nely,nelz])
x.create_geometry()
x.create_mesh()
# x.visualize()
# x.finalize()
x.add_forcebc()
x.getNodesForPhysicalGroup(dimTag=(1,2))
x.nodal_forces()
x.add_fixedbc()
x.getNodesForPhysicalGroup(dimTag=(1,3))
z=x.get_node_coord()
x.free_dof()
x.densities()
x.simp_formula()
x.elementtype()
x.elemprop()
x.element_nodes()
#r1=x.my_coord[4][1],x.my_coord[0][21],x.my_coord[1][21],x.my_coord[2][21],x.my_coord[3][21]
x.gauss_points()
x.len_quad_pts()
x.basis_function()
x.split_shapefunc()
x.Jacob()
x.Bmat()
x.constitutive_matrix()
x.density_filter()
# x.ref_element()
x.element_stiffness_matrix()
x.globalstiffness_matrix()
x.nodal_displacements()
x.densityestimation()
#r=x.U.reshape(30,3)
#r2=r[8][1],r[24][1],r[25][1],r[26][1],r[10][1] 
#plt.plot(r1,r2) # plotting by taking x coordinates of nodes along the x- aixs(8,24,25,26,10) w.r.t disp along y 
#plt.show()

#print(comp)
# print(x.centroid)
print(x.pyh_dens)


class my_neural_network:
    def __init__(self,X_train, Y_train,layer_dims ):
        self.X_train = X_train
        self.Y_train=Y_train   
        self.layer_dims=layer_dims
        parameters = self.initialize_parameters_deep()
        print(parameters)
 
    def initialize_parameters_deep(self):
        """
        Arguments:
        layer_dims -- python array (list) containing the dimensions of each layer in the network
        
        Returns:
        parameters -- python dictionary containing  parameters "W1", "b1", ..., "WL", "bL":
                        Wl -- weight matrix of shape (layer_dims[l], layer_dims[l-1])
                        bl -- bias vector of shape (layer_dims[l], 1)
        """
        
        parameters = {}
        L = len(self.layer_dims) # number of layers in the network

        for l in range(1, L):
            parameters["W" + str(l)] = np.random.uniform(-1/np.sqrt(self.layer_dims[l - 1]),1/np.sqrt(self.layer_dims[l - 1]),(self.layer_dims[l], self.layer_dims[l-1]))
            parameters["b" + str(l)] = np.zeros((self.layer_dims[l], 1))
            parameters["gamma" + str(l)] = np.ones((self.layer_dims[l], 1))
            parameters["beta" + str(l)] = np.zeros((self.layer_dims[l], 1))
            assert(parameters['W' + str(l)].shape == (self.layer_dims[l], self.layer_dims[l - 1]))
            assert(parameters['b' + str(l)].shape == (self.layer_dims[l], 1))
            assert(parameters['gamma' + str(l)].shape == (self.layer_dims[l], 1))
            assert(parameters['beta' + str(l)].shape == (self.layer_dims[l], 1))
        return parameters

    def linear_forward(self,X_train, W, b):
        """
        Implement the linear part of a layer's forward propagation.

        Arguments:
        A -- activations from previous layer (or input data): (size of previous layer, number of examples)
        W -- weights matrix: numpy array of shape (size of current layer, size of previous layer)
        b -- bias vector, numpy array of shape (size of the current layer, 1)

        Returns:
        Z -- the input of the activation function, also called pre-activation parameter 
        cache -- a python tuple containing "A", "W" and "b" ; stored for computing the backward pass efficiently
        """
        
      
        Z=np.dot(W,X_train)+b
        
     
        cache = (X_train, W, b)
        
        return Z, cache


    def batchnorm_forward(self,Z, gamma, beta, eps=1e-5):
        N, D = Z.shape
        
        sample_mean = Z.mean(axis=0)
        sample_var = Z.var(axis=0)
        
        std = np.sqrt(sample_var + eps)
        x_centered = Z - sample_mean
        x_norm = x_centered / std
        out = gamma * x_norm + beta
        
        cache = (x_norm, x_centered, std, gamma)

        return out, cache

    def sigmoid(self,Z):
        """
        Compute the sigmoid of z
        """

        A=1/(1+np.exp(-Z))
        cache = Z
        return A, cache

    def leaky_relu(self,Z):
        A=np.where(Z>0,Z,Z*0.01)
        cache=Z
        return A, cache

    def relu(self,Z):
        A=np.maximum(0,Z)
        assert(A.shape == Z.shape)
        cache = Z 
        return A, cache

    def soft_max(self,Z):
        A=exp(Z)/np.sum(exp(Z))
        cache=Z
        return A, cache

    def linear_activation_batch_forward(self,X_train_prev, W, b,gamma,beta, activation):
        """
        Implement the forward propagation for the LINEAR->ACTIVATION layer

        Arguments:
        x_train_prev -- activations from previous layer (or input data): (size of previous layer, number of examples)
        W -- weights matrix: numpy array of shape (size of current layer, size of previous layer)
        b -- bias vector, numpy array of shape (size of the current layer, 1)
        activation -- the activation to be used in this layer, stored as a text string: "sigmoid" or "relu"

        Returns:
        A -- the output of the activation function, also called the post-activation value 
        cache -- a python tuple containing "linear_cache" and "activation_cache";
                stored for computing the backward pass efficiently
        """
        
        if activation == "sigmoid":
            Z, linear_cache = self.linear_forward(X_train_prev, W, b)
            Z_hat, batchnorm_cache = self.batchnorm_forward(Z, gamma, beta, eps=1e-5)
            A, activation_cache = self.sigmoid(Z_hat)
        
        elif activation == "relu":
            Z, linear_cache = self.linear_forward(X_train_prev, W, b)
            Z_hat, batchnorm_cache = self.batchnorm_forward(Z, gamma, beta, eps=1e-5)
            A, activation_cache = self.relu(Z_hat)
            
        elif activation == "leaky_relu":
            Z, linear_cache = self.linear_forward(X_train_prev, W, b)
            Z_hat, batchnorm_cache = self.batchnorm_forward(Z, gamma, beta, eps=1e-5)
            A, activation_cache = self.leaky_relu(Z_hat)
        elif activation == "soft_max":
            Z, linear_cache = self.linear_forward(X_train_prev, W, b)
            ###################################################################################################################siggu
            Z_hat, batchnorm_cache = self.batchnorm_forward(Z, gamma, beta, eps=1e-5)
            A, activation_cache = self.soft_max(Z)
        cache = (linear_cache, batchnorm_cache, activation_cache)

        return A, cache

    def L_model_forward(self,X_train, parameters):
        """
        Implement forward propagation for the [LINEAR->RELU]*(L-1)->LINEAR->SIGMOID computation
        
        Arguments:
        X -- data, numpy array of shape (input size, number of examples)
        parameters -- output of initialize_parameters_deep()
        
        Returns:
        AL -- activation value from the output (last) layer
        caches -- list of caches containing:
                    every cache of linear_activation_forward() (there are L of them, indexed from 0 to L-1)
        """

        caches = []
        A = X_train
        L = len(parameters) // 2                  # number of layers in the neural network
        
        #  [LINEAR -> RELU]*(L-1). Add "cache" to the "caches" list.
        # The for loop starts at 1 because layer 0 is the input
        for l in range(1, L):
           X_train_prev = A 
           W,b= parameters["W" + str(l)],parameters["b" + str(l)]
            
           A, cache = self.linear_activation_batch_forward(X_train_prev,W,b,"leaky_relu")
           caches.append(cache)
        
        #  LINEAR -> SIGMOID. Add "cache" to the "caches" list.
        W,b= parameters["W" + str(L)], parameters["b" + str(L)] 
            
        AL, cache =self.linear_activation_batch_forward(A,W,b, "soft_max")
        caches.append(cache)

            
        return AL, caches


    def 






network= my_neural_network(np.transpose(x.centroid),x.pyh_dens,[3,20,20,20,20,20,1])

