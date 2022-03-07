
import numpy as np
from common import gmsh_helper
from new import Rectangle_beam
from parameters import Parameters


class FE_solver:
    def __init__(self, geometry: Rectangle_beam, params:Parameters) -> None:
        self.geometry= geometry
        self.params= params
        self.helper= gmsh_helper()



        ################################# Initialisations ################################

        self.integration_points, self.weights= self.helper.gauss_points(self.params.integration_type)

        #For now no need for shape functions so omitted here
        _, self.shapefunc_dertv= self.helper.basis_function(self.integration_points)
        self.shapefunc_dertv= self.shapefunc_dertv.reshape(len(self.weights), -1)

        self.determinants= self.helper.Jacob(self.integration_points)
        self.determinants= self.determinants.reshape(-1, len(self.weights))

        self.init_density()
        self.nodetags, self.centroids= self.helper.element_nodes()

    
    
    def Bmat(self):
        '''
        input : shape function derivatives with respect to x, y ,z axis
        output: all stacked strain displacmenet matrix (B) together i.e. No. of B matrix is realted to  number of gauss points 
        '''
        
        self.B_mat = []                        # it consists of all the B matrices(27) stacked together 
        for row in self.shapefunc_dertv:
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
            self.B_mat.append(np.hstack(B))
        self.B_mat=np.array(self.B_mat)  

    

    def constitutive_matrix(self):
        '''
        it defines the threee dimensional constitutive matrix for an isotropic element
        '''
    
        C=np.zeros((6,6))
        C[0,0]=1-self.params.nu
        C[0,1]=self.params.nu
        C[0,2]=C[0,1]
        C[1,0]=C[0,1]
        C[1,1]=C[0,0]
        C[1,2]=C[0,1]
        C[3,3]=(1-2*self.params.nu)/2
        C[2,2]=C[0,0]
        C[4,4]=C[3,3]
        C[5,5]=C[3,3]
        C[2,0]=C[0,1]
        C[2,1]=C[0,1] 
        self.C= 1/((1+self.params.nu)*(1-2*self.params.nu))*C
    

    
    def element_stiffness_matrix(self):
        '''
        obtaining individual element stiffeness matrices from Bmat and consitutive matrix functions
        '''

        self.ke = []
        for j in range(self.params.num_elems):
            k=[]
            for i in range(len(self.B_mat)):
                k.append(self.determinants[j][i]*self.weights[i]*(np.matmul(np.transpose(self.B_mat[i]),np.matmul(self.C, self.B_mat[i]))))
            k= np.sum(k,axis = 0) 
            self.ke.append(k)
        self.ke = np.array(self.ke)

    
    def init_density(self):
        '''
        initialisation of design and physical variables (densities)
        '''

        self.des_dens=np.ones((self.params.num_elems,1))*self.params.volfrac
        self.pyh_dens=self.des_dens                  # physical densities are assigned a constant and unifrom values(initially)

    def simp_formula(self):
        '''
        defining the formula for the modified SIMP method(relation betwen density and youngs modulus)
        '''

        E0, Emin, p= self.params.E0, self.params.Emin, self.params.p
        nelx, nely, nelz= self.params.nelx, self.params.nely, self.params.nelz
        
        msimp=Emin+(np.transpose(self.pyh_dens.flatten())**p)*(E0-Emin)

        return msimp
    
    def globalstiffness_matrix(self):
        '''
        forming of global stiffness matrix by assembling all the element matrices which is obtained from the element_stiffness_matrix function
        '''

        self.kg=np.zeros((self.params.tdof,self.params.tdof))
        msimp= self.simp_formula()
        
        for i in range(self.params.num_elems):
            nodes=np.vstack((self.nodetags[i]*3, self.nodetags[i]*3+1, self.nodetags[i]*3+2)).flatten('F')
            x,y=np.meshgrid(nodes, nodes)
            self.kg[y,x]+= msimp[i]*self.ke[i]
    
    def nodal_forces(self):
        '''
        defining the nodal forces and nodal displacements
        '''
        self.F=np.zeros(self.params.tdof)
        
        # Todo check for multiple entities in a physical group
        forcedofyc= self.helper.getNodesForPhysicalGroup(dimTag=(1,2))[1]
        self.F[forcedofyc]= self.params.force
        
    def nodal_displacements(self):
        fixeddofc= self.helper.getNodesForPhysicalGroup(dimTag=(2,3))
        freedof= self.helper.free_dof() 
        
        #solving for nodal displacements at free dofs
        x,y=np.meshgrid(freedof,freedof)
        self.U=np.zeros((self.params.tdof,1))
        self.U[freedof]=np.linalg.solve(self.kg[y,x],self.F[freedof])
    
    def solve(self):
        '''Call all the individual functions for an automatic solve'''
        self.Bmat()
        self.constitutive_matrix()
        self.element_stiffness_matrix()
        self.globalstiffness_matrix()
        self.nodal_forces()
        self.nodal_displacements()


params= Parameters()
geometry= Rectangle_beam()
solver= FE_solver(geometry, params)
solver.solve()
print(solver.U)

