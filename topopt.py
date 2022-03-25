
from matplotlib import pyplot as plt
import numpy as np
from common import GMSH_helper
from parameters import Parameters
from geometry import Rectangle_beam
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import  cho_factor, LinAlgError
class FE_solver:
    ''' Finite Element Solver class '''

    def __init__(self, params:Parameters) -> None:

        self.params= params
        self.geometry= Rectangle_beam(params)
        self.geometry.geom_automatic()
        self.helper= GMSH_helper(params)



        ################################# Initialisations ################################
        self.integration_points, self.weights= self.helper.gauss_points()
        _, self.shapefunc_dertv= self.helper.basis_function(self.integration_points)
        self.shapefunc_dertv= self.shapefunc_dertv.reshape(len(self.weights), -1)

        self.determinants= self.helper.Jacob(self.integration_points)
        self.determinants= self.determinants.reshape(-1, len(self.weights))

        self.nodetags, self.element_dofs, self.centroids= self.helper.element_nodes()
        
        self.init_densities()
        self.prepare_system_of_equations()

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
    
    def init_densities(self):
        '''
        Initialisation of physical densities
        '''
        # physical densities are assigned a constant and unifrom values(initially)
        self.phy_dens= np.ones(self.params.num_elems)*self.params.volfrac   

    def simp_formula(self):
        '''
        defining the formula for the modified SIMP method(relation betwen density and youngs modulus)
        '''

        E0, Emin, p= self.params.E0, self.params.Emin, self.params.p
        
        msimp=Emin+(0.01+self.phy_dens)**p*(E0-Emin)

        return msimp

    def element_stiffness_matrix(self):
        '''
        Obtain individual element stiffeness matrices from Bmat and consitutive matrix functions
        '''
        k_ref= []
        for i in range(len(self.B_mat)):
            k_ref.append(self.weights[i]*np.matmul(np.transpose(self.B_mat[i]),np.matmul(self.C, self.B_mat[i])))
        k_ref= np.array(k_ref)
        
        # Multiply the jacobian determinants to get the element stiffness matrices
        self.ke= (self.determinants[...,np.newaxis, np.newaxis]*k_ref[np.newaxis]).sum(axis=1)

    def globalstiffness_matrix(self):
        '''
        Form global stiffness matrix by assembling all the element matrices which is obtained from the element_stiffness_matrix function
        '''
        msimp= self.simp_formula()
        ke2= self.ke.reshape(self.params.num_elems, -1)

        iK= np.repeat(self.element_dofs, self.params.node_per_ele*self.params.n_dim, axis=0).ravel()
        jK= np.repeat(self.element_dofs, self.params.node_per_ele*self.params.n_dim, axis=1).ravel()
        vK= (ke2*msimp[:, np.newaxis]).ravel()

        self.kg = coo_matrix((vK, (iK, jK)), shape= (self.params.tdof, self.params.tdof)).tocsc()

    def nodal_forces(self):

        self.F=np.zeros((self.params.tdof, self.params.num_load))

        forcedof= self.helper.getDofsForNodeTags(self.geometry.forceNodeTags)
        for i in range(self.params.num_load):
            self.F[forcedof[i], i]= self.params.force[i]

    def nodal_displacements(self):
        '''Solving for nodal displacements at free dofs'''

        fixeddof= self.helper.getDofsForNodeTags(self.geometry.fixedNodeTags)
        freedof= self.helper.free_dof(fixeddof) 

        self.U=np.zeros((self.params.tdof, self.params.num_load))
        for i in range(self.params.num_load):
            self.U[freedof, i]= spsolve(self.kg[freedof,:][:,freedof], self.F[freedof, i])

    def plot_disp(self):
        center= self.helper.getDofsForNodeTags(self.geometry.centerNodeTags)[1]
        nodes_c= np.arange(0,self.params.nelx+1)
        fig,ax=plt.subplots(figsize = (8,6))

        for i in range(self.params.num_load):
            ax.plot(nodes_c, self.U[center, i], marker = '*', color = 'blue')

        ax.set_title("Nodes along x-axis vs Displacement along y-axis",size = 14)
        ax.set_xlabel("x coordinates",size = 12)
        ax.set_ylabel("U_y",size = 12)
        # fig.savefig("Rectangle_beam")  
        plt.show()

        
    def elemental_compliance(self):
        """
        Compute element compliance and its derivative
        """

        self.Jelem= np.zeros(self.params.num_elems)
        for i in range(self.params.num_elems):
            for j in range(self.params.num_load):   
                # this order of loops is chosen so that compiler can do optimizations for simd
                self.Jelem[i]+= np.dot(np.dot(self.U[self.element_dofs[i], j], self.ke[i]), self.U[self.element_dofs[i], j])
        
        E0, Emin, p= self.params.E0, self.params.Emin, self.params.p     #just using these as local variables
        self.d_Jelem= np.zeros(self.params.num_elems)
        for i in range(self.params.num_load):    
            self.d_Jelem-= p*(E0-Emin)*self.phy_dens**(p-1)* self.Jelem

        #for some reasons not known this vectorized implementation is slower than the iterative one
        # Jelem= np.matmul(self.ke, self.U[self.element_dofs][:,:,np.newaxis]).squeeze(-1)
        # Jelem= (Jelem * self.U[self.element_dofs]).sum(axis=1)

    def prepare_system_of_equations(self):
        self.Bmat()
        self.constitutive_matrix()
        self.element_stiffness_matrix()
        self.nodal_forces()

    def solve(self, new_density= None):
        '''
        solve for nodal displacements and element compliance for given elemental densities 
        returns 
            nodal_displacements 
            derivatives of element compliances
        '''
        if new_density is not None:
            self.phy_dens= new_density

        self.globalstiffness_matrix()
        self.nodal_displacements()
        self.elemental_compliance()
        return self.U, self.d_Jelem


if __name__ == '__main__':
    # import cProfile, pstats
    # profiler = cProfile.Profile()
    # profiler.enable()

    params= Parameters()
    solver= FE_solver(params)
    density= solver.phy_dens
    solver.solve(density)
    solver.plot_disp()
    # solver.geometry.visualize()
    # profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats('tottime')
    # stats.print_stats()   


