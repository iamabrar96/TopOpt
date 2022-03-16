
from matplotlib import pyplot as plt
import numpy as np
from common import GMSH_helper
from parameters import Parameters
from geometry import Rectangle_beam
import gmsh # todo remove this later
class FE_solver:
    def __init__(self, params:Parameters) -> None:
        self.params= params
        self.helper= GMSH_helper(params)



        ################################# Initialisations ################################
        self.integration_points, self.weights= self.helper.gauss_points()
        # self.weights/=2
        #For now no need for shape functions so omitted here
        _, self.shapefunc_dertv= self.helper.basis_function(self.integration_points)
        self.shapefunc_dertv= self.shapefunc_dertv.reshape(len(self.weights), -1)
        self.determinants= self.helper.Jacob(self.integration_points)
        self.determinants= self.determinants.reshape(-1, len(self.weights))

        self.nodetags, self.element_dofs, self.centroids= self.helper.element_nodes()
        self.elem_nodes_arrangement= np.array([7,4,0,3,6,5,1,2])
        offset= np.tile(np.arange(self.params.n_dim, dtype='int'), self.params.node_per_ele)
        self.elem_nodes_arrangement= np.repeat(self.elem_nodes_arrangement*self.params.n_dim, self.params.n_dim)+ offset
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
        initialisation of design and physical variables (densities)
        '''

        self.phy_dens= np.ones(self.params.num_elems)*self.params.volfrac   # physical densities are assigned a constant and unifrom values(initially)

    def simp_formula(self):
        '''
        defining the formula for the modified SIMP method(relation betwen density and youngs modulus)
        '''

        E0, Emin, p= self.params.E0, self.params.Emin, self.params.p
        
        msimp=Emin+(self.phy_dens**p)*(E0-Emin)

        return msimp

    def element_stiffness_matrix(self):
        '''
        obtaining individual element stiffeness matrices from Bmat and consitutive matrix functions
        '''
        self.ke = []
        for j in range(self.params.num_elems):
            k=[]
            for i in range(len(self.B_mat)):
                k.append(self.determinants[j][i]*self.weights[i]*np.matmul(np.transpose(self.B_mat[i]),np.matmul(self.C, self.B_mat[i])))
            k= np.sum(k,axis = 0)
            self.ke.append(k)
            
        self.ke = np.array(self.ke)

    def globalstiffness_matrix(self):
        '''
        forming of global stiffness matrix by assembling all the element matrices which is obtained from the element_stiffness_matrix function
        '''
        msimp= self.simp_formula()

        self.kg=np.zeros((self.params.tdof,self.params.tdof))
        u,v=np.meshgrid(self.elem_nodes_arrangement, self.elem_nodes_arrangement)

        for i in range(self.params.num_elems):
            x,y=np.meshgrid(self.element_dofs[i], self.element_dofs[i])
            self.kg[y,x]+= self.ke[i][v,u] * msimp[i]
        
    
    def nodal_forces(self):
        '''
        defining the nodal forces and nodal displacements
        '''
        self.F=np.zeros(self.params.tdof)
        
        # Todo check for multiple entities in a physical group
        forcedofy= self.helper.getNodesForPhysicalGroup(dimTag=(1,2))[1]
        self.F[forcedofy]= self.params.force
        
    def nodal_displacements(self):
        fixeddof= self.helper.getNodesForPhysicalGroup(dimTag=(2,3))
        freedof= self.helper.free_dof(fixeddof) 
        #solving for nodal displacements at free dofs
        x,y=np.meshgrid(freedof,freedof)
        self.U=np.zeros(self.params.tdof)
        self.U[freedof]=np.linalg.solve(self.kg[y,x],self.F[freedof])

    def plot_disp(self):
        center= self.helper.getNodesForPhysicalGroup(dimTag=(1,4))[1]
        nodes, nodes_c, _= gmsh.model.mesh.getNodes(includeBoundary=True, returnParametricCoord=False)
        nodes_c= np.arange(0,self.params.nelx+1)
        plt.plot(nodes_c, self.U[center])
        plt.show()
        
    def elemental_compliance(self):
        self.Jelem= []
        u,v=np.meshgrid(self.elem_nodes_arrangement, self.elem_nodes_arrangement)
        for i in range(self.params.num_elems):
            self.Jelem.append(np.dot(np.dot(self.U[self.element_dofs[i]], self.ke[i][v,u]), self.U[self.element_dofs[i]]))
        self.Jelem= np.array(self.Jelem)        
    
    def prepare_system_of_equations(self):
        self.Bmat()
        self.constitutive_matrix()
        self.element_stiffness_matrix()
        self.nodal_forces()

    def solve(self, new_density):
        '''
        solve for nodal displacements and element compliance for given elemental densities 
        returns 
            nodal_displacements 
            element compliances
        '''

        self.phy_dens= new_density
        self.globalstiffness_matrix()
        self.nodal_displacements()
        self.elemental_compliance()
        return self.U, self.Jelem


if __name__ == '__main__':
    params= Parameters()
    geometry= Rectangle_beam(params)
    geometry.geom_automatic()
    solver= FE_solver(params)
    density= solver.phy_dens
    solver.solve(density)
    solver.plot_disp()
