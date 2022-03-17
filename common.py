import gmsh
from matplotlib.pyplot import axis
import numpy as np

from parameters import Parameters
class GMSH_helper():
    def __init__(self, params: Parameters) -> None:
        self.params= params

    def getNodesForPhysicalGroup(self,dimTag=(1,2)):
        p=gmsh.model.mesh.getNodesForPhysicalGroup(*dimTag)[0]
        p[1:]= np.roll(p[1:],-1)
        groupdof= (p-1)*self.params.n_dim                                         
        groupdof=np.vstack([groupdof+i for i in range(self.params.n_dim)])
        return groupdof

    def free_dof(self, fixeddof):
        'Elimination approach'

        dof=np.arange(0,self.params.tdof)
        freedof= np.setdiff1d(dof, fixeddof.flatten())   
        return freedof

    @property
    def element_type(self):
        #Todo get element type based on the geometry not hardcoded
        my_element= 5 #gmsh.model.mesh.getElementType(familyName="Hexahedron", order=1, serendip = False)
        return my_element

    def elemprop(self):
        prop=gmsh.model.mesh.getElementProperties(elementType=self.my_element)
        return prop

    def element_nodes(self):
        '''
        this function gives the nodes, dofs and coordinates of each hexahydron element
        '''

        nodetags, coord, _ =gmsh.model.mesh.getNodesByElementType(elementType= self.element_type,tag = -1, returnParametricCoord = False)
        nodetags= nodetags.reshape(self.params.num_elems, -1).astype('int')-1 # since in python indices start with zero
        
        offset= np.tile(np.arange(self.params.n_dim, dtype='int'), self.params.node_per_ele)  #[0,1,2, 0,1,2, 0,1,...]
        element_dofs= np.repeat(self.params.n_dim * nodetags, self.params.n_dim, axis=1) + offset  #[3*nodetag, 3*nodetag+1, 3*nodetag+2 ...]

        coord= coord.reshape(self.params.num_elems, -1)
        centroids=[]
        for j in range(self.params.num_elems):
            centroids.append((np.sum(coord[j].reshape(-1,3),axis=0)/self.params.node_per_ele))
        centroids=np.vstack(centroids)
        return nodetags, element_dofs, centroids
    
    def gauss_points(self):
        '''
        this function takes the output of elementtype function as input  and fourth order gauss quadrature rule is applied
        '''
        points, weights=gmsh.model.mesh.getIntegrationPoints(elementType= self.element_type, integrationType= self.params.integration_type)
        return points, weights
 
    def basis_function(self, integration_points):
        _, shapefunc,_ = gmsh.model.mesh.getBasisFunctions(elementType= self.element_type,localCoord=integration_points,functionSpaceType="Lagrange")
        _, shapefunc_dertv,_= gmsh.model.mesh.getBasisFunctions(elementType= self.element_type,localCoord= integration_points,functionSpaceType="GradLagrange")

        return shapefunc, shapefunc_dertv

                            # conversion of list to array of shapefunction derivative 

    def Jacob(self, integration_points):
        _, determinants, _=gmsh.model.mesh.getJacobians(elementType=self.element_type,localCoord= integration_points,tag = -1, task = 0, numTasks = 1)
        return determinants
    
    def get_entities_for_physical_group(self):
        physical_group=gmsh.model.getPhysicalGroups()
    
    
    
    def finalize(self):
        gmsh.finalize()
    

class Topology_viz:
    def __init__(self, params: Parameters):
        self.params= params
        self.step = 0
        self.t = gmsh.view.add("Topology Visualization")

    def add_view(self, temp):
        densities= temp.copy()
        densities[densities<self.params.density_cutoff] = 0.0
        self.step+=1
        ele_tag, _= gmsh.model.mesh.getElementsByType(5)

        gmsh.view.addModelData(
                self.t,
                self.step,
                self.params.geometry_type,
                "ElementData",
                ele_tag,  # tags of all 3d elements
                densities[:, None])  # data, per element should be of shape (n, 1)
    
    def visualize(self):
        gmsh.fltk.run()