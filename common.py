import gmsh
import numpy as np

from parameters import Parameters
class gmsh_helper():
    def __init__(self, params: Parameters) -> None:
        self.params= params

    def getNodesForPhysicalGroup(self,dimTag=(1,2)):
        p=gmsh.model.mesh.getNodesForPhysicalGroup(*dimTag)
        fixeddof= (p[0]-1)*self.params.n_dim                                            # node numbers of fixddof
        # node coordinates of fixeddof 
        fixeddofc=np.vstack([fixeddof,fixeddof+1,fixeddof+2])          # node coordinates of fixeddof 
        return fixeddofc

    def free_dof(self, fixeddof):
        'Elimination approach'

        dof=np.arange(0,self.params.tdof)
        freedof= np.setdiff1d(dof, fixeddof.flatten())                                                               # it contains all the dofs except the fixeddof
        return freedof

    @property
    def elementtype(self):
        #Todo get element type based on the geometry not hardcoded
        my_element= 5 #gmsh.model.mesh.getElementType(familyName="Hexahedron", order=1, serendip = False)
        return my_element

    def elemprop(self):
        prop=gmsh.model.mesh.getElementProperties(elementType=self.my_element)
        return prop

    def element_nodes(self):
        '''
        this function gives the nodes and coordinates of each hexahydron element
        '''

        nodetags, coord, _ =gmsh.model.mesh.getNodesByElementType(elementType= self.element_type,tag = -1, returnParametricCoord = True)
        
        nodetags= nodetags.reshape(self.params.num_elems, -1) -1 # since in python indices start with zero
        coord= coord.reshape(self.params.num_elems, -1)

        centroids=[]
        for j in range(self.params.num_elems):
            centroids.append((np.sum(coord[j].reshape(-1,3),axis=0)/self.params.node_per_ele))
        centroids=np.vstack(centroids)

        return nodetags, centroids
    
    def gauss_points(self, integration_type="CompositeGauss4"):
        '''
        this function takes the output of elementtype function as input  and fourth order gauss quadrature rule is applied
        '''
        points, weights=gmsh.model.mesh.getIntegrationPoints(elementType= self.element_type, integrationType= integration_type)
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
