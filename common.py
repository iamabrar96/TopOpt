from platform import node
import gmsh
from matplotlib.pyplot import axis
import numpy as np
import torch

from parameters import Parameters
class GMSH_helper():
    def __init__(self, params: Parameters) -> None:
        self.params= params

    def getDofsForNodeTags(self, tags):
        assert isinstance(tags, list)
        nodeDof= []
        for tag in tags:
            temp= (tag-1)*self.params.n_dim                                         
            nodeDof.append(np.vstack([temp+i for i in range(self.params.n_dim)]).T)
        return nodeDof

    def free_dof(self, fixeddof):
        'Elimination approach'

        dof=np.arange(0,self.params.tdof)
        freedof= np.setdiff1d(dof, fixeddof)   
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
        nodetags[:]= nodetags[:,[2,6,7,3,1,5,4,0]] #transform according to the reference element
        
        offset= np.tile(np.arange(self.params.n_dim, dtype='int'), self.params.node_per_ele)  #[0,1,2, 0,1,2, 0,1,...]
        element_dofs= np.repeat(self.params.n_dim * nodetags, self.params.n_dim, axis=1) + offset  #[3*nodetag, 3*nodetag+1, 3*nodetag+2 ...]
        coord= coord.reshape(self.params.num_elems, -1)

        return nodetags, element_dofs
    
    def get_element_centers(self):
        return gmsh.model.mesh.getBarycenters(5,-1, False, True).reshape(-1,3)
    
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
    
    def finalize(self):
        gmsh.finalize()
    

class Topology_viz:
    def __init__(self, params: Parameters):
        self.params= params
        self.step = 0
        self.t = gmsh.view.add("Topology Visualization")

    def add_view(self, densities):
        filter= densities>self.params.density_cutoff
        self.step+=1
        ele_tag, _= gmsh.model.mesh.getElementsByType(5)

        if len(filter)==0:
            raise Exception("There are zero elements that satisfy the density cutoff. Please reduce the cutoff")

        gmsh.view.addModelData(
                self.t,
                self.step,
                self.params.geometry_type,
                "ElementData",
                ele_tag[filter],  # tags of all 3d elements
                densities[filter][:, None])  # data, per element should be of shape (n, 1)
    
    def visualize(self):
        gmsh.option.setNumber('Geometry.Curves',0)
        gmsh.option.setNumber('Geometry.Points',0)
        gmsh.option.setNumber('Mesh.SurfaceEdges',0)
        gmsh.option.setNumber('Mesh.VolumeEdges',0)

        gmsh.fltk.run()

def device_common():
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        print('device is assigned as CUDA')
        torch.cuda.empty_cache()

    else:
        device = torch.device("cpu")   
        print('device is assigned as CPU')

    return device

def loss_fn(density_predicted, Jelem, volfrac, p, obj0):
    loss = torch.sum(torch.div(Jelem,density_predicted**penal))/obj0; 
    volConstraint =((torch.mean(density_predicted)/desiredVolumeFraction) - 1.0);
    return loss, volConstraint
