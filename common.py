

class gmsh_helper():
    def __init__(self) -> None:
        pass 
    
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
