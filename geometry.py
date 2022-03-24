import gmsh
import numpy as np
from parameters import Parameters


class Rectangle_beam:
    def __new__(cls, params:Parameters):
        if params.geometry_type == "Rectangle_beam":
            geom = super(Rectangle_beam, Rectangle_beam).__new__(Rectangle_beam)
        elif params.geometry_type == "Michell_beam":
            geom = super(Rectangle_beam, Michell_beam).__new__(Michell_beam)
        elif params.geometry_type == "multiple__load_case":
            geom = super(Rectangle_beam, multiple__load_case).__new__(multiple__load_case)
        elif params.geometry_type == "Mid_cantilever":
            geom = super(Rectangle_beam, Mid_cantilever).__new__(Mid_cantilever)
        else:
            raise Exception("Invalid geometry type selected")

        return geom

    def __init__(self, params: Parameters):
        self.params= params
        self.point1=[0,0,0]
        self.point2= [params.nelx, params.nely, params.nelz]
        gmsh.initialize()
        gmsh.model.add(type(self).__name__) #important to use the same name as the geometry type
        gmsh.option.setNumber("Mesh.RecombineAll",1)
    
    def create_geometry(self):
        '''
        This function takes dimensions(x,y,z) and the number of divisions on each axis as input 
        '''

        self.box =   gmsh.model.occ.addBox(*self.point1,*self.point2)
        gmsh.model.occ.synchronize()
        for i in [9,10,11,12] :
            gmsh.model.mesh.setTransfiniteCurve(i,self.params.nelx+1)
        for i in [1,3,5,7] : 
            gmsh.model.mesh.setTransfiniteCurve(i,self.params.nelz+1)
        for i in [2,4,6,8]:
            gmsh.model.mesh.setTransfiniteCurve(i,self.params.nely+1)
        for i in range(1,7):
            gmsh.model.mesh.setTransfiniteSurface(i)

    
    def create_mesh(self):
        gmsh.model.mesh.setTransfiniteVolume(self.box)
        gmsh.model.mesh.generate(3)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.recombine() 
        gmsh.model.mesh.renumberNodes()
    
    def visualize(self):
        gmsh.option.setNumber('Geometry.CurveLabels',1)
        gmsh.option.setNumber('Geometry.PointLabels',1)
        gmsh.option.setNumber('Mesh.NodeLabels',1)
        gmsh.fltk.run()
    
    def get_node_coord(self, dimTag=(-1,-1)): 
        '''
        This function gives the nodes and coordinates of the entire geometry
        '''

        return gmsh.model.mesh.getNodes(*dimTag, includeBoundary=True,returnParametricCoord=True)
    
    
    

    def add_forcebc(self):
        '''
        It defines the nodal points where load  is applied
        load applied on the right end of the beam
        '''

        self.forceNodeTags= np.sort(gmsh.model.mesh.getNodes(1, 5, includeBoundary=True)[0])
     
   
    def add_fixedbc(self):
        '''
        It defines the nodal points where there is no displacement .i.e. the boundary condition are fixed .
        nodes on the left surface of the beam are fixed
        '''

        self.fixedNodeTags= np.sort(gmsh.model.mesh.getNodes(2,1, includeBoundary=True)[0])

    
    def add_center(self):
        self.centerNodeTags= gmsh.model.mesh.getNodes(1, 9, includeBoundary=True)[0]

        #boundary nodes come in fornt like so [bn1, bn2, ..in_n..] so transform it such that [bn1, ..in_n.., bn2]
        self.centerNodeTags[-1], self.centerNodeTags[-2]= self.centerNodeTags[-2], self.centerNodeTags[-1]
        self.centerNodeTags[:]= np.roll(self.centerNodeTags[:],1) 

    def geom_automatic(self):
        #give this sucker a good name

        self.create_geometry()
        self.create_mesh()
        self.add_forcebc()
        self.add_fixedbc()
        self.add_center()
        


'''
different types of boundary conditions
'''
class Michell_beam(Rectangle_beam):
    def add_forcebc(self):
        # load applied on the whole top surface excluding boundary nodes
        self.forceNodeTags= np.sort(gmsh.model.mesh.getNodes(2, 4, includeBoundary=False)[0])
    def add_fixedbc(self):
        # fixed  boundary nodes of the bottom surface
        self.fixedNodeTags= np.sort(gmsh.model.mesh.getNodes(2,3, includeBoundary=True)[0])
        self.fixedNodeTags=self.fixedNodeTags[:4]       
        
class Mid_cantilever(Rectangle_beam):
    def add_forcebc(self):
        # load applied on the mid point of the right end of the beam
        self.forceNodeTags= (gmsh.model.mesh.getNodes(1, 5, includeBoundary=False)[0])
        print(self.forceNodeTags)
        # when nelz is even
        self.forceNodeTags= self.forceNodeTags[int((self.params.nelz/2)-1)] 
        # when nelz is odd
        #self.forceNodeTags= self.forceNodeTags[int(self.params.nelz/2)]
        self.forceNodeTags=self.forceNodeTags.astype(int)
        print(self.forceNodeTags)
    def add_fixedbc(self):
        self.fixedNodeTags= np.sort(gmsh.model.mesh.getNodes(2,1, includeBoundary=True)[0])

class multiple__load_case(Rectangle_beam):
    def add_forcebc(self):
        # mid point load applied opposite to each other in two different direction(up, down)
        # when nelz is even
        self.forceNodeTags= (gmsh.model.mesh.getNodes(1, 5, includeBoundary=False)[0])
        self.forceNodeTags=np.atleast_1d(self.forceNodeTags[int((self.params.nelz/2)-1)])
        self.forceNodeTags2= (gmsh.model.mesh.getNodes(1, 7, includeBoundary=False)[0])
        self.forceNodeTags1=np.atleast_1d(self.forceNodeTags2[int((self.params.nelz/2)-1)])
        self.forceNodeTags=np.concatenate((self.forceNodeTags,self.forceNodeTags1))
        print(self.forceNodeTags)
        # when nelz is odd
        # self.forceNodeTags= (gmsh.model.mesh.getNodes(1, 5, includeBoundary=False)[0])
        # self.forceNodeTags=np.atleast_1d(self.forceNodeTags[int((self.params.nelz/2))])
        # self.forceNodeTags2= (gmsh.model.mesh.getNodes(1, 7, includeBoundary=False)[0])
        # self.forceNodeTags1=np.atleast_1d(self.forceNodeTags2[int((self.params.nelz/2))])
        # self.forceNodeTags=np.concatenate((self.forceNodeTags,self.forceNodeTags1))
    def add_fixed(self):
        #  nodes on the left surface of the beam are fixed
        self.fixedNodeTags= np.sort(gmsh.model.mesh.getNodes(2,1, includeBoundary=True)[0])

if __name__ == '__main__':
    params= Parameters()
    params.geometry_type = "Mid_cantilever"
    geom= Rectangle_beam(params)
    geom.create_geometry()
    geom.create_mesh()
    gmsh.fltk.run()
    assert(isinstance(geom, Mid_cantilever))