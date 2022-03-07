import gmsh
import numpy as np
from parameters import Parameters


class Rectangle_beam:
   

    def __init__(self, params: Parameters):
        self.params= params
        self.point1=[0,0,0]
        self.point2= [params.nelx, params.nely, params.nelz]
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
        self.fixed_bc= gmsh.model.addPhysicalGroup(2,[1],3)
        #gmsh.model.setPhysicalName(1, self.fixed_bc, "Fixed boundary condition")

    def visualize(self):
        gmsh.fltk.run()
    
    def finalize(self):
        gmsh.finalize()

params= Parameters()
rec_geo=Rectangle_beam(params)
rec_geo.create_geometry()
rec_geo.create_mesh()
rec_geo.visualize()
rec_geo.finalize()

'''
different types of boundary conditions
'''
class Michell_beam(Rectangle_beam):
    def add_forcebc(self):
        self.force_bc=gmsh.model.addPhysicalGroup(0,[14,26,17],2)
    def add_fixedbc(self):
         self.fixed_bc= gmsh.model.addPhysicalGroup(1,[1,5],3)
class Distributed_load_beam(Rectangle_beam):
    def add_forcebc(self):
        self.force_bc=gmsh.model.addPhysicalGroup(1,[1,5],2)
    def add_fixedbc(self):
        self.fixed_bc= gmsh.model.addPhysicalGroup(2,[4],3)
class mid_cantilever(Rectangle_beam):
    def add_forcebc(self):
        self.force_bc=gmsh.model.addPhysicalGroup(0,[11,12],2)
    def add_fixedbc(self):
        self.fixed_bc= gmsh.model.addPhysicalGroup(2,[1],3)