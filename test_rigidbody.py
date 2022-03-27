import numpy as np
import gmsh
from parameters import Parameters
from geometry import Rectangle_beam
from topopt import FE_solver
from scipy.sparse.linalg import spsolve

class Displacement_constraint_on_exterior_nodes(Rectangle_beam):

    def add_fixedbc(self):
        self.fixedNodeTags= []
        for i in range(1, 7):
            self.fixedNodeTags.append(gmsh.model.mesh.getNodes(2,i, includeBoundary=True)[0].astype('int'))

    def add_forcebc(self):
        self.forceNodeTags= []

def test__rigid_body_translation():
    params= Parameters()
    params.nelx, params.nely, params.nelz= 3,3,3 
    params.geometry_type= Displacement_constraint_on_exterior_nodes
    params.force= np.array([[0,0,0]])
    params.disp= np.array([[0,0.5,0]])

    solver= FE_solver(params)
    solver.globalstiffness_matrix()
    solver.reaction_forces_at_fixeddof()
    solver.nodal_displacements_at_freedof()
    U= solver.U
    solver.F= solver.F.reshape(-1,3)[0]
    print(F)
    for i in range(params.num_load):
        U[:, i]= spsolve(solver.kg, F[:, i])
    
    print(solver.U)
if __name__ == '__main__':
    test__rigid_body_translation()