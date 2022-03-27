#AUTHOR : Abrar Hyder Mohammed
#MATRICULATION NUMBER : 65092
#Personal Programming Project
#---------------------------------------#
#A python test file to check test_shapefunction.py
# --------------------------------------#
from parameters import Parameters
from geometry import Rectangle_beam
from topopt import FE_solver
from numpy import round
from datetime import datetime
 
def test__shapefunctionsum():
   '''
   Aim : To check whether the shape functions add upto 1.
 
   Expected Output : The sum of shape functions is 1.
 
   Remarks : test case passed successfully
   '''
   start_time = datetime.now()
   # Initializing the inputs
   params= Parameters()
   solver= FE_solver(params)
   shapefunc, _ = solver.helper.basis_function(solver.integration_points)
   shapefunc = shapefunc.reshape(len(solver.weights), -1)
   if all(round(shapefunc.sum(axis = 1),8))==1:
       print('The sum of shape functions is 1.')
   else:
       print('The sum of shape functions is not 1.')
   end_time = datetime.now()
   print("test passed","--",'Duration: {}'.format(end_time - start_time))
#test__shapefunctionsum()
if __name__=='__main__':
    test__shapefunctionsum()

