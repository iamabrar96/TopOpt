import math
import numpy as np
from common import Topology_viz
from parameters import Parameters
from topopt import FE_solver
from scipy.sparse import coo_matrix
from tqdm import tqdm
class SimpOptimizer:
    def __init__(self, params:Parameters) -> None:
        self.params = params
        self.solver= FE_solver(params)
        self.topology_visualizer= Topology_viz(params)

        ################################# Initialisations ################################
        self.density_filter2()

    def densityestimation(self):
        
        phy_dens= self.solver.phy_dens
        old_dens= phy_dens
        difference=1

        self.compliances = []
        msimp= self.solver.simp_formula()

        for loop in tqdm(range(self.params.max_loop)):

            if not difference>self.params.tol:
                break 

            
            ''''minimum compliance objective function and sensitivity analysis'''
            
            _, Jelem, d_Jelem= self.solver.solve(phy_dens) 
            self.compliances.append(Jelem.dot(msimp))
            
            d_vol=np.ones(self.params.num_elems)

            if self.params.filter==1:   #case: density filter
                d_Jelem= self.H.dot(d_Jelem/self.HS)
                d_vol= self.H.dot(d_vol/self.HS)
            
            elif self.params.filter==2: #case: sensitivity filter
                d_Jelem= np.dot(self.H ,old_dens*d_Jelem)/(self.HS*np.maximum(1e-3, old_dens))
                
           
            ''' Optimality criteria update scheme'''
            ##### implementing the bisection algorithm to predict the lambda value ######
            l1=0 
            l2=1e9
            forward=0.2
            dens_backward= old_dens-forward
            dens_forward= old_dens+forward
            sqrt_d_Jelem= np.sqrt(-d_Jelem/ d_vol)
            while (l2-l1)/(l1+l2)>1e-3:
                lmid=0.5*(l2+l1)
                new_dens= np.maximum(0.0,np.maximum(dens_backward, np.minimum(1.0,np.minimum(dens_forward, old_dens*sqrt_d_Jelem/np.sqrt(lmid)))))
                
                if self.params.filter==1:
                    phy_dens= self.H.dot(new_dens/self.HS)
                
                elif self.params.filter==2:
                    phy_dens=new_dens     
                            
                if np.sum(phy_dens)>self.params.volfrac*self.params.num_elems:
                    l1=lmid
                else :
                    l2=lmid
            self.topology_visualizer.add_view(phy_dens)

            difference=np.max(abs(new_dens-old_dens))
            old_dens=new_dens
        print(Jelem)
        return phy_dens
    
    def density_filter2(self):
        # ToDo add comment here
        centroids = self.solver.helper.get_element_centers()
        self.H= np.zeros((self.params.num_elems, self.params.num_elems))
        self.HS= np.empty(self.params.num_elems)

        for i in range(self.params.num_elems-1):
            curr_centroid = centroids[i]
            dist= np.linalg.norm(curr_centroid - centroids[i+1 :], axis=1)
            j= np.where(dist< self.params.rmin)[0]
            val= self.params.rmin - dist[j]
            self.H[i, i+j+1]= val 

        self.H+= self.H.T - np.diag(np.diag(self.H))
        self.HS = np.sum(self.H, axis=0)

        

    def density_filter(self):
        '''
        preparing the filter function
        '''
        rmin= self.params.rmin
        nelx, nely, nelz= self.params.nelx, self.params.nely, self.params.nelz

        self.ih=[1]*self.params.num_elems*(2*(int(np.ceil(rmin)-1))+1)**2
        self.jh=[1]*self.params.num_elems*(2*(int(np.ceil(rmin)-1))+1)**2
        self.sh=[0]*len(self.ih)
        self.cn=0                  # counter 
        for k1 in range(1,nelz+1):
            for i1 in range(1,nelx+1):
                for j1 in range(1,nely+1):
                    self.e1=(k1-1)*nelx*nely + (i1-1)*nely + (j1 -1)
                    for k2 in range  (np.maximum(k1-math.floor(rmin),1),np.minimum(k1+math.floor(rmin),nelz)+1):
                        for i2 in range(np.maximum(i1-math.floor(rmin),1),np.minimum(i1+math.floor(rmin),nelx)+1):
                            for j2 in range(np.maximum(j1-math.floor(rmin),1),np.minimum(j1+math.floor(rmin),nely)+1):
                               self.e2=((k2-1)*nelx*nely) +(i2-1)*nely+ (j2 -1)
                               if self.cn<self.params.num_elems*(2*(int(np.ceil(rmin)-1))+1)**2:  
                                    self.ih[self.cn]=self.e1                        # row indexing
                                    self.jh[self.cn]=self.e2                        # column indexing
                                    self.sh[self.cn]=np.maximum(0,rmin-math.sqrt((i1-i2)**2 + (j1-j2)**2 +(k1-k2)**2))
                               else:
                                    self.ih.append(self.e1)
                                    self.jh.append(self.e2)
                                    self.sh.append(np.maximum(0,rmin-math.sqrt((i1-i2)**2 + (j1-j2)**2 +(k1-k2)**2)))

                               self.cn=self.cn+1

        self.row=np.array(self.ih)
        self.column=np.array(self.jh)
        self.val=np.array(self.sh)                       
        self.H=coo_matrix((self.val,(self.row,self.column)),shape=(nelx*nely*nelz,nelx*nely*nelz))
        self.HS=np.squeeze(np.asarray(self.H.sum(axis=1)))

    def display_topology(self):
        self.topology_visualizer.visualize()



if __name__=='__main__':
    
    params = Parameters()
    op = SimpOptimizer(params)
    dens= op.densityestimation()
    op.display_topology()