

from common import gmsh_helper
from new import Rectangle_beam


class solver:
    def __init__(self, geometry: Rectangle_beam) -> None:
        self.geometry= geometry
        self.helper= gmsh_helper()
    
    '''
    input : shape function derivatives with respect to x, y ,z axis
    output: all stacked strain displacmenet matrix (B) together i.e. No. of B matrix is realted to  number of gauss points 
    '''
    def Bmat(self):
        self.a = []                        # it consists of all the B matrices(27) stacked together 
        for row in self.shapefunc_dertv.reshape(-1, 24):
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
            self.a.append(np.hstack(B))
        self.a=np.array(self.a)  

    '''
    it defines the threee dimensional constitutive matrix for an isotropic element
    '''

    def constitutive_matrix(self):
        
        C=np.zeros((6,6))
        C[0,0]=1-nu
        C[0,1]=nu
        C[0,2]=C[0,1]
        C[1,0]=C[0,1]
        C[1,1]=C[0,0]
        C[1,2]=C[0,1]
        C[3,3]=(1-2*nu)/2
        C[2,2]=C[0,0]
        C[4,4]=C[3,3]
        C[5,5]=C[3,3]
        C[2,0]=C[0,1]
        C[2,1]=C[0,1] 
        self.C= 1/((1+nu)*(1-2*nu))*C
    
    # def ref_element(self):
    #     k= np.zeros((24,24))
    #     for i,weight in enumerate(self.weights):
    #         k+= weight*np.matmul(self.a[i].T, np.matmul(self.C, self.a[i]))
    #     print(k[10,9])
    #     exit()
    '''
    obtaining individual element stiffeness matrices from Bmat and consitutive matrix functions
    '''

    def element_stiffness_matrix(self):
        self.ke = []
        for j in range(self.determinants.shape[0]):
            self.k=[]
            for i in range (self.a.shape[0]):
                self.k.append(self.determinants[j][i]*self.weights[i]*(np.matmul(np.transpose(self.a[i]),np.matmul(self.C, self.a[i]))))
            self.k= np.sum(self.k,axis = 0) 
            #self.k=np.array(self.k)
            self.ke.append(self.k)
        self.ke = np.array(self.ke)

    '''
    forming of global stiffness matrix by assembling all the element matrices which is obtained from the element_stiffness_matrix function
    '''
    def globalstiffness_matrix(self):
        self.kg=np.zeros((tdof,tdof))
        for i in range(self.ke.shape[0]):
            self.Nodes=np.vstack((self.my_nodes[i]*3-3,self.my_nodes[i]*3-2,self.my_nodes[i]*3-1)).flatten('F')
            self.Nodes=np.array(self.Nodes)
            x,y=np.meshgrid(self.Nodes,self.Nodes)
            self.kg[y,x]+= self.msimp[:,i]*self.ke[i]
    
    '''
    preparing the filter function
    '''
    def density_filter(self):
        self.ih=[1]*no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2
        self.jh=[1]*no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2
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
                               if self.cn<no_of_ele*(2*(int(np.ceil(rmin)-1))+1)**2:  
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
        self.H=coo_matrix((self.val,(self.row,self.column)),shape=(nelx*nely*nelz,nelx*nely*nelz)).tocsc()
        self.HS=np.sum(self.H,axis=1)
    '''
    initialisation of design and physical variables (densities)
    '''
    def densities(self):
        self.des_dens=np.ones((no_of_ele,1))*volfrac
        self.pyh_dens=self.des_dens                  # physical densities are assigned a constant and unifrom values(initially)
    '''
    defining the formula for the modified SIMP method(relation betwen density and youngs modulus)
    '''
    def simp_formula(self):
        self.msimp=Emin+(np.transpose(self.pyh_dens.reshape(nelx*nely*nelz,1))**p)*(E0-Emin)

    '''
    defining the nodal forces and nodal displacements
    '''
    def nodal_forces(self):
        self.F=np.zeros((tdof,1))
        self.F[self.forcedofyc[0],0]=-1
        self.F[self.forcedofyc[1],0]=-1
        self.F[self.forcedofyc[2],0]=-1

    def nodal_displacements(self):
        self.U=np.zeros((tdof,1))
        x,y=np.meshgrid(self.freedof,self.freedof)
        self.U[self.freedof]=np.linalg.solve(self.kg[y,x],self.F[self.freedof])
    
##############################################################################
        ''' 
            iteration  process

        '''
############################################################################

    def densityestimation(self):
 
        loop=0
        difference=1
        #while difference>tol and loop < max_loop:
        if difference>tol:
            if loop<max_loop:

                # update global stiffness using x.globalstiffness()
                # solve for nodal disp using x.nodaldisp()
                # complicance
                # flter
                # Update densities
                ''''minimum compliance objective function and sensitivity analysis'''
                comp = []

                for i in range(no_of_ele):
                    temp= np.dot(self.U[self.my_coord_edof[i]].T,np.dot(self.ke[i],self.U[self.my_coord_edof[i]]))
                    comp.append(temp)
                comp= np.array(comp).reshape(no_of_ele,1)
                #print(comp)
                der_comp=(comp*( -p*(E0-Emin)*self.pyh_dens**(p-1) ))
                der_vol=np.ones((no_of_ele,1))
                ''' use of filtering function to improve the sensitivity analysis  '''
                der_comp=self.H*( der_comp/self.HS)
                der_vol= self.H*( der_vol/self.HS)
                ''' Optimality criteria update scheme'''
                ##### implementing the bisection algorithm to predict the lambda value ######
                l1=0 
                l2=1e9
                forward=0.2
                while (l2-l1)/(l1+l2)>1e-3:
                    lmid=0.5*(l2+l1)
                    new_density= np.maximum(0.0,np.maximum(self.des_dens-forward,np.minimum(1.0,np.minimum(self.des_dens+forward,np.multiply(self.des_dens,np.sqrt(-(der_comp)/ der_vol/lmid))))))
                    self.pyh_dens= (self.H*new_density)/self.HS
                    if np.sum(self.pyh_dens)>volfrac*no_of_ele:
                        l1=lmid
                    else :
                        l2=lmid
                difference=abs(new_density-self.des_dens)
                self.des_dens=new_density
                loop=loop+1

