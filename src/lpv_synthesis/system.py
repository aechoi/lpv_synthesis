import sympy as sp
import numpy as np
from collections import namedtuple
import control as c

class System:
    def __init__(self,A,B,C,D,nv,ne,ny,nw,nd,nu):
        """A,B,C,D are sympy matrices, ns are int"""
        """Creates the initial system object"""
        
        self.A=A
        self.B=B
        self.C=C
        self.D=D
        self.nv = int(nv)
        self.ne = int(ne)
        self.ny = int(ny)
        self.nw = int(nw)
        self.nd = int(nd)
        self.nu = int(nu)
        self.B1=B[:,:nw]
        self.B2=B[:,nw:nw+nd]
        self.B3=B[:,nw+nd:]
        self.C1=C[:nv,:]
        self.C2=C[nv:nv+ne,:]
        self.C3=C[nv+ne:,:]
        self.D11=D[:nv,:nw]
        self.D12=D[:nv,nw:nw+nd]
        self.D13=D[:nv,nw+nd:]
        self.D21=D[nv:nv+ne,:nw]
        self.D22=D[nv:nv+ne,nw:nw+nd]
        self.D23=D[nv:nv+ne,nw+nd:]
        self.D31=D[nv+ne:,:nw]
        self.D32=D[nv+ne:,nw:nw+nd]
        self.D33=D[nv+ne:,nw+nd:]
        
        

    def generate_Grho_tilde(self,Psi11,Psi22):
        """Psi11 and Psi22 are control library transfer function and state space objects"""
        """Assigns Grho_tilde as an attribute which is a named tuple with its own attributes: A,B,C,D"""

        nv=self.nv
        ne=self.ne
        ny=self.ny
        nw=self.nw
        nd=self.nd
        nu=self.nu

        #generate the matrix that left multiplies Grho in equation 18
        Psi11ss=c.ss(Psi11)
        Psi11ss=Psi11ss #.minreal()

        Psi11ssA=sp.Matrix(Psi11ss.A)
        Psi11ssB=sp.Matrix(Psi11ss.B)
        Psi11ssC=sp.Matrix(Psi11ss.C)
        Psi11ssD=sp.Matrix(Psi11ss.D)

        nx_psi11=sp.shape(Psi11ssA)[0]

        Psi11_B_expanded=sp.BlockMatrix([Psi11ssB,sp.zeros(nx_psi11,nv)])
        Psi11_C_expanded=sp.BlockMatrix([[Psi11ssC],
                                         [sp.zeros(ne+ny,nx_psi11)]])
        Psi11_D_expanded=sp.BlockMatrix([[Psi11ssD,sp.zeros(nv,ne+ny)],
                                         [sp.zeros(ne+ny,nv),sp.eye(ne+ny)]])

        #generate the matrix that left multiplies Grho in equation 18
        Psi22ss=c.ss(Psi22)
        Psi22ss=Psi22ss #.minreal()
        #invert the state space
        Psi22Ainv=Psi22ss.A-Psi22ss.B@np.linalg.inv(Psi22ss.D)@Psi22ss.C
        Psi22Binv=Psi22ss.B@np.linalg.inv(Psi22ss.D)
        Psi22Cinv=-1*np.linalg.inv(Psi22ss.D)@Psi22ss.C
        Psi22Dinv=np.linalg.inv(Psi22ss.D)

        Psi22Ainv=sp.Matrix(Psi22Ainv)
        Psi22Binv=sp.Matrix(Psi22Binv)
        Psi22Cinv=sp.Matrix(Psi22Cinv)
        Psi22Dinv=sp.Matrix(Psi22Dinv)
        nx_psi22inv=sp.shape(Psi22Ainv)[0]

        Psi22inv_B_expanded=sp.BlockMatrix([Psi22Binv,sp.zeros(nx_psi22inv,nw)])
        
        Psi22inv_C_expanded=sp.BlockMatrix([[Psi22Cinv],
                                          [sp.zeros(nd+nu,nx_psi22inv)]])
        Psi22inv_D_expanded=sp.BlockMatrix([[Psi22Dinv,sp.zeros(nw,nd+nu)],
                                         [sp.zeros(nd+nu,nw),sp.eye(nd+nu)]])

        #Implement equation  18
        Grho_psi22inv_product=self.sys1_tosys2_seriesconnect(Psi22Ainv,Psi22inv_B_expanded,Psi22inv_C_expanded,Psi22inv_D_expanded,self.A,self.B,self.C,self.D)
        Grhotilde=self.sys1_tosys2_seriesconnect(Grho_psi22inv_product.A,Grho_psi22inv_product.B,Grho_psi22inv_product.C,Grho_psi22inv_product.D,Psi11ss.A,Psi11_B_expanded,Psi11_C_expanded,Psi11_D_expanded)

        return System(Grho_tilde.A,Grho_tilde.B,Grho_tilde.C,Grho_tilde.D,nv,ne,ny,nw,nd,nu)

    
    def sys1_tosys2_seriesconnect(self,A1,B1,C1,D1,A2,B2,C2,D2): # series connect 2 LTI/LPV systems together
        """This is a helper function for generate_Grho_tilde"""
        print(B1, sp.shape(B1))
        print(C1, sp.shape(C1))
        print(D1, sp.shape(D1))
        print(B2, sp.shape(B2))
        Aseries=sp.BlockMatrix([[A1, sp.zeros(np.shape(A1)[0],sp.shape(A2)[1])],
                         [B2@C1,A2]])
        Bseries=sp.BlockMatrix([[B1],
                                [B2@D1]])
        Cseries=sp.BlockMatrix([D2@C1,C2])
        Dseries=D2@D1
        
        series_system = namedtuple('series_system', ['A', 'B','C','D'])
        series_system_tuple= series_system(Aseries,Bseries,Cseries,Dseries)
        
        return series_system_tuple

if __name__=="__main__":
    rho = sp.symbols('rho')
    mysys=System(sp.Matrix([rho]),sp.Matrix([[1,0,1]]),sp.Matrix([0]),sp.eye(3),1,1,1,1,1,1)
    mysys.generate_Grho_tilde(c.ss(0,0,0,1),c.ss(0,0,0,1))
    
        
