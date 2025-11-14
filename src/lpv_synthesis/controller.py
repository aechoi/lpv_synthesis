import numpy as np
import cvxpy as cp
import sympy as sp
from collections import namedtuple
import control as c
class Controller:
    def __init__(Grho,Psi11,Psi22,gamma,rho_lower_lim,rho_upper_lim,rho_samples):    

        # Employ Lemma 1 so that you can use sup_rho norm(Fl(Grho_tilde,Krho))<=gamma
        # in order to show that sup_delta norm(Fu(Fl(Grho,Krho),Delta))<=gamma holds.

        Grho_tilde=generate_Grho_tilde(Grho,Psi11,Psi22,Grho.ne,Grho.ny,Grho.nd,Grho.nu,Grho.nv,Grho.nw)# Equation 18

        ne1=Grho.nv
        ne2=Grho.ne
        ny=Grho.ny
        nd1=Grho.nw
        nd2=Grho.nd
        nu=Grho.nu
    
        A=Grho_tilde.A
        nx=np.shape(Grho_tilde.A)[0]

        B=Grho_tilde.B
        B11=B[:,:nd1]
        B12=B[:,nd1:nd1+nd2]
        B1=B[:,:nd1+nd2]
        B2=B[:,nd1+nd2:]

        C=Grho_tilde.C
        C11=C[:ne1,:]
        C12=C[ne1:ne1+ne2,:]
        C1=C[:ne1+ne2,:]
        C2=C[ne1+ne2:,:]


        # Now, solve theorem 1 for P and Q (or X and Y using ref 4 notation) using SDP
        #But all matrices in the LMI are a non-convex function of rho. So I guess we
        # must sample for rho in the set of possible rho and solve the LMIs for all
        #possible rhos.

        #Define optimization problem accordingly:
        X=cp.Variable((nx,nx), symmetric=True)
        Y=cp.Variable((nx,nx), symmetric=True)

        LMI_11c=cp.bmat([
            [X,np.eye(nx)],
            [np.eye(nx),Y]])
        
        constraints=[LMI_11c>>0]

        rho_vec=np.linspace(rho_lower_lim,rho_upper_lim,rho_samples)

        for rho_val in rho_vec: 

            a=A.subs(rho_val)
            
            b1=B1.subs(rho_val)
            b11=B11.subs(rho_val)
            b12=B12.subs(rho_val)
            b2=B2.subs(rho_val)

            c1=C1.subs(rho_val)
            c11=C11.subs(rho_val)
            c12=C12.subs(rho_val)
            c2=C2.subs(rho_val)
                       
            Ahat=a-b2@c12
            Atilde=a-b12@c2
            
            LMI_11a=cp.bmat([
            [Y@Ahat.T+Ahat@Y-gamma*b2@b2.T, Y@c11.T           , b1],
            [c11@Y                        , -gamma*np.eye(ne1), np.zeros((ne1,nd))], 
            [b1.T                         , np.zeros((nd,ne1)), -gamma*np.eye(nd)]])
            
            LMI_11b=cp.bmat([
            [Atilde.T@X+X@Atilde-gamma*c2.T@c2, X@b11             , c1.T],
            [b11.T@X                          , -gamma*np.eye(nd1), np.zeros((nd1,ne))], 
            [c1                               , np.zeros((ne,nd1)), -gamma*np.eye(ne)]])
            
            constraints.append(LMI_11a<<0)
            constraints.append(LMI_11b<<0)


        prob = cp.Problem(cp.Minimize(0.0),constraints)  # Depending on this, the result may be different to that of the paper, play around with this.
        prob.solve()

        print("Problem status:", prob.status)
        print("Solution X is")
        print(X.value)
        print("Solution Y is")
        print(Y.value)

        #Based on reference 4, use X and Y to calcaulte Krho in the form Ak,Bk,Ck,Dk
        Yinv=np.linalv.inv(Y)
        Xinv=np.linalv.inv(X)
        Q=X-Yinv
        Qinv=np.linalg.inv(Q)

        F=-(gamma*B2.T@Yinv+C12)
        L=-(gamma*Xinv@C2.T+B12)
        Af=A+B2@F
        Cf=C1.T+F.T
        H=-(Af.T@Yinv+Yinv@Af+1/gamma*Cf.T@Cf+1/gamma*Yinv@B1@B1.T@Yinv)

        Ak=A+1/gamma*(Qinv@X@L@B12.T+B1@B1.t)@Yinv+B2@F+Qinv@X@L@C2-Qinv@H
        Bk=-Qinv@X@L
        Ck=F
        Dk=sp.zeros(nu,ny)

        #Return controller object
        self.Ak=Ak
        self.Bk=Bk
        self.Ck=Ck
        self.Dk=Dk

    def sys1_tosys2_seriesconnect(A1,B1,C1,D1,A2,B2,C2,D2): # series connect 2 LTI/LPV systems together

        Aseries=sp.BlockMatrix([[A1, sp.zeros(np.shape(A1)[0],np.shape(A2)[1])],
                                [B2@C1,A2]])
        Bseries=sp.BlockMatrix([[B1],
                                [B2@D1]])
        Cseries=sp.BlockMatrix([D2@C1,C2])
        Dseries=D2@D1
        
        series_system = namedtuple('series_system', ['A', 'B','C','D'])
        series_system= namedtuple(Aseries,Bseries,Cseries,Dseries)
        
        return series_system

    def generate_Grho_tilde(Grho,Psi11,Psi22,ne,ny,nd,nu,nv,nw):

        #generate the matrix that left multiplies Grho in equation 18
        Psi11ss=c.ss(Psi11)
        Psi11ss=Psi11ss.minreal()
        nx_psi11=np.shape(Psi11ss.A)[0]

        Psi11_B_expanded=sp.BlockMatrix([Psi11ss.B,sp.zeros(nx_psi11,nv)])
        Psi11_C_expanded=sp.BlockMatrix([[Psi11ss.C],
                                         [sp.zeros(ne+ny,nx_psi11)]])
        Psi11_D_expanded=sp.BlockMatrix([[Psi11ss.D,sp.zeros(nv,ne+ny)],
                                         [sp.zeros(ne+ny,nv),np.eye(ne+ny)]])

        #generate the matrix that left multiplies Grho in equation 18
        Psi22ss=c.c.ss(Psi22).minreal()
        Psi22ss=c.minreal(Psi22ss)
        #invert the state space
        Psi22Ainv=Psi22ss.A-Psi22ss.B@np.linalg.inv(Psi22ss.D)-Psi22ss.C
        Psi22Binv=Psi22ss.B@np.linalg.inv(Psi22ss.D)
        Psi22Cinv=-1*np.linalg.inv(Psi22ss.D)@Psi22ss.C
        Psi22Dinv=np.linalg.inv(Psi22ss.D)
        
        Pssi22invss=c.ss(Psi22Ainv,Psi22Binv,Psi22Cinv,Psi22Dinv)
        
        Psi22invss=c.minreal(Psi22invss)
        
        nx_psi22inv=np.shape(Psi22invss.A)[0]

        Psi22inv_B_expanded=sp.BlockMatrix([Psi22invss.B,sp.zeros(nx_psi22inv,nw)])
        Psi22inv_C_expanded=sp.BlockMatrix([[Psi22invss.C],
                                          [sp.zeros(nd+nu,nx_psi22inv)]])
        Psi22inv_D_expanded=sp.BlockMatrix([[Psi22invss.D,sp.zeros(nw,nd+nu)],
                                         [sp.zeros(nd+nu,nw),np.eye(nd+nu)]])

        #Implement equation  18
        Grho_psi22inv_product=sys1_tosys2_seriesconnect(Pssi22invss.A,Psi22inv_B_expanded,Psi22inv_C_expanded,Psi22inv_D_expanded,Grho.A,Grho.B,Grho.C,Grho.D)
        Grhotilde=sys1_tosys2_seriesconnect(Grho_psi22inv_product.A,Grho_psi22inv_product.B,Grho_psi22inv_product.C,Grho_psi22inv_product.D,Psi11ss.A,Psi11_B_expanded,Psi11_C_expanded,Psi11_D_expanded)

        return Grhotilde
        

    
