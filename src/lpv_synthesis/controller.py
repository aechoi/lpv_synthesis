import numpy as np
import cvxpy as cp
import sympy as sp
from collections import namedtuple
import control as c
class Controller:
    def __init__(self,Grho_tilde_modified,gamma):    
        """Obtain K for a given K(rho) as a sympy matrix given a RP level gamma"""
        """ The input system is Grho_tilde after Alex's change of basis"""
        
        # Employ Lemma 1 so that you can use sup_rho norm(Fl(Grho_tilde,Krho))<=gamma
        # in order to show that sup_delta norm(Fu(Fl(Grho,Krho),Delta))<=gamma holds.        
        rho_vec=list(Grho_tilde_modified.rho_arr)

        ne1=Grho_tilde_modified.nv
        ne2=Grho_tilde_modified.ne
        ny=Grho_tilde_modified.ny
        nd1=Grho_tilde_modified.nw
        nd2=Grho_tilde_modified.nd
        nu=Grho_tilde_modified.nu

        ne=ne1+ne2
        nd=nd1+nd2
    
        A=Grho_tilde_modified.A
        nx=sp.shape(A)[0]

        B=Grho_tilde_modified.B
        B11=B[:,:nd1]
        B12=B[:,nd1:nd1+nd2]
        B1=B[:,:nd1+nd2]
        B2=B[:,nd1+nd2:]
        
        C=Grho_tilde_modified.C
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


        for rho_val in rho_vec: 

            rho_val={rho: rho_val}

            a=np.array(A.subs(rho_val))
            
            b1=np.array(B1.subs(rho_val))
            b11=np.array(B11.subs(rho_val))
            b12=np.array(B12.subs(rho_val))
            b2=np.array(B2.subs(rho_val))

            c1=np.array(C1.subs(rho_val))
            c11=np.array(C11.subs(rho_val))
            c12=np.array(C12.subs(rho_val))
            c2=np.array(C2.subs(rho_val))
                       
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
        X=sp.Matrix(X)
        Y=sp.Matrix(Y)
        
        Yinv=Y.inv()
        Xinv=X.inv()
        Q=X-Yinv
        Qinv=Q.inv()

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

if __name__ == "__main__":
    from system import System
    rho = sp.symbols('rho')
    mysys=System(sp.Matrix([rho]),sp.Matrix([[1,0,1]]),sp.Matrix([[0],[1],[0]]),sp.eye(3),1,1,1,1,1,1,np.linspace(0.2,2,3))
    controller=Controller(mysys,0.8)
    
    

    

    
