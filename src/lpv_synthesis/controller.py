import numpy as np
import cvxpy as cp
import sympy as sp

class Controller:
    def __init__(Grho,Psi,gamma,nv,nw,ne1,ne2,ny,nd1,nd2,nu,nx,rho_lower_lim,rho_upper_lim,rho_samples,A,B1,B2,C1,C2):    

        nd=dn1+nd2
        ne=ne1+ne2
        
        Psi11=Psi[:nv,:nv]
        Psi22=Psi[nv:,nv:]
        Psi22inv=np.linalg.inv(Psi22)
        
        mat1=np.block([[Psi11,np.zeros((nv,ne+ny))],
                      [np.zeros((ne+ny,nv)),np.eye(ne+ny)]])
        mat2=np.block([[Psi22inv,np.zeros((nw,nd+nu))],
                      [np.zeros((nd+nu,nw)),np.eye(nd+nu)]])

        Grho_tilde=mat1@Grho@mat2    # Equation 18

        # Employ Lemma 1 so that you can use sup_rho norm(Fl(Grho_tilde,Krho))<=gamma
        # in order to show that sup_delta norm(Fu(Fl(Grho,Krho),Delta))<=gamma holds.

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

        for rho in rho_vec:   #should be for Grho_tilde
            a=A(rho)
            b1=B1(rho)
            b2=B2(rho)
            c1=C1(rho)
            c2=C2(rho)

            c11=c1[:ne1,:]
            c12=c1[ne1:,:]
            b11=c1[:,:nd1]
            b12=c1[:nd1,:]
            
            
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
