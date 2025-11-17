"""Main loop for control synthesis of an LPV system with IQC defined uncertainty.

Details taken from the paper "Robust Synthesis for Linear Parameter Varying Systems
Using Integral Quadratic Constraints" by S. Wang, H. Pfifer, and P. Seiler.
"""

import control as ct
#from lpv_synthesis import iqc
import system
import numpy as np
import sympy as sp


def robust_synthesis(
    plant_dynamics: ct.StateSpace,
    convergence_tolerance: float = 1e-3,
    max_iterations=50,
) -> None:
    """Main loop for control synthesis

    Args:
        plant_dynamics: The state-space representation of the plant dynamics. Should be
            in some form that allows for extraction of i/o dimensions of interconnect
            components.
        convergence_tolerance: The tolerance for convergence of the RP level.
        max_iterations: The maximum number of iterations to perform before giving up.
    """

    last_rp_level = float("inf")
    iqc_filters = []
    rp_levels = []

    iqc_filter = iqc.get_initializing_filter(plant_dynamics)

    for _ in range(max_iterations):
        controller = synthesize_controller(iqc_filter)
        iqc_filter, rp_level = analyze(controller)
        iqc_filters.append(iqc_filter)
        rp_levels.append(rp_level)

        if abs(rp_level - last_rp_level) < convergence_tolerance:
            break
        last_rp_level = rp_level
    else:
        print(
            f"Warning: Synthesis did not converge after {max_iterations} iterations. Final RP level: {rp_level}. Final RP difference: {abs(rp_level - last_rp_level)}"
        )
    return controller, (iqc_filters, rp_levels)


if __name__ == "__main__":

    #Create Grho as described in the paper:
    #e=We*y
    #v=Wu*u
    #y=G*(u+v)+d

    # This leads to the tf table:

    #    ##  w   #   d   #   u   #
    ##############################
    # v  ##  0   #   0   #  Wu   #
    # e  ## We*G #  We   # We*G  #
    # y  ##  G   #   1   #  G    #

    s=ct.tf([1,0],[1])

    Wu=(s+1)/(s+100)
    We=(0.1*s+10)/(s+1)
    Wu=ct.ss(Wu) #.minreal
    We=ct.ss(We) #.#minreal

    WuA=sp.Matrix(Wu.A)
    WuB=sp.Matrix(Wu.B)
    WuC=sp.Matrix(Wu.C)
    WuD=sp.Matrix(Wu.D)
    nxWu=sp.shape(WuA)[0]

    WeA=sp.Matrix(We.A)
    WeB=sp.Matrix(We.B)
    WeC=sp.Matrix(We.C)
    WeD=sp.Matrix(We.D)
    nxWe=sp.shape(WeA)[0]

    #Define G as in the paper - this is not Grho
    rho = sp.symbols('rho')
    GA=sp.Matrix([-1/(sp.sqrt(133.6-16.8*rho))])
    GB=sp.Matrix([1/(sp.sqrt(133.6-16.8*rho))])
    GC=sp.Matrix([sp.sqrt(4.8*rho-8.6)])
    GD=sp.Matrix([0.0])
    nxG=1

    #Define Guwd system (the 3rd row of the transfer function matrix above)
    GuwdA=GA
    GuwdB=sp.Matrix([[GB,sp.Matrix([0]),GB]])
    GuwdC=sp.Matrix([GC])
    GuwdD=sp.Matrix([[GD,sp.Matrix([1]),GD]])

    #We require the state equtions for Wu, Guwd but also interconnect We@Guwd
    
    WeG=system.System.sys1_tosys2_seriesconnect(GuwdA,GuwdB,GuwdC,GuwdD,WeA,WeB,WeC,WeD)
    nxWeG=sp.shape(WeG.A)[0]
    
    #Define Grho
    rho_arr=np.linspace(2,7,11)

    ne=1
    nv=1
    ny=1
    nw=1
    nd=1
    nu=1
    
    GrhoA=sp.BlockMatrix([[WuA,sp.zeros(nxWu,nxWe+nxG)],
                          [sp.zeros(2,1),WeG.A]])
    GrhoB=sp.BlockMatrix([[WuB,sp.zeros(nxWu,nd+nw)],
                          [WeG.B[:,0],WeG.B[:,1:3]]]) # had to split this up to put in block matrix form, similar in D
    GrhoC=sp.BlockMatrix([[WuC,sp.Matrix([0]),sp.Matrix([0])],
                          [sp.Matrix([0]),WeG.C[:,0],WeG.C[:,1]],
                          [sp.Matrix([0]),sp.Matrix([0]),GuwdC]])
    GrhoD=sp.BlockMatrix([[WuD,sp.zeros(nv,nd+nw)],
                          [WeG.D[:,0],WeG.D[:,1:3]],
                          [GuwdD[:,0],GuwdD[:,1:3]]]) 

    Grho=system.System(GrhoA,GrhoB,GrhoC,GrhoD,nv,ne,ny,nw,nd,nu,rho_arr)

    #Define Psi11 and Psi22 to generate Grhotilde
    Psi11=ct.ss([0],[0],[0],[1])
    Psi22=ct.ss([0],[0],[0],[1])

    Grho_tilde=Grho.generate_Grho_tilde(Psi11,Psi22)
    
    plant = ct.ss()
    controller, _ = robust_synthesis(plant)
