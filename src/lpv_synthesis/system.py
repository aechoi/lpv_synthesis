# """Module for system-level struct"""

# from collections import namedtuple

# import numpy as np


# class InterconnectSystem:
#     """features needed

#     args for input and output dims
#     addressable matrices (sys.A(rho), sys.B1, etc)
#     callable (pass in inputs, get outputs)
#     transform into other forms (diff number of inputs/outputs)

#     """

#     def __init__(self, lpv_state_matrices, state_dim, input_dims, output_dims):
#         pass


# def simplify_system(system):
#     """Normalize the D matrices of the lpv system.

#     Assume you have some A, B1, B2, C1, C2, D11, D12, D21, D22 matrices which
#     are functions of a parameter vector rho. We want to create another system
#     object which is transformed from the original such that D11=D22=0 and D12
#     is [0, I_nu] and D21 is [0, I_ny]. It should still be parameterized by rho.

#     Args:
#         system: The LPV system to be simplified (plant).

#     Returns:
#         The simplified LPV system and a function for retrieving the original
#         variables.
#     """
#     # Systems are the interstitial systems between each change of variables
#     # the variable transforms take the vectors, and transform them back to the
#     # original system variables. This is dependent on rho.
#     system_1, var_transform_1 = cov_1(system)
#     system_2, var_transform_2 = cov_2(system_1)
#     system_3, var_transform_3 = cov_3(system_2)
#     system_4, var_transform_4 = cov_4(system_3)
#     system_5, var_transform_5 = cov_5(system_4)
#     system_6, var_transform_6 = cov_6(system_5)

#     def var_transform(io, rho):
#         io_5 = var_transform_6(io, rho)
#         io_4 = var_transform_5(io_5, rho)
#         io_3 = var_transform_4(io_4, rho)
#         io_2 = var_transform_3(io_3, rho)
#         io_1 = var_transform_2(io_2, rho)
#         io_orig = var_transform_1(io_1, rho)
#         return io_orig

#     return system_6, var_transform


# def cov_1(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims

#     D22_1 = np.zeros((ny, nu))

#     system.D22 = lambda rho: D22_1

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         io[3] = D22(rho) @ io[1] + io[3]
#         return io

#     return system, io_transform


# def cov_2(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims
#     nx = system.state_dim

#     def svd_gen(rho):
#         U12, d12, Vh12 = np.linalg.svd(D12(rho))
#         dim_diff12 = np.abs(len(U12) - len(Vh12))
#         Sigma12 = np.diag(d12)
#         u_reind12 = np.r_[np.arange(-dim_diff12, 0), np.arange(-len(U12), -dim_diff12)]
#         U12 = U12[:, u_reind12]

#         U21, d21, Vh21 = np.linalg.svd(D21)
#         dim_diff21 = np.abs(len(U21) - len(Vh21))
#         Sigma21 = np.diag(d21)
#         v_reind21 = np.r_[np.arange(-dim_diff21, 0), np.arange(-len(Vh21), -dim_diff21)]
#         Vh21 = Vh21[v_reind21, :]

#         Su1 = Vh21.T
#         Su2 = Vh12.T @ (1 / Sigma12)
#         Sy1 = U12.T
#         Sy2 = (1 / Sigma21) @ U21.T
#         return Su1, Su2, Sy1, Sy2

#     Su1 = lambda rho: svd_gen(rho)[0]
#     Su2 = lambda rho: svd_gen(rho)[1]
#     Sy1 = lambda rho: svd_gen(rho)[2]
#     Sy2 = lambda rho: svd_gen(rho)[3]

#     B1_2 = lambda rho: B1(rho) @ Su1(rho)
#     B2_2 = lambda rho: B2(rho) @ Su2(rho)
#     C1_2 = lambda rho: Sy1(rho) @ C1(rho)
#     C2_2 = lambda rho: Sy2(rho) @ C2(rho)
#     D11_2 = lambda rho: Sy1(rho) @ D11(rho) @ Su1(rho)
#     D12_2 = lambda rho: Sy1(rho) @ D12(rho) @ Su2(rho)
#     D21_2 = lambda rho: Sy2(rho) @ D21(rho) @ Su1(rho)
#     D22_2 = lambda rho: Sy2(rho) @ D22(rho) @ Su2(rho)

#     system_1 = InterconnectSystem(
#         lpv_state_matrices=(A, B1_2, B2_2, C1_2, C2_2, D11_2, D12_2, D21_2, D22_2),
#         state_dim=nx,
#         input_dims=(nd, nu),
#         output_dims=(ne, ny),
#     )

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         io[0] = Su1(rho) @ io[0]
#         io[1] = Su2(rho) @ io[1]
#         io[2] = Sy1(rho) @ io[2]
#         io[0] = Sy2(rho) @ io[3]

#     return system_1, io_transform


# def cov_3(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims

#     # Calcs
#     # dim_u2 =

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         return io

#     return system, io_transform


# def cov_4(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims

#     # Calcs

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         return io

#     return system, io_transform


# def cov_5(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims

#     # Calcs

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         return io

#     return system, io_transform


# def cov_6(system):
#     """The first change of variables for the LPV system simplification."""
#     A, B1, B2, C1, C2, D11, D12, D21, D22 = system.matrices
#     nd, nu = system.input_dims
#     ne, ny = system.output_dims

#     # Calcs

#     # [u1; u2; y1; y2] = [d; u; e; y]
#     def io_transform(io, rho):
#         return io

#     return system, io_transform
import sympy as sp
import numpy as np
from collections import namedtuple
import control as c


class System:
    def __init__(self, A, B, C, D, nv, ne, ny, nw, nd, nu):
        """A,B,C,D are sympy matrices, ns are int"""
        """Creates the initial system object"""

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.nv = int(nv)
        self.ne = int(ne)
        self.ny = int(ny)
        self.nw = int(nw)
        self.nd = int(nd)
        self.nu = int(nu)

    def generate_Grho_tilde(self, Psi11, Psi22):
        """Psi11 and Psi22 are control library transfer function and state space objects"""
        """Assigns Grho_tilde as an attribute which is a named tuple with its own attributes: A,B,C,D"""

        nv = self.nv
        ne = self.ne
        ny = self.ny
        nw = self.nw
        nd = self.nd
        nu = self.nu

        # generate the matrix that left multiplies Grho in equation 18
        Psi11ss = c.ss(Psi11)
        Psi11ss = Psi11ss  # .minreal()

        Psi11ssA = sp.Matrix(Psi11ss.A)
        Psi11ssB = sp.Matrix(Psi11ss.B)
        Psi11ssC = sp.Matrix(Psi11ss.C)
        Psi11ssD = sp.Matrix(Psi11ss.D)

        nx_psi11 = sp.shape(Psi11ssA)[0]

        Psi11_B_expanded = sp.BlockMatrix([Psi11ssB, sp.zeros(nx_psi11, nv)])
        Psi11_C_expanded = sp.BlockMatrix([[Psi11ssC], [sp.zeros(ne + ny, nx_psi11)]])
        Psi11_D_expanded = sp.BlockMatrix(
            [
                [Psi11ssD, sp.zeros(nv, ne + ny)],
                [sp.zeros(ne + ny, nv), sp.eye(ne + ny)],
            ]
        )

        # generate the matrix that left multiplies Grho in equation 18
        Psi22ss = c.ss(Psi22)
        Psi22ss = Psi22ss  # .minreal()
        # invert the state space
        Psi22Ainv = Psi22ss.A - Psi22ss.B @ np.linalg.inv(Psi22ss.D) @ Psi22ss.C
        Psi22Binv = Psi22ss.B @ np.linalg.inv(Psi22ss.D)
        Psi22Cinv = -1 * np.linalg.inv(Psi22ss.D) @ Psi22ss.C
        Psi22Dinv = np.linalg.inv(Psi22ss.D)

        Psi22Ainv = sp.Matrix(Psi22Ainv)
        Psi22Binv = sp.Matrix(Psi22Binv)
        Psi22Cinv = sp.Matrix(Psi22Cinv)
        Psi22Dinv = sp.Matrix(Psi22Dinv)
        nx_psi22inv = sp.shape(Psi22Ainv)[0]

        Psi22inv_B_expanded = sp.BlockMatrix([Psi22Binv, sp.zeros(nx_psi22inv, nw)])

        Psi22inv_C_expanded = sp.BlockMatrix(
            [[Psi22Cinv], [sp.zeros(nd + nu, nx_psi22inv)]]
        )
        Psi22inv_D_expanded = sp.BlockMatrix(
            [
                [Psi22Dinv, sp.zeros(nw, nd + nu)],
                [sp.zeros(nd + nu, nw), sp.eye(nd + nu)],
            ]
        )

        # Implement equation  18
        Grho_psi22inv_product = self.sys1_tosys2_seriesconnect(
            Psi22Ainv,
            Psi22inv_B_expanded,
            Psi22inv_C_expanded,
            Psi22inv_D_expanded,
            self.A,
            self.B,
            self.C,
            self.D,
        )
        Grho_tilde = self.sys1_tosys2_seriesconnect(
            Grho_psi22inv_product.A,
            Grho_psi22inv_product.B,
            Grho_psi22inv_product.C,
            Grho_psi22inv_product.D,
            Psi11ss.A,
            Psi11_B_expanded,
            Psi11_C_expanded,
            Psi11_D_expanded,
        )

        self.Grho_tilde = Grho_tilde

    def sys1_tosys2_seriesconnect(
        self, A1, B1, C1, D1, A2, B2, C2, D2
    ):  # series connect 2 LTI/LPV systems together
        """This is a helper function for generate_Grho_tilde"""
        print(B1, sp.shape(B1))
        print(C1, sp.shape(C1))
        print(D1, sp.shape(D1))
        print(B2, sp.shape(B2))
        Aseries = sp.BlockMatrix(
            [[A1, sp.zeros(np.shape(A1)[0], sp.shape(A2)[1])], [B2 @ C1, A2]]
        )
        Bseries = sp.BlockMatrix([[B1], [B2 @ D1]])
        Cseries = sp.BlockMatrix([D2 @ C1, C2])
        Dseries = D2 @ D1

        series_system = namedtuple("series_system", ["A", "B", "C", "D"])
        series_system_tuple = series_system(Aseries, Bseries, Cseries, Dseries)

        return series_system_tuple


if __name__ == "__main__":
    rho = sp.symbols("rho")
    mysys = System(
        sp.Matrix([rho]),
        sp.Matrix([[1, 0, 1]]),
        sp.Matrix([0]),
        sp.eye(3),
        1,
        1,
        1,
        1,
        1,
        1,
    )
    mysys.generate_Grho_tilde(c.ss(0, 0, 0, 1), c.ss(0, 0, 0, 1))
