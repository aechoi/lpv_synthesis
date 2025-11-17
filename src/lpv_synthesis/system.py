from collections import namedtuple
import copy

import sympy as sp
import numpy as np
import control as c


def make_block_property(matrix_name, r_slice, c_slice):
    def getter(self):
        return getattr(self, matrix_name)[r_slice, c_slice]

    def setter(self, value):
        getattr(self, matrix_name)[r_slice, c_slice] = value

    return property(getter, setter)


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

        self.params = [A, B, C, D, nv, ne, ny, nw, nd, nu]

        self._B_cols = [0, nw, nw + nd, nw + nd + nu]
        self._C_rows = [0, nv, nv + ne, nv + ne + ny]
        self._make_block_properties()

    def __str__(self):
        A_str = sp.srepr(self.A) if self.A is None else sp.pretty(self.A)
        B_str = sp.srepr(self.B) if self.B is None else sp.pretty(self.B)
        C_str = sp.srepr(self.C) if self.C is None else sp.pretty(self.C)
        D_str = sp.srepr(self.D) if self.D is None else sp.pretty(self.D)

        # Split into lines
        A_lines = A_str.split("\n")
        B_lines = B_str.split("\n")
        C_lines = C_str.split("\n")
        D_lines = D_str.split("\n")

        # Determine consistent row heights
        top_height = max(len(A_lines), len(B_lines))
        bot_height = max(len(C_lines), len(D_lines))

        # Determine column widths
        left_width = max(
            max(len(line) for line in A_lines), max(len(line) for line in C_lines)
        )
        right_width = max(
            max(len(line) for line in B_lines), max(len(line) for line in D_lines)
        )

        # Helper pad
        pad = lambda s, w: s + " " * (w - len(s))

        # Build top and bottom rows
        top = [
            pad(A_lines[i] if i < len(A_lines) else "", left_width)
            + " │ "
            + pad(B_lines[i] if i < len(B_lines) else "", right_width)
            for i in range(top_height)
        ]

        sep = ["─" * (left_width - 1) + "─┼─" + "─" * (right_width - 1)]

        bottom = [
            pad(C_lines[i] if i < len(C_lines) else "", left_width)
            + " │ "
            + pad(D_lines[i] if i < len(D_lines) else "", right_width)
            for i in range(bot_height)
        ]

        return "\n".join(top + sep + bottom + ["\n"])

    def _make_block_properties(self):
        # B (1 x 3 blocks): [B1 B2 B3]
        for j, name in enumerate(("B1", "B2", "B3")):
            c0, c1 = self._B_cols[j], self._B_cols[j + 1]
            setattr(
                self.__class__,
                name,
                make_block_property("B", slice(None), slice(c0, c1)),
            )

        # C (3 x 1 blocks): [C1; C2; C3]
        for i, name in enumerate(("C1", "C2", "C3")):
            r0, r1 = self._C_rows[i], self._C_rows[i + 1]
            setattr(
                self.__class__,
                name,
                make_block_property("C", slice(r0, r1), slice(None)),
            )

        # D (3 x 3 blocks)
        for i, rname in enumerate(("1", "2", "3")):
            r0, r1 = self._C_rows[i], self._C_rows[i + 1]
            for j, cname in enumerate(("1", "2", "3")):
                c0, c1 = self._B_cols[j], self._B_cols[j + 1]
                propname = f"D{rname}{cname}"
                setattr(
                    self.__class__,
                    propname,
                    make_block_property("D", slice(r0, r1), slice(c0, c1)),
                )

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

        Psi11_B_expanded = sp.BlockMatrix([Psi11ssB, sp.zeros(nx_psi11, ne + ny)])
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

        Psi22inv_B_expanded = sp.BlockMatrix(
            [Psi22Binv, sp.zeros(nx_psi22inv, nd + nu)]
        )

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
        Grhotilde = self.sys1_tosys2_seriesconnect(
            Grho_psi22inv_product.A,
            Grho_psi22inv_product.B,
            Grho_psi22inv_product.C,
            Grho_psi22inv_product.D,
            Psi11ssA,
            Psi11_B_expanded,
            Psi11_C_expanded,
            Psi11_D_expanded,
        )
        return System(
            Grhotilde.A, Grhotilde.B, Grhotilde.C, Grhotilde.D, nv, ne, ny, nw, nd, nu
        )

    def sys1_tosys2_seriesconnect(
        self, A1, B1, C1, D1, A2, B2, C2, D2
    ):  # series connect 2 LTI/LPV systems together
        """This is a helper function for generate_Grho_tilde"""
        print(A1, sp.shape(A1))
        print(B1, sp.shape(B1))
        print(C1, sp.shape(C1))
        print(D1, sp.shape(D1))
        print(B2, sp.shape(B2))

        Aseries = sp.BlockMatrix(
            [[A1, sp.zeros(sp.shape(A1)[0], sp.shape(A2)[1])], [B2 @ C1, A2]]
        )
        Bseries = sp.BlockMatrix([[B1], [B2 @ D1]])
        Cseries = sp.BlockMatrix([D2 @ C1, C2])
        Dseries = D2 @ D1

        series_system = namedtuple("series_system", ["A", "B", "C", "D"])
        series_system_tuple = series_system(Aseries, Bseries, Cseries, Dseries)

        return series_system_tuple


def full_svd(A: sp.Matrix) -> tuple:
    """Return the full SVD of a sympy matrix A with zeros listed first.

    Sympy only provides the condensed SVD. The full SVD is needed for the
    preprocessing steps. This function computes a full SVD by padding with
    orthogonal matrices via QR decomposition. The SVD is formatted in such a
    way that the singular value matrix, S, has the zeros padded either above
    or to the left of the non-zero singular values (as needed in preprocessing)

    Args:
        A: The sympy matrix to decompose.

    Returns:
        A tuple containing the full U, S, and V matrices of the SVD."""
    U, S, V = A.singular_value_decomposition()  # condensed
    # U: m x r, S: r x r diagonal, V: n x r
    m, r, n = U.rows, S.rows, V.rows

    # augment V to n x n
    N = V.T.nullspace()
    if N:
        # V_full = V.row_join(sp.Matrix.hstack(*N))
        V_full = sp.Matrix.hstack(*N).row_join(V)
        V_full = V_full.QRdecomposition()[0]
    else:
        V_full = V

    # augment U to m x m
    N = U.T.nullspace()
    if N:
        U_full = sp.Matrix.hstack(*N).row_join(U)
        # U.row_join(sp.Matrix.hstack(*N))
        U_full = U_full.QRdecomposition()[0]
    else:
        U_full = U

    # # pad S with zeros to (m x n)
    # S_full = sp.zeros(m, n)
    # S_full[m - r :, n - r :] = S

    # the zero part of S is never used, so don't augment
    return U_full, S, V_full


def simplify_system(system):
    """Normalize the D matrices of the lpv system.

    Assume you have some A, B1, B2, C1, C2, D11, D12, D21, D22 matrices which
    are functions of a parameter vector rho. We want to create another system
    object which is transformed from the original such that D11=D22=0 and D12
    is [0, I_nu] and D21 is [0, I_ny]. It should still be parameterized by rho.

    Args:
        system: The LPV system to be simplified. inputs x, d, u; outputs dx, e, y

    Returns:
        The simplified LPV system and a function for retrieving the original
        variables.
    """
    # Systems are the interstitial systems between each change of variables
    # the variable transforms take the vectors, and transform them back to the
    # original system variables. This is dependent on rho.
    system_1, var_transform_1 = cov_1(system)
    system_2, var_transform_2 = cov_2(system_1)
    system_3, var_transform_3 = cov_3(system_2)
    system_4, var_transform_4 = cov_4(system_3)
    system_5, var_transform_5 = cov_5(system_4)
    system_6, var_transform_6 = cov_6(system_5)

    def var_transform(io):
        io_5 = var_transform_6(io)
        io_4 = var_transform_5(io_5)
        io_3 = var_transform_4(io_4)
        io_2 = var_transform_3(io_3)
        io_1 = var_transform_2(io_2)
        io_orig = var_transform_1(io_1)
        return io_orig

    return system_6, var_transform


def cov_1(sys: System):
    """The first change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    D22_1 = sp.zeros(
        sys.ne, sys.nd
    )  # maybe need to update this. Need to figure out what the 2i/o plant format is
    new_sys.D22 = D22_1

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        io[3] = sys.D22 @ io[1] + io[3]
        return io

    return new_sys, io_transform


def cov_2(sys: System):
    """The second change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    U12, S12, V12 = full_svd(sys.D12)
    U21, S21, V21 = full_svd(sys.D21)

    Su1 = V21
    Su2 = V12 @ S12.inv()
    Sy1 = U12.T
    Sy2 = S21.inv() @ U21.T

    new_sys.B1 = sys.B1 @ Su1
    new_sys.B2 = sys.B2 @ Su2
    new_sys.C1 = Sy1 @ sys.C1
    new_sys.C2 = Sy2 @ sys.C2
    new_sys.D11 = Sy1 @ sys.D11 @ Su1
    new_sys.D12 = Sy1 @ sys.D12 @ Su2
    new_sys.D21 = Sy2 @ sys.D21 @ Su1
    new_sys.D22 = Sy2 @ sys.D22 @ Su2

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        io[0] = Su1 @ io[0]
        io[1] = Su2 @ io[1]
        io[2] = Sy1.T @ io[2]
        io[3] = Sy2.inv() @ io[3]
        return io

    return new_sys, io_transform


def cov_3(sys: System):
    """The first change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    dim_u = sys.D12.rank()
    dim_y = sys.D21.rank()

    r, c = sys.D11.shape
    D1111 = sys.D11[: r - dim_u, : c - dim_y]
    D1112 = sys.D11[: r - dim_u, -dim_y:]
    D1121 = sys.D11[-dim_u:, : c - dim_y]
    D1122 = sys.D11[-dim_u:, -dim_y:]

    K_inf = -(
        D1122 + D1121 @ (sp.eye(c - dim_y) - D1111.T @ D1111).inv() @ D1111.T @ D1112
    )

    D11 = sys.D11
    D11[-dim_u:, -dim_y:] = D1122 + K_inf
    new_sys.D11 = D11

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        io[1] = io[1] + K_inf @ io[3]
        return io

    return new_sys, io_transform


def cov_4(sys):
    """The fourth change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    X = sys.D11

    new_sys.A = sys.A + sys.B1 @ X.T @ (sp.eye(X.shape[0]) - X @ X.T).inv() @ sys.C1
    new_sys.B1 = sys.B1 @ (sp.eye(X.shape[0]) - X.T @ X) ** (-sp.S(1) / 2)
    new_sys.B2 = (
        sys.B2
        + sys.B1 @ X.T @ (sp.eye(X.shape[0]) - X @ X.T) ** (-sp.S(1) / 2) @ sys.D12
    )
    new_sys.C1 = (sp.eye(X.shape[0]) - X @ X.T) ** (-sp.S(1) / 2) @ sys.C1
    new_sys.C2 = sys.C2 + sys.D21 @ X.T @ (sp.eye(X.shape[0]) - X @ X.T).inv() @ sys.C1
    new_sys.D11 = sp.zeros(*sys.D11.shape)
    new_sys.D12 = (sp.eye(X.shape[0]) - X @ X.T) ** (-sp.S(1) / 2) @ sys.D12
    new_sys.D21 = sys.D21 @ (sp.eye(X.shape[0]) - X.T @ X) ** (-sp.S(1) / 2)
    new_sys.D22 = sys.D21 @ X.T @ (sp.eye(X.shape[0]) - X @ X.T).inv() @ sys.D12

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        new_io = copy.deepcopy(io)

        new_io[0] = (
            X.T @ (sp.eye(X.shape[0]) - X @ X.T) ** (-sp.S(1) / 2) @ io[2]
            + (sp.eye(X.shape[1]) - X.T @ X) ** (-sp.S(1) / 2) @ io[0]
        )
        new_io[2] = (sp.eye(X.shape[0]) - X @ X.T) ** (-sp.S(1) / 2) @ io[2] + (
            sp.eye(X.shape[0]) - X @ X.T
        ) ** (-sp.S(1) / 2) @ X @ io[0]
        return new_io

    return new_sys, io_transform


def cov_5(sys):
    """The fifth change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    new_sys.D22 = sp.zeros(sys.ne, sys.nd)  # same issue as cov_1

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        io[3] = sys.D22 @ io[1] + io[3]
        return io

    return new_sys, io_transform


def cov_6(sys):
    """The last change of variables for the LPV system simplification."""
    new_sys = copy.deepcopy(sys)

    U12, S12, V12 = full_svd(sys.D12)
    U21, S21, V21 = full_svd(sys.D21)

    Su1 = V21
    Su2 = V12 @ S12.inv()
    Sy1 = U12.T
    Sy2 = S21.inv() @ U21.T

    new_sys.B1 = sys.B1 @ Su1
    new_sys.B2 = sys.B2 @ Su2
    new_sys.C1 = Sy1 @ sys.C1
    new_sys.C2 = Sy2 @ sys.C2
    new_sys.D11 = Sy1 @ sys.D11 @ Su1
    new_sys.D12 = Sy1 @ sys.D12 @ Su2
    new_sys.D21 = Sy2 @ sys.D21 @ Su1
    new_sys.D22 = Sy2 @ sys.D22 @ Su2

    # [u1; u2; y1; y2] = [d; u; e; y]
    def io_transform(io):
        io[0] = Su1 @ io[0]
        io[1] = Su2 @ io[1]
        io[2] = Sy1.T @ io[2]
        io[3] = Sy2.inv() @ io[3]
        return io

    return new_sys, io_transform


###

if __name__ == "__main__":
    rho = sp.symbols("rho")
    mysys = System(
        sp.Matrix([rho]),
        sp.Matrix([[1, 0, 1]]),
        sp.Matrix([[0], [1], [0]]),
        sp.eye(3),
        1,
        1,
        1,
        1,
        1,
        1,
    )
    # Grho_tilde = mysys.generate_Grho_tilde(c.ss(0, 0, 0, 1), c.ss(0, 0, 0, 1))
    print(mysys.A, mysys.B, mysys.C, mysys.D)

    new_sys, var_transform = simplify_system(mysys)
