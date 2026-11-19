"""Hydrogen-like orbital math shared by the table generator and the tests.

Stdlib only. Lengths in Bohr radii.

Scaling
-------
Every radial function is tabulated for its *reference* charge Z0 = n, so in
the natural variable u = Z r / n all orbitals decay like exp(-2u). For an
orbital with effective charge Z let kappa = Z / n. Then

    R_Z(r)      = kappa^(3/2) R_ref(kappa r)
    rho_Z(x)    = kappa^3     rho_ref(kappa x)
    col_Z(p, d) = kappa^2     col_ref(kappa p, d)

Angular parts (real spherical harmonics), with x the position and r = |x|:
    s : Y^2 = 1 / (4 pi)
    p : Y^2 = (e . x)^2 / r^2         e = sqrt(3 / 4pi) * unit axis
    d : Y^2 = (x^T Q x)^2 / r^4        Q symmetric, normalisation folded in

Column densities from 1-D moment tables
---------------------------------------
For a ray x = p + t d (|d| = 1, p . d = 0, b = |p|) define

    I_k(b) = integral R(r)^2 t^k / r^(2l) dt,    r = sqrt(b^2 + t^2)

Odd powers of t integrate to zero, so (reference units):
    s : col = I_0 / (4 pi)
    p : col = (e.p)^2 I_0 + (e.d)^2 I_2
    d : col = a^2 I_0 + (4 beta^2 + 2 a c) I_2 + c^2 I_4
        a = p^T Q p,  beta = p^T Q d,  c = d^T Q d
"""
import math

# Table rows: one per radial function, in this order.
RADIAL = [(1, 0), (2, 0), (3, 0), (4, 0), (2, 1), (3, 1), (4, 1), (3, 2)]
ROW_NAMES = ["1s", "2s", "3s", "4s", "2p", "3p", "4p", "3d"]
SUBSHELL_LETTER = "spdf"


# ------------------------------------------------------------ polynomials
def _padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0.0) + (b[i] if i < len(b) else 0.0) for i in range(n)]


def _pscale(a, s):
    return [c * s for c in a]


def _pshift(a):  # multiply by x
    return [0.0] + list(a)


def laguerre(k, alpha):
    """Coefficients (ascending) of the generalised Laguerre L_k^alpha(x)."""
    prev, cur = [1.0], [1.0 + alpha, -1.0]
    if k == 0:
        return prev
    for j in range(1, k):
        # L_{j+1} = ((2j+1+alpha-x) L_j - (j+alpha) L_{j-1}) / (j+1)
        nxt = _padd(_padd(_pscale(cur, 2 * j + 1 + alpha), _pscale(_pshift(cur), -1.0)),
                    _pscale(prev, -(j + alpha)))
        prev, cur = cur, _pscale(nxt, 1.0 / (j + 1))
    return cur


def peval(c, x):
    acc = 0.0
    for coef in reversed(c):
        acc = acc * x + coef
    return acc


class Radial:
    """Reference radial function (Z = n) of orbital (n, l), in u = r.

    R_ref(u) = N exp(-u) (2u)^l L_{n-l-1}^{2l+1}(2u)
    We keep q(u) = R_ref(u) / u^l = N exp(-u) 2^l L(2u) to avoid 0/0.
    """

    def __init__(self, n, l):
        self.n, self.l = n, l
        self.norm = math.sqrt(8.0 * math.factorial(n - l - 1) / (2 * n * math.factorial(n + l)))
        lag = laguerre(n - l - 1, 2 * l + 1)
        # coefficients of 2^l L(2u) in u
        self.poly = [c * (2.0 ** k) * (2.0 ** l) for k, c in enumerate(lag)]

    def q(self, u):
        return self.norm * math.exp(-u) * peval(self.poly, u)

    def R(self, u):
        return self.q(u) * u ** self.l

    def R2_over_r2l(self, u):
        """R_ref(u)^2 / u^(2l), finite at u = 0."""
        v = self.q(u)
        return v * v


RADIALS = [Radial(n, l) for n, l in RADIAL]


def row_of(n, l):
    return RADIAL.index((n, l))


# ------------------------------------------------------------ angular shapes
S15 = math.sqrt(15.0 / (4 * math.pi))
D_SHAPES = {
    "xy":    [[0, S15 / 2, 0], [S15 / 2, 0, 0], [0, 0, 0]],
    "yz":    [[0, 0, 0], [0, 0, S15 / 2], [0, S15 / 2, 0]],
    "xz":    [[0, 0, S15 / 2], [0, 0, 0], [S15 / 2, 0, 0]],
    "x2-y2": [[S15 / 2, 0, 0], [0, -S15 / 2, 0], [0, 0, 0]],
    "z2":    [[-math.sqrt(5 / (16 * math.pi)), 0, 0],
              [0, -math.sqrt(5 / (16 * math.pi)), 0],
              [0, 0, 2 * math.sqrt(5 / (16 * math.pi))]],
}
P_AXES = {"x": 0, "y": 1, "z": 2}
P_NORM = math.sqrt(3.0 / (4 * math.pi))


class Shape:
    """One real orbital: radial row plus angular part."""

    def __init__(self, n, l, name):
        self.n, self.l, self.name = n, l, name
        self.row = row_of(n, l)
        self.radial = RADIALS[self.row]
        self.e = None
        self.Q = None
        if l == 1:
            self.e = [0.0, 0.0, 0.0]
            self.e[P_AXES[name]] = P_NORM
        elif l == 2:
            self.Q = D_SHAPES[name]

    @property
    def label(self):
        return f"{self.n}{SUBSHELL_LETTER[self.l]}" + (self.name if self.l else "")

    def angular2(self, x):
        """Y^2 at direction x (need not be normalised)."""
        r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
        if self.l == 0:
            return 1.0 / (4 * math.pi)
        if self.l == 1:
            v = dot(self.e, x)
            return v * v / r2
        v = qform(self.Q, x, x)
        return v * v / (r2 * r2)

    def angular2_max(self):
        if self.l == 0:
            return 1.0 / (4 * math.pi)
        if self.l == 1:
            return P_NORM ** 2
        return max(abs(ev) for ev in sym_eigvals(self.Q)) ** 2

    def density_ref(self, x):
        """|psi|^2 in reference units (Z = n) at point x."""
        r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
        u = math.sqrt(r2)
        q2 = self.radial.R2_over_r2l(u)
        if self.l == 0:
            return q2 / (4 * math.pi)
        if self.l == 1:
            v = dot(self.e, x)
            return q2 * v * v          # (R^2/r^2) (e.x)^2
        v = qform(self.Q, x, x)
        return q2 * v * v              # (R^2/r^4) (xQx)^2

    def density(self, z, x):
        k = z / self.n
        return k ** 3 * self.density_ref([k * c for c in x])

    def column_ref(self, moments, p, d):
        """Column density in reference units given moments(k, b) -> I_k(b)."""
        b = math.sqrt(dot(p, p))
        if self.l == 0:
            return moments(0, b) / (4 * math.pi)
        if self.l == 1:
            ep, ed = dot(self.e, p), dot(self.e, d)
            return ep * ep * moments(0, b) + ed * ed * moments(2, b)
        a = qform(self.Q, p, p)
        beta = qform(self.Q, p, d)
        c = qform(self.Q, d, d)
        return (a * a * moments(0, b) + (4 * beta * beta + 2 * a * c) * moments(2, b)
                + c * c * moments(4, b))

    def column(self, moments, z, p, d):
        k = z / self.n
        return k * k * self.column_ref(moments, [k * c for c in p], d)


def all_shapes():
    out = []
    for n, l in RADIAL:
        if l == 0:
            out.append(Shape(n, l, ""))
        elif l == 1:
            out += [Shape(n, l, a) for a in ("x", "y", "z")]
        else:
            out += [Shape(n, l, s) for s in D_SHAPES]
    return out


# ------------------------------------------------------------ small linear algebra
def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def qform(Q, a, b):
    return sum(a[i] * Q[i][j] * b[j] for i in range(3) for j in range(3))


def sym_eigvals(Q):
    """Eigenvalues of a symmetric 3x3 matrix (closed form)."""
    p1 = Q[0][1] ** 2 + Q[0][2] ** 2 + Q[1][2] ** 2
    q = (Q[0][0] + Q[1][1] + Q[2][2]) / 3
    if p1 == 0:
        return [Q[0][0], Q[1][1], Q[2][2]]
    p2 = sum((Q[i][i] - q) ** 2 for i in range(3)) + 2 * p1
    p = math.sqrt(p2 / 6)
    B = [[(Q[i][j] - (q if i == j else 0)) / p for j in range(3)] for i in range(3)]
    detB = (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
            - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
            + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]))
    phi = math.acos(max(-1.0, min(1.0, detB / 2))) / 3
    e1 = q + 2 * p * math.cos(phi)
    e3 = q + 2 * p * math.cos(phi + 2 * math.pi / 3)
    return [e1, 3 * q - e1 - e3, e3]


def mean_r2(n, l, z):
    return n * n / (2 * z * z) * (5 * n * n + 1 - 3 * l * (l + 1))
