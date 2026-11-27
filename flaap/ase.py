"""Aeroservoelasticity (ASE): actuators, sensors, controllers and closed-loop analysis.

The aeroelastic plant is the RFA state space of
:class:`flaap.timedomain.AeroelasticSystem` with control-surface modes.
Its inputs are ``[w_g, w_g', delta, delta', delta'']`` per control
surface, and its outputs are ``eta``, ``eta'``, ``eta''`` (labels
``eta:x``, ``etad:x``, ``etadd:x``) plus any registered output groups,
such as loads from :func:`flaap.loads.station_output_group`.  This module
closes the loop around it::

           w_g, w_g'                                     loads ...
              |                                             ^
   u_cmd  +----------+  delta, delta',  +-----------+  eta, eta', eta''  +---------+  y
 -------->| actuator |----------------->| aeroelastic|------------------>| sensors |------+
   ^      +----------+  delta''         +-----------+                    +---------+      |
   |                                                                                      |
   +-------------------------------------- controller K(s) <------------------------------+

Signals and blocks
------------------
Every block is a :class:`flaap.timedomain.StateSpace` with named inputs
and outputs.  :func:`connect` joins any set of blocks by signal name, in
the manner of MATLAB's ``connect``: an input is fed by the output of the
same name.  Algebraic loops through feed-through terms are solved
exactly: with block outputs ``y = C x + D u`` and interconnection
``u = F y + E w``, ``y = (I - D F)^{-1} (C x + D E w)``.  Such loops
arise from acceleration sensors, which see ``delta''``, combined with a
proper controller.

Actuator
--------
:class:`Actuator` is a strictly proper transfer function
``delta / u_cmd = N(s) / D(s)``, realised in controllable companion form
``x' = A x + B u``, ``delta = C x``.  Its three outputs are consistent
derivatives of the actuator state::

    delta   = C x
    delta'  = C A x + C B u           (C B = 0 for relative degree >= 2)
    delta'' = C A^2 x + C A B u       (exact for relative degree >= 2)

Predefined forms (unit DC gain):

* first order ``1 / (tau s + 1)``. Relative degree 1, so ``delta''``
  would contain ``u'``. That impulsive term is dropped, and a
  second-order or higher actuator is recommended when ``delta''``
  matters;
* second order ``w^2 / (s^2 + 2 zeta w s + w^2)``;
* third order, the second-order form times ``1 / (tau s + 1)``.

Sensors
-------
:class:`Sensor` gives the vertical displacement, rate or acceleration
of a structural point: ``y = phi . eta (+ phi_c . delta)``, and the
same with ``eta'``/``eta''`` and ``delta'``/``delta''``.
:func:`section_sensor` covers the typical section and
:func:`beam_sensor` a point on a beam-modal wing.  A sensor can also
read a single generalised coordinate (``phi = e_j``).

Controllers
-----------
:func:`gain`, :func:`tf`, :func:`pid` (with a derivative filter),
:func:`lead_lag`, :func:`lead_compensator`, :func:`lowpass`,
:func:`washout`, :func:`notch`, :func:`sumblk`, :func:`series`,
:func:`parallel`, :func:`append` and :func:`feedback`.

Sign convention: the controller output is the actuator command
``u_cmd = K(s) y``, so the sign of the feedback is part of ``K``.  The
loop transfer function, broken at the actuator command, is
``L(s) = -K(s) G(s)``, where ``G`` is the transfer from ``u_cmd`` to
``y``.  The closed loop is then the standard negative-feedback loop
``1 + L``.

Analyses (:class:`ASEModel`)
----------------------------
* ``flutter(velocities, controller)``: eigenvalues of the closed-loop
  state matrix versus ``V`` (p-method tracking of the structural
  roots).  ``controller=None`` gives the open loop, with the actuator
  holding ``delta = 0``.
* ``loop(V, controller)`` / ``margins(V, controller)``: ``L(i w)``, gain
  and phase margins, and the Nyquist encirclement count ``N`` compared
  with the number ``P`` of unstable open-loop poles (stable iff
  ``N_ccw = P``).
* ``gust_response(V, H, Uds, controller)``: 1-cos gust time history.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla

from .timedomain import StateSpace

__all__ = [
    "ss", "gain", "tf", "pid", "lead_lag", "lead_compensator", "lowpass", "washout", "notch",
    "sumblk", "connect", "series", "parallel", "append", "feedback", "rename",
    "Actuator", "Sensor", "section_sensor", "beam_sensor", "ASEModel", "Margins",
    "loop_margins", "nyquist_encirclements", "bode_data",
]


# ======================================================================
# elementary LTI blocks
# ======================================================================
def _names(x):
    return [x] if isinstance(x, str) else list(x)


def ss(A, B, C, D, inputs, outputs, states=None, prefix: str = "x") -> StateSpace:
    """:class:`StateSpace` with named inputs/outputs (arrays promoted to 2-D)."""
    inputs, outputs = _names(inputs), _names(outputs)
    A = np.atleast_2d(np.asarray(A, float)) if np.size(A) else np.zeros((0, 0))
    n = A.shape[0]
    B = np.asarray(B, float).reshape(n, len(inputs))
    C = np.asarray(C, float).reshape(len(outputs), n)
    D = np.asarray(D, float).reshape(len(outputs), len(inputs))
    states = list(states) if states is not None else [f"{prefix}{i + 1}" for i in range(n)]
    return StateSpace(A, B, C, D, states, inputs, outputs)


def gain(K, inputs, outputs) -> StateSpace:
    """Static gain ``y = K u`` (``K`` scalar or (nout, nin))."""
    inputs, outputs = _names(inputs), _names(outputs)
    K = np.broadcast_to(np.asarray(K, float), (len(outputs), len(inputs))) if np.ndim(K) == 0 \
        else np.asarray(K, float).reshape(len(outputs), len(inputs))
    return ss(np.zeros((0, 0)), np.zeros((0, len(inputs))), np.zeros((len(outputs), 0)), K, inputs, outputs)


def tf(num, den, input: str = "u", output: str = "y", prefix: str | None = None) -> StateSpace:
    """SISO transfer function ``num(s)/den(s)`` (coefficients in descending
    powers; proper) in controllable canonical form."""
    num = np.atleast_1d(np.asarray(num, float))
    den = np.atleast_1d(np.asarray(den, float))
    den = np.trim_zeros(den, "f")
    num = np.trim_zeros(num, "f") if np.any(num) else np.zeros(1)
    if num.size > den.size:
        raise ValueError("improper transfer function")
    a0 = den[0]
    den = den / a0
    num = np.r_[np.zeros(den.size - num.size), num / a0]
    n = den.size - 1
    D = num[0]
    r = num[1:] - D * den[1:]           # strictly proper remainder (descending, length n)
    if n == 0:
        return gain(D, input, output)
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0
    C = r[None, :]
    pre = prefix or f"{output}.x"
    return ss(A, B, C, [[D]], input, output, [f"{pre}{i + 1}" for i in range(n)])


def tf_eval(num, den, s):
    """Evaluate ``num(s)/den(s)``."""
    return np.polyval(num, s) / np.polyval(den, s)


def pid(kp: float, ki: float = 0.0, kd: float = 0.0, tau_d: float | None = None,
        input: str = "e", output: str = "u") -> StateSpace:
    """``kp + ki/s + kd s/(tau_d s + 1)``. The derivative filter
    ``tau_d`` defaults to ``kd / kp / 10``, or 1e-3 s."""
    terms_num = [np.array([kp])]
    terms_den = [np.array([1.0])]
    if ki:
        terms_num.append(np.array([ki]))
        terms_den.append(np.array([1.0, 0.0]))
    if kd:
        td = tau_d if tau_d is not None else (abs(kd / kp) / 10 if kp else 1e-3)
        terms_num.append(np.array([kd, 0.0]))
        terms_den.append(np.array([td, 1.0]))
    num, den = np.array([0.0]), np.array([1.0])
    for n_, d_ in zip(terms_num, terms_den):
        num = np.polyadd(np.polymul(num, d_), np.polymul(n_, den))
        den = np.polymul(den, d_)
    return tf(num, den, input, output)


def lead_lag(k: float, w_zero: float, w_pole: float, input: str = "u", output: str = "y") -> StateSpace:
    """``k (s/w_zero + 1) / (s/w_pole + 1)``. This is a lead for
    ``w_pole > w_zero`` and a lag otherwise."""
    return tf([k / w_zero, k], [1.0 / w_pole, 1.0], input, output)


def lead_compensator(phase_deg: float, w_center: float, k: float = 1.0, input: str = "u",
                     output: str = "y") -> StateSpace:
    """Lead (``phase_deg > 0``) or lag (``< 0``) network with its maximum
    phase ``phase_deg`` at ``w_center`` and gain ``k`` at that frequency."""
    s = np.sin(np.radians(abs(phase_deg)))
    a = (1 + s) / (1 - s)
    wz, wp = (w_center / np.sqrt(a), w_center * np.sqrt(a))
    if phase_deg < 0:
        wz, wp = wp, wz
    g = abs(tf_eval([1 / wz, 1], [1 / wp, 1], 1j * w_center))
    return lead_lag(k / g, wz, wp, input, output)


def lowpass(w: float, order: int = 2, zeta: float = 0.7071, input: str = "u", output: str = "y") -> StateSpace:
    """Unit-DC-gain low-pass filter: first order ``w/(s+w)`` or second order
    ``w^2/(s^2 + 2 zeta w s + w^2)``."""
    if order == 1:
        return tf([w], [1.0, w], input, output)
    if order == 2:
        return tf([w * w], [1.0, 2 * zeta * w, w * w], input, output)
    raise ValueError("order must be 1 or 2")


def washout(w: float, input: str = "u", output: str = "y") -> StateSpace:
    """High-pass (washout) filter ``s / (s + w)``."""
    return tf([1.0, 0.0], [1.0, w], input, output)


def notch(w0: float, zeta_num: float = 0.02, zeta_den: float = 0.5, input: str = "u",
          output: str = "y") -> StateSpace:
    """Notch ``(s^2 + 2 zn w0 s + w0^2) / (s^2 + 2 zd w0 s + w0^2)``."""
    return tf([1.0, 2 * zeta_num * w0, w0 * w0], [1.0, 2 * zeta_den * w0, w0 * w0], input, output)


def sumblk(output: str, inputs, signs=None) -> StateSpace:
    """Summing junction ``output = sum(sign_i * input_i)``."""
    inputs = _names(inputs)
    signs = np.ones(len(inputs)) if signs is None else np.asarray(signs, float)
    return gain(signs[None, :], inputs, output)


def rename(sys: StateSpace, inputs: dict | None = None, outputs: dict | None = None,
           prefix: str | None = None) -> StateSpace:
    """Copy of ``sys`` with renamed signals (and optionally prefixed states)."""
    ins = [(inputs or {}).get(n, n) for n in sys.inputs]
    outs = [(outputs or {}).get(n, n) for n in sys.outputs]
    st = [f"{prefix}{s}" for s in sys.states] if prefix else list(sys.states)
    return StateSpace(sys.A.copy(), sys.B.copy(), sys.C.copy(), sys.D.copy(), st, ins, outs)


def connect(blocks, inputs, outputs, check: bool = True) -> StateSpace:
    """Interconnect named blocks.

    Every block input whose name matches an output of some block is fed
    by that signal.  Any remaining input must be one of the external
    ``inputs`` (otherwise it is set to zero, and ``check`` raises).
    ``outputs`` may list any internal signal or external input.
    Algebraic loops are solved exactly (an error is raised if
    ``I - D F`` is singular, i.e. the loop is ill-posed).
    """
    blocks = list(blocks)
    inputs, outputs = _names(inputs), _names(outputs)
    A = sla.block_diag(*[b.A for b in blocks]) if blocks else np.zeros((0, 0))
    A = np.atleast_2d(A) if A.size else np.zeros((sum(b.A.shape[0] for b in blocks),) * 2)
    n = A.shape[0]
    nu = sum(len(b.inputs) for b in blocks)
    ny = sum(len(b.outputs) for b in blocks)
    B = np.zeros((n, nu))
    C = np.zeros((ny, n))
    D = np.zeros((ny, nu))
    states, uin, yout = [], [], []
    i0 = j0 = k0 = 0
    for b in blocks:
        nb, mb, pb = b.A.shape[0], len(b.inputs), len(b.outputs)
        B[i0:i0 + nb, j0:j0 + mb] = b.B
        C[k0:k0 + pb, i0:i0 + nb] = b.C
        D[k0:k0 + pb, j0:j0 + mb] = b.D
        states += list(b.states)
        uin += list(b.inputs)
        yout += list(b.outputs)
        i0, j0, k0 = i0 + nb, j0 + mb, k0 + pb
    if len(set(yout)) != len(yout):
        dup = sorted({y for y in yout if yout.count(y) > 1})
        raise ValueError(f"duplicate output signals: {dup}")
    yidx = {y: i for i, y in enumerate(yout)}
    F = np.zeros((nu, ny))
    E = np.zeros((nu, len(inputs)))
    missing = []
    for j, name in enumerate(uin):
        if name in yidx:
            F[j, yidx[name]] = 1.0
        elif name in inputs:
            E[j, inputs.index(name)] = 1.0
        else:
            missing.append(name)
    if check and missing:
        raise ValueError(f"unconnected inputs: {sorted(set(missing))}")
    IDF = np.eye(ny) - D @ F
    Mi = np.linalg.solve(IDF, np.hstack([C, D @ E]))
    Cy, Dy = Mi[:, :n], Mi[:, n:]
    Acl = A + B @ F @ Cy
    Bcl = B @ (F @ Dy + E)
    rows_C, rows_D = [], []
    for o in outputs:
        if o in yidx:
            rows_C.append(Cy[yidx[o]])
            rows_D.append(Dy[yidx[o]])
        elif o in inputs:
            rows_C.append(np.zeros(n))
            e = np.zeros(len(inputs))
            e[inputs.index(o)] = 1.0
            rows_D.append(e)
        else:
            raise ValueError(f"unknown output signal '{o}'")
    Cc = np.array(rows_C).reshape(len(outputs), n)
    Dc = np.array(rows_D).reshape(len(outputs), len(inputs))
    return StateSpace(Acl, Bcl, Cc, Dc, states, inputs, outputs)


def series(*blocks) -> StateSpace:
    """Cascade ``blocks[0] -> blocks[1] -> ...``.  Signals are matched by
    name where possible; otherwise outputs of block ``i`` feed the
    inputs of block ``i+1`` in order (the inputs are renamed)."""
    blocks = list(blocks)
    out = [blocks[0]]
    for b in blocks[1:]:
        prev = out[-1]
        if not set(b.inputs) & set(prev.outputs):
            if len(b.inputs) != len(prev.outputs):
                raise ValueError("cannot match signals by order")
            b = rename(b, inputs=dict(zip(b.inputs, prev.outputs)))
        out.append(b)
    ext = [u for u in blocks[0].inputs]
    return connect(out, ext, out[-1].outputs)


def parallel(*blocks, output: str | None = None) -> StateSpace:
    """Sum of blocks with identical inputs and single (or same) outputs."""
    out_name = output or blocks[0].outputs[0]
    rb = [rename(b, outputs={b.outputs[0]: f"__par{i}"}, prefix=f"p{i}.") for i, b in enumerate(blocks)]
    s = sumblk(out_name, [f"__par{i}" for i in range(len(blocks))])
    return connect(rb + [s], list(blocks[0].inputs), [out_name])


def append(*blocks) -> StateSpace:
    """Block-diagonal stacking (no connections)."""
    ins = sum((list(b.inputs) for b in blocks), [])
    outs = sum((list(b.outputs) for b in blocks), [])
    return connect(blocks, ins, outs, check=False)


def feedback(plant: StateSpace, controller: StateSpace, inputs=None, outputs=None) -> StateSpace:
    """Close ``controller`` around ``plant`` by signal names: controller
    inputs are plant outputs and controller outputs are plant inputs.
    The sign is in the controller (positive interconnection)."""
    ext = [u for u in plant.inputs if u not in controller.outputs] if inputs is None else _names(inputs)
    outs = list(plant.outputs) if outputs is None else _names(outputs)
    return connect([plant, controller], ext, outs)


# ======================================================================
# actuator and sensors
# ======================================================================
@dataclass
class Actuator:
    """Control-surface actuator ``delta / u_cmd = num(s) / den(s)`` (strictly proper).

    ``name`` is the control-surface name of the aeroelastic plant (its
    inputs are ``name``, ``name_dot`` and ``name_ddot``).  The command
    input is ``name_cmd``.
    """

    name: str
    num: np.ndarray
    den: np.ndarray

    @classmethod
    def first_order(cls, name: str, tau: float) -> "Actuator":
        return cls(name, np.array([1.0]), np.array([tau, 1.0]))

    @classmethod
    def second_order(cls, name: str, wn: float, zeta: float = 0.7) -> "Actuator":
        return cls(name, np.array([wn * wn]), np.array([1.0, 2 * zeta * wn, wn * wn]))

    @classmethod
    def third_order(cls, name: str, wn: float, zeta: float = 0.7, tau: float = 0.0) -> "Actuator":
        den = np.polymul([1.0, 2 * zeta * wn, wn * wn], [tau, 1.0])
        return cls(name, np.array([wn * wn]), den)

    @property
    def order(self) -> int:
        return int(np.trim_zeros(np.asarray(self.den, float), "f").size - 1)

    @property
    def relative_degree(self) -> int:
        return self.order - int(np.trim_zeros(np.asarray(self.num, float), "f").size - 1)

    @property
    def cmd(self) -> str:
        return f"{self.name}_cmd"

    def response(self, w):
        return tf_eval(self.num, self.den, 1j * np.asarray(w, float))

    def state_space(self) -> StateSpace:
        """``u_cmd -> [delta, delta', delta'']`` (see module docstring)."""
        if self.relative_degree < 1:
            raise ValueError("actuator transfer function must be strictly proper")
        base = tf(self.num, self.den, self.cmd, self.name, prefix=f"act:{self.name}.")
        A, B, C = base.A, base.B, base.C
        Cf = np.vstack([C, C @ A, C @ A @ A])
        Df = np.vstack([np.zeros((1, 1)), C @ B, C @ A @ B])
        if self.relative_degree >= 2:
            Df[1] = 0.0
        return StateSpace(A, B, Cf, Df, list(base.states), [self.cmd],
                          [self.name, f"{self.name}_dot", f"{self.name}_ddot"])


@dataclass
class Sensor:
    """Point (or generalised-coordinate) sensor.

    ``y = phi . eta + phi_c . delta`` for ``kind='disp'``, and the first
    or second derivatives for ``'rate'`` and ``'acc'``.  ``phi`` has
    length ``ns`` (structural DOFs) and ``phi_c`` length ``n_ctrl``.
    ``gain`` scales the output (sensor sensitivity).
    """

    name: str
    phi: np.ndarray
    kind: str = "acc"
    phi_c: np.ndarray | None = None
    gain: float = 1.0

    def block(self, struct_names, ctrl_names) -> StateSpace:
        tag = {"disp": ("eta", ""), "rate": ("etad", "_dot"), "acc": ("etadd", "_ddot")}[self.kind]
        ins = [f"{tag[0]}:{s}" for s in struct_names]
        K = [np.asarray(self.phi, float) * self.gain]
        if ctrl_names and self.phi_c is not None:
            ins += [f"{c}{tag[1]}" for c in ctrl_names]
            K.append(np.asarray(self.phi_c, float) * self.gain)
        return gain(np.concatenate(K)[None, :], ins, self.name)


def section_sensor(sec, x: float, kind: str = "acc", name: str | None = None,
                   flap_dof: bool | str = "ctrl") -> Sensor:
    """Vertical motion (up) of chord point ``x`` (m aft of mid-chord) of a
    :class:`~flaap.models.typical_section.TypicalSection`::

        z = h - (x - a b) alpha - (x - c b) beta [x > c b]

    ``flap_dof='ctrl'``: the flap is a control input (2 structural DOFs);
    ``'struct'``: the flap is the third structural DOF.
    """
    b = sec.b
    xa, xc = sec.a * b, sec.c * b
    pf = -(x - xc) if x > xc else 0.0
    name = name or f"{kind}@{x:+.3f}"
    if flap_dof == "struct":
        return Sensor(name, np.array([1.0, -(x - xa), pf]), kind)
    return Sensor(name, np.array([1.0, -(x - xa)]), kind, np.array([pf]))


def beam_sensor(modal, x: float, y: float, x_ea: float, kind: str = "acc", name: str | None = None,
                flap: tuple | None = None, nstruct: int | None = None) -> Sensor:
    """Vertical motion of point ``(x, y)`` of a beam-modal wing,
    ``z = w(y) - (x - x_ea) ry(y)`` (``ry`` nose up), interpolated
    linearly between the beam nodes.  ``flap = (x_h, y0, y1)`` adds
    the flap rotation ``-(x - x_h)`` when the point lies on the flap.
    """
    ny = modal.nodes[:, 1]
    o = np.argsort(ny)
    phi = np.array([np.interp(y, ny[o], modal.shapes[o, 2, i]) - (x - x_ea) * np.interp(y, ny[o], modal.shapes[o, 4, i])
                    for i in range(modal.nmodes)])
    if nstruct is not None:
        phi = np.r_[phi, np.zeros(nstruct - phi.size)]
    pc = None
    if flap is not None:
        x_h, y0, y1 = flap
        pc = np.array([-(x - x_h) if (x > x_h and y0 <= y <= y1) else 0.0])
    return Sensor(name or f"{kind}@({x:.2f},{y:.2f})", phi, kind, pc)


# ======================================================================
# loop analysis helpers
# ======================================================================
@dataclass
class Margins:
    """Stability margins of a SISO loop ``L(i w)`` (negative feedback ``1 + L``)."""

    gain_margins: list            # (factor, dB, w) at phase crossovers (L real negative)
    phase_margins: list           # (deg, w) at gain crossovers (|L| = 1)
    gm_upper: float               # smallest factor > 1 that destabilises (inf if none)
    gm_lower: float               # largest factor < 1 that destabilises (0 if none)
    pm: float                     # smallest |phase margin| [deg] (inf if no crossover)
    w_pm: float
    w_gm: float                   # frequency of the upper gain margin
    P: int                        # unstable open-loop poles
    N_ccw: int                    # counter-clockwise encirclements of -1
    nyquist_stable: bool
    closed_loop_stable: bool | None = None

    @property
    def gm_upper_db(self):
        return 20 * np.log10(self.gm_upper) if np.isfinite(self.gm_upper) else np.inf

    @property
    def gm_lower_db(self):
        return 20 * np.log10(self.gm_lower) if self.gm_lower > 0 else -np.inf

    def summary(self) -> str:
        s = [f"P = {self.P} unstable open-loop poles, N_ccw = {self.N_ccw} -> "
             f"{'stable' if self.nyquist_stable else 'UNSTABLE'} (Nyquist)"]
        if self.closed_loop_stable is not None:
            s[-1] += f", eigenvalues: {'stable' if self.closed_loop_stable else 'UNSTABLE'}"
        s.append(f"gain margins: upper {self.gm_upper_db:+.2f} dB @ {self.w_gm:.2f} rad/s, "
                 f"lower {self.gm_lower_db:+.2f} dB")
        s.append(f"phase margin: {self.pm:.2f} deg @ {self.w_pm:.2f} rad/s")
        return "\n".join(s)


def frequency_grid(A, wmin: float, wmax: float, n: int = 3000, dense: int = 400) -> np.ndarray:
    """Log grid refined around the lightly damped poles of ``A``."""
    w = list(np.geomspace(wmin, wmax, n))
    for lam in np.linalg.eigvals(A):
        if lam.imag > wmin and lam.imag < wmax:
            bw = max(abs(lam.real), 1e-4 * lam.imag)
            w += list(lam.imag + bw * np.sinh(np.linspace(-6, 6, dense)) / 10)
    w = np.array(w)
    return np.unique(w[(w >= wmin) & (w <= wmax)])


def nyquist_encirclements(L) -> int:
    """Counter-clockwise encirclements of -1 by ``L(i w)`` over
    ``w in (-inf, inf)``, given ``L`` on ``w >= 0`` (ascending, ending where
    ``|L| << 1``).  Uses ``Delta arg(1 + L)`` and conjugate symmetry."""
    f = 1.0 + np.asarray(L)
    th = np.unwrap(np.angle(f))
    return int(np.round(2 * (th[-1] - th[0]) / (2 * np.pi)))


def loop_margins(L, w, P: int = 0) -> Margins:
    """Gain/phase margins and the Nyquist count from samples ``L(i w)``."""
    L = np.asarray(L)
    w = np.asarray(w, float)
    gms = []
    im = L.imag
    for j in np.nonzero(np.sign(im[:-1]) * np.sign(im[1:]) < 0)[0]:
        f = im[j] / (im[j] - im[j + 1])
        Lr = L[j].real + f * (L[j + 1].real - L[j].real)
        wj = w[j] + f * (w[j + 1] - w[j])
        if Lr < 0:
            g = -1.0 / Lr
            gms.append((g, 20 * np.log10(g), wj))
    pms = []
    mag = np.abs(L)
    for j in np.nonzero(np.sign(mag[:-1] - 1) * np.sign(mag[1:] - 1) < 0)[0]:
        f = (1 - mag[j]) / (mag[j + 1] - mag[j])
        Lc = L[j] + f * (L[j + 1] - L[j])
        wj = w[j] + f * (w[j + 1] - w[j])
        pm = (np.degrees(np.angle(Lc)) + 180.0 + 180.0) % 360.0 - 180.0
        pms.append((pm, wj))
    up = [g for g in gms if g[0] > 1]
    lo = [g for g in gms if g[0] < 1]
    gmu = min(up, key=lambda g: g[0]) if up else (np.inf, np.inf, np.nan)
    gml = max(lo, key=lambda g: g[0]) if lo else (0.0, -np.inf, np.nan)
    pm = min(pms, key=lambda p: abs(p[0])) if pms else (np.inf, np.nan)
    N = nyquist_encirclements(L)
    return Margins(gms, pms, gmu[0], gml[0], abs(pm[0]), pm[1], gmu[2], P, N, N == P)


def bode_data(L, w) -> dict:
    """Magnitude (dB) and unwrapped phase (deg) of ``L(i w)``."""
    L = np.asarray(L)
    return {"w": np.asarray(w), "mag_db": 20 * np.log10(np.abs(L) + 1e-300),
            "phase_deg": np.degrees(np.unwrap(np.angle(L)))}


# ======================================================================
# ASE model
# ======================================================================
class ASEModel:
    """Closed-loop aeroservoelastic model.

    Parameters
    ----------
    system : AeroelasticSystem
        Plant with ``n_ctrl`` control surfaces.
    actuators : list of Actuator
        One per control surface (``Actuator.name`` must match the
        control-mode name).
    sensors : list of Sensor
    rho : float
        Air density.
    """

    def __init__(self, system, actuators, sensors, rho: float):
        self.system = system
        self.actuators = list(actuators)
        self.sensors = list(sensors)
        self.rho = float(rho)
        ns = system.ns
        self.struct_names = list(system.names[:ns])
        self.ctrl_names = list(system.names[ns:])
        for a in self.actuators:
            if a.name not in self.ctrl_names:
                raise ValueError(f"actuator '{a.name}' is not a control mode {self.ctrl_names}")

    @property
    def commands(self):
        return [a.cmd for a in self.actuators]

    @property
    def sensor_names(self):
        return [s.name for s in self.sensors]

    # ------------------------------------------------------------------
    def plant(self, V: float, gust: bool = True, extra_outputs=None) -> StateSpace:
        """Open-loop plant (actuators + aeroelastic system + sensors).

        Inputs: ``[w_g, w_g'] + commands``.  Outputs: sensor signals,
        control deflections ``delta, delta', delta''``, then
        ``extra_outputs`` (default: all aeroelastic outputs).
        """
        sysm = self.system
        P = sysm.state_space(V, self.rho, gust=gust, outputs=True)
        gust = "wg" in P.inputs
        blocks = [P] + [a.state_space() for a in self.actuators]
        blocks += [s.block(self.struct_names, self.ctrl_names) for s in self.sensors]
        ext = (["wg", "wg_dot"] if gust else []) + self.commands
        defl = sum(([a.name, f"{a.name}_dot", f"{a.name}_ddot"] for a in self.actuators), [])
        extra = list(P.outputs) if extra_outputs is None else _names(extra_outputs)
        return connect(blocks, ext, self.sensor_names + defl + extra)

    def closed_loop(self, V: float, controller: StateSpace | None = None, gain: float = 1.0,
                    gust: bool = True, extra_outputs=None, external_commands: bool = False) -> StateSpace:
        """Closed loop with ``u_cmd = gain * K(s) y``.

        ``external_commands=True`` adds an external input
        ``<cmd>_ext`` summed into each command (for loop breaking or
        reference commands).  ``controller=None`` or ``gain = 0`` gives
        the open loop.
        """
        Pl = self.plant(V, gust, extra_outputs)
        ext = [u for u in Pl.inputs if u not in self.commands]
        blocks = [Pl]
        if controller is not None and gain != 0.0:
            Kc = rename(controller, outputs={c: f"{c}__k" for c in self.commands}, prefix="K:")
            Kc = StateSpace(Kc.A, Kc.B, gain * Kc.C, gain * Kc.D, Kc.states, Kc.inputs, Kc.outputs)
            blocks.append(Kc)
            for c in self.commands:
                if f"{c}__k" in Kc.outputs:
                    sig = [f"{c}__k"] + ([f"{c}_ext"] if external_commands else [])
                    blocks.append(sumblk(c, sig))
                elif external_commands:
                    blocks.append(sumblk(c, [f"{c}_ext"]))
                else:
                    blocks.append(gain_zero(c))
        else:
            for c in self.commands:
                blocks.append(sumblk(c, [f"{c}_ext"]) if external_commands else gain_zero(c))
        if external_commands:
            ext = ext + [f"{c}_ext" for c in self.commands]
        return connect(blocks, ext, list(Pl.outputs))

    # ------------------------------------------------------------------
    def flutter(self, velocities, controller: StateSpace | None = None, gain: float = 1.0,
                w_initial=None):
        """Closed-loop flutter analysis (p-method on the closed-loop state matrix)."""
        from .flutter import p_method
        sysm = self.system
        if w_initial is None:
            w = np.sqrt(np.maximum(np.real(np.linalg.eigvals(np.linalg.solve(sysm.M, sysm.K))), 0))
            w_initial = np.sort(w[w > 1e-6])
        build = lambda V: self.closed_loop(V, controller, gain, gust=False, extra_outputs=[]).A
        return p_method(build, velocities, self.rho, len(w_initial), sysm.bref,
                        names=self.struct_names, w_initial=w_initial,
                        k_valid=getattr(sysm.rfa, "kmax", np.inf))

    def flutter_speed(self, V_range, controller=None, gain: float = 1.0, n: int = 80) -> float:
        """Lowest speed with an unstable closed-loop eigenvalue (scan plus
        bisection), or ``nan``."""
        Vs = np.linspace(V_range[0], V_range[1], n)
        f = lambda V: np.max(np.linalg.eigvals(self.closed_loop(V, controller, gain, gust=False,
                                                                 extra_outputs=[]).A).real)
        s = np.array([f(V) for V in Vs])
        j = np.nonzero((s[:-1] < 0) & (s[1:] >= 0))[0]
        if not j.size:
            return np.nan
        from scipy.optimize import brentq
        return brentq(f, Vs[j[0]], Vs[j[0] + 1], xtol=1e-6 * Vs[j[0] + 1])

    def loop_system(self, V: float, controller: StateSpace) -> StateSpace:
        """Open loop broken at the (single) actuator command:
        ``<cmd>_ext -> controller output``, so that ``L = -T``."""
        if len(self.commands) != 1:
            raise ValueError("loop analysis needs a single control loop")
        c = self.commands[0]
        Pl = self.plant(V, gust=False, extra_outputs=[])
        Kc = rename(controller, outputs={c: "__loop_out"}, prefix="K:")
        Pl = rename(Pl, inputs={c: f"{c}_ext"})
        return connect([Pl, Kc], [f"{c}_ext"], ["__loop_out"])

    def loop(self, V: float, controller: StateSpace, w=None, gain: float = 1.0):
        """``(w, L(i w), T)``: loop transfer ``L = -gain K G`` on ``w``."""
        T = self.loop_system(V, controller)
        if w is None:
            wm = np.sqrt(np.abs(np.linalg.eigvals(np.linalg.solve(self.system.M, self.system.K))))
            w = frequency_grid(T.A, 0.02 * wm.min(), 200 * wm.max())
        L = -gain * T.freq_response(w)[:, 0, 0]
        return np.asarray(w), L, T

    def margins(self, V: float, controller: StateSpace, w=None, gain: float = 1.0) -> Margins:
        """Gain/phase margins and the Nyquist check at speed ``V``.  ``P``
        counts the unstable poles of the open-loop system (plant,
        actuator and controller)."""
        w, L, T = self.loop(V, controller, w, gain)
        P = int(np.sum(np.linalg.eigvals(T.A).real > 1e-9))
        m = loop_margins(L, w, P)
        cl = self.closed_loop(V, controller, gain, gust=False, extra_outputs=[])
        m.closed_loop_stable = bool(np.all(np.linalg.eigvals(cl.A).real < 1e-9))
        return m

    # ------------------------------------------------------------------
    def gust_response(self, V: float, H: float, Uds: float, controller: StateSpace | None = None,
                      gain: float = 1.0, t_end: float = 2.0, dt: float = 1e-3, outputs=None) -> dict:
        """1-cos gust response (gust reference at the leading edge of the
        RFA gust columns).  Returns a dict with ``t``, ``y`` (nt, nout),
        ``outputs``, ``wg`` and ``peak`` (max |y| per output)."""
        from .gust import one_minus_cosine
        cl = self.closed_loop(V, controller, gain, gust=True, extra_outputs=outputs)
        t = np.arange(0.0, t_end + 0.5 * dt, dt)
        wg, wgd = one_minus_cosine(t, V, H, Uds)
        y, _ = cl.simulate(t, np.c_[wg, wgd])
        return {"t": t, "y": y, "outputs": list(cl.outputs), "wg": wg,
                "peak": dict(zip(cl.outputs, np.max(np.abs(y), axis=0)))}


def gain_zero(output: str) -> StateSpace:
    """A zero source for an unused command signal."""
    return ss(np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((1, 0)), np.zeros((1, 0)), [], [output])
