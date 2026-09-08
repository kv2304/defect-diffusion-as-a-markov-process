"""
Trap-limited defect diffusion on a 2D lattice as an absorbing Markov chain.

[E] = eV
[t] = s
"""

import os
from typing import Any

import numpy as np
from scipy.sparse import csr_array, diags_array, eye_array
from scipy.sparse.linalg import splu

NU = 1.0e13     # attempt frequency, s^-1
E_MIG = 0.5     # migration barrier, eV
E_B = 0.3       # trap binding energy, eV


def chain(L, traps, kT, E_mig=E_MIG, E_b=E_B):
    """
    Assemble the jump chain: Q (interior transitions), absorption weights,
    and trap-dependent waiting times tau on an L x L lattice, row 0 absorbing.
    """
    if traps.shape != (L, L):
        raise ValueError("traps must mask the full L x L lattice.")

    row, col = np.divmod(np.arange(L * L), L)
    z = 4 - (row == 0) - (row == L - 1) - (col == 0) - (col == L - 1)
    row, col, z = row[L:], col[L:], z[L:]

    rn = row + np.array([[-1], [1], [0], [0]])
    cn = col + np.array([[0], [0], [-1], [1]])
    on = (rn >= 0) & (rn < L) & (cn >= 0) & (cn < L)
    inside = on & (rn >= 1)

    site = np.nonzero(inside)[1]
    Q = csr_array((1.0 / z[site], (site, rn[inside] * L + cn[inside] - L)), shape=(z.size, z.size))
    absorbed = ((on & ~inside) / z).sum(axis=0)

    tau = np.exp((E_mig + E_b * traps[1:].ravel()) / kT) / (z * NU)

    return Q, absorbed, tau


def escape(Q, tau, lu: Any = None):
    """
    Solve (I - Q) t = tau for mean escape times, reusing the LU
    factorisation lu across calls.
    """
    if lu is None:
        lu = splu((eye_array(Q.shape[0]) - Q).tocsc())

    return lu.solve(tau), lu


def disorder(L, cs, betas, n_real, seed=0, E_mig=E_MIG, E_b=E_B):
    """
    Average escape time over trap disorder, density, and temperature,
    using linearity in tau to avoid refactorising per configuration.
    """
    rng = np.random.default_rng(seed)
    none, full = np.zeros((L, L), bool), np.ones((L, L), bool)

    Q = chain(L, none, 1.0, E_mig, E_b)[0]
    tau_free = np.array([chain(L, none, 1 / b, E_mig, E_b)[2] for b in betas])
    tau_trap = np.array([chain(L, full, 1 / b, E_mig, E_b)[2] for b in betas])
    ramp = tau_free[:, 0] / tau_free[0, 0]
    rho = tau_trap[:, 0] / tau_free[:, 0]

    G, lu = escape(Q, tau_free.T)
    G = G.mean(axis=0)

    t, sd = np.empty((2, betas.size, cs.size))
    for j, c in enumerate(cs):
        n = round(c * L * L)
        score = rng.random((n_real, L * L))
        traps = score <= np.partition(score, n - 1, axis=1)[:, n - 1, None]
        share = escape(Q, (tau_free[0] * traps[:, L:]).T, lu)[0].mean(axis=0)
        t[:, j] = G + ramp * (rho - 1) * share.mean()
        sd[:, j] = ramp * (rho - 1) * share.std(ddof=1)

    return t, sd, G

def measure(L, cs, betas, n_real, E_mig=E_MIG, E_b=E_B) -> dict:
    """
    Verify the escape-time model against closed-form limits and extract
    the effective activation energy crossover with density.
    """
    L_CONV = 70         # sites, second lattice for the size check
    M_MAX = 20          # sites, 1D chain for the closed-form check
    C_FIG, X_FIG = 0.05, 3.0    # the density and E_b/kT the figures use
    xs = betas[:, betas.shape[1] // 2] * E_b
    ex = np.exp(xs)
    none, full = np.zeros((L, L), bool), np.ones((L, L), bool)

    def slopes(y):
        """Arrhenius slope d ln<t>/d beta in each window, one per density."""
        return np.array([np.polyfit(betas[a], np.log(y[a]), 1)[0] for a in range(xs.size)])

    def crossover(c, f):
        """Density at which f reaches 1/2. logit(f) = logit(c) + E_b/kT exactly,
        so the grid is read in that coordinate rather than in f itself."""
        lg = np.log(c / (1 - c))
        return 1 / (1 + np.exp(-np.interp(0.0, np.log(f / (1 - f)), lg)))

    cs_conv = np.round(cs * L_CONV ** 2) / L_CONV ** 2
    t, sd, G = disorder(L, cs, betas.ravel(), n_real, 0, E_mig, E_b)
    t_conv = disorder(L_CONV, cs_conv, betas.ravel(), n_real, 0, E_mig, E_b)[0]

    # reuse one factorisation below
    j_fig = int(np.abs(cs - C_FIG).argmin())
    a_fig = int(np.abs(xs - X_FIG).argmin())
    c_fig, kT_fig = cs[j_fig], E_b / xs[a_fig]
    Q, absorbed, tau_free = chain(L, none, kT_fig, E_mig, E_b)
    G_fig, lu = escape(Q, tau_free)

    # float64 bound on the solve, eps * cond_1(I - Q), off the same factorisation.
    IQ = eye_array(Q.shape[0]) - Q
    cond = float(np.abs(IQ).sum(axis=0).max() * lu.solve(np.ones(tau_free.size), "T").max())
    floor = np.finfo(float).eps * cond

    # check row sums, Q trap-independence, coordination
    rng = np.random.default_rng(0)
    n_fig = round(c_fig * L * L)
    pair = [np.isin(np.arange(L * L), rng.choice(L * L, n_fig, False))
            .reshape(L, L) for _ in range(3)]
    q_same = all(np.array_equal(chain(L, m, 0.07, E_mig, E_b)[0].data, Q.data) for m in pair)
    z = np.rint(1.0 / Q.data[Q.indptr[:-1]]).astype(int)

    # 1D chain: t, N known in closed form
    e_chain = e_fund = 0.0
    for M in range(3, M_MAX + 1):
        i = np.arange(1, M)
        Qc = diags_array([np.full(M - 2, 0.5)] * 2, offsets=[-1, 1], format="csr")  # type: ignore
        tc, luc = escape(Qc, np.ones(M - 1))
        Nc = luc.solve(np.eye(M - 1))
        Nex = (2 * np.minimum(i[:, None], i) * (M - np.maximum(i[:, None], i)) / M)
        e_chain = max(e_chain, np.abs(tc / (i * (M - i)) - 1).max(),
                      np.abs(Nc.sum(axis=1) / (i * (M - i)) - 1).max())
        e_fund = max(e_fund, np.abs(Nc / Nex - 1).max())

    # single trap: closed form via N_im
    m, T = M_MAX // 3, 40.0
    im = np.arange(1, M_MAX)
    Qm = diags_array([np.full(M_MAX - 2, 0.5)] * 2, offsets=[-1, 1], format="csr")  # type: ignore
    tau_c = np.ones(M_MAX - 1)
    tau_c[m - 1] = T
    t_trap = im * (M_MAX - im) + (T - 1) * 2 * np.minimum(im, m) * (M_MAX - np.maximum(im, m)) / M_MAX
    e_trap = np.abs(escape(Qm, tau_c)[0] / t_trap - 1).max()

    # trap-free lattice: Poisson eq. in row index, quadratic
    tau_0 = tau_free * z
    row = np.arange(z.size) // L + 1
    t_geom = tau_0 / 2 * row * (2 * L - 1 - row)
    e_geom = np.abs(G_fig / t_geom - 1).max()

    # sanity check: absorption probability = 1
    e_sure = np.abs(escape(Q, absorbed, lu)[0] - 1.0).max()

    tau_fig = chain(L, pair[0], kT_fig, E_mig, E_b)[2]
    t_fig = escape(Q, tau_fig, lu)[0]

    # check <t> matches sampling error
    trapped = (1 - cs) + cs * np.exp(betas.ravel()[:, None] * E_b)
    res = np.abs(t - trapped * G[:, None])
    fac_z = np.sqrt(((res * np.sqrt(n_real) / sd) ** 2).mean())

    # check G independent of c, E_b
    G_move = 0.0
    for c in (0.0, 0.001, 0.05, 0.3, 0.9, 1.0):
        for eb in (0.1, E_b, 0.7, 1.2):
            tau_mean = (1 - c) * tau_free + c * chain(L, full, kT_fig, E_mig, eb)[2]
            G_est = (escape(Q, tau_mean, lu)[0].mean() / ((1 - c) + c * np.exp(eb / kT_fig)))
            G_move = max(G_move, abs(G_est / G_fig.mean() - 1))
    tau_anchor = 0.95 * tau_free + 0.05 * chain(L, full, kT_fig, E_mig, E_b)[2]
    t_anchor = escape(Q, tau_anchor, lu)[0][-L:].mean()

    # check f's low/high-density limits
    eps = 1.0e-6
    lim_lo = eps / (ex * (1 - eps) + eps)
    lim_hi = eps * ex / (eps * ex + 1 - eps)
    c_lim = np.concatenate([lim_lo, 1 - lim_hi])
    f_lim = c_lim * np.tile(ex, 2) / ((1 - c_lim) + c_lim * np.tile(ex, 2))
    lim_err = np.abs(np.concatenate([f_lim[:xs.size], 1 - f_lim[xs.size:]]) / eps - 1).max()

    # effective activation energy from Arrhenius slope
    nb = betas.shape[1]
    f_an = cs * ex[:, None] / ((1 - cs) + cs * ex[:, None])
    Ea, Ea_an = slopes(t.reshape(xs.size, nb, -1)), E_mig + E_b * f_an
    f_meas = (Ea - E_mig) / E_b

    # check finite-size stability at matched c
    cstar = np.array([crossover(cs, f) for f in f_meas])
    cstar_an = 1 / (1 + ex)
    f_conv = (slopes(t_conv.reshape(xs.size, nb, -1)) - E_mig) / E_b
    cstar_conv = np.array([crossover(cs_conv, f) for f in f_conv])
    keep = (cs > cs_conv[0]) & (cs < cs_conv[-1])
    lg, lg_conv = np.log(cs / (1 - cs)), np.log(cs_conv / (1 - cs_conv))
    conv_f = max(np.abs(1 / (1 + np.exp(-np.interp(lg, lg_conv,
                 np.log(fc / (1 - fc)))))[keep] - f_meas[a][keep]).max()
                 for a, fc in enumerate(f_conv))

    # check scatter falls as 1/sqrt(n)
    score = np.random.default_rng(1).random((n_real, L * L))
    many = score <= np.partition(score, n_fig - 1, axis=1)[:, n_fig - 1, None]
    share = escape(Q, (tau_free * many[:, L:]).T, lu)[0]
    lift = (np.exp(E_b / kT_fig) - 1.0) / G_fig

    conv_n = np.array([1, 2, 4, 5, 8, 10, 16, 20, 25, 40, 50, 80, 100])
    batch = [np.array([1.0 + lift * share[:, j * n:(j + 1) * n].mean(axis=1)
                       for j in range(n_real // n)]) for n in conv_n]
    conv_sd = np.array([b.std(axis=1).mean() / b.mean() for b in batch])
    conv_slope = float(np.polyfit(np.log(conv_n), np.log(conv_sd), 1)[0])

    # check factorisation across full sweep
    mid = np.arange(xs.size) * nb + nb // 2
    fac_ratio = t[mid] / G[mid][:, None]
    fac_an = (1 - cs) + cs * ex[:, None]
    arr = [int(np.abs(cs - v).argmin()) for v in (0.0016, 0.017, 0.086, 0.44)]

    ends = np.abs(Ea[:, [0, -1]] / Ea_an[:, [0, -1]] - 1).max(axis=0)

    return dict(
        L=L, L_CONV=L_CONV, cs=cs, xs=xs, betas=betas, n_real=n_real,
        q_same=q_same, e_chain=e_chain, e_fund=e_fund, e_trap=e_trap,
        e_geom=e_geom, e_sure=e_sure,
        fac_res=(res / t).mean(), fac_z=fac_z, G_move=G_move, lim_err=lim_err,
        Ea=Ea, Ea_an=Ea_an, ends=ends, f_meas=f_meas, f_an=f_an,
        cstar=cstar, cstar_an=cstar_an, cstar_conv=cstar_conv, conv_f=conv_f,
        tau_bulk=tau_free[(L // 2 - 1) * L + L // 2],
        t_far_free=G_fig[-L:].mean(), t_anchor=t_anchor,
        c_fig=c_fig, x_fig=xs[a_fig], kT_fig=kT_fig, a_fig=a_fig,
        trap_fig=pair[0], t_map=np.r_[np.zeros(L), t_fig].reshape(L, L),
        conv_n=conv_n, conv_sd=conv_sd, conv_slope=conv_slope,
        cond=cond, floor=floor,
        fac_ratio=fac_ratio, fac_an=fac_an,
        fig_trapped=(1 - c_fig) + c_fig * ex[a_fig],
        arr=arr, arr_Ea=Ea[:, arr],
        arr_t=t.reshape(xs.size, nb, -1)[:, :, arr])

def figures(results, outdir="figures"):
    """
    Generate figures using results dict from measure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.transforms import Bbox

    os.makedirs(outdir, exist_ok=True)
    plt.rcParams.update({
        "figure.figsize": (3.5, 2.6), "savefig.dpi": 300,
        "savefig.bbox": "tight", "font.family": "serif",
        "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
        "axes.linewidth": 0.8, "lines.linewidth": 1.4,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "legend.frameon": False, "axes.titlesize": 7.5, "axes.titlepad": 4.0,
    })
    COLOURS = ["#009E73", "#D55E00", "#CC79A7", "#0072B2"]  
    PURPLES = plt.get_cmap("Purples")(np.linspace(0.42, 0.95, results["xs"].size))
    PAIR = (4.9, 3.2)       
    L, cs, xs = results["L"], results["cs"], results["xs"]
    c_fig, x_fig = results["c_fig"], results["x_fig"]
    kT_fig = results["kT_fig"]
    tr = np.nonzero(results["trap_fig"])

    # lattice.svg: random traps, one absorbing edge, three reflecting
    edge = np.arange(-0.5, L, 1.0)
    fig, ax = plt.subplots(figsize=(4.3, 3.4))
    ax.vlines(edge, -0.5, L - 0.5, color="0.92", lw=0.25, zorder=0)
    ax.hlines(edge, -0.5, L - 0.5, color="0.92", lw=0.25, zorder=0)
    ax.axhspan(-0.5, 0.5, color=COLOURS[1], alpha=0.16, zorder=1)
    ax.plot(tr[1], tr[0], "s", ms=1.9, mec="none", color=COLOURS[0], zorder=3)
    ax.plot([-0.5, L - 0.5], [-0.5, -0.5], color=COLOURS[1], lw=2.2, solid_capstyle="butt", zorder=4)
    for xy in ([[-0.5, -0.5], [-0.5, L - 0.5]], [[L - 0.5, L - 0.5], [-0.5, L - 0.5]], [[-0.5, L - 0.5], [L - 0.5, L - 0.5]]):
        ax.plot(xy[0], xy[1], color="0.45", lw=1.8, solid_capstyle="butt", zorder=4)
    ax.annotate("absorbing edge, row 0", (L / 2 - 0.5, -2.6),
                ha="center", fontsize=7, color=COLOURS[1])
    ax.annotate("reflecting", (-3.1, L / 2), rotation=90, va="center",
                ha="center", fontsize=7, color="0.35")
    ax.annotate("reflecting", (L / 2 - 0.5, L + 3.1), ha="center",
                va="center", fontsize=7, color="0.35")
    ax.set_xlim(-5.2, L + 0.5)
    ax.set_ylim(L + 5.0, -5.0)
    ax.set_aspect("equal")
    ax.set_xlabel("column (sites)")
    ax.set_ylabel("row (sites), distance from the edge")
    ax.set_title(rf"$L={L}$: $L^2={L * L}$ sites, $L(L-1)={L * (L - 1)}$" rf" interior" "\n" rf"$c={c_fig:.4f}$, {tr[0].size} traps" r" (green)")

    fig.savefig(f"{outdir}/lattice.svg")

    # escape_map.svg: mean escape time per site, one configuration, traps on top
    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    im = ax.imshow(results["t_map"] * 1e9, cmap="magma", interpolation="nearest")
    ax.plot(tr[1], tr[0], "s", ms=1.7, mec="none", color=COLOURS[0])
    cb = fig.colorbar(im, ax=ax, pad=0.03)
    cb.set_label(r"$t$, mean escape time (ns)")
    ax.set_xlabel("column (sites)")
    ax.set_ylabel("row (sites)")
    ax.set_title(rf"one placement, $c={c_fig:.4f}$, traps in green"
                 "\n" rf"$E_b/kT={x_fig:.0f}$ ($kT={kT_fig:.3f}$ eV)")
    fig.savefig(f"{outdir}/escape_map.svg")

    # factorisation.svg: L the trapping factor against density at each kT
    fr, fa = results["fac_ratio"], results["fac_an"]
    cn, csd = results["conv_n"], results["conv_sd"]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    for k in range(xs.size):
        ax0.loglog(cs, fa[k], lw=0.9, color=PURPLES[k], zorder=1)
        ax0.loglog(cs, fr[k], "o", ms=3.0, mec="none", color=PURPLES[k], zorder=2)
        ax0.annotate(rf"${xs[k]:.0f}$", (cs[-1] * 1.15, fa[k, -1]), fontsize=7,
                     va="center", color="0.35")
    ax0.set_xlabel("trap density $c$")
    ax0.set_ylabel(r"$\langle t \rangle / G$")
    ax0.set_title(rf"points solved, lines $(1-c)+ce^{{E_b/kT}}$;"
                  rf" worst {np.abs(fr / fa - 1).max():.1e}")
    ax0.annotate(r"labelled by $E_b/kT$", (0.03, 0.93), xycoords="axes fraction", fontsize=7.5, color="0.35")

    ax1.loglog(cn, csd * 100, "o", ms=3.5, mec="none", color=PURPLES[results["a_fig"]])
    ax1.loglog(cn, csd[0] * 100 / np.sqrt(cn), "--", lw=0.9, color="0.35",
               label=rf"$n^{{-1/2}}$, fitted ${results['conv_slope']:.3f}$")
    ax1.set_xlabel("placements averaged, $n$")
    ax1.set_ylabel(r"scatter over sites in $t/G$ (%)")
    ax1.set_title("the scatter it leaves is sampling error")
    ax1.legend(loc="lower left", handlelength=1.4, handletextpad=0.5)

    fig.tight_layout()
    fig.savefig(f"{outdir}/factorisation.svg", bbox_inches=Bbox([[0, 0], (7.0, 2.9)]))

    # arrhenius.svg: ln<t> against beta, straight, and the slope is E_a_eff
    betas, arr = results["betas"], results["arr"]
    at, aE = results["arr_t"], results["arr_Ea"]
    blues = plt.get_cmap("Blues")(np.linspace(0.42, 0.95, len(arr)))
    fig, ax = plt.subplots(figsize=PAIR)
    for d in range(len(arr)):
        y = np.log(at[:, :, d])
        ax.plot(betas.ravel(), y.ravel(), "o", ms=2.2, mec="none", color=blues[d], label="measured" if d == 0 else None)
        for a in range(xs.size):
            ax.plot(betas[a], y[a].mean() + aE[a, d] * (betas[a] - betas[a].mean()), ls="--", lw=0.8, color="black", zorder=4, label="fit in each window" if a == d == 0 else None)
        ax.annotate(rf"$c={cs[arr[d]]:.4f}$", (betas[-1, -1] + 0.25, y[-1, -1]), fontsize=6.5, va="center", color="0.25")
    ax.set_xlim(5.2, 21.2)
    ax.set_xlabel(r"$\beta = 1/kT$ (eV$^{-1}$)")
    ax.set_ylabel(r"$\ln(\langle t \rangle\,/\,\mathrm{s})$")
    ax.set_title("slope fitted in each window (dashed);" "\n" + r"at $E_b/kT=3$:  $E_a^{\mathrm{eff}}=$ " + ", ".join(f"{v:.3f}" for v in aE[results["a_fig"]]) + " eV")
    ax.legend(loc="upper left", handlelength=1.6, handletextpad=0.5)
    fig.tight_layout()
    fig.savefig(f"{outdir}/arrhenius.svg", bbox_inches=Bbox([[0, 0], PAIR]))

    # crossover.svg: f against density, against the analytic curve, c* 
    f_meas, f_an = results["f_meas"], results["f_an"]
    cstar, cstar_an = results["cstar"], results["cstar_an"]
    cf = np.logspace(np.log10(cs[0]), np.log10(cs[-1]), 400)
    worst = np.abs(cstar / cstar_an - 1).max()
    fig, ax = plt.subplots(figsize=PAIR)
    ax.axhline(0.5, color="0.85", lw=0.6, ls=":", zorder=0)
    for a in range(xs.size):
        ax.plot(cf, cf * np.exp(xs[a]) / ((1 - cf) + cf * np.exp(xs[a])), ls="--", lw=0.8, color="black", zorder=2, label="analytic" if a == 0 else None)
        ax.plot(cs, f_meas[a], "o", ms=2.6, mec="none", color=PURPLES[a], zorder=3, label="measured" if a == 0 else None)
        ax.plot([cstar[a]] * 2, [0.0, 0.5], lw=0.8, color=PURPLES[a], zorder=1)
        ax.annotate(rf"${xs[a]:.0f}$", (cstar[a], 0.545), fontsize=7, ha="center", color="0.25")
    ax.set_xscale("log")
    ax.set_xlim(cs[0] * 0.75, cs[-1] * 1.5)
    ax.set_ylim(-0.06, 1.15)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xlabel(r"$c$, trap density (fraction of the $L^2$ sites)")
    ax.set_ylabel(r"$f$, fraction of time spent trapped")
    ax.set_title(r"$c^*=$ " + ", ".join(f"{v:.4f}" for v in cstar) + r" at $E_b/kT=2,3,4,5$" + "\n" + r"against $1/(1+e^{E_b/kT})$:" f"  worst {worst:.3%}")
    ax.legend(loc="upper left", handlelength=1.4, handletextpad=0.5)
    fig.tight_layout()
    fig.savefig(f"{outdir}/crossover.svg", bbox_inches=Bbox([[0, 0], PAIR]))


if __name__ == "__main__":
    L, N_REAL = 50, 400
    XS = np.array([2.0, 3.0, 4.0, 5.0])     # E_b/kT 
    N_BETA, SPAN = 9, 0.10                  # fit window, +-10% in beta
    cs = np.round(np.logspace(-3, np.log10(0.98), 18) * L * L) / L ** 2
    betas = XS[:, None] / E_B * (1 + np.linspace(-SPAN, SPAN, N_BETA))

    results = measure(L, cs, betas, N_REAL)

    print(f"lattice  L = {L}, {L * (L - 1)} interior sites, {N_REAL} placements\n"
          f"\n{'check':<16}{'quantity':<27}measured | expected")
    print(f"jump chain      Q over 3 placements        {'bit-identical' if results['q_same'] else 'THEY DIFFER'} | trap-independent")
    print(f"closed forms    1D t and N, 2D t, absorb   {max(results['e_chain'], results['e_fund'], results['e_trap'], results['e_geom'], results['e_sure']):.1e} | < 1e-12")
    print(f"FACTORISATION   <t> vs [(1-c)+ce^x] G      {results['fac_res']:.2e} | < 1e-3, sampling")
    print(f"                residual in its own s.e.   {results['fac_z']:.2f} | 1")
    print(f"geometry factor G over 6 c x 4 E_b         {results['G_move']:.1e} | < 1e-12")
    print(f"f limits        round trip to f = 1e-6     {results['lim_err']:.1e} | < 1e-6")
    print(f"ACTIVATION E    E_a vs E_mig + E_b f       {np.abs(results['Ea'] / results['Ea_an'] - 1).max():.3%} | < 0.5%")
    print(f"                at the dilute, dense ends  {results['ends'][0]:.4%}, {results['ends'][1]:.4%} | < 0.5%")
    print(f"CROSSOVER       c* vs 1/(1+e^x), worst     {np.abs(results['cstar'] / results['cstar_an'] - 1).max():.3%} | < 2%")
    print(f"convergence     L {L} -> {results['L_CONV']}: c*, then f     {np.abs(results['cstar_conv'] / results['cstar'] - 1).max():.3%}, {results['conv_f']:.4f} | < 2%, < 0.01")
    print(f"precision       worst check | eps cond   {max(results['e_chain'], results['e_fund'], results['e_trap'], results['e_geom'], results['e_sure']):.1e} | {results['floor']:.1e}, cond {results['cond']:.0f}")
    print(f"\nanchors  tau bulk {results['tau_bulk'] * 1e12:.2f} ps;  far row"
          f" {results['t_far_free'] * 1e9:.1f} ns untrapped,"
          f" {results['t_anchor'] * 1e9:.1f} ns at c = 0.05, E_b/kT = 3")

    figures(results)
