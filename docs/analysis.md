# Analysis

Figures produced by `python diffusion.py`. 

## Escape map: $t$ per site

<p align="center">
<img src="../figures/escape_map.svg">
</p>

The axes are row and column in sites; colour is the mean escape time from that site in ns, for
one placement of the 144 traps, marked in green. Row 0 is the absorbing edge, drawn along the
top, and the colour bands run parallel to it: distance from that edge sets the escape time, and
almost nothing else does.

With no traps the field is exactly $\tfrac{1}{2}\tau_0\,r(2L-1-r)$, reproduced at all 2450 sites
to 9e-14: at $E_b/kT = 3$ a bulk site waits 3.71 ps between jumps and the far row sits at 18.2
ns. Traps lift this bodily by $(1-c) + c\,e^{E_b/kT}$, 2.10 at this density, and the row
averages of this one placement hold that lifted form to 2 %. What is left drifts from $-0.4$ %
at the absorbing edge to $+2$ % at the far row, rather than scattering.

The reason is that $t_i = \sum_j N_{ij}\tau_j$ weights every site a defect visits before it
escapes, not just its own. A defect starting far from the edge crosses most of the lattice
first, so a single trap lifts a whole region rather than the site it sits on. That is why the
green squares carry no colour of their own, and why this arrangement, which happens to run
trap-rich away from the edge, drifts upward with distance instead of scattering site to site.

## Factorisation: $\langle t \rangle / G$ against $c$

<p align="center">
<img src="../figures/factorisation.svg">
</p>

Both panels divide an escape time by $G$, the trap-free time at the same site and temperature.
If traps only rescale the clock, that ratio depends on $c$ and $E_b/kT$ alone.

The left panel tests that across the whole sweep, on log axes. Points are solved values of
$\langle t \rangle / G$ averaged over 400 placements; solid lines are $(1-c) + c\,e^{E_b/kT}$
with nothing fitted; the four pairs are labelled $2, 3, 4, 5$ at their right ends by $E_b/kT$,
light to dark. The ratio runs from 1 to 145 over three decades in $c$ and four temperatures, and
the worst departure anywhere on it is 4.2e-03. The entire effect of the traps is that one
bracket. Only the marginal probability enters it, never the joint distribution over placements:
the arrangement is free to be correlated. These traps are negatively correlated, a fixed number
drawn without replacement; clustered, positively-correlated traps are expected to fit just as
well.

The right panel is what that 4.2e-03 is made of. At a single placement the ratio still varies
site to site, by 1.55 % at $c = 0.058$ and $E_b/kT = 3$: 144 traps are a small sample of the
region a defect visits before it leaves. Averaging $n$ placements before taking the ratio
shrinks that scatter as $n^{-1/2}$, the dashed reference line, fitted $-0.483$ over two decades
and reaching 0.170 % at a hundred. This is sampling error about an exact answer, not a departure
from it. The batching is into disjoint groups, which is why it stops at $n = 100$, where four
independent groups still remain per point.

The route itself never moves. $Q$ comes out bit-identical for every placement, and for a crystal
with no traps at all: trapping only ever rescales a clock that geometry has already set.

## Arrhenius: $E_a^{\mathrm{eff}}$ and $c^{\ast}$

<p align="center">
<img src="../figures/arrhenius.svg" width="49%">
<img src="../figures/crossover.svg" width="49%">
</p>

The left panel is the Arrhenius plot. $\ln\langle t \rangle$ against $\beta$ is straight within
each temperature window, one cluster of points per window, and the dashed line through each
cluster is the fit whose slope is $E_a^{\mathrm{eff}}$. Four trap densities are drawn, light to
dark and labelled by $c$ at their right ends. Raising the density lifts the slope from
$E_{mig}$ towards $E_{mig} + E_b$: one migration mechanism reports different activation
energies in crystals differing only in how many traps they hold.

The right panel plots the measured $f$, extracted as $(E_a^{\mathrm{eff}} - E_{mig})/E_b$ from
those slopes, against $c$ on a log axis. Points are measured, the dashed curves are the
analytic $f$, each labelled $2, 3, 4, 5$ by $E_b/kT$ above the vertical line marking its own
$c^{\ast}$. The points lie on the curves at every density and every temperature, each crossing
$f = 1/2$ where that line stands.

The crossover moves left as binding strengthens relative to $kT$: a deeper trap holds the
defect longer, so fewer are needed to dominate the time budget. The dilute limit is out of
reach at the largest $E_b/kT$: the sparsest lattice swept holds two traps in 2500 sites,
$c = 8 \times 10^{-4}$, and those two already carry 10.6 % of the diffusion time.

At $E_b/kT = 5$ the crossover is $c^{\ast} = 0.0067$: a trap detains the defect
$e^{E_b/kT} = 148$ times longer than an ordinary site, so one trap in $1 + 148$ sites
accumulates as much waiting as the 148 around it. Below that density an Arrhenius plot returns
$E_{mig} = 0.5$ eV; above it the slope climbs towards $E_{mig} + E_b = 0.8$ eV, and the
diffusion time rises by a factor of 130 across the densities swept.
