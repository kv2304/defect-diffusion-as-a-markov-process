# Method

The state after $n$ jumps is the site the defect occupies, and where it goes next depends on
that site alone: the barrier is set by its binding energy, the directions by its bonds.

## Given

<p align="center">
<img src="../figures/lattice.svg">
</p>

Defects hop on an $L \times L$ lattice, $L = 50$. Row 0 is an absorbing edge and the other three
sides reflect, so the far row is a mirror and the lattice is half of a symmetric system $2L-1$
rows across. Each of the $L(L-1) = 2450$ interior sites is a state $S_i$, each jump an arrow to
a neighbour, and row 0 is the single absorbing state $S_0$.

A fraction $c$ of sites are traps, drawn without replacement so each site has marginal
probability $c$, swept over 18 logarithmic points from $8 \times 10^{-4}$ to $0.98$. Leaving an
ordinary site costs $E_{mig} = 0.5$ eV and leaving a trap costs $E_{mig} + E_b$ with
$E_b = 0.3$ eV, both at attempt frequency $\nu = 10^{13}$ s⁻¹. Temperature enters only as
$E_b/kT$, swept over $\lbrace 2, 3, 4, 5 \rbrace$, so a trap holds a defect $e^{E_b/kT}$ times
longer than an ordinary site: 7 times at one end of the sweep, 148 at the other.

## Step 1 — build the transition matrix  (`chain`)

State the Markov property. With $X_n$ the site occupied after $n$ jumps,

```math
\Pr\big(X_{n+1} = S_j \mid X_n = S_i,\, X_{n-1}, \ldots, X_0\big)
 = \Pr\big(X_{n+1} = S_j \mid X_n = S_i\big) = P_{ij}
```

so the process is carried by the single matrix $P$. Collect the occupation probabilities into a
row vector $\mathbf{p}_n$ with $(\mathbf{p}_n)_i = \Pr(X_n = S_i)$, and list the interior states
first:

```math
\mathbf{p}_{n+1} = \mathbf{p}_n P, \qquad
P = \begin{pmatrix} Q & \mathbf{a} \\ 0 & 1 \end{pmatrix}
```

where $Q$ is the interior-to-interior block and $\mathbf{a}$ the column of probabilities of
reaching $S_0$ on the next jump. Absorption is permanent, so
$\Pr(X_{n+1} = S_0 \mid X_n = S_0) = 1$ is the bottom row and the block beneath $Q$ is zero.

Fill $Q$. Escape from a site is isotropic, so with $z_i$ the coordination of site $i$ — 4 in the
bulk, 3 on an edge, 2 in a corner —

```math
Q_{ij} = \frac{1}{z_i} \;\text{ for } j \text{ a neighbour of } i, \qquad Q_{ij} = 0 \;\text{ otherwise}
```

Bonds running into row 0 leave $Q$ and become entries of $\mathbf{a}$, so every row of $P$
closes on 1.

Worked instance. On a $3 \times 3$ lattice the interior states $S_1, \ldots, S_6$ are the sites
$(1,0), (1,1), (1,2), (2,0), (2,1), (2,2)$ with $z = 3, 4, 3, 2, 3, 2$:

```math
Q=\begin{pmatrix}
0 & \tfrac{1}{3} & 0 & \tfrac{1}{3} & 0 & 0 \\
\tfrac{1}{4} & 0 & \tfrac{1}{4} & 0 & \tfrac{1}{4} & 0 \\
0 & \tfrac{1}{3} & 0 & 0 & 0 & \tfrac{1}{3} \\
\tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0 \\
0 & \tfrac{1}{3} & 0 & \tfrac{1}{3} & 0 & \tfrac{1}{3} \\
0 & 0 & \tfrac{1}{2} & 0 & \tfrac{1}{2} & 0
\end{pmatrix},
\qquad
\mathbf{a}=\begin{pmatrix}
\tfrac{1}{3} \\ \tfrac{1}{4} \\ \tfrac{1}{3} \\ 0 \\ 0 \\ 0
\end{pmatrix}
```

Only $S_1, S_2, S_3$ neighbour row 0, hence three non-zero entries in $\mathbf{a}$.

Fill $\tau$. Each of the $z_i$ bonds is attempted at frequency $\nu$ and cleared with
probability $e^{-E_{a,i}/kT}$, so the escape rate from site $i$ is $z_i \nu e^{-E_{a,i}/kT}$ and
its reciprocal is the mean wait:

```math
\tau_i = \frac{e^{E_{a,i}/kT}}{z_i \nu},
\qquad E_{a,i} = \begin{cases} E_{mig} & \text{ordinary site} \\ E_{mig} + E_b & \text{trap} \end{cases}
```

Return $Q$, $\mathbf{a}$ and $\tau$. Note which carries what: only $\tau$ carries a trap, and
$Q$ is assembled from coordination alone.

## Step 2 — solve for the mean escape time  (`escape`)

Let $A = \min\lbrace n : X_n = S_0 \rbrace$ be the jump on which the defect is absorbed. The
escape time is the sum of the waits along the path, and the quantity wanted is its expectation
from each starting state:

```math
T = \sum_{k=0}^{A-1} \tau_{X_k}, \qquad t_i = \mathbb{E}\big[T \mid X_0 = S_i\big]
```

Condition on the first jump. The defect waits $\tau_i$, then moves to some $S_j$:

```math
t_i = \tau_i + \sum_j Q_{ij}\, t_j
```

the sum running over interior states only, since a jump to $S_0$ ends the walk. Gather into a
vector and solve:

```math
(I - Q)\,t = \tau, \qquad t = N\tau, \qquad N = (I - Q)^{-1}
```

Justify the inverse. Absorption is certain, because every interior state has a path out, so
$Q^n \to 0$ and $N = \sum_{n \ge 0} Q^n$ converges. Its entries are then

```math
N_{ij} = \mathbb{E}\Big[\textstyle\sum_{n \lt A} \mathbf{1}\lbrace X_n = S_j \rbrace \;\Big|\; X_0 = S_i\Big]
```

where $\mathbf{1}\lbrace \cdot \rbrace$ is 1 when its condition holds and 0 otherwise, so the
sum counts visits to $S_j$ before absorption and $N_{ij}$ is the mean number of them. Hence
$t = N\tau$ is time-per-visit summed over visits, and only the mean of each wait enters.

Implementation. $N$ is dense and is never formed. $(I - Q)$ is sparse and is factorised once
with `scipy.sparse.linalg.splu`, and since $Q$ holds no trap and no temperature, that one
factorisation serves every density, temperature and arrangement.

Closed form for a check. With a uniform activation energy, $\tau_i = \tau_0 / z_i$ where
$\tau_0 = e^{E_{mig}/kT}/\nu$, so $z_i \tau_i$ is constant. Multiply the recursion by $z_i$ to
get $z_i t_i = \tau_0 + \sum_{j \sim i} t_j$, the sum running over the neighbours of $i$, and
since $z_i t_i$ counts $t_i$ once per bond,

```math
\sum_{j \sim i} (t_j - t_i) = -\tau_0 \qquad \Longrightarrow \qquad t(r) = \tfrac{1}{2}\,\tau_0\, r\,(2L - 1 - r)
```

This is a discrete Poisson equation on the lattice graph, with $t = 0$ on row 0 as the absorbing
boundary and a missing bond on a reflecting side imposing zero flux. Its solution depends on the
row index $r$ alone.

## Step 3 — average over trap placements  (`disorder`)

The traps are randomly placed, so the average is taken over both the possible walks and the
different trap arrangements, with the walks averaged at each fixed arrangement:

```math
\langle t_i \rangle = \mathbb{E}_{\text{traps}}\,\mathbb{E}_{\text{paths}}\big[T \mid X_0 = S_i\big]
```

The two expectations commute because $t$ is linear in $\tau$, so the outer one passes through
$N$ and lands on $\tau$ alone. Site $j$ is a trap with marginal probability $c$, so

```math
\langle \tau_j \rangle = \Big[(1-c) + c\,e^{E_b/kT}\Big]\,\frac{e^{E_{mig}/kT}}{z_j \nu}
```

That bracket is the same at every site, so it leaves the sum:

```math
\langle t_i \rangle = \Big[(1-c) + c\,e^{E_b/kT}\Big]
                      \sum_j N_{ij}\, \frac{e^{E_{mig}/kT}}{z_j \nu}
```

The remaining sum is the geometry factor $G_i$: one value per starting site, free of $c$ and
$E_b$, and equal to the trap-free escape time at the same $kT$.

## Step 4 — extract the effective activation energy  (`measure`)

Take logs, with $\beta = 1/kT$ and $G_i$ depending on $\beta$ only through $e^{\beta E_{mig}}$:

```math
\ln\langle t_i \rangle = \beta E_{mig}
                     + \ln\!\Big[(1-c) + c\,e^{\beta E_b}\Big] + \mathrm{const}_i,
\qquad \mathrm{const}_i = \ln \sum_j \frac{N_{ij}}{z_j \nu}
```

$\mathrm{const}_i$ carries no $\beta$, so every site gives the same slope, and so does their
average. Differentiate:

```math
E_a^{\mathrm{eff}} = \frac{\mathrm{d}\ln\langle t \rangle}{\mathrm{d}\beta} = E_{mig} + E_b f, \qquad
f = \frac{c\,e^{\beta E_b}}{(1-c) + c\,e^{\beta E_b}}
```

where $f$ is the fraction of time spent trapped. Set $f = 1/2$ for the crossover density:

```math
c^{\ast} = \frac{1}{1 + e^{E_b/kT}}
```

Fitting. Fit the slope over a window of $\pm 10\%$ in $\beta$ about each $E_b/kT$, 9 points per
window: narrow enough that the fit returns the derivative at its centre. Read $c^{\ast}$ in the
coordinate where the relation is exactly straight,

```math
\ln\!\frac{f}{1-f} = \ln\!\frac{c}{1-c} + \frac{E_b}{kT}
```

Results of these fits, and of the closed-form and stability checks, are tabulated in the
[README](../README.md).
