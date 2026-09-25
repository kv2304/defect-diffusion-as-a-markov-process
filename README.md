# Defect diffusion as a Markov Process

Absorbing Markov chain model of trap-limited defect diffusion on a 2D lattice, with varying trap
density and temperature. Computes the mean escape time from every site and recovers
the crossover trap density $c^{\ast}$ from the collapse of the Arrhenius slope onto the
trap-limited barrier, agreeing with the analytical solution
$c^{\ast} = 1/(1 + e^{E_b/kT})$ to 0.342 %.

## Physical concepts

A defect is a missing or misplaced atom in a crystal lattice. It migrates by hopping to a
neighbouring site, each hop clearing an activation barrier at a rate set by temperature. Sites
that bind it more tightly are traps, and the trap density $c$ is the fraction of sites that are
traps.

Hopping is memoryless: the site $X_{n+1}$ occupied after $n+1$ hops depends on the site $X_n$
occupied after $n$, and on no earlier one.

```math
\Pr\big(X_{n+1} = S_j \mid X_n = S_i,\, X_{n-1}, \ldots, X_0\big)
 = \Pr\big(X_{n+1} = S_j \mid X_n = S_i\big) = P_{ij}
```

Every transition is then an entry of the one matrix $P$. Listing the interior sites before the
absorbing edge splits $P$ into the interior block $Q$ and the column of absorption
probabilities, and $N = (I-Q)^{-1}$ holds the mean number of visits to each site before
absorption, so $N$ acting on the site waits $\tau$ is the mean escape time from every site.

```math
t = N\tau, \qquad N = (I-Q)^{-1}, \qquad \tau_i = \frac{e^{E_{a,i}/kT}}{z_i\nu}
```

Raising $c$ raises that time without altering the route taken, so the activation energy an
Arrhenius plot returns lies between the migration barrier $E_{mig}$ and the trap-limited
$E_{mig} + E_b$, crossing between them at the density where traps hold the defect half the time.

```math
c^{\ast} = \frac{1}{1 + e^{E_b/kT}}
```

## Results

Complete method and analysis can be read in [docs/](docs/).

| | measured | independent | agreement |
|---|---|---|---|
| crossover $c^{\ast}$, four temperatures | 0.1193, 0.0474, 0.0179, 0.0067 | 0.1192, 0.0474, 0.0180, 0.0067 from $1/(1 + e^{E_b/kT})$ | **worst 0.342 %** |
| $E_a^{\mathrm{eff}}$, 4 temperatures x 18 densities | 72 fitted slopes | $E_{mig} + E_b f$ at each | worst 0.195 % |
| factorisation of $\langle t \rangle$ | $\langle t \rangle$ over 400 placements | $[(1-c) + c\,e^{E_b/kT}]\,G$ | 6.34e-04 |
| geometry factor $G$ | 24 values, 6 densities x 4 binding energies | the trap-free $G$ at the same $kT$ | worst 2.2e-15, or 10 $\varepsilon$ |
| jump chain $Q$ | $Q$ for three placements | assembled from $z_i$ alone | bit-identical |
| five closed forms | the sparse solve | 1D chain $t$, $N_{ij}$, $\tfrac{1}{2}\tau_0 r(2L{-}1{-}r)$, one trap, absorption certain | worst 1.1e-13 |
| $f$ inverted for $c$ | $f$ re-evaluated at the recovered $c$ | the targets, $f = 10^{-6}$ and $1 - 10^{-6}$ | 8.2e-11 |
| $c^{\ast}$ stability, two lattice sizes | range 0.444 % | 2450 and 4830 interior sites | $f$ moves 0.0028 |
| scatter in $t/G$ as placements are averaged | exponent $-0.483$ | $-1/2$ for a sampling law | 3.4 % |
| worst error against the closed forms | 1.1e-13 | float64 bound $\varepsilon\,\mathrm{cond}_1(I-Q)$ = 2.4e-12 | 22x inside the bound |

## Run

```
git clone https://github.com/kv2304/defect-diffusion-as-a-markov-process.git
cd defect-diffusion-as-a-markov-process
pip install -r requirements.txt
python diffusion.py     # ~4s, redraws figures/
```
