# multimode-parametric-network

MATLAB (also Octave ≥ 8) tools for the Gaussian steady state of a network of damped
bosonic modes coupled by two-mode squeezing and beam splitters. A network is given
as two lists of tuples. The conditional steady state under continuous heterodyne
monitoring is also supported.

```matlab
addpath('path/to/multimode-parametric-network')

squeezings    = [1 2 0.25;            % (i, j, g)
                 2 3 0.10 pi/2];      % (i, j, g, phase)
beamsplitters = [1 3 0.40];           % (i, j, h) or (i, j, h, phase)

sim = simulate_network(3, squeezings, beamsplitters);            % unconditional
sim = simulate_network(3, squeezings, beamsplitters, 'Eta', 1);  % heterodyne-monitored

logneg(sim.sigma, 1, 2)       % log-negativity between modes 1 and 2
logneg(sim.sigma, 1, [2 3])   % ... between mode 1 and modes {2,3}
purity(sim.sigma, [1 2])      % purity of the reduced state of modes 1,2
```

## Input formats

Each coupling list can be
- a numeric matrix, one row per coupling: `[i j amp]` or `[i j amp phase]`
- a cell array of tuples: `{[1 2 0.25], {2, 3, 0.1, pi/2}}`
- `[]` for none.

There are three ways to set a phase, and they add up: complex `amp`, a 4th tuple
entry (radians), or the `'SqueezingPhases'` / `'BeamsplitterPhases'` options
(a vector with one phase per tuple). Use the options to sweep or randomize phases
without rebuilding the tuple lists.

Other options: `'Kappa'` (damping rate, scalar or per mode, default 1),
`'Eta'` (heterodyne efficiency, scalar or per mode, default 0 = unmonitored),
`'Detunings'` (mode frequencies, default 0).

## Model

ħ = 1, rotating frame, every mode damped at rate κⱼ into vacuum:

| tuple | Hamiltonian term | Langevin drift |
|---|---|---|
| squeezing `(i,j,g)` | `i(g aᵢ†aⱼ† − g* aᵢaⱼ)` | `ȧᵢ += g aⱼ†`, `ȧⱼ += g aᵢ†` |
| beam splitter `(i,j,h)` | `h aᵢ†aⱼ + h* aⱼ†aᵢ` | `ȧᵢ += −i h aⱼ`, `ȧⱼ += −i h* aᵢ` |

Quadratures are `x = a + a†` and `p = −i(a − a†)`, ordered `(x1,p1,x2,p2,…)`, so the
vacuum covariance matrix is the identity and purity is `1/sqrt(det σ)`.

- `Eta = 0`: `Aσ + σAᵀ + D = 0` (Lyapunov equation; it is an error above threshold).
- `Eta > 0`: `Aσ + σAᵀ + D − (σ − I) W (σ − I) = 0` with `W = diag(ηⱼκⱼ/2) ⊗ I₂`.
  This is the conditional state under heterodyne monitoring of every output. It is
  solved with an ordered Schur decomposition, so no Control System Toolbox is needed.

The squeezing convention matches `graphtests.m` (`Adj(i,j) = r·e^{iθ}` ↔ `g`). Its
heterodyne-noise parameter `s = sec²b + n·tan²b` corresponds to `η = 2/(1+s)`.

## Files

| file | purpose |
|---|---|
| `simulate_network.m` | main solver |
| `logneg.m`, `purity.m`, `reduced_cm.m`, `symplectic_eigs.m` | entanglement / purity of (sub)systems |
| `apply_loss.m` | pure-loss channel (`GC` in `graphtestsasym.m`) |
| `analysis/run_comb_analysis.m` | the `graphtests.m` + `graphtestsasym.m` analysis rebuilt with `simulate_network` |
| `analysis/comb_pumps.m`, `analysis/comb_index.m` | tuple lists for the pumped 44-mode frequency comb |
| `tests/test_network.m` | two small examples with checks (2-mode squeezer; squeezer + 3-mode beam-splitter ring with flux) |

## Comb analysis (`analysis/run_comb_analysis.m`)

The comb has 44 modes with labels ±1…±22 and 1–15 pumps, each with r = 0.25 κ.
For modes (−1, 1) the script computes the heterodyne-conditioned E_N and purity in six cases:
- sym / asym pump geometry (asym: from pump 3 on, each pump couples every comb mode
  to its own new idler mode, so 15 pumps means 616 modes)
- equal phases / random phases / random phases with noise (η = 1/2)

Outputs go to `analysis/results/`: the PNGs `symasymcomplogneg.png` and
`symasymcomppur.png`, plus `comb_results.{csv,mat}` (a full run takes about 10 min in Octave because the asym
system grows to 616 modes). The plots overlay the values
hardcoded in `graphtests.m` as circles. The equal-phase curves reproduce those values
exactly. The random-phase curves depend on the phase draw (`seed`), and the original
draws were unseeded.

## Tests

```matlab
run('tests/test_network.m')
```
