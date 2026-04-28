# PeriodicFluidModel

MATLAB code accompanying the paper:

> Margolius, B. and O'Reilly, O. (2026). *Asymptotic periodic analysis of cyclic stochastic fluid flows.* [Journal reference to be added]

This package computes the asymptotic periodic distribution of a cyclic stochastic fluid model with a time-varying periodic generator $T(t) = T(t+1)$. Given a background Markov chain with generator $T(t)$ and fluid rates $C = \text{diag}(c_1, \ldots, c_m)$, the code computes $P_j(t, x) = P(X(t) \leq x, J(t) = j)$ in the asymptotic periodic limit.

---

## Requirements

- MATLAB R2021a or later (uses `groupsummary`, `parfor`, complex-step derivatives)
- Wegert's `PhasePlot` package for phase portrait visualisation (optional but recommended)
  - Available from: [https://www.mathworks.com/matlabcentral/fileexchange/44145](https://www.mathworks.com/matlabcentral/fileexchange/44145)
  - **Important:** In `pltphase.m`, comment out the line `axis off` to enable axis labels. The version of `pltphase.m` included with this package already has this modification applied.
- `mysimc.m` (included) for Monte Carlo simulation

---

## Files

| File | Description |
|------|-------------|
| `PeriodicFluidModel.m` | Main class. Constructs the model and runs the full pipeline. |
| `ExampleDefinitions.m` | Defines all named example models from the paper. |
| `runExample.m` | Runs a named example end-to-end with optional save/load. |
| `saveModel.m` | Saves computed results to a `.mat` file. |
| `loadModel.m` | Loads saved results and reconstructs the model. |
| `plotAllForModel.m` | Generates all diagnostic plots for a named example. |
| `mysimc.m` | Monte Carlo simulation of the fluid model (Ross rejection method). |
| `wegert/` | Modified Wegert phase portrait code (`axis off` commented out). |

---

## Quick Start

```matlab
% Add package to path
addpath('path/to/PeriodicFluidModel');
addpath('path/to/PeriodicFluidModel/wegert');

% Run Model 2 from the paper
model = runExample('model2', 'n', 400, 'imag_max', 20*pi, 'save', true);

% Generate all plots (saves to figures/model2/)
plotAllForModel('model2', 'save_figs', true, 'close_figs', true);
```

---

## Model Setup

Define your own model by providing a periodic generator function and fluid rate vector:

```matlab
% Example: 3-phase model with sinusoidal rates
a = @(t) 1 + 0.8*sin(2*pi*t);
T = @(t) [-2*a(t)   a(t)    a(t);
           a(t)   -3*a(t)   2*a(t);
           a(t)    a(t)   -2*a(t)];
C = [1 -2 0];   % must be in order [+...+ -...- 0...0]

model = PeriodicFluidModel(T, C, 'n', 400, 'imag_max', 20*pi);
```

**Important:** `C` must be in canonical order — all positive rates first, then negative rates, then zero rates. The rows and columns of `T(t)` must be ordered consistently.

---

## Pipeline

The computation proceeds in five steps, each building on the previous. Steps are triggered automatically when needed:

```matlab
model.computeGamma()         % Step 1: stationary distribution of J(t)
model.findZeros()            % Step 2: zeros of det(I - Wf(0,1,lambda))
model.computeP0()            % Step 3: boundary mass P_j(t,0)
model.computeEigenfunctions() % Step 4: eigenfunction matrix E
```

After Step 4, use the query methods:

```matlab
EX  = model.meanFluid()          % E[X(t)], 1 x n
VX  = model.varFluid()           % Var[X(t)], 1 x n
Ptx = model.cdf(x_vals)         % P(X(t)<=x, J=j), Nx x (m*n)
ptx = model.density(x_vals)     % d/dx P(X(t)<=x, J=j), Nx x (m*n)
```

---

## Constructor Options

```matlab
model = PeriodicFluidModel(T_fun, C_diag, Name, Value)
```

| Option | Default | Description |
|--------|---------|-------------|
| `'n'` | 100 | Time discretisation points per period |
| `'imag_max'` | $20\pi$ | Upper imaginary bound for eigenvalue search |
| `'K_rect'` | 3 | Controls rectangle height $B = (K\_rect + 0.5) \times \max(\text{im\_step})$ |
| `'WfMethod'` | `'expm'` | Method for computing $W_f$: `'expm'` or `'ode'` |
| `'ResMethod'` | `'adjugate'` | Residue method: `'adjugate'` or `'eigvec'` |

**Choosing `n`:** Use `n = 100` for quick checks; `n = 400` for publication-quality results. The eigenfunction computation cost scales as $O(n^2 m^3)$ per eigenvalue.

**Choosing `imag_max`:** Should not be close to a multiple of $2\pi/c_j$ for any positive rate $c_j$. Use $20\pi$ as a starting point; increase if the CDF or density shows truncation artifacts (horizontal striping in density plots).

---

## Plot Methods

```matlab
model.plotGamma()                          % gamma_j(t): stationary phase dist
model.plotDrift()                          % drift E[dX/dt | J(t)]
model.plotP0()                             % P_j(t,0): mass at fluid level 0
model.plotMeanFluid()                      % E[X(t)] +/- sigma(t)
model.plotCDF(phase, x_levels)            % CDF for one phase
model.plotDensity(x_vals, phases)         % density contourf plots
model.plotPhasePlot(realRange, imagRange, gridSize)  % phase portrait
model.plotCDFwithSim(pdistcum, phase, x_levels)     % CDF with simulation
```

All plot methods auto-compute prerequisites as needed.

---

## Simulation

```matlab
% Run simulation (requires mysimc.m on path)
pdistcum = model.simulate(50000, 20, 5, 100);
%                          ^      ^   ^   ^
%                      events  runs max  mesh
%                           per run fluid

% Overlay with analytic CDF
model.plotCDFwithSim(pdistcum, 4, [0.2 0.5 1.0 2.0 5.0]);
```

Simulation results can be saved and reloaded:
```matlab
save('saved_models/mymodel_sim.mat', 'pdistcum');
load('saved_models/mymodel_sim.mat', 'pdistcum');
```

`plotAllForModel` automatically runs and saves the simulation if no saved file is found.

---

## Save and Load

```matlab
% Save after computation
saveModel(model, 'mymodel', 'saved_models');

% Load later (reconstructs T_fun from ExampleDefinitions)
model = loadModel('mymodel', 'saved_models');

% Run a named example with save/load
model = runExample('model2', 'n', 400, 'save', true);
model = runExample('model2', 'load', true);   % skip recomputation
```

**Note:** Function handles cannot be saved to `.mat` files. `loadModel` reconstructs `T_fun` from `ExampleDefinitions.m` using the stored `example_name`. For custom models, either add your model to `ExampleDefinitions.m` or reconstruct `T_fun` manually after loading.

---

## Example Suite

The following examples from the paper are available via `ExampleDefinitions`:

```matlab
ExampleDefinitions()   % list all available examples
```

| Name | Description |
|------|-------------|
| `model1` | Paper Model 1: piecewise-constant $T$, $C=[4\ 2\ 1\ {-12}\ 0\ 0]$ |
| `model2` | Paper Model 2: sinusoidal $T$, $C=[1\ 2\ 3\ {-4}\ {-3}\ 0]$ |
| `const2` | Constant-rate version of model2 |
| `ex1`–`ex5` | Sinusoidal and constant-rate 6-phase examples |
| `ex6`–`ex7` | Sinusoidal variants with larger negative rates |
| `ex8`–`ex9` | 5-phase models without zero-rate phase |
| `ex10` | Small transition rates (scaled by 0.01) |
| `ex11`–`ex14` | Variants with different $C$ structures and piecewise-constant $T$ |
| `ex15`–`ex16` | 5-phase models |
| `ex17`–`ex18` | 6-phase sinusoidal models |

See `ExampleDefinitions.m` for full details and the mapping from original v17 case numbers.

---

## Known Limitations

- **Large transition rates:** When rates are large relative to the fluid rates, `expm` accuracy degrades. Use a finer mesh (`n`) or the `'ode'` `WfMethod` option.
- **Density striping:** Truncating the eigenfunction series at `imag_max` produces Gibbs-phenomenon oscillations in the density. Increase `imag_max` or apply a window function to reduce this.
- **Coincident asymptotic families:** Models where two $S_+$ phases have the same asymptotic real part $\bar{T}_{jj}/c_j$ require the coincident family logic in `findZeros`. These are handled automatically but may require tighter `tol_dedup` (default `1e-9`).
- **`imag_max` near zero heights:** Choose `imag_max` away from multiples of $2\pi/c_j$ to avoid boundary effects in the eigenvalue count.
- **Very small transition rates (ex10):** When transition rates are very small relative to fluid rates (e.g. scaled by 0.01), the eigenvalues cluster near the origin with very small negative real parts and closely-spaced imaginary parts. The simulation fails due to duplicate timestamps in the event history, and the analytic computation requires a much larger `imag_max` and finer mesh to resolve the closely-spaced eigenvalue families. This example is included for completeness but reliable numerical results require special treatment.

---

## Algorithm Overview

The algorithm follows Margolius & O'Reilly (2026). The key steps are:

1. **Gamma:** Solve the ODE $\dot{\gamma}(t) = \gamma(t) T(t)$ over one period to get the stationary phase distribution.
2. **Zeros:** Find zeros of $\det(I - W_f(0,1,\lambda))$ in $\{\text{Re}(\lambda) \leq 0\}$ using the argument principle (rectangle region) and asymptotic seeds (above the rectangle).
3. **P0:** Build the $(mn \times mn)$ matrix $H = \sum_k A_k$ and solve for $P(t,0)$.
4. **Eigenfunctions:** Compute $\phi(t, \lambda_k) = P(t,0) \cdot A_k$ for each eigenvalue.
5. **Distribution:** Evaluate $P(t,x) = \sum_k e^{\lambda_k x} \phi(t,\lambda_k) + \gamma(t)$.

---

## Citation

If you use this code, please cite:

```
Margolius, B. and O'Reilly, O. (2026). Asymptotic periodic analysis of 
cyclic stochastic fluid flows. [Journal], [Volume], [Pages].
```

---

## License

[To be specified — CC BY 4.0 recommended for Mendeley Data]
