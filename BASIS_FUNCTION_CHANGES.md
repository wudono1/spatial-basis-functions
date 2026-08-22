# Basis-Function Extension of the CORGI Downscaling Model

**A complete record of changes made to `corgi-examples` to replace the Gaussian-process
covariance of Zheng et al. (2025) with a basis-function representation.**

---

## 0. What changed, in one page

### The original model (MEE paper, Equation 1)

The fine-resolution covariate field is a Gaussian process over Basic Areal Units (BAUs):

$$X_s = \mathbf{1}_{N_D}\gamma_0 + R\gamma + \omega(\mathbb{S}), \qquad
\omega(\cdot) \sim GP\big(0,\ \sigma^2 \exp(-\|s-s'\|/\psi)\big)$$

so the fine-scale covariance is an $N_D \times N_D$ matrix built directly from pairwise
BAU distances:

$$C_S = \sigma^2 K_\psi^{(\text{BAU})}, \qquad K_\psi^{(\text{BAU})}[i,j] = \exp(-d_{ij}/\psi)$$

### The new model (basis functions, Cressie)

The spatial process is instead expanded in $r$ fixed basis functions with random coefficients:

$$\omega(s) = \phi(s)^\top \alpha + \delta(s), \qquad \alpha \sim N(0, \sigma^2 K_\psi)$$

where $\phi(s)$ is the $r$-vector of basis functions evaluated at $s$, and $K_\psi$ is now an
$r \times r$ exponential kernel over the **basis-function centroids (knots)**, not over BAUs.
Assuming $\delta(\mathbb{S}) = 0$ (imposed by the PI), the fine-scale covariance becomes

$$C_S = \Phi\,(\sigma^2 K_\psi)\,\Phi^\top$$

where $\Phi$ is the $N_D \times r$ matrix of basis functions at the BAU centroids
(the paper's $\phi(\mathbb{S})^\top$).

### The consequence for the block level

Aggregation is unchanged: $X_B = W_B X_s$ with $W_B$ the $M_D \times N_D$ averaging matrix.
Therefore

$$K_B(\psi) = W_B \Phi\, K_\psi\, \Phi^\top W_B^\top = \Psi_B K_\psi \Psi_B^\top,
\qquad \Psi_B \equiv W_B\Phi \ \ (M_D \times r)$$

**This single substitution is the entire mathematical change.** The block-level law keeps the
same form, $[X_B \mid \tilde\gamma, \sigma^2, \psi] = N(R_B\tilde\gamma,\ \sigma^2 K_B(\psi))$,
so every downstream formula — the collapsed likelihood, the Metropolis step, the Gibbs
updates — carries over untouched. Only *how* $K_B$ is constructed differs.

> **The key insight that made this tractable:** because $K_B$ retains the same shape
> ($M_D \times M_D$) and meaning, the ~150 lines of MCMC machinery built on top of it are
> blind to the change. Only two lines in the sampler had to be rewritten.

### The two constraints this introduces

1. **Rank.** $\operatorname{rank}(K_B) \le \min(M_D, r)$. With $\delta = 0$ there is no nugget
   to restore full rank, so **$r \ge M_D$ is a well-posedness requirement**, not a tuning
   preference. Violating it makes $K_B$ singular and the likelihood undefined.
2. **Speed.** The payoff. The old code rebuilt an $N_D \times N_D$ matrix every MCMC
   iteration; the new code rebuilds an $r \times r$ one. With $N_D = 3000$, $r = 120$ that is a
   ~600× reduction in the dominant per-iteration cost.

### Dimensions in this project

| Symbol | Meaning | Value | Source |
|---|---|---|---|
| $N_D$ | fine BAUs | **3000** | `50 × 60` grid, `simulate_covariate_data.R` |
| $r$ | basis functions | **120** | `10 × 12` knots |
| $M_D$ | coarse blocks | **30** | `kk_x=5 × kk_y=6` |
| $p+1$ | design columns | **1** | intercept only (`coar_covar = rep(1, Md)`) |

$M_D \le r \ll N_D$ holds comfortably: $30 \le 120 \ll 3000$. Both constraints are satisfied
with room to spare.

---

## 1. Notation map: paper ↔ code

The original code used variable names that conflict with the papers. This was a persistent
source of confusion, so the new functions were renamed to match the papers.

| Paper | Code (old) | Code (new) | Meaning |
|---|---|---|---|
| $\psi$ | `phi` | `psi` | range parameter (inside the exponential) |
| $\sigma^2$ | `sigmasq` | `sigmasq` | variance parameter |
| $\phi \equiv (\psi,\sigma^2)$ | — | — | the *pair*; never named in code |
| $X_B$ | `y` | `xB` | coarse response |
| $R_B$ | `X` | `RB` | **coarse** design ($M_D \times p$) |
| $\tilde R$ | `X` (in `pred*`) | `Rtilde` | **fine** design ($N_D \times p$) |
| $\tilde\gamma=(\gamma_0,\gamma)$ | `bb` | `gam` | regression coefficients |
| $\tilde\mu_{\tilde\gamma}$ | `bb_hat` | `gam_hat` | GLS estimate / posterior mean |
| $K_\psi$ | `bau_corrmat` | `Kpsi` | correlation matrix |
| $K_B(\psi)$ | `blk_corrmat` | `KB` | block-level correlation |
| $S^2(\psi)$ | `aa` | `S2` | generalised residual sum of squares |
| $\Psi_B = W_B\phi(\mathbb{S})^\top$ | — | `PsiB` | *new*; not named in the papers |
| $\phi(\mathbb{S})^\top$ | — | `Phi` | basis functions at BAU centroids |

> **Two traps worth remembering.**
> **(a)** The code's `phi` is the paper's $\psi$, *not* the paper's $\phi$.
> **(b)** `X` means different things in different functions: the **coarse** design $R_B$ in
> `fitGPBAU_new`, but the **fine** design $\tilde R$ in the prediction functions. They are
> related by $R_B = W_B\tilde R$. This is why the prediction argument is named `Rtilde`, not `RB`.
>
> **`Ψ_B` is my notation, not the paper's.** Equation (14) writes
> $K_B(\psi) = W_B\phi(\mathbb{S})^\top K_\psi \phi(\mathbb{S})W_B^\top$. The left and right
> factor pairs are transposes of each other, so naming the left pair $\Psi_B$ collapses this to
> $\Psi_B K_\psi \Psi_B^\top$. Identical matrix — the name just lets us compute the
> $\psi$-independent part once instead of every iteration. The paper does the same thing for
> the design matrix when it defines $R_B \equiv W_B\tilde R$.

---

## 2. Chronological record of changes

### Step 0 — Data generation (`data-examples/simulate_covariate_data.R`)

**Done before this log begins; the starting point.**

The original generated each field as a GP over the full fine grid:

```r
xx <- t(rgp(1, fine_grid, cov.mod = "powexp", nugget = 0, smooth = 1,
            sill = sill_vals[j], range = range_vals[i]))
```

This was replaced with a basis-function construction: draw coefficients $\alpha$ from a GP over
the **knots**, then map up to the fine grid.

```r
x_seq_knots <- seq(0, 1, length = 10)
y_seq_knots <- seq(0, 1, length = 12)
fine_grid_knots <- as.matrix(expand.grid(x_seq_knots, y_seq_knots))   # 120 knots

D   <- rdist(fine_grid, fine_grid_knots)          # 3000 x 120
Phi <- exp(-D^2 / (2 * length_scale^2))           # Gaussian RBF basis

alpha <- t(rgp(1, fine_grid_knots, cov.mod = "powexp", nugget = 0, smooth = 1,
               sill = sill_vals[j], range = range_vals[i]))            # cov = sill * exp(-d/range)
xx <- Phi %*% alpha
```

**Why this establishes the simulation truth.** `rgp` with `cov.mod = "powexp", smooth = 1`
gives covariance $\text{sill}\cdot\exp(-d/\text{range})$. Matching against
$\operatorname{cov}(\alpha) = \sigma^2 K_\psi$ with $K_\psi = \exp(-d/\psi)$ gives an exact
correspondence, which is what makes Step D's recovery check possible:

$$\psi_{\text{true}} = \texttt{range\_vals[i]}, \qquad
\sigma^2_{\text{true}} = \texttt{sill\_vals[j]}, \qquad
\gamma_{0,\text{true}} = 0$$

The intercept is zero because $\alpha$ has mean zero and no mean term is added.

Note the basis kernel and the coefficient kernel are *different*: the basis functions are
Gaussian RBFs $\exp(-d^2/2\ell^2)$, while $K_\psi$ is exponential $\exp(-d/\psi)$. They play
different roles and there is no requirement that they match.

---

### Step 1 — The new sampler (`spds/src/mcmc_sd.cpp` → `fitGPBAU_new`)

`fitGPBAU` was copied to `fitGPBAU_new`, leaving the original untouched for comparison.

**1a. Renaming.** All local variables were renamed per the table in §1 (code unchanged,
comments rewritten). The R-facing list keys were deliberately **left alone** —
`prior["u_phi"]`, `starting["phi"]`, `tuning["se_phi"]`, and the output keys
`"phi"/"sigmasq"/"beta"` — so no R code had to change. Only the C++ locals differ:

```cpp
double u_psi = prior["u_phi"];   // key unchanged, local renamed
```

**1b. New inputs**, read from `mod_param` exactly as `aggmat` is:

```cpp
arma::mat Phi = mod_param["Phi"];                           // Nd x r
arma::mat centroid_distmat = mod_param["centroid_distmat"]; // r x r
arma::mat PsiB = aggmat * Phi;                              // Md x r -- computed ONCE
```

`PsiB` is hoisted out of the MCMC loop because it does not depend on $\psi$. This is where the
speedup comes from.

**1c. The actual change — two lines, twice.** In the initial evaluation:

```cpp
// OLD:  Kpsi = expCorr(bau_distmat, psi);        // Nd x Nd
//       KB   = aggmat * Kpsi * aggmat.t();
Kpsi = expCorr(centroid_distmat, psi);            // r x r
KB   = PsiB * Kpsi * PsiB.t();                    // Md x Md
```

and identically in the Metropolis candidate block (`cand_Kpsi`, `cand_KB`).

**Everything else in the sampler is unchanged** — `KB_log_det`, `quad_corrmat`, `quad_log_det`,
`proj_xB`, `gam_hat`, `ee`, `S2`, the collapsed log-likelihood, the logit Jacobian, the
acceptance test, and both Gibbs draws. They depend only on `KB`, `RB`, `xB`, and `KB` kept its
shape and meaning.

**1d. Bug fix.** The original had `cand_S2` computed from `ee` (current residuals) instead of
`cand_ee`. Corrected in the new function to match Eq. (14):
$S^2(\psi') = (X_B - R_B\tilde\mu_{\tilde\gamma}(\psi'))^\top K_B(\psi')^{-1}(\cdots)$.

**1e. Why the degrees-of-freedom term did not change.** `pp = RB.n_cols` still counts only the
*regression* coefficients. The basis coefficients $\alpha$ are **integrated out** into $K_B$ as
a random effect, not estimated as fixed parameters, so $r$ never enters
$\tilde u_{\sigma^2} = u_{\sigma^2} + (M_D - p - 1)/2$.

---

### Step 2 — Smoke test and the rank discovery (`sandbox/smoke_test_fitgpbau_new.R`, new file)

A throwaway script, deliberately outside `spds/` so it is not part of the package. It calls
`fitGPBAU_new` **directly** with raw C++ list keys, on a tiny problem, and reports
`rank(PsiB)` and the minimum eigenvalue of `KB` before invoking the sampler.

**What it revealed.** Setting $r < M_D$ produces a Cholesky failure. The reason:

$$\operatorname{rank}(C_S) \le r \ \xrightarrow{\ W_B(\cdot)W_B^\top\ }\
\operatorname{rank}(K_B) \le \min(M_D, r)$$

If $r < M_D$ then $K_B$ ($M_D \times M_D$) has at least $M_D - r$ zero eigenvalues. It is
positive **semi**-definite, not positive definite, so `log_det_sympd` ($\log|K_B| = -\infty$),
`solve` (no inverse), and `mvnrnd` (no factorisation) all fail.

**Why this is not fixable by better linear algebra.** With $\delta(\mathbb{S})=0$ and
$r < M_D$, the model claims some contrast across blocks has *exactly zero variance*. $X_B$ then
has no density with respect to Lebesgue measure on $\mathbb{R}^{M_D}$ — the likelihood does not
exist. The Woodbury identity does not help either: it requires the form $A + UCV^\top$ with $A$
invertible, and pure $\Psi_B K_\psi \Psi_B^\top$ has no such $A$.

**The two possible remedies:**
- **Path A** — require $r \ge M_D$. Chosen, since $120 \ge 30$ holds here.
- **Path B** — reinstate $\delta(\mathbb{S})$, which adds $\tau^2 W_B W_B^\top$ (diagonal, full
  rank) to $K_B$ and removes the constraint entirely. **Ruled out by the PI**, who requires
  $\delta = 0$.

Test cases run: `r = 8 > Md = 6` passes; `r = 6 = Md` passes (boundary); `r = 3 < Md` fails as
predicted.

**Caveat retained for the record.** The smoke test builds `Phi` from `rnorm()`, which is
full-rank almost surely. Real Gaussian RBFs are far more collinear, and block-averaging makes
it worse. So $r \ge M_D$ is **necessary but not sufficient** — what actually matters is
$\operatorname{rank}(\Psi_B) = M_D$, which is why the wrapper checks it directly.

---

### Step 3 — The R wrapper (`spds/R/gpbau.R` → `gpbau_new`)

`gpbau()` was copied to `gpbau_new()`, differing in only two respects.

**3a. A validation block**, placed immediately after `Md`/`Nd` are known:

- (a) presence of `Phi` and `centroid_distmat`
- (c) conformability: `nrow(Phi) == Nd`, `centroid_distmat` is `r × r`
- (d) `r >= Md` — the necessary condition, with an explicit error message
- (e) `qr(aggmat %*% Phi)$rank == Md` — the *sufficient* condition, which is the one that
  actually catches a collinear basis collapsing under aggregation

Plus a conditioning check after the starting values are set: it computes the condition number
of $K_B$ at the starting $\psi$ and **warns** (not errors) if it exceeds $10^{12}$. A warning
rather than an error because conditioning varies with $\psi$ over the chain.

**Why validate in R rather than C++:** readable messages naming the offending dimensions,
raised before any expensive computation, and before `foreach` scatters work to workers whose
errors come back mangled.

**3b. Class `"gpbau_new"`, not `"gpbau"`.** Deliberate. `predict()` dispatches on class; if the
new objects were tagged `"gpbau"`, `predict()` would route them to `predGPBAU2/3`, which still
build the *old* covariance — returning silently wrong fields rather than erroring. Different
class ⇒ fail loud. (`thinning()` is a plain function, not S3, so it works on either class.)

---

### Step 4 — Starting values (`spds/R/utils.R` → `findinit_new`, `gpbauLoglik_new`)

**What `findinit` does.** MCMC needs a starting point. `findinit` grid-searches for the
approximate maximum-likelihood $(\psi, \sigma^2)$ and hands it to the sampler as `starting`.
The grid is 50 $\psi$ values × 5 $\sigma^2$ values = 250 likelihood evaluations. The $\sigma^2$
grid is built *adaptively per $\psi$* by a moment-matching heuristic: read
$\texttt{diag}(K_B)[1]$ (the block-1 variance at $\sigma^2 = 1$), then bracket
$[\text{var}(X_B)/\text{ratio},\ 3\,\text{var}(X_B)/\text{ratio}]$.

**Why `gpbauLoglik` exists separately.** It is the objective function evaluated at each of the
250 grid points. Importantly it is **not** the same likelihood as the C++: the C++ uses the
*collapsed* likelihood (with $\tilde\gamma, \sigma^2$ integrated out), while `gpbauLoglik` is
the plain multivariate-normal density at a specific $(\psi, \sigma^2, \tilde\gamma)$ triple
with $\tilde\gamma$ plugged in as $\overline{X_B}$. That is adequate for locating a rough mode.

**The change.** Identical grid logic; only the covariance construction was swapped:

```r
# old: bau_covmat  <- sigmasq * exp(-distmat / phi)        # Nd x Nd, rebuilt 50 times
#      coar_covmat <- ww %*% bau_covmat %*% t(ww)
PsiB <- mod_param$aggmat %*% mod_param$Phi                 # hoisted out of the loop
Kpsi <- sigmasq * exp(-data$centroid_distmat / psi)        # r x r
coar_covmat <- PsiB %*% Kpsi %*% t(PsiB)
```

**Why this mattered.** Left alone, chains would have started from values optimised under the
*old* model. Not fatal (MCMC recovers) but it wastes burn-in and invalidates runtime
comparisons. Side benefit: the old version rebuilt a $3000\times3000$ matrix for each of 50
$\psi$ values; the new one builds $120\times120$, so the grid search is dramatically faster.

---

### Step 5 — Wiring the experiment (`data-examples/`)

**5a. Persist the basis at generation time** (`simulate_covariate_data.R`). `Phi` was being
computed but not saved — only its ingredients. Now `Phi` and `centroid_distmat` are saved
directly:

```r
centroid_distmat <- rdist(fine_grid_knots, fine_grid_knots)   # r x r, added
save(fine_grid, coar_grid, range_vals, sill_vals, length_scale, fine_grid_knots,
     knot_spacing, fine_xx, coar_xx, Phi, centroid_distmat,
     file = paste0(outpath, "covar_dat.rdata"))
```

**Why this is the highest-value integrity change:** the model you *fit* now uses the identical
$\Phi$ the data was *generated* from, with no possibility of silent drift from a rebuilt basis.
The lowercase `phi` was also renamed to `Phi` to avoid colliding with the range parameter's
name once loaded into the global environment.

Note `D` is not saved (a pure intermediate), and `fine_grid` is mutated with
`x_cut`/`y_cut`/`par_label` columns before saving — so anything rebuilding $\Phi$ downstream
must subset to `c("lon","lat")`.

**5b. `downscale_gpbau.R`.** `mod_param` now carries the basis inputs; `bau_distmat` was
removed once nothing read it (an $N_D\times N_D$ matrix ≈ 72 MB, copied to every
`multisession` worker). An up-front dimension/rank check runs before the `foreach` loops:

```r
mod_param <- list(aggmat = aggmat, Phi = Phi, centroid_distmat = centroid_distmat)

r_basis  <- ncol(Phi); PsiB_chk <- aggmat %*% Phi
stopifnot(nrow(Phi) == Nd, nrow(centroid_distmat) == r_basis,
          r_basis >= Md, qr(PsiB_chk)$rank == Md)
```

Calls were switched in place: `findinit` → `findinit_new`, `gpbau` → `gpbau_new`. A comment
notes that swapping the single name `gpbau_new` back to `gpbau` reverts to the old model.

---

### Step 6 — Parameter recovery (`sim1-downscaling/downscale_param_recovery.R`, new file)

**Purpose.** Verify the sampler recovers the *known* simulation truth before trusting anything
downstream. For each of the 40 fits it extracts the thinned chains, computes posterior means
and 95% intervals, and checks coverage against $\psi_{\text{true}} = \texttt{range\_vals[i]}$,
$\sigma^2_{\text{true}} = \texttt{sill\_vals[j]}$, $\gamma_0 = 0$. Also reports relative bias,
ESS, acceptance rates, and writes caterpillar/scatter/trace plots.

#### 6a. First run — a bug found

```
coverage:  gamma0 0.850   psi 0.925   sigmasq 0.900
acceptance rate:  Min. 0.0000
non-covering fits:
  1 2  psi   0.1  0.03333333  0.03333333  0.03333333   <- zero-width interval
  7 4  psi   0.7  1.00000000  1.00000000  1.00000000
  8 3  psi   0.8  1.00000000  1.00000000  1.00000000
```

`lower95 == upper95 == post_mean` means **the chain never moved**, and the stuck values are
exactly $u_\psi = 0.1/3$ and $v_\psi = 1$ — the prior bounds.

**Mechanism.** `findinit_new`'s grid `seq(0.1/3, 1, length = 50)` had endpoints *equal to* the
prior bounds. When the profile likelihood is maximised at an endpoint, `starting$phi` lands
exactly on a bound. Then

$$\text{logit}(\psi, u, v) = \log\frac{\psi - u}{v - \psi} = \log 0 = -\infty
\ \Longrightarrow\ \texttt{cur\_logpost} = -\infty$$

and the acceptance test evaluates $-\infty - (-\infty) = \text{NaN}$, which compares `FALSE`
forever. Because $\sigma^2$ and $\tilde\gamma$ are updated **only inside the acceptance block**,
they froze too — producing a chain that looks converged with a zero-width credible interval.

**Three fixes applied:**

| Fix | File | What |
|---|---|---|
| A | `downscale_gpbau.R` | grid made strictly interior: `eps <- 1e-2*diff(prior$phi_unif)`, `phi_vals <- seq(u+eps, v-eps, length=50)`. Also derives from `prior$phi_unif` so the two can't drift apart. |
| B | `gpbau.R` (`gpbau_new`) | clamps `starting$phi` to 1% inside the bounds **and warns**, so boundary-seeking fits are visible rather than masked. Protects all callers. |
| C | `mcmc_sd.cpp` (`fitGPBAU_new`) | `if (!std::isfinite(cur_logpost)) Rcpp::stop(...)` — converts a silent freeze into a loud error. Not an alternative to A/B; a net for unknown triggers. |

A and B are alternatives for the root cause; C is complementary. All three were applied.

#### 6b. Second run — clean

```
coverage:  gamma0 0.925   psi 1.000   sigmasq 0.975
rel. bias: psi +7.8%   sigmasq -7.4%
acceptance: Min. 0.2956, Median 0.4019, Max 0.4688
ESS:       gamma0 522   psi 464   sigmasq 488   (of 500 draws)
non-covering: 1 fit -- (7,3) sigmasq, truth 3, interval [0.706, 2.98]
```

**Interpretation.** With 40 fits the Monte Carlo SE on a true 0.95 coverage is
$\sqrt{0.95 \times 0.05/40} \approx 0.034$, so all three are within ~1 SE. The single miss is
off by 0.02; at 95% nominal you expect ~2 misses in 40. ESS ≈ 464–522 of 500 means the thinned
draws are nearly independent. **The sampler is validated.**

The mirror-image biases ($\psi$ **+7.8%**, $\sigma^2$ **−7.4%**) reflect the classic
range/variance confounding: a longer range with smaller variance produces nearly the same
block-level covariance as a shorter range with larger variance. $\psi$ over-covering at 40/40
is the same story — with only 30 blocks and heavy aggregation smoothing, $\psi$ is weakly
identified and its intervals are wide. This is a property of the design, not a defect.

---

### Step 7 — Prediction (`mcmc_sd.cpp` → `predGPBAU2_new`, `predGPBAU3_new`; `gpbau.R` → `predict.gpbau_new`)

**Purpose in the project.** Fitting only produces parameters; nobody downstream wants a range
parameter. Prediction produces the actual deliverable — $K$ draws of the fine-resolution field
$X_s$. The draws (not a point estimate) are what make the Berkson uncertainty propagation in
Stage 2 possible.

#### 7a. `predGPBAU2_new` — joint, and why it is a restructure

The naive port fails. Under the basis model $C_S = \sigma^2\Phi K_\psi\Phi^\top$ is
$3000\times3000$ with rank only 120, so the conditional covariance is singular and `mvnrnd`
cannot factorise it — and forming it costs ~72 MB per draw.

**The solution: condition in $\alpha$-space** using *conditional simulation by kriging*
(Journel 1974; Chilès & Delfiner §7.3). **This technique is not in your papers** — it is a
standard geostatistics method that produces the same conditional law Algorithm 1 specifies,
with better numerics.

Because $\delta(\mathbb{S}) = 0$, the relation $X_B = R_B\tilde\gamma + \Psi_B\alpha$ holds
*exactly*. Conditioning on $X_B$ is therefore a **linear constraint** on the 120-dimensional
$\alpha$, not a noisy observation. With $V = \sigma^2 K_\psi$ and $c = X_B - R_B\tilde\gamma$:

$$\alpha_{\text{sim}} \sim N(0, V), \qquad
\alpha_{\text{cond}} = \alpha_{\text{sim}} + L\,(c - \Psi_B\alpha_{\text{sim}}), \qquad
L = V\Psi_B^\top(\Psi_B V\Psi_B^\top)^{-1}$$

$$X_s = \tilde R\tilde\gamma + \Phi\,\alpha_{\text{cond}}$$

**Proof it gives the exact target.** Rewrite $\alpha_{\text{cond}} = (I - L\Psi_B)\alpha_{\text{sim}} + Lc$,
an affine map of a Gaussian. Mean $= Lc$. For the covariance, with $S = L\Psi_B V$ (symmetric)
and using $L\Psi_B V\Psi_B^\top L^\top = S$:

$$\operatorname{Var} = (I - L\Psi_B)V(I - \Psi_B^\top L^\top) = V - S - S + S = V - S$$

which is exactly $V - L\Psi_B V$, the standard Gaussian-conditioning covariance. $\blacksquare$

**Why it dodges the rank problem.** $L\Psi_B$ is idempotent, so $A = I - L\Psi_B$ projects onto
$\ker(\Psi_B)$ of dimension $r - M_D = 90$. Hence $\operatorname{rank}(V - S) \le 90 < 120$ —
**the target covariance really is singular**, because knowing $X_B$ pins down 30 of the 120
directions exactly. But the trick never forms $V - S$. It factorises only $V = \sigma^2 K_\psi$
(full rank $r$) and solves with $\Psi_B V\Psi_B^\top = \sigma^2 K_B$ (full rank $M_D$ when
$r \ge M_D$). The singular distribution is realised *implicitly* as the image of a full-rank
draw under a projection.

*Analogy:* to sample on a line in $\mathbb{R}^2$, don't Cholesky a singular $2\times2$
covariance — draw a scalar and map it onto the line.

In the code, $\sigma^2$ cancels: $L = \sigma^2 K_\psi\Psi_B^\top(\sigma^2 K_B)^{-1} = K_\psi\Psi_B^\top K_B^{-1}$.

**Built-in correctness check.** Since $\Psi_B L = I$, every joint draw satisfies
$W_B X_s = X_B$ **exactly**. `downscale_gpbau.R` now prints this residual.
**Measured: `4.330e-14`** — machine precision (double eps ≈ 2.2e-16). This is the single
strongest evidence the conditioning is correct.

> **Historical note.** The original `predGPBAU2` was already mildly degenerate:
> $\operatorname{rank}(\Sigma_{1|2}) = N_D - M_D = 2970$ of 3000. It survived because
> round-off on 30 near-zero eigenvalues usually let Cholesky through. The basis model drops
> that to 90 of 3000. Same pathology, three orders of magnitude worse.

#### 7b. `predGPBAU3_new` — independent, a lighter edit

The per-location loop survives; only the per-row quantities change, computed without ever
forming $C_S$:

$$C_S[i,i] = \sigma^2\phi_i^\top K_\psi\phi_i, \qquad
C_S[i,\cdot]W_B^\top = \sigma^2\phi_i^\top K_\psi\Psi_B^\top$$

**The one subtle trap.** The old code wrote `cond_var = sigmasq - ...` because `expCorr` has 1
on its diagonal. Under the basis model the marginal variance is
$\sigma^2\phi_i^\top K_\psi\phi_i$, which is **not 1** and **varies by location** depending on
how much basis support sits at $s_i$. Using `sigmasq` would give uniform, wrong predictive
variances — an error visible only as mis-calibrated coverage. `Phi * Kpsi` is hoisted out of
the location loop; a `cond_var < 0` guard handles round-off.

#### 7c. `predict.gpbau_new` (`gpbau.R`)

Needed for two reasons, neither a bare passthrough:

1. **S3 dispatch.** `predict()` looks up `predict.<class>`; without this method, `"gpbau_new"`
   objects error.
2. **Adapter work.** Chooses `post_sams` vs `post_sams_thin`; **transposes `beta`** in the
   unthinned branch (the C++ stores `gam_save` as $p\times$niter but `pred*` reads
   `gam_save.row(iter)`, expecting niter$\times p$; `thinning()` already transposes);
   defaults `nreport`; switches on `type`.

**Current setting:** `downscale_gpbau.R` uses the default `type = "independent"`.
See §4 for why switching to `"joint"` is recommended.

---

### Step 8 — Downstream (`downscale_comp.R`, `sim2-sdm/`)

**Design decision that made this nearly free:** the output contract was preserved —
`xx_pred[[i]][[j]]` is still an $N_D \times K$ matrix saved to `sim_gpbau_pred.RData`. All
three consumers (`downscale_comp.R:14`, `pois_replicate_fit.R:24`, `pois_eval.R:19`) index it
identically, so no edits were required.

**What each downstream script is for:**

| Script | Role |
|---|---|
| `downscale_comp.R` | Section 4.1: field-level coverage (`gpbau_pred_ci`) vs true `fine_xx`; compares GPB against thin-plate splines |
| `sim2-sdm/simulate_pois_data.R` | generates Poisson count responses from the *true* field; `nrep = 100` replicates |
| `sim2-sdm/pois_replicate_fit.R` | **consumes the predictions**; fits GLM-oracle / GLM-plugin / GLM-Berkson |
| `sim2-sdm/pois_eval.R` | scores those fits (bias, coverage, RMSE) |
| `sim2-sdm/comp_figures.R` | renders Section 4.2 figures |

The scientific point is that **Berkson beats plugin** — propagating downscaling uncertainty is
better calibrated than treating the downscaled map as exact.

**Fixes required:** `pois_replicate_fit.R:33` and `pois_eval.R:23` contained the placeholder
`source("directory/to/source_functions.R")`. Corrected to `source("sim2-sdm/source_functions.R")`
— **not** `"source_functions.R"`, because `source()` resolves relative to the working
directory (`data-examples/`), not the script's own location.

**Runtime warning.** `save()` is *after* the replication loop in both scripts, so interrupting
partway writes **nothing**. For plumbing verification, reduce `nrep` at the source
(`simulate_pois_data.R:29`) — both downstream scripts derive `nrep <- length(rep_data)`, so one
edit propagates consistently.

---

### Housekeeping

| Change | File | Why |
|---|---|---|
| Removed `bau_distmat` from C++ and `mod_param` | `mcmc_sd.cpp`, `gpbau.R`, `downscale_gpbau.R` | 72 MB × every parallel worker, for a matrix nothing read |
| Added recovery script to pipeline | `run_all_rscripts.R` | validation becomes part of the reproducible run |
| Documented the variant and `r ≥ M_D` | `spds/DESCRIPTION` | the constraint is a modelling requirement, not a tuning note |
| Modelling constraint comment | `simulate_covariate_data.R` | only place a reader sees both $r$ and $M_D$ |

---

## 3. Complete file manifest

| File | Status | What |
|---|---|---|
| `data-examples/simulate_covariate_data.R` | modified | basis-function generation; saves `Phi`, `centroid_distmat`; `phi`→`Phi` |
| `spds/src/mcmc_sd.cpp` | modified | `fitGPBAU_new`, `predGPBAU2_new`, `predGPBAU3_new` added; originals untouched |
| `spds/R/gpbau.R` | modified | `gpbau_new`, `predict.gpbau_new` added; originals untouched |
| `spds/R/utils.R` | modified | `findinit_new`, `gpbauLoglik_new` added; originals untouched |
| `data-examples/sim1-downscaling/downscale_gpbau.R` | modified | basis `mod_param`, rank check, `*_new` calls, interior grid, aggregation check |
| `data-examples/sim1-downscaling/downscale_param_recovery.R` | **new** | Step D validation |
| `sandbox/smoke_test_fitgpbau_new.R` | **new** | rank/PD smoke test |
| `data-examples/run_all_rscripts.R` | modified | added recovery script |
| `data-examples/sim2-sdm/pois_replicate_fit.R` | modified | `source()` path fix |
| `data-examples/sim2-sdm/pois_eval.R` | modified | `source()` path fix |
| `spds/DESCRIPTION` | modified | documents variant + constraint |

**Nothing in the original code path was deleted.** `fitGPBAU`, `predGPBAU2`, `predGPBAU3`,
`gpbau`, `predict.gpbau`, `findinit`, `gpbauLoglik` are all intact, so the old model remains
runnable for comparison.

---

## 4. Open items

1. **`type = "joint"` in sim1.** Currently `"independent"`. Joint and independent have
   *identical marginals*, so pointwise coverage is unaffected — but only joint satisfies the
   aggregation invariant, and the Berkson propagation in `sim2-sdm` fits a GLM across many
   locations at once, where spatially independent covariate error averages out far faster than
   correlated error and would **understate** Stage-2 uncertainty. The original likely chose
   independent because the old joint predictor was $O(N_D^3)$; the $\alpha$-space version is not.
   **Recommended: switch.**

2. **The framing issue (most important for the writeup).** `Phi` in
   `simulate_covariate_data.R` is *literally the same matrix* passed to `mod_param`. The true
   field is exactly $\Phi\alpha$, so there is **zero model misspecification** — you generated
   the answer key with the same formula you're grading with. Real fields are not built from 120
   Gaussian RBFs; in practice the basis is an *approximation*. A reviewer will notice.
   **Robustness check:** generate from the original `rgp` over the full fine grid (still
   commented out in the script) and fit with the basis model. Keep the matched case as the
   "sampler is correct" check; add the mismatched case as the "method is useful" check.

3. **Full `sim2-sdm` run** at `nrep = 100`. Also worth checking whether
   `pois_replicate_fit.R`/`pois_eval.R` set `plan(multisession)` — they use `%dofuture%` but
   may be running the 40-cell grid sequentially.

4. **Weak identification of $\psi$** at large ranges. Documented, not a defect; worth
   reporting rather than tuning away.

---

## 5. Rebuild and run

```r
# from corgi-examples/
Rcpp::compileAttributes("spds")   # REQUIRED whenever C++ exports change
devtools::document("spds")        # regenerates NAMESPACE from @export tags
devtools::install("spds")         # workers need library(spds), not load_all
# restart R

setwd("~/Documents/GitHub/corgi-examples/data-examples")
source("simulate_covariate_data.R")
source("sim1-downscaling/downscale_gpbau.R")          # prints aggregation check
source("sim1-downscaling/downscale_param_recovery.R") # prints coverage
source("sim1-downscaling/downscale_comp.R")
source("sim2-sdm/simulate_pois_data.R")
source("sim2-sdm/pois_replicate_fit.R")
source("sim2-sdm/pois_eval.R")
source("sim2-sdm/comp_figures.R")
```

**Gotchas learned the hard way:**

- **Working directory must be `data-examples/`.** All scripts use `../outputs/`; running from
  the repo root silently writes outside the repo (`dir.create` suppresses warnings).
- **`Rcpp::compileAttributes` is not a C++ parser.** Nothing may sit between the export
  attribute and the function signature — it scans that span *as* the signature and splices
  comment text into the generated wrapper. Multi-line comments inside the parameter list break
  it too. Keep notes above the attribute.
- **Disable "Restore .RData at startup"** in RStudio. Stale globals can mask a failed `load()`
  and let checks pass on leftover objects.

**Health indicators:**

| Check | Expected | Last measured |
|---|---|---|
| `rank(W_B Phi)` | $= M_D = 30$ | 30 |
| aggregation error (joint) | $\approx 0$ | **4.33e-14** |
| MH acceptance | 0.20–0.45 | 0.30–0.47 |
| coverage $\psi$ / $\sigma^2$ / $\gamma_0$ | $\approx 0.95$ | 1.000 / 0.975 / 0.925 |
| ESS (of 500 thinned) | high | 464–522 |
| field coverage `gpbau_pred_ci` | $\approx 0.95$ | *to verify* |
