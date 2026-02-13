from numpy.linalg import solve, LinAlgError
import pandas as pd
import numpy as np


# -----------------------------
# 1) Entropy balancing MAIC core
# -----------------------------
class EntropyBalancer:
    """
    Entropy balancing via exponential tilting:
        w_i ∝ q_i * exp(A_i^T λ)
    subject to E_w[A] = b and sum(w)=1.

    A: (n x m) constraint matrix (each column is a moment/indicator to match)
    b: (m,) target moments (proportions/means)
    q: baseline weights (defaults uniform)
    """
    def __init__(self, A, b, q=None, max_iter=200, tol=1e-10, line_search=True):
        A = np.asarray(A, dtype=float)
        b = np.asarray(b, dtype=float).ravel()
        assert A.ndim == 2
        assert b.ndim == 1 and A.shape[1] == b.shape[0], "A (n×m) and b (m,) mismatch"

        self.A = A
        self.b = b
        self.n, self.m = A.shape
        self.q = np.ones(self.n) / self.n if q is None else np.asarray(q, dtype=float) / np.sum(q)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.line_search = bool(line_search)

        self.lam_ = np.zeros(self.m)
        self.w_ = np.ones(self.n) / self.n
        self.converged_ = False
        self.n_iter_ = 0

    def _weights_from_lambda(self, lam):
        s = self.A @ lam
        ex = np.exp(s) * self.q
        w = ex / ex.sum()
        return w

    def _moment_and_hessian(self, lam):
        w = self._weights_from_lambda(lam)
        m = self.A.T @ w  # E_w[A]
        C = self.A - m    # center rows by current moments
        H = (C * w[:, None]).T @ C  # Cov_w(A)
        return w, m, H


    def fit(self, lam0=None):
        lam = np.zeros(self.m) if lam0 is None else np.asarray(lam0, dtype=float)
        for it in range(1, self.max_iter + 1):
            w, m, H = self._moment_and_hessian(lam)
            g = m - self.b
            if np.linalg.norm(g, ord=np.inf) < self.tol:
                self.converged_ = True
                self.lam_, self.w_, self.n_iter_ = lam, w, it
                return self

            # Newton step: H step = g
            try:
                step = solve(H, g)
            except LinAlgError:
                step = solve(H + 1e-10 * np.eye(H.shape[0]), g)

            # Backtracking line search for stability
            if self.line_search:
                norm_g = np.linalg.norm(g, ord=np.inf)
                step_scale = 1.0
                for _ in range(25):
                    lam_try = lam - step_scale * step
                    _, m_try, _ = self._moment_and_hessian(lam_try)
                    if np.linalg.norm(m_try - self.b, ord=np.inf) < norm_g:
                        lam = lam_try
                        break
                    step_scale *= 0.5
                else:
                    lam = lam - 0.1 * step
            else:
                lam = lam - step

        self.converged_ = False
        self.lam_, self.w_, self.n_iter_ = lam, self._weights_from_lambda(lam), self.max_iter
        return self

    @staticmethod
    def ess(w):
        w = np.asarray(w, dtype=float)
        return 1.0 / np.sum(w**2)

    def bootstrap(
        self,
        B,
        estimand_fn=None,
        *,
        method="multiplier",
        seed=None,
        warm_start=True,
        max_iter=None,
        tol=None,
        line_search=None,
    ):
        """
        Returns a DataFrame with full bootstrap distributions for:
          - ESS
          - optional estimand(s) via estimand_fn(w, idx=None|idx)
          - convergence + max imbalance diagnostics
        """
        rng = np.random.default_rng(seed)
        A0, b0 = self.A, self.b
        n, m = A0.shape

        lam0 = getattr(self, "lam_", np.zeros(m))
        if not warm_start:
            lam0 = np.zeros(m)

        max_iter_eff = self.max_iter if max_iter is None else int(max_iter)
        tol_eff = self.tol if tol is None else float(tol)
        line_search_eff = self.line_search if line_search is None else bool(line_search)

        rows = []
        for draw in range(B):
            if method == "multiplier":
                # Bayesian / multiplier bootstrap
                q_b = rng.exponential(1.0, size=n)
                q_b /= q_b.sum()
                A_b = A0
                idx = None
            elif method == "resample":
                idx = rng.integers(0, n, size=n)
                A_b = A0[idx, :]
                q_b = np.ones(n) / n
            else:
                raise ValueError("method must be 'multiplier' or 'resample'")

            eb_b = EntropyBalancer(A_b, b0, q=q_b, max_iter=max_iter_eff, tol=tol_eff, line_search=line_search_eff)
            eb_b.fit(lam0=lam0)

            w_b = eb_b.w_
            max_imb = float(np.max(np.abs(A_b.T @ w_b - b0)))

            out = {
                "draw": draw,
                "method": method,
                "converged": bool(eb_b.converged_),
                "n_iter": int(eb_b.n_iter_),
                "max_imbalance": max_imb,
                "ess": float(self.ess(w_b)),
            }

            if estimand_fn is not None:
                val = estimand_fn(w_b, idx=idx)
                if isinstance(val, dict):
                    out.update(val)
                else:
                    out["estimand"] = float(val)

            rows.append(out)

            # keep warm start close
            if warm_start and eb_b.converged_:
                lam0 = eb_b.lam_

        return pd.DataFrame(rows)

# ----------------------------------------------------
# 2) Helper: build constraints for marginal frequencies
# ----------------------------------------------------
def build_marginal_frequency_constraints(df, specs, target_totals=None):
    """
    specs[col] describes how to balance each covariate.

    For categorical/binary: create indicator constraints for all-but-one level (reference dropped).
      Constraint column: 1(x==level)
      Target b: target proportion for that level
    For continuous: include raw column as constraint (target is a mean).

    target_totals (optional): {col: N_target} if targets provided as COUNTS.
    """
    A_cols, b_vals, names = [], [], []

    for col, meta in specs.items():
        typ = meta["type"]

        if typ == "continuous":
            A_cols.append(df[col].to_numpy(dtype=float))
            b_vals.append(float(meta["target"]))
            names.append(f"mean({col})")
            continue

        if typ not in ("categorical", "binary"):
            raise ValueError(f"Unknown type={typ} for {col}")

        x = df[col]
        levels = list(meta.get("levels") or pd.Series(x).dropna().unique())
        levels.sort(key=lambda v: str(v))
        drop = meta.get("drop", levels[-1])
        keep_levels = [lv for lv in levels if lv != drop]

        tgt = meta["target"]  # dict level -> count or proportion
        vals = np.array([tgt[lv] for lv in levels], dtype=float)
        s = vals.sum()

        # convert counts -> proportions if target_totals provided (or normalize)
        if target_totals and col in target_totals:
            Nstar = float(target_totals[col])
            if np.isclose(s, Nstar):
                props = vals / Nstar
            else:
                props = vals / s
        else:
            props = vals / s  # if already proportions summing to 1, this is unchanged

        prop_map = {lv: p for lv, p in zip(levels, props)}

        for lv in keep_levels:
            ind = (x == lv).to_numpy(dtype=float)
            A_cols.append(ind)
            b_vals.append(float(prop_map[lv]))
            names.append(f"Pr({col}={lv})")

    A = np.column_stack(A_cols).astype(float)
    b = np.array(b_vals, dtype=float)
    return A, b, names
