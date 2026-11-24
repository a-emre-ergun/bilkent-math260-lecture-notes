import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
from scipy import stats

plt.rcParams.update({
    "font.size": 10, "mathtext.fontset": "cm", "font.family": "serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "path",
})
INK = "#222222"; FILL = "#9ecae1"; FILL2 = "#fdae6b"; LINE2 = "#d95f02"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"); os.makedirs(OUT, exist_ok=True)

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)

def bare(ax, xlabel=None):
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])
    ax.set_ylim(bottom=0)
    if xlabel: ax.set_xlabel(xlabel)

x = np.linspace(-3.6, 3.6, 600)
phi = stats.norm.pdf(x)

def normal_panel(ax, regions, ticks, labels, notes=()):
    ax.plot(x, phi, color=INK, lw=1.4)
    for lo, hi in regions:
        m = (x >= lo) & (x <= hi)
        ax.fill_between(x[m], phi[m], color=FILL, ec="none")
    ax.set_xticks(ticks); ax.set_xticklabels(labels)
    for (txt, xy, xytext) in notes:
        ax.annotate(txt, xy=xy, xytext=xytext, ha="center",
                    arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
    bare(ax)
    ax.set_ylim(0, 0.48)

# 1-3. Rejection regions for the z test --------------------------------------
za, za2 = 1.645, 1.96
for name, regions, ticks, labels, notes in [
    ("rr_upper", [(za, 9)], [0, za], ["$0$", r"$z_\alpha$"],
     [(r"area $\alpha$" "\n(rejection region)", (2.15, 0.02), (2.7, 0.22))]),
    ("rr_lower", [(-9, -za)], [-za, 0], [r"$-z_\alpha$", "$0$"],
     [(r"area $\alpha$" "\n(rejection region)", (-2.15, 0.02), (-2.7, 0.22))]),
    ("rr_two", [(-9, -za2), (za2, 9)], [-za2, 0, za2],
     [r"$-z_{\alpha/2}$", "$0$", r"$z_{\alpha/2}$"],
     [(r"area $\alpha/2$", (-2.35, 0.012), (-2.9, 0.18)),
      (r"area $\alpha/2$", (2.35, 0.012), (2.9, 0.18))]),
]:
    fig, ax = plt.subplots(figsize=(3.6, 1.7))
    normal_panel(ax, regions, ticks, labels, notes)
    save(fig, name)

# 4. Power of the lower-tailed t test: central vs noncentral t ---------------
nu, delta = 12, -2.2
tt = np.linspace(-6.5, 4, 700)
f0 = stats.t.pdf(tt, nu); f1 = stats.nct.pdf(tt, nu, delta)
tc = -stats.t.ppf(0.95, nu)
fig, ax = plt.subplots(figsize=(4.6, 2.0))
m = tt <= tc
ax.fill_between(tt[m], f1[m], color=FILL2, alpha=0.7, ec="none", label=r"power $=1-\beta(\mu')$")
ax.fill_between(tt[m], f0[m], color=FILL, ec="none", label=r"$\alpha$")
ax.plot(tt, f0, color=INK, lw=1.4)
ax.plot(tt, f1, color=LINE2, lw=1.4, ls="--")
ax.text(1.5, 0.30, "$t_{n-1}$\n($H_0$ true)", ha="left", fontsize=9)
ax.text(-4.4, 0.33, "noncentral $t$\n" r"($\mu=\mu'$)", ha="center", fontsize=9, color=LINE2)
ax.axvline(tc, color=INK, lw=0.7, ls=":")
ax.set_xticks([tc, 0]); ax.set_xticklabels([r"$-t_{\alpha,n-1}$", "$0$"])
ax.legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(1.18, 0.3))
bare(ax); ax.set_ylim(0, 0.45)
save(fig, "noncentral_t_power")

# 5. P-value vs rejection region, alpha = .01 --------------------------------
fig, axes = plt.subplots(1, 2, figsize=(6.4, 1.9))
for ax, z, title in [(axes[0], 3.0, r"Case I: $z\geq z_{.01}=2.33$"),
                     (axes[1], 1.4, r"Case II: $z<z_{.01}=2.33$")]:
    ax.plot(x, phi, color=INK, lw=1.4)
    m = x >= z
    ax.fill_between(x[m], phi[m], color=FILL, ec="none")
    ax.axvline(2.33, color=LINE2, lw=0.9, ls="--")
    ax.set_xticks([0, 2.33, z]); ax.set_xticklabels(["$0$", "$2.33$" if z < 2.33 else "$2.33\\;$", "$\\;z$" if z > 2.33 else "$z\\;$"])
    ax.annotate(r"$P$-value" + (r" $\leq .01$" if z > 2.33 else r" $>.01$"),
                xy=(z + 0.35, stats.norm.pdf(z + 0.35) / 2), xytext=(1.2, 0.30) if z > 2.33 else (2.45, 0.30), ha="left",
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_title(title, fontsize=9); bare(ax); ax.set_ylim(0, 0.45)
fig.tight_layout(); save(fig, "pvalue_cases")

# 6. P-value for the z test ---------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(7.2, 1.9))
specs = [
    (axes[0], [(1.2, 9)], [0, 1.2], ["$0$", "$z$"], r"upper-tailed: $1-\Phi(z)$"),
    (axes[1], [(-9, -1.2)], [-1.2, 0], ["$z$", "$0$"], r"lower-tailed: $\Phi(z)$"),
    (axes[2], [(-9, -1.2), (1.2, 9)], [-1.2, 0, 1.2], [r"$-|z|$", "$0$", r"$|z|$"],
     r"two-tailed: $2[1-\Phi(|z|)]$"),
]
for ax, reg, ticks, labels, title in specs:
    normal_panel(ax, reg, ticks, labels)
    ax.set_title(title, fontsize=9)
fig.tight_layout(); save(fig, "pvalue_z")

# 7. Two ways a contingency table arises -------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6))
ax = axes[0]
for i, cx in enumerate([1, 2.6, 4.2]):
    ax.add_patch(Ellipse((cx, 1.6), 1.1, 0.8, fc="#deebf7", ec=INK))
    ax.text(cx, 1.6, f"Pop {i+1}", ha="center", va="center")
    ax.annotate("", xy=(cx, 0.75), xytext=(cx, 1.18), arrowprops=dict(arrowstyle="->"))
    ax.add_patch(Ellipse((cx, 0.5), 0.6, 0.38, fc="white", ec=INK))
    ax.text(cx, 0.5, "sample", ha="center", va="center", fontsize=7)
ax.text(5.0, 1.6, r"$\cdots$", ha="center", va="center")
ax.text(2.6, 2.3, r"$I$ populations, $J$ categories in each", ha="center")
ax.text(2.6, -0.05, "stratified sampling (one sample per population)", ha="center", fontsize=8)
ax.set_title("1. Multiple populations, one classification", fontsize=9)
ax = axes[1]
ax.add_patch(Ellipse((2.6, 1.6), 3.2, 1.0, fc="#deebf7", ec=INK))
ax.text(2.6, 1.6, "one population", ha="center", va="center")
ax.annotate("", xy=(2.6, 0.75), xytext=(2.6, 1.1), arrowprops=dict(arrowstyle="->"))
ax.add_patch(Ellipse((2.6, 0.5), 1.4, 0.4, fc="white", ec=INK))
ax.text(2.6, 0.5, "one sample", ha="center", va="center", fontsize=7)
ax.text(0.9, 2.3, "factor 1: $I$ categories", ha="center", fontsize=9)
ax.text(4.3, 2.3, "factor 2: $J$ categories", ha="center", fontsize=9)
ax.text(2.6, -0.05, r"simple random sampling; each individual falls in a cell $(i,j)$," "\n"
        "e.g. (female, blue), (male, yellow)", ha="center", fontsize=8, va="top")
ax.set_title("2. One population, two classifications", fontsize=9)
for ax in axes:
    ax.set_xlim(0, 5.3); ax.set_ylim(-0.45, 2.55); ax.axis("off")
fig.tight_layout(); save(fig, "contingency_designs")

# 8. True regression line with deviations ------------------------------------
rng = np.random.default_rng(3)
b0, b1, s = 1.0, 0.8, 0.7
xs = np.sort(rng.uniform(0.5, 6, 14)); ys = b0 + b1 * xs + rng.normal(0, s, xs.size)
fig, ax = plt.subplots(figsize=(4.4, 2.6))
xx = np.linspace(0, 6.5, 2)
ax.plot(xx, b0 + b1 * xx, color=INK, lw=1.4)
ax.scatter(xs, ys, s=14, color=LINE2, zorder=3)
for k in [3, 7, 10]:
    ax.plot([xs[k]] * 2, [b0 + b1 * xs[k], ys[k]], color="#555", lw=0.8)
    ax.text(xs[k] + 0.08, (b0 + b1 * xs[k] + ys[k]) / 2, r"$\varepsilon_i$", fontsize=9, va="center")
ax.text(4.6, b0 + b1 * 6.4 + 0.1, "true regression line\n" r"$E(Y)=\beta_0+\beta_1x$", fontsize=9, ha="center")
ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlim(0, 6.8); ax.set_ylim(0, 8.3)
save(fig, "true_regression_line")

# 9. Conditional distribution of Y at x = x1 ---------------------------------
fig, ax = plt.subplots(figsize=(4.4, 2.8))
ax.plot(xx, b0 + b1 * xx, color=INK, lw=1.4)
for xk, lab in [(1.8, "$x_1$"), (4.5, "$x_2$")]:
    mu = b0 + b1 * xk
    yy = np.linspace(mu - 3 * s, mu + 3 * s, 200)
    ax.plot(xk + 1.6 * stats.norm.pdf(yy, mu, s) * s, yy, color=LINE2, lw=1.2)
    ax.plot([xk, xk], [0, mu + 3 * s], color="#888", lw=0.6)
    ax.plot([0, xk], [mu, mu], color="#888", lw=0.6, ls="--")
    ax.scatter([xk], [mu], s=12, color=INK, zorder=3)
mu1 = b0 + b1 * 1.8
ax.set_xticks([1.8, 4.5]); ax.set_xticklabels(["$x_1$", "$x_2$"])
ax.set_yticks([mu1]); ax.set_yticklabels([r"$\beta_0+\beta_1x_1$"])
ax.annotate(r"$N(\beta_0+\beta_1x_1,\ \sigma^2)$", xy=(2.45, mu1), xytext=(2.9, 1.0),
            arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
ax.text(5.9, b0 + b1 * 6.3, r"$E(Y\mid x)=\beta_0+\beta_1x$", fontsize=9, ha="right")
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.set_xlim(0, 6.6); ax.set_ylim(0, 7.6)
save(fig, "conditional_normal")

# 10. Candidate lines through the data ---------------------------------------
fig, ax = plt.subplots(figsize=(4.4, 2.6))
ax.scatter(xs, ys, s=14, color=LINE2, zorder=3)
bb1, bb0 = np.polyfit(xs, ys, 1)
ax.plot(xx, bb0 + bb1 * xx, color=INK, lw=1.4)
ax.plot(xx, 0.6 + 0.25 * xx, color="#888", lw=1.2, ls="--")
ax.text(6.55, bb0 + bb1 * 6.5, "$y=b_0+b_1x$", va="center", fontsize=9)
ax.text(6.55, 0.6 + 0.25 * 6.5, "$y=a_0+a_1x$", va="center", fontsize=9, color="#555")
ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlim(0, 8); ax.set_ylim(0, 8.3)
save(fig, "candidate_lines")

# 11. Three scatter plots for R^2 --------------------------------------------
rng = np.random.default_rng(7)
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.1), sharey=True)
xr = np.linspace(0.5, 5.5, 12)
for ax, noise, title in [(axes[0], 0.0, "SSE $=0$"), (axes[1], 0.35, "SSE small"),
                         (axes[2], 1.6, "SSE large")]:
    yr = 1 + 0.9 * xr + (rng.normal(0, noise, xr.size) if noise else 0)
    if noise > 1: yr = 4 + rng.normal(0, noise, xr.size)
    ax.scatter(xr, yr, s=12, color=LINE2)
    p = np.polyfit(xr, yr, 1); ax.plot([0, 6], np.polyval(p, [0, 6]), color=INK, lw=1)
    ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("$x$")
axes[0].set_ylabel("$y$")
fig.tight_layout(); save(fig, "r2_scatterplots")
print("ok", sorted(os.listdir(OUT)))
