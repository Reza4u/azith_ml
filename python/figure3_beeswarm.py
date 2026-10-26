"""Figure 3: SHAP beeswarm for the final Random Forest.

Reads the exact Shapley values written by R/08_final_model_shap.R
(output/private/shap_exact_S.csv and shap_exact_X.csv; participant-level, never
share them) and writes output/figures/Figure_3.png and .tiff.

Run from the repository root:  python3 python/figure3_beeswarm.py
Needs numpy, pandas and matplotlib.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager

OUT = os.environ.get("AZITH_OUT_DIR", "output")
S = pd.read_csv(os.path.join(OUT, "private", "shap_exact_S.csv"))
X = pd.read_csv(os.path.join(OUT, "private", "shap_exact_X.csv"), keep_default_na=False)

lab = {"macro_numeric": "Faecal macrophages (2–5 per HPF)", "wbc_numeric": "Stool white-cell count",
       "resident_area": "Residence (urban)", "abx_cat": "Antibiotic before admission†",
       "age": "Age (months)", "diarrhea_dur": "Diarrhoea duration (days)",
       "straining": "Straining/tenesmus", "mucus_RE": "Mucus in stool",
       "blood_RE": "Blood in stool", "pcell_numeric": "Stool pus-cell count"}
# category order used for colouring (low -> high)
order_levels = {"resident_area": ["Rural Area", "Urban Area"], "straining": ["No", "Yes"],
                "mucus_RE": ["Small", "Moderate", "Plenty"],
                "blood_RE": ["Absent", "Small", "Moderate", "Plenty"],
                "abx_cat": ["None", "Ciprofloxacin", "Metronidazole", "Cipro_Metro", "Other"]}
# top -> bottom, as in Supplementary Table S10 (blood and mucus tie at 0.024)
feats = ["macro_numeric", "wbc_numeric", "resident_area", "abx_cat", "age",
         "diarrhea_dur", "straining", "mucus_RE", "blood_RE", "pcell_numeric"]
assert set(feats) == set(S.columns), (feats, list(S.columns))
feats = feats[::-1]


def colour_value(col):
    v = X[col]
    if col in order_levels:
        lv = order_levels[col]
        t = v.astype(str).str.strip().str.lower().map({l.lower(): i for i, l in enumerate(lv)})
        assert t.notna().all(), (col, v.unique())
        v = t.astype(float)
    v = v.astype(float)
    lo, hi = np.nanpercentile(v, 5), np.nanpercentile(v, 95)
    if hi <= lo:
        lo, hi = v.min(), v.max()
    return np.clip((v - lo) / (hi - lo) if hi > lo else v * 0 + 0.5, 0, 1)


for fam in ["Liberation Sans", "Arial", "DejaVu Sans"]:
    if any(fam in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = fam
        break
cmap = LinearSegmentedColormap.from_list("bp", ["#1A85FF", "#D41159"])
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
row_h = 0.4
for yi, col in enumerate(feats):
    s = S[col].to_numpy()
    c = colour_value(col).to_numpy()
    # beeswarm offsets (same idea as shap.summary_plot)
    nbins = 100
    quant = np.round(nbins * (s - s.min()) / (s.max() - s.min() + 1e-12))
    inds = np.argsort(quant + np.random.RandomState(0).randn(len(s)) * 1e-6)
    layer, last_bin, ys = 0, -1, np.zeros(len(s))
    for ind in inds:
        if quant[ind] != last_bin:
            layer = 0
        ys[ind] = np.ceil(layer / 2) * ((layer % 2) * 2 - 1)
        layer += 1
        last_bin = quant[ind]
    ys *= 0.9 * (row_h / (np.max(ys + 1) if np.max(ys) > 0 else 1))
    order = np.random.RandomState(1).permutation(len(s))
    ax.scatter(s[order], yi + ys[order], c=c[order], cmap=cmap, vmin=0, vmax=1, s=9, alpha=0.9,
               linewidths=0, rasterized=False, zorder=3)
ax.axvline(0, color="#666666", lw=0.8, zorder=2)
ax.set_yticks(range(len(feats)))
ax.set_yticklabels([lab[f] for f in feats], fontsize=10.5)
ax.set_ylim(-0.6, len(feats) - 0.4)
for yi in range(len(feats)):
    ax.axhline(yi, color="#E5E5E5", lw=0.6, ls=(0, (1, 2)), zorder=1)
ax.set_xlabel("SHAP value (change in predicted probability of failure)", fontsize=10.5,
              fontweight="bold")
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_linewidth(0.8)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", labelsize=9.5)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, ticks=[0, 1], aspect=30, fraction=0.03, pad=0.02)
cb.ax.set_yticklabels(["Low", "High"], fontsize=9.5)
cb.outline.set_visible(False)
cb.set_label("Feature value", rotation=270, labelpad=10, fontsize=10.5, fontweight="bold")
fig.tight_layout()
os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
fig.savefig(os.path.join(OUT, "figures", "Figure_3.png"), dpi=300, facecolor="white")
fig.savefig(os.path.join(OUT, "figures", "Figure_3.tiff"), dpi=300, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Figure 3 written; mean |SHAP|:")
print(S.abs().mean().sort_values(ascending=False).round(3).to_string())
