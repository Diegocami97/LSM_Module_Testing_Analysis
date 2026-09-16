#!/usr/bin/env python3
# Diego Venegas-Vargas
# DAMIC-M Collaboration
# Combines ROOT cluster trees by grid point/extension, computes front/back
# double-Gaussian Fe-55 energy peaks, and outputs CSV + comparison plots.
"""
Parameter scan energy peaks: read ROOT cluster trees, combine files by grid point
and extension, compute front/back double-Gaussian energy peaks (like MakeComposite_Fe55_Clusters_Extensions.cc)
and output CSV + comparison plots (like parameter_scan_comparison.py).

File naming convention: ..._<image_id>_ext<ext_num>.root
  - Images with same base (grid point) and same ext, e.g. ..._0_ext1.root and ..._1_ext1.root,
    are combined before fitting.
"""

import re
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Optional: uproot for reading ROOT files
try:
    import uproot
except ImportError:
    uproot = None

# ------------------------- Config (match .cc) -------------------------
ENE_MIN = 4.5          # keV
ENE_MAX = 7.0          # keV
SIGMA_FRONT_CUT = 1.0  # pixels; front = sigma_xy < this
G1_LO, G1_HI = 5.3, 6.1   # first peak (Kα) fit range [keV]
NBINS = 50
E_KA_NOMINAL = 5.898   # Mn Kα
E_KB_NOMINAL = 6.490   # Mn Kβ


def parse_file_list(filelist_path, pattern_two_images=True):
    """
    Parse list of ROOT paths. Each line: path to ..._<image_id>_ext<ext>.root
    Returns: dict (grid_point, ext_num) -> list of full paths to combine.
    grid_point = basename minus _<image_id>_ext<ext>.root (so 0 and 1 share same key).
    """
    path_by_key = defaultdict(list)
    # Basename pattern: anything, then _<digits>_ext<1-4>.root
    pat = re.compile(r"^(.+)_(\d+)_ext([1-4])\.root$", re.IGNORECASE)
    with open(filelist_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line)
            if not p.suffix.lower() == ".root":
                continue
            base = p.name
            m = pat.match(base)
            if not m:
                continue
            grid_point, image_id, ext_num = m.group(1), m.group(2), int(m.group(3))
            path_by_key[(grid_point, ext_num)].append(str(p))
    return dict(path_by_key), list(path_by_key.keys())


def s_cut_variable(rmsx, rmsy, image_number):
    """Match .cc: sCutVar = sxy normally; sx for image 6/9; sy for 7/10."""
    sx = np.asarray(rmsx, dtype=float)
    sy = np.asarray(rmsy, dtype=float)
    sxy = np.sqrt(0.5 * (sx**2 + sy**2))
    if image_number in (6, 9):
        return sx
    if image_number in (7, 10):
        return sy
    return sxy


def load_energy_arrays(files, tree_name="clustersRec", image_number=5):
    """Load Energy, wSTD_XY, fSTD_X, fSTD_Y from tree; combine all files. Return arrays."""
    if uproot is None:
        raise RuntimeError("uproot is required. Install with: pip install uproot")
    energies = []
    rmsxy_list = []
    rmsx_list = []
    rmsy_list = []
    for fp in files:
        with uproot.open(fp) as f:
            if tree_name not in f:
                continue
            t = f[tree_name]
            branches = ["Energy", "wSTD_XY", "fSTD_X", "fSTD_Y"]
            if not all(b in t.keys() for b in branches):
                continue
            arr = t.arrays(branches, library="np")
            energies.append(arr["Energy"])
            rmsxy_list.append(arr["wSTD_XY"])
            rmsx_list.append(arr["fSTD_X"])
            rmsy_list.append(arr["fSTD_Y"])
    if not energies:
        return None, None, None, None
    E = np.concatenate([np.ravel(e) for e in energies])
    rmsx = np.concatenate([np.ravel(r) for r in rmsx_list])
    rmsy = np.concatenate([np.ravel(r) for r in rmsy_list])
    s_cut = s_cut_variable(rmsx, rmsy, image_number)
    return E, s_cut, rmsx, rmsy


def fill_histograms(E, s_cut, bins_edges):
    """Fill front (s_cut < SIGMA_FRONT_CUT) and back (s_cut >= SIGMA_FRONT_CUT) energy histograms in [ENE_MIN, ENE_MAX]."""
    mask_ew = (E >= ENE_MIN) & (E <= ENE_MAX)
    E_ew = E[mask_ew]
    s_ew = s_cut[mask_ew]
    front_mask = s_ew < SIGMA_FRONT_CUT
    back_mask = ~front_mask
    hist_front, _ = np.histogram(E_ew[front_mask], bins=bins_edges)
    hist_back, _ = np.histogram(E_ew[back_mask], bins=bins_edges)
    return hist_front, hist_back, bins_edges


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def double_gauss(x, A1, mu1, s1, A2, mu2, s2):
    return gauss(x, A1, mu1, s1) + gauss(x, A2, mu2, s2)


def fit_double_gaussian_peaks(counts, bin_edges):
    """
    Match .cc: fit first Gaussian in [G1_LO, G1_HI], then second in dynamic window.
    Returns (mu1, mu1_err, mu2, mu2_err) for Kα and Kβ, or (nan, 0, nan, 0) on failure.
    """
    from scipy.optimize import curve_fit
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    n = np.sum(counts)
    if n < 10:
        return np.nan, 0.0, np.nan, 0.0
    # First peak in [G1_LO, G1_HI]
    mask1 = (bin_centers >= G1_LO) & (bin_centers <= G1_HI)
    if np.sum(counts[mask1]) < 5:
        return np.nan, 0.0, np.nan, 0.0
    x1 = bin_centers[mask1]
    y1 = counts[mask1].astype(float)
    p0_1 = (np.max(y1), 5.9, 0.12)
    try:
        popt1, pcov1 = curve_fit(gauss, x1, y1, p0=p0_1, bounds=([0, 4.5, 0.03], [1e6, 6.5, 0.5]))
        mu1, sigma1 = popt1[1], popt1[2]
        mu1_err = np.sqrt(pcov1[1, 1]) if pcov1 is not None else 0.0
    except Exception:
        return np.nan, 0.0, np.nan, 0.0
    # Second peak window: mu1 + 2*sigma1, width 1 keV (match .cc)
    g2_left = max(ENE_MIN, mu1 + 2 * sigma1)
    g2_right = min(ENE_MAX, g2_left + 1.0)
    if g2_right <= g2_left:
        g2_right = min(ENE_MAX, g2_left + 0.5)
    mask2 = (bin_centers >= g2_left) & (bin_centers <= g2_right)
    if np.sum(counts[mask2]) < 3:
        return mu1, mu1_err, np.nan, 0.0
    x2 = bin_centers[mask2]
    y2 = counts[mask2].astype(float)
    p0_2 = (np.max(y2), (g2_left + g2_right) / 2, 0.12)
    try:
        popt2, pcov2 = curve_fit(gauss, x2, y2, p0=p0_2, bounds=([0, 5.5, 0.03], [1e6, 7.5, 0.5]))
        mu2 = popt2[1]
        mu2_err = np.sqrt(pcov2[1, 1]) if pcov2 is not None else 0.0
    except Exception:
        return mu1, mu1_err, np.nan, 0.0
    return mu1, mu1_err, mu2, mu2_err


def process_one_group(files, image_number, bins_edges):
    """Load files, fill front/back histos, fit; return (mu_KA_f, mu_KB_f, mu_KA_b, mu_KB_b) and errors."""
    E, s_cut, _, _ = load_energy_arrays(files, image_number=image_number)
    if E is None or len(E) == 0:
        return (np.nan, np.nan, np.nan, np.nan), (0.0, 0.0, 0.0, 0.0)
    hist_front, hist_back, _ = fill_histograms(E, s_cut, bins_edges)
    mu1_f, mu1_fe, mu2_f, mu2_fe = fit_double_gaussian_peaks(hist_front, bins_edges)
    mu1_b, mu1_be, mu2_b, mu2_be = fit_double_gaussian_peaks(hist_back, bins_edges)
    mus = (mu1_f, mu2_f, mu1_b, mu2_b)  # Kα front, Kβ front, Kα back, Kβ back
    errs = (mu1_fe, mu2_fe, mu1_be, mu2_be)
    return mus, errs


def run_scan(filelist_path, image_number=5, module_label="DM02", output_dir=".", plot=True):
    """Parse list, process each (grid_point, ext), build results and optional plots."""
    path_groups, keys = parse_file_list(filelist_path)
    if not path_groups:
        print("No valid ROOT files found in list.")
        return
    bins_edges = np.linspace(ENE_MIN, ENE_MAX, NBINS + 1)
    # Unique grid points in stable order (use first occurrence order)
    seen_gp = []
    for (gp, ext) in keys:
        if gp not in seen_gp:
            seen_gp.append(gp)
    scan_ids = [f"S{i+1}" for i in range(len(seen_gp))]
    grid_point_to_scan_id = {gp: scan_ids[i] for i, gp in enumerate(seen_gp)}
    exts = ["Ext1", "Ext2", "Ext3", "Ext4"]
    # Results: [scan_id][module][ext] -> value (we have single module here)
    mu_KA_front = defaultdict(lambda: defaultdict(dict))
    mu_KB_front = defaultdict(lambda: defaultdict(dict))
    mu_KA_back = defaultdict(lambda: defaultdict(dict))
    mu_KB_back = defaultdict(lambda: defaultdict(dict))
    for (grid_point, ext_num) in keys:
        sid = grid_point_to_scan_id[grid_point]
        ext_name = exts[ext_num - 1]
        files = path_groups[(grid_point, ext_num)]
        mus, errs = process_one_group(files, image_number, bins_edges)
        mu_KA_front[sid][module_label][ext_name] = mus[0]
        mu_KB_front[sid][module_label][ext_name] = mus[1]
        mu_KA_back[sid][module_label][ext_name] = mus[2]
        mu_KB_back[sid][module_label][ext_name] = mus[3]
        print(f"  {sid} {module_label} {ext_name}: Kα_f={mus[0]:.3f} Kβ_f={mus[1]:.3f} Kα_b={mus[2]:.3f} Kβ_b={mus[3]:.3f}")
    # Fill missing (scan_id, module, ext) with nan so plots don't break
    for sid in scan_ids:
        for e in exts:
            if e not in mu_KA_front[sid][module_label]:
                mu_KA_front[sid][module_label][e] = np.nan
                mu_KB_front[sid][module_label][e] = np.nan
                mu_KA_back[sid][module_label][e] = np.nan
                mu_KB_back[sid][module_label][e] = np.nan
    # CSV
    outpath = Path(output_dir)
    outpath.mkdir(parents=True, exist_ok=True)
    csv_path = outpath / "energy_peaks_scan.csv"
    with open(csv_path, "w") as f:
        f.write("scan_id,module,extension,mu_KA_front,mu_KB_front,mu_KA_back,mu_KB_back\n")
        for sid in scan_ids:
            for e in exts:
                f.write(f"{sid},{module_label},{e},"
                        f"{mu_KA_front[sid][module_label][e]:.4f},"
                        f"{mu_KB_front[sid][module_label][e]:.4f},"
                        f"{mu_KA_back[sid][module_label][e]:.4f},"
                        f"{mu_KB_back[sid][module_label][e]:.4f}\n")
    print(f"Wrote {csv_path}")
    # Plots (same style as parameter_scan_comparison)
    if plot:
        _plot_scan_summary(mu_KA_front, mu_KB_front, "Front events", module_label, exts, scan_ids, outpath)
        _plot_scan_summary(mu_KA_back, mu_KB_back, "Back events", module_label, exts, scan_ids, outpath)
    return mu_KA_front, mu_KB_front, mu_KA_back, mu_KB_back


def _plot_scan_summary(mu_KA, mu_KB, title, module_label, exts, scan_ids, output_dir):
    """One figure: Kβ top, Kα bottom, vs (Module, Extension), one point per scan_id."""
    modules = [module_label]
    x_positions = {}
    x_all, xtick_labels = [], []
    for mi, m in enumerate(modules):
        for ei, e in enumerate(exts):
            x = mi * len(exts) + ei
            x_positions[(m, e)] = x
            x_all.append(x)
            xtick_labels.append(f"{m}\n{e}")
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"][:len(scan_ids)]
    jitter = np.linspace(-0.15, 0.15, len(scan_ids))
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True, constrained_layout=True)
    for si, sid in enumerate(scan_ids):
        xs, ys = [], []
        for m in modules:
            for e in exts:
                xs.append(x_positions[(m, e)] + jitter[si])
                ys.append(mu_KB[sid][m][e])
        ax_top.scatter(xs, ys, marker=markers[si], label=sid)
    ax_top.axhline(E_KB_NOMINAL, linestyle="--", color="gray")
    ax_top.set_ylabel(r"$\mu$ [keV]")
    ax_top.set_title("Kβ")
    for si, sid in enumerate(scan_ids):
        xs, ys = [], []
        for m in modules:
            for e in exts:
                xs.append(x_positions[(m, e)] + jitter[si])
                ys.append(mu_KA[sid][m][e])
        ax_bot.scatter(xs, ys, marker=markers[si], label=sid)
    ax_bot.axhline(E_KA_NOMINAL, linestyle="--", color="gray")
    ax_bot.set_ylabel(r"$\mu$ [keV]")
    ax_bot.set_title("Kα")
    ax_bot.set_xticks(x_all)
    ax_bot.set_xticklabels(xtick_labels)
    ax_top.legend(ncol=min(len(scan_ids), 8), loc="upper center", bbox_to_anchor=(0.5, 1.15), frameon=True)
    fig.suptitle(f"{title} — peak mean vs Extension ({module_label})", y=1.02)
    out_file = output_dir / f"scan_energy_peaks_{title.lower().replace(' ', '_')}.png"
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_file}")


def main():
    ap = argparse.ArgumentParser(description="Parameter scan: compute Fe-55 energy peaks per grid point and extension from ROOT cluster trees.")
    ap.add_argument("filelist", help="Path to file list (one ROOT path per line; names ..._<image_id>_ext<1-4>.root)")
    ap.add_argument("--image-number", type=int, default=5, help="Image number for front/back cut (5=sxy, 6/9=sx, 7/10=sy)")
    ap.add_argument("--module", default="DM02", help="Module label for output (e.g. DM02, 109)")
    ap.add_argument("--output-dir", "-o", default=".", help="Output directory for CSV and plots")
    ap.add_argument("--no-plot", action="store_true", help="Skip generating plots")
    args = ap.parse_args()
    run_scan(args.filelist, image_number=args.image_number, module_label=args.module,
             output_dir=args.output_dir, plot=not args.no_plot)


if __name__ == "__main__":
    main()
