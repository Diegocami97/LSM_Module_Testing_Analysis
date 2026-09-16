#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Diego Venegas-Vargas
# DAMIC-M Collaboration
# Fits Gaussian peaks to CCD charge spectra from FITS images or ROOT clusters.
"""
Fit Gaussian peaks to CCD charge spectra (from FITS images or ROOT clusters).

Usage examples:
  # From FITS (ADU -> electrons via --gain if --inadu is set)
  python fit_peaks.py -f '/path/to/*.fits' --inadu true -k 3.74 --ne 0 50 --de 0.05 --ext 1 -o out

  # From ROOT tree (expects TTree "clustersRec" with pixels_E in keV)
  python fit_peaks.py -r clusters.root -k 3.74 --ne 0 50 --de 0.05 -o out

  # Per-extension / per-file gains (e.g. each ROOT file has a different gain)
  python fit_peaks.py -r ext4.root --gain-map gains.json -k 3.74 --ne 0 50 -o out
  # gains.json: {"ext4.root": 3.72, "ext5.root": 3.75}  (keys: basename or no extension)

  # Combine multiple ROOT files; gain from each file's "info" tree (ext{N}_gainext); -k optional
  python fit_peaks.py -r 'ext4_*.root' --root-ext 4 --ne 0 50 -o out

  # Multiple ROOT files (options first, then files; works when shell expands a glob)
  python fit_peaks.py --root-ext 3 --ne 200 210 -o out path/to/DM02/*ext3.root
"""

import json
import os
import sys
import click

# To avoid ROOT help printout unless explicitly asking for help
if ("-h" not in sys.argv) and ("--help" not in sys.argv):
    import ROOT
    from array import array
    import numpy as np
    from glob import glob
    from astropy.io import fits
    # SciPy is optional; only needed for the global multi-peak fit mode
    try:
        from scipy.optimize import curve_fit  # type: ignore
    except Exception:
        curve_fit = None
else:
    import ROOT  # type: ignore
    from array import array  # type: ignore
    import numpy as np  # type: ignore
    from glob import glob  # type: ignore
    from astropy.io import fits  # type: ignore
    try:
        from scipy.optimize import curve_fit  # type: ignore
    except Exception:
        curve_fit = None


def global_multi_peak_fit(hist_e, ne_ini, ne_max):
    """
    Global multi-Gaussian fit to extract peak centroids in electrons.

    Model:
      y(q) = sum_k A_k * exp(-(q - (i_k + d_k))^2 / (2*sigma^2)) + b0

    where:
      i_k runs from ne_ini..ne_max (integer peak number),
      d_k is a small offset for each peak,
      sigma is shared, b0 is a constant background.

    Returns:
      x_num_electron, centers_e, e_centers_e, sigma_list, e_sigma_list
    """
    if curve_fit is None:
        raise RuntimeError(
            "Global multi-peak fit requested but SciPy is not available. "
            "Install scipy or avoid using --fit-mode global."
        )

    nbins = hist_e.GetNbinsX()
    bin_centers = np.array([hist_e.GetBinCenter(i) for i in range(1, nbins + 1)], dtype=float)
    counts = np.array([hist_e.GetBinContent(i) for i in range(1, nbins + 1)], dtype=float)

    # Restrict to a slightly expanded region around [ne_ini, ne_max]
    mask = (bin_centers >= (ne_ini - 2.0)) & (bin_centers <= (ne_max + 2.0))
    x = bin_centers[mask]
    y = counts[mask]

    if x.size == 0 or np.all(y <= 0):
        raise RuntimeError("Histogram is empty in requested fit range; cannot perform global fit.")

    n_peaks = ne_max - ne_ini + 1

    def multi_gauss(xvals, *params):
        n_loc = (len(params) - 2) // 2
        sigma = params[-2]
        b0 = params[-1]
        y_model = np.zeros_like(xvals, dtype=float)
        for k in range(n_loc):
            A = params[2 * k]
            d = params[2 * k + 1]
            i_peak = ne_ini + k
            mu = i_peak + d
            y_model += A * np.exp(-0.5 * ((xvals - mu) / sigma) ** 2)
        y_model += b0
        return y_model

    # Initial guesses and bounds
    A0 = max((y.max() - np.median(y)) / max(n_peaks, 1), 1.0)
    sigma0 = 0.3
    b0 = max(np.median(y), 0.0)

    p0 = []
    lo = []
    hi = []
    for _ in range(n_peaks):
        p0.extend([A0, 0.0])     # A_k, d_k
        lo.extend([0.0, -0.5])   # amplitudes >=0, offsets in [-0.5, 0.5]
        hi.extend([np.inf, 0.5])
    p0.append(sigma0)
    lo.append(0.05)
    hi.append(2.0)
    p0.append(b0)
    lo.append(0.0)
    hi.append(np.inf)

    popt, pcov = curve_fit(
        lambda xx, *pp: multi_gauss(xx, *pp),
        x,
        y,
        p0=p0,
        bounds=(lo, hi),
        maxfev=40000,
    )

    centers_e = []
    e_centers_e = []
    # one shared sigma; propagate its uncertainty equally to all
    sigma_val = popt[-2]
    sigma_err = float(np.sqrt(pcov[-2, -2])) if pcov.size else float("nan")

    for k in range(n_peaks):
        d_idx = 2 * k + 1
        d_hat = popt[d_idx]
        i_peak = ne_ini + k
        mu_hat = i_peak + d_hat
        centers_e.append(mu_hat)
        err = float(np.sqrt(pcov[d_idx, d_idx])) if pcov.size else float("nan")
        e_centers_e.append(err)

    x_num_electron = list(range(ne_ini, ne_max + 1))
    sigma_list = [sigma_val] * len(x_num_electron)
    e_sigma_list = [sigma_err] * len(x_num_electron)

    return x_num_electron, centers_e, e_centers_e, sigma_list, e_sigma_list


def read_gain_from_info_tree(filepath, ext_num, in_adu_per_e=True):
    """
    Read gain from ROOT file's TTree "info", branch ext{N}_gainext.
    The branch is stored in ADU/e- (same as -k). Returns that value, or None if not found.
    """
    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        return None
    tinfo = f.Get("info")
    if not tinfo or tinfo.GetEntries() == 0:
        f.Close()
        return None
    branch_name = f"ext{ext_num}_gainext"
    if not tinfo.FindBranch(branch_name):
        f.Close()
        return None
    gainext = np.array([0.0], dtype=np.float64)
    tinfo.SetBranchAddress(branch_name, gainext)
    tinfo.GetEntry(0)
    val = float(gainext[0])
    f.Close()
    if val <= 0:
        return None
    return val if in_adu_per_e else (1.0 / val)


E_PER_ELECTRON = 3.74  # eV per electron (keV -> e-: energy_keV * 1000 / E_PER_ELECTRON)

# Fe-55 X-ray line energies [keV] -> electron equivalent for residual extrapolation
FE55_KA_KEV = 5.898
FE55_KB_KEV = 6.490


def keV_to_electrons(kev):
    """Convert energy in keV to electrons using E_PER_ELECTRON."""
    return kev * 1000.0 / E_PER_ELECTRON


def fill_hist_from_clusters_tree(tree, hist_e):
    """
    Fill hist_e (in electrons) from TTree 'clustersRec'.
    Branch pixels_E is vector<vector<float>>: per-cluster list of pixel energies in keV.
    Converts each pixel energy to electrons (keV * 1000 / 3.74) and fills the histogram.
    """
    n_ent = tree.GetEntries()
    for i in range(n_ent):
        tree.GetEntry(i)
        try:
            pe = tree.pixels_E
        except AttributeError:
            pe = None
        if pe is None:
            continue
        for cluster in pe:
            try:
                for e_kev in cluster:
                    hist_e.Fill(float(e_kev) * 1000.0 / E_PER_ELECTRON)
            except TypeError:
                hist_e.Fill(float(cluster) * 1000.0 / E_PER_ELECTRON)
    return hist_e


def build_energy_sigma_xy(tree, nbins_e=50, e_min=4.5, e_max=7.0, nbins_s=50, s_min=0.0, s_max=2.0):
    """
    Build a TH2F of cluster Energy (keV) vs sigma_xy (pixels) from a 'clustersRec' tree.
    Expects branches: Energy (cluster keV), rmsx, rmsy (pixels).
    """
    h2 = ROOT.TH2F("hE_sigma_xy", "Cluster Energy vs #sigma_{xy}",
                   nbins_e, e_min, e_max,
                   nbins_s, s_min, s_max)
    n_ent = tree.GetEntries()
    for i in range(n_ent):
        tree.GetEntry(i)
        if not (hasattr(tree, "Energy") and hasattr(tree, "rmsx") and hasattr(tree, "rmsy")):
            continue
        try:
            energies = tree.Energy
            sx_list = tree.rmsx
            sy_list = tree.rmsy
        except AttributeError:
            continue
        nC = len(energies)
        for n in range(nC):
            E_kev = float(energies[n])
            sx = float(sx_list[n])
            sy = float(sy_list[n])
            sxy = (0.5 * (sx * sx + sy * sy)) ** 0.5
            if e_min <= E_kev <= e_max:
                h2.Fill(E_kev, sxy)
    h2.SetDirectory(0)
    return h2


def get_histe(infile, hist_e, inADU, gain, cols=(8, -1), rows=(0, -1), ext=0):
    """Fill a ROOT TH1 from FITS images (optionally converting ADU -> e-)."""
    lof = glob(infile)
    CF = (1.0 / gain) if inADU else 1.0
    for f in lof:
        print(f"[READ] {f}")
        d = fits.getdata(f, ext).astype(float)[rows[0]:rows[1], cols[0]:cols[1]]
        for qij in d.ravel():
            hist_e.Fill(qij * CF)
    return hist_e


@click.command(help="Script to fit gaussian peaks to the reconstructed image")
@click.option("-f", "--infile", type=str,
              help="Input FITS pattern (e.g., '/path/*.fits') with calibrated data.")
@click.option("--inadu", type=bool, default=False,
              help="If input data are in ADU (FITS path). Converts ADU->e- using --gain.")
@click.option("-r", "--root", type=str, multiple=True,
              help="ROOT file(s): repeat -r for each file, or pass one quoted glob (e.g. -r 'path/*.root'). Quote the glob so the shell does not expand it.")
@click.option("--root-ext", type=int, default=0,
              help="Extension number for reading gain from 'info' tree (branch ext{N}_gainext). Default 0.")
@click.option("-k", "--gain", type=float, default=None,
              help="Gain (ADU per e-). Optional when ROOT 'info' tree has ext{N}_gainext or --gain-map has the file. Required for FITS unless --gain-map matches.")
@click.option("--ne", type=(int, int), default=(0, 1350),
              help="Minimum and maximum number of electrons to be fitted (inclusive).")
@click.option("--de", type=float, default=0.05,
              help="Bin width for the spectrum in electrons.")
@click.option("--display", type=bool, default=False,
              help="Display canvases instead of saving PNGs.")
@click.option("--cols", type=(int, int), default=(8, -1),
              help="Column range [start, end) for --infile.")
@click.option("--rows", type=(int, int), default=(0, -1),
              help="Row range [start, end) for --infile.")
@click.option("--erange", type=(float, float), default=(0.35, 0.5),
              help="Half-window around each peak center for fit: [left, right] => mu±(left,right).")
@click.option("--fixed-window", is_flag=True, default=False,
              help="Center fit window at integer peak number i (no chain drift). Use for stable residuals and extrapolation to high e-.")
@click.option("--fit-mode", type=click.Choice(["chained", "fixed", "global"]), default="chained",
              help="Peak fitting mode: 'chained' (original running mean), 'fixed' (integer-centered windows), or 'global' (global multi-Gaussian fit in Python).")
@click.option("--ext", type=int, default=1, help="FITS extension for data.")
@click.option("-o", "--output", type=str, default='.', help="Output directory.")
@click.option("--residual-fit-max", type=float, default=500.0,
              help="Max reconstructed charge (e-) to use for parabola fit in residual plot.")
@click.option("--residual-reliable-max", type=float, default=800.0,
              help="Reconstructed charge (e-) above which data are annotated as unreliable.")
@click.option("--extrapolate-max", type=float, default=0.0,
              help="Extend residual plot x-axis and parabola to this reconstructed charge (e-). E.g. 2000 to extrapolate to 2000 e-.")
@click.option("--residual-xrange", type=(float, float), default=(0.0, 0.0),
              help="Optional x-axis range [xmin, xmax] (e-) for the non-linearity residual plot.")
@click.option("--residual-yrange", type=(float, float), default=(0.0, 0.0),
              help="Optional y-axis range [ymin, ymax] (e-) for the non-linearity residual plot.")
@click.option("--residual-fit-model", type=click.Choice(["pol2", "scaledpol2", "pol3", "scaledpol3", "rational"]), default="pol2",
              help="Model used to fit residuals vs reconstructed charge. 'scaled*' uses t=x/1000 for numerical stability; 'rational' is (a+b t+c t^2)/(1+d t).")
@click.option("--no-residual-fit", is_flag=True, default=False,
              help="Do not draw the parabola fit on the residual plot (data points only).")
@click.option("--charge-xrange", type=(float, float), default=(0.0, 0.0),
              help="Optional x-axis range [xmin, xmax] in electrons for the charge distribution plot.")
@click.option("--charge-plot-title", type=str, default=None,
              help="Optional title for charge distribution plot (e.g. 'MASK Selected Rows - Channel 3').")
@click.option("--gain-map", type=click.Path(exists=True), default=None,
              help="JSON file mapping filename -> gain (ADU/e-) for per-file gains. Keys: basename (e.g. 'file.root') or no extension (e.g. 'file'). -k used if key missing.")
@click.argument("root_files", nargs=-1, type=click.Path(exists=False))
def main(infile, inadu, root, root_ext, gain, ne, de, display, cols, rows, erange, fixed_window, fit_mode, ext, output,
         residual_fit_max, residual_reliable_max, extrapolate_max, residual_xrange, residual_yrange, residual_fit_model, no_residual_fit,
         charge_xrange, charge_plot_title, gain_map, root_files):
    # Resolve ROOT input: -r (tuple when multiple=True) and/or positional root_files (e.g. shell-expanded glob)
    root_spec = None
    if root:
        if len(root) == 1 and "*" in root[0]:
            root_spec = root[0]
        else:
            root_spec = list(root)
    if root_files:
        if len(root_files) == 1 and "*" in root_files[0] and root_spec is None:
            root_spec = root_files[0]
        else:
            extra = list(root_files)
            if root_spec is not None and isinstance(root_spec, list):
                root_spec = root_spec + extra
            elif root_spec is None:
                root_spec = extra
    # Basic checks
    if (infile is None) and (root_spec is None):
        raise click.UsageError("Provide either --infile (FITS) or ROOT files (e.g. -r 'path/*.root' or list files after options).")
    if (infile is not None) and (root_spec is not None):
        raise click.UsageError("Provide only one of --infile or ROOT files, not both.")

    # Resolve gain: use --gain-map if provided and key found, else -k
    gain_map_dict = {}
    if gain_map:
        with open(gain_map) as f:
            gain_map_dict = json.load(f)

    def _effective_gain(input_path, base_tag_from_glob=None):
        """Resolve gain from map (by basename or basename no ext) or fallback to -k. Returns (gain, from_map). gain may be None if -k not set and not in map."""
        if gain_map_dict:
            base = os.path.basename(input_path) if input_path else (base_tag_from_glob or "")
            base_no_ext = os.path.splitext(base)[0] if base else ""
            for key in (base, base_no_ext):
                if key in gain_map_dict:
                    return float(gain_map_dict[key]), True
        return (float(gain) if gain is not None else None), False

    os.makedirs(output, exist_ok=True)

    # Batch mode handling
    ROOT.gROOT.SetBatch(not display)

    ne_ini, ne_max = ne

    # Build histogram (either from ROOT TTree draw or manual fill from FITS)
    hist_name = "hist_e"
    hist_e = None
    energy_sigma_tree = None
    base_tag = "spectra"

    if infile is None:
        # Resolve ROOT file list from root_spec (string = glob, list = explicit files)
        if isinstance(root_spec, str):
            root_paths = sorted(glob(root_spec))
        else:
            root_paths = root_spec
        if not root_paths:
            raise RuntimeError(f"No ROOT files found for: {root_spec}")

        nbins = int((ne_max - ne_ini) / de)
        hist_e = ROOT.TH1D(hist_name, hist_name, nbins, ne_ini - 3, ne_max + 3)
        hist_e.Sumw2()

        if len(root_paths) == 1:
            f = ROOT.TFile.Open(root_paths[0])
            if not f or f.IsZombie():
                raise RuntimeError(f"Failed to open ROOT file: {root_paths[0]}")
            clusters = f.Get("clustersRec")
            if clusters is None:
                raise RuntimeError("Could not find TTree 'clustersRec' in ROOT file.")
            fill_hist_from_clusters_tree(clusters, hist_e)
            hist_e.SetDirectory(0)
            energy_sigma_tree = clusters
            base_tag = os.path.splitext(os.path.basename(root_paths[0]))[0]
            input_for_gain = root_paths[0]
        else:
            # Combine multiple ROOT files with TChain
            chain = ROOT.TChain("clustersRec")
            for rp in root_paths:
                chain.Add(rp)
            if chain.GetEntries() == 0:
                raise RuntimeError("TChain has no entries.")
            fill_hist_from_clusters_tree(chain, hist_e)
            hist_e.SetDirectory(0)
            energy_sigma_tree = chain
            base_tag = "combined"
            input_for_gain = root_paths[0]
    else:
        # From FITS images: convert to electrons if --inadu
        nbins = int((ne_max - ne_ini) / de)
        hist_e = ROOT.TH1D(hist_name, hist_name, nbins, ne_ini - 3, ne_max + 3)
        hist_e.Sumw2()
        globbed = glob(infile)
        input_for_gain = globbed[0] if globbed else None
        gain, _ = _effective_gain(input_for_gain or "")
        if gain is None or gain <= 0:
            raise click.UsageError(
                "FITS input requires gain. Provide -k (e.g. -k 3.74) or --gain-map with an entry for this file."
            )
        get_histe(infile, hist_e, inadu, gain, cols=cols, rows=rows, ext=ext)
        base_tag = os.path.splitext(os.path.basename(globbed[0]))[0] if globbed else "fits_input"

    # Resolve gain for ROOT: info tree (per file or mean over files) > gain-map > -k
    if infile is None:
        # root_paths already set above when infile is None
        gains_adu = []
        for rp in root_paths:
            g = read_gain_from_info_tree(rp, root_ext, in_adu_per_e=True)
            if g is not None:
                gains_adu.append(g)
        if gains_adu:
            gain = float(np.mean(gains_adu))
            if len(gains_adu) > 1:
                print(f"[gain] Using mean gain from info tree: {gain:.6f} ADU/e- (over {len(gains_adu)} files)")
            else:
                print(f"[gain] Using gain from info tree: {gain:.6f} ADU/e- for {os.path.basename(root_paths[0])}")
        else:
            gain, gain_from_map = _effective_gain(input_for_gain)
            if gain_from_map and input_for_gain:
                print(f"[gain] Using gain {gain} for {os.path.basename(input_for_gain)} (from --gain-map)")
    else:
        gain, gain_from_map = _effective_gain(input_for_gain)
        if gain_from_map and input_for_gain:
            print(f"[gain] Using gain {gain} for {os.path.basename(input_for_gain)} (from --gain-map)")

    if gain is None or gain <= 0:
        raise click.UsageError(
            "No gain available. For ROOT: ensure files have TTree 'info' with branch ext{N}_gainext (use --root-ext N), "
            "or provide -k or --gain-map. For FITS: provide -k or --gain-map."
        )

    if not hist_e:
        raise RuntimeError("Failed to create/fill histogram.")

    # Global ROOT style: remove stats/fit boxes from all plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    # Style
    hist_e.GetXaxis().SetTitle("n_{e^{-}}")
    hist_e.GetYaxis().SetTitle("Counts")
    hist_e.SetDirectory(0)
    hist_e.SetStats(0)

    # --- Fit individual Gaussian peaks ---
    gauss = []
    centers_e, e_centers_e = [], []   # centroid and error in electrons (fit is on spectrum in e-)
    sigma, e_sigma = [], []          # in electrons
    x_num_electron = []

    # Legacy / ROOT-based methods: chained and fixed-window
    if fit_mode in ("chained", "fixed"):
        mu_running = ne_ini - 1
        fit_opts = "LEMRQ0"
        use_fixed = (fit_mode == "fixed") or fixed_window

        for i in range(ne_ini, ne_max + 1):
            if use_fixed:
                center = float(i)
            else:
                center = mu_running + 1
            x_min = center - erange[0]
            x_max = center + erange[1]

            fgi = ROOT.TF1(f"ne_{i}", "gaus", x_min, x_max)
            fgi.SetParLimits(1, x_min, x_max)   # mean
            fgi.SetParLimits(2, 0.01, 1.0)      # sigma in electrons
            fgi.SetLineColor(ROOT.kRed + 1)

            hist_e.Fit(fgi, fit_opts)
            # Retry once if the fit latched onto a wrong region
            if use_fixed:
                if abs(fgi.GetParameter(1) - i) > 2:
                    fgi.SetParLimits(1, x_min, x_max)
                    hist_e.Fit(fgi, fit_opts)
            else:
                if (fgi.GetParameter(1) - mu_running) > 2:
                    fgi.SetParLimits(1, x_min, x_max)
                    hist_e.Fit(fgi, fit_opts)

            mu_running = fgi.GetParameter(1)

            # Store: centroid and error in electrons; sigma in electrons
            centers_e.append(fgi.GetParameter(1))
            e_centers_e.append(fgi.GetParError(1))
            sigma.append(fgi.GetParameter(2))
            e_sigma.append(fgi.GetParError(2))
            x_num_electron.append(i)
            gauss.append(fgi)
    elif fit_mode == "global":
        # Global multi-Gaussian fit in Python (SciPy) on the electron spectrum
        x_num_electron, centers_e, e_centers_e, sigma, e_sigma = global_multi_peak_fit(
            hist_e, ne_ini, ne_max
        )
    else:
        raise RuntimeError(f"Unknown fit mode: {fit_mode}")

    # ADU only where needed (calibration, normalized, delta mu)
    centers = [c * gain for c in centers_e]
    e_centers = [e * gain for e in e_centers_e]

    # Print centroid vs peak number (line value) for comparison
    print("")
    print("  n_e (line)   centroid (e-)   ± err       residual (e-)")
    print("  " + "-" * 55)
    for i, (n_e, ce) in enumerate(zip(x_num_electron, centers_e)):
        res = ce - n_e
        ee = e_centers_e[i] if i < len(e_centers_e) else 0.0
        print(f"  {n_e:4d}         {ce:12.6f}   ± {ee:.6f}   {res:12.6f}")
    print("")

    # Canvas: spectrum with fitted Gaussians
    c_spec = ROOT.TCanvas("c_spec", "Spectrum with Gaussian fits", 900, 700)
    hist_e.Draw("HIST")
    for fgi in gauss:
        fgi.Draw("SAME")
    c_spec.SetLogy()
    if display:
        c_spec.Update()
    else:
        out_png = os.path.join(output, f"{base_tag}_FittedGauss.pdf")
        c_spec.SaveAs(out_png)
        print(f"[SAVE] {out_png}")

    # --- Charge distribution with peak-number markers (non-linearity visualization) ---
    c_charge = ROOT.TCanvas("c_charge", "Charge distribution with peak numbers", 900, 700)
    c_charge.SetLogy(1)
    hist_charge = hist_e.Clone("hist_charge")
    hist_charge.SetTitle("Pixel Charge distribution" + (" - " + charge_plot_title if charge_plot_title else ""))
    hist_charge.GetXaxis().SetTitle("Reconstructed pixel charge [e^{-}]")
    hist_charge.GetYaxis().SetTitle("counts")
    hist_charge.SetLineWidth(2)
    hist_charge.SetLineColor(ROOT.kBlue)
    hist_charge.SetFillColorAlpha(ROOT.kBlue, 0.35)
    hist_charge.SetFillStyle(1000)
    hist_charge.SetStats(0)
    # Optional user-specified x-axis range for the charge distribution plot
    if charge_xrange and len(charge_xrange) == 2:
        x_min, x_max = charge_xrange
        if x_max > x_min:
            hist_charge.GetXaxis().SetRangeUser(x_min, x_max)
    hist_charge.Draw("HIST")
    ymax = hist_charge.GetMaximum()
    ymin = hist_charge.GetMinimum()
    if ymin <= 0:
        ymin = hist_charge.GetBinContent(hist_charge.GetMaximumBin()) * 1e-5
    # Use histogram data range for line length: min positive bin content to max
    ymin_data = ymax
    for b in range(1, hist_charge.GetNbinsX() + 1):
        c = hist_charge.GetBinContent(b)
        if c > 0:
            ymin_data = min(ymin_data, c)
    if ymin_data >= ymax:
        ymin_data = ymin
    ymin_data = max(ymin_data, 1e-10)
    peak_lines = []
    peak_labels = []
    # Draw dashed lines at the reconstructed centroids (in e-) and label them with the peak number n_e
    # This matches the reference style where x-axis is reconstructed charge but line labels are integer peak numbers.
    for n_e, mu_e in zip(x_num_electron, centers_e):
        line = ROOT.TLine(mu_e, ymin_data, mu_e, ymax)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(7)
        line.SetLineWidth(1)
        line.Draw("SAME")
        peak_lines.append(line)
        txt = ROOT.TLatex()
        txt.SetTextSize(0.024)
        txt.SetTextColor(ROOT.kRed)
        txt.DrawLatex(mu_e + 0.15, ymax * 0.92, str(n_e))
        peak_labels.append(txt)
    c_charge.Modified()
    c_charge.Update()
    if display:
        c_charge.Update()
    else:
        out_charge = os.path.join(output, f"{base_tag}_charge_distribution_peak_numbers.pdf")
        c_charge.SaveAs(out_charge)
        print(f"[SAVE] {out_charge}")

    # --- Cluster Energy vs sigma_xy (for ROOT inputs only, when branches exist) ---
    if infile is None and energy_sigma_tree is not None:
        try:
            hE_sig = build_energy_sigma_xy(energy_sigma_tree)
        except Exception:
            hE_sig = None
        if hE_sig:
            c_es = ROOT.TCanvas("c_es", "Energy vs sigma_xy", 900, 700)
            hE_sig.SetTitle("Cluster Energy vs #sigma_{xy}")
            hE_sig.GetXaxis().SetTitle("Energy [keV]")
            hE_sig.GetYaxis().SetTitle("#sigma_{xy} [pix]")
            hE_sig.Draw("COLZ")
            c_es.Modified()
            c_es.Update()
            if display:
                c_es.Update()
            else:
                out_es = os.path.join(output, f"{base_tag}_energy_vs_sigma_xy.pdf")
                c_es.SaveAs(out_es)
                print(f"[SAVE] {out_es}")

    # --- Calibration: peak centers vs electron number (ADU vs e-) ---
    c_cent = ROOT.TCanvas("c_cent", "Centers vs n_e", 900, 700)

    x = x_num_electron
    y = centers           # ADU
    ey = e_centers        # ADU
    ex = [0.0] * len(x)

    print("x (e-)    y (ADU)    ey (ADU)")
    print( "---------------------------")
    for xi, yi, eyi in zip(x, y, ey):
        print(f"{xi:5d}  {yi:8.2f}  {eyi:8.2f}")    

    g_cent = ROOT.TGraphErrors(
        len(y),
        array('d', [float(xx) for xx in x]),
        array('d', [float(yy) for yy in y]),
        array('d', ex),
        array('d', [float(e) for e in ey]),
    )
    g_cent.SetMarkerStyle(21)
    g_cent.SetMarkerSize(0.8)
    g_cent.Draw("AP E1")             # draw before GetHistogram()
    c_cent.Update()

    hframe = g_cent.GetHistogram()
    if hframe:
        hframe.GetXaxis().SetTitle("n_{e^{-}}")
        hframe.GetYaxis().SetTitle("#mu_{n_{e^{-}}} [ADU/e^{-}]")

    pol1 = ROOT.TF1("fitfunc_cent", "pol1", ne_ini, ne_max)
    pol1.SetLineColor(ROOT.kRed + 1)
    pol1.SetLineStyle(9)
    pol1.SetParameter(0, 0)          # allow small offset; don't hard-fix
    pol1.SetParLimits(1, 0, 2 * gain)
    pol1.SetParameter(1, gain)

    g_cent.Fit(pol1, "EMR0")
    pol1.Draw("SAME")
    c_cent.Update()

    # # # ---------- Extend calibration to Fe-55 and annotate predicted ADU ----------
    # # e_per_electron = 3.74  # eV / e-
    # # E_targets_eV = {"K#alpha": 5895.0, "K#beta": 6490.0}
    # # ne_targets = {name: (E / e_per_electron) for name, E in E_targets_eV.items()}
    # # mu_pred_ADU = {name: pol1.Eval(ne) for name, ne in ne_targets.items()}

    # # xmax_needed = 1.05 * max(ne_targets.values())
    # # hframe = g_cent.GetHistogram()
    # # if hframe:
    # #     hframe.GetXaxis().SetRangeUser(ne_ini, max(xmax_needed, hframe.GetXaxis().GetXmax()))
    # #     ymax_needed = 1.10 * max(hframe.GetMaximum(), max(mu_pred_ADU.values()))
    # #     hframe.SetMaximum(ymax_needed)

    # # # Draw the fitted line out to xmax_needed
    # # pol1.SetRange(ne_ini, xmax_needed)
    # # pol1.Draw("SAME")

    # # col_map = {"K#alpha": ROOT.kBlue + 1, "K#beta": ROOT.kGreen + 2}
    # # for name, ne_val in ne_targets.items():
    # #     mu_val = mu_pred_ADU[name]
    # #     col = col_map.get(name, ROOT.kBlue + 1)

    # #     v = ROOT.TLine(ne_val, 0, ne_val, mu_val);  v.SetLineStyle(7); v.SetLineColor(col); v.Draw()
    # #     h = ROOT.TLine(ne_ini, mu_val, ne_val, mu_val); h.SetLineStyle(7); h.SetLineColor(col); h.Draw()
    # #     m = ROOT.TMarker(ne_val, mu_val, 20); m.SetMarkerColor(col); m.SetMarkerSize(1.0); m.Draw()

    # #     txt = ROOT.TLatex()
    # #     txt.SetTextSize(0.030)
    # #     txt.DrawLatex(ne_val * 0.70, mu_val * 1.02,
    # #                   f"{name}: n_{{e}}={ne_val:.1f},  #mu={mu_val:.1f} ADU")

    # # c_cent.Modified(); c_cent.Update()

    # # # ---------- Fe-55 guides, labels, and axis extension (robust for TGraph) ----------
    # # e_per_electron = 3.74  # eV / e-
    # # E_targets_eV = {"K#alpha": 5895.0, "K#beta": 6490.0}

    # # # n_e for each target and predicted ADU from the linear fit
    # # ne_targets  = {name: (E / e_per_electron) for name, E in E_targets_eV.items()}
    # # mu_pred_ADU = {name: pol1.Eval(ne)        for name, ne in ne_targets.items()}

    # # # Determine needed axis limits
    # # xmax_needed = 1.05 * max(ne_targets.values())
    # # # Keep current ymin; expand ymax to include predicted points
    # # # (If you prefer to keep the lower part tight, just don't set SetMinimum.)
    # # current_ymin = g_cent.GetMinimum()
    # # current_ymax = g_cent.GetMaximum()
    # # ymax_needed  = 1.10 * max(current_ymax, max(mu_pred_ADU.values()))

    # # # >>> IMPORTANT: use TGraph axis setters, not the frame's SetRangeUser <<<
    # # g_cent.GetXaxis().SetLimits(ne_ini, xmax_needed)  # extend X
    # # g_cent.SetMaximum(ymax_needed)                    # extend Y up
    # # # g_cent.SetMinimum(current_ymin)                # (optional) keep same bottom
    # # c_cent.Modified(); c_cent.Update()

    # # # Draw solid fit over the measured range [ne_ini, ne_max]
    # # pol1.SetRange(ne_ini, ne_max)
    # # pol1.SetLineStyle(1)   # solid
    # # pol1.Draw("SAME")

    # # # Draw dashed fit for the extrapolated segment (ne_max -> xmax_needed)
    # # pol1_ex = ROOT.TF1("fitfunc_cent_ex", "pol1", ne_max, xmax_needed)
    # # pol1_ex.SetParameters(pol1.GetParameter(0), pol1.GetParameter(1))
    # # pol1_ex.SetLineColor(pol1.GetLineColor())
    # # pol1_ex.SetLineStyle(2)  # dashed
    # # pol1_ex.Draw("SAME")

    # # # Guides, markers, labels
    # # col_map = {"K#alpha": ROOT.kBlue + 1, "K#beta": ROOT.kGreen + 2}
    # # for name, ne_val in ne_targets.items():
    # #     mu_val = mu_pred_ADU[name]
    # #     col = col_map.get(name, ROOT.kBlue + 1)

    # #     v = ROOT.TLine(ne_val, current_ymin if current_ymin is not None else 0.0, ne_val, mu_val)
    # #     v.SetLineStyle(7); v.SetLineColor(col); v.Draw()

    # #     h = ROOT.TLine(ne_ini, mu_val, ne_val, mu_val)
    # #     h.SetLineStyle(7); h.SetLineColor(col); h.Draw()

    # #     m = ROOT.TMarker(ne_val, mu_val, 20)
    # #     m.SetMarkerColor(col); m.SetMarkerSize(1.0); m.Draw()

    # #     txt = ROOT.TLatex()
    # #     txt.SetTextSize(0.030)
    # #     txt.SetTextColor(col)
    # #     txt.DrawLatex(ne_val * 0.70, mu_val * 1.02,
    # #                   f"{name}: n_{{e}}={ne_val:.1f},  #mu={mu_val:.1f} ADU/e-")

    # # c_cent.Modified(); c_cent.Update()



    # # # Extract linear fit parameters and print with correct unit conversions
    # # a = pol1.GetParameter(0)   # ADU
    # # b = pol1.GetParameter(1)   # ADU/e-
    # # print(f"[poly1] Calibration: offset a = {a:.3f} ADU, gain b = {b:.6f} ADU/e-")
    # # adu_per_keV = b * (1000.0 / e_per_electron)    # ADU/keV
    # # eV_per_ADU  = e_per_electron / b               # eV/ADU
    # # print(f"[poly1] => {adu_per_keV:.3f} ADU/keV, {eV_per_ADU:.3f} eV/ADU")

    # # Move stats box if present
    # try:
    #     stat = c_cent.GetPrimitive("stats")
    #     if stat:
    #         stat.SetX1NDC(0.18); stat.SetY1NDC(0.60)
    #         stat.SetX2NDC(0.50); stat.SetY2NDC(0.80)
    #         c_cent.Modified(); c_cent.Update()
    # except Exception:
    #     pass

    # if display:
    #     c_cent.Update()
    # else:
    #     out_png = os.path.join(output, f"{base_tag}_calibration_ne.pdf")
    #     c_cent.SaveAs(out_png)
    #     print(f"[SAVE] {out_png}")

    # # --- Fit individual Gaussian peaks window-by-window (robust) ---
    # gauss = []
    # centers, e_centers = [], []   # ADU
    # sigma, e_sigma = [], []       # electrons
    # x_num_electron = []

    # fit_opts = "QNR"  # quiet, no draw; enforce fit Range

    # def fit_one_peak(i, hist, erange):
    #     """Return (ok, fgi) where ok=False if the fit is not reliable."""
    #     x_min = i - erange[0]
    #     x_max = i + erange[1]

    #     # Seed from local maximum in the window
    #     ax   = hist.GetXaxis()
    #     bmin = max(1, ax.FindBin(x_min))
    #     bmax = min(hist.GetNbinsX(), ax.FindBin(x_max))
    #     maxbin, maxval = bmin, 0.0
    #     for b in range(bmin, bmax + 1):
    #         v = hist.GetBinContent(b)
    #         if v > maxval:
    #             maxval, maxbin = v, b
    #     mu_seed  = ax.GetBinCenter(maxbin)
    #     amp_seed = maxval
    #     sig_seed = 0.15  # electrons (tune if needed)

    #     # If there are basically no counts, bail early
    #     if amp_seed < 5:
    #         return False, None

    #     fgi = ROOT.TF1(f"ne_{i}", "gaus", x_min, x_max)
    #     fgi.SetParameters(amp_seed, mu_seed, sig_seed)
    #     fgi.SetParLimits(1, x_min, x_max)   # mean constrained to the window
    #     fgi.SetParLimits(2, 0.02, 0.60)     # sigma (e-)
    #     fgi.SetLineColor(ROOT.kRed + 1)

    #     hist.Fit(fgi, fit_opts)

    #     # Sanity checks: mean near intended i, reasonable errors & width
    #     mu_fit  = fgi.GetParameter(1)
    #     emu_fit = fgi.GetParError(1)
    #     amp_fit = fgi.GetParameter(0)
    #     sig_fit = fgi.GetParameter(2)

    #     ok = (
    #         amp_fit > 0 and
    #         0.015 <= sig_fit <= 0.8 and
    #         abs(mu_fit - i) <= 0.45 and
    #         emu_fit <= 0.30
    #     )
    #     return ok, fgi

    # for i in range(ne_ini, ne_max + 1):
    #     ok, fgi = fit_one_peak(i, hist_e, erange)
    #     if not ok:
    #         # stop at first failure to avoid cascading garbage (or 'continue' if you prefer)
    #         break
    #     centers.append(gain * fgi.GetParameter(1))     # ADU
    #     e_centers.append(gain * fgi.GetParError(1))    # ADU
    #     sigma.append(fgi.GetParameter(2))              # e-
    #     e_sigma.append(fgi.GetParError(2))             # e-
    #     x_num_electron.append(i)
    #     gauss.append(fgi)

    # # Canvas: spectrum with fitted Gaussians
    # c_spec = ROOT.TCanvas("c_spec", "Spectrum with Gaussian fits", 900, 700)
    # hist_e.Draw("HIST")
    # for fgi in gauss:
    #     fgi.Draw("SAME")
    # c_spec.SetLogy()
    # if display:
    #     c_spec.Update()
    # else:
    #     out_png = os.path.join(output, f"{base_tag}_FittedGauss.pdf")
    #     c_spec.SaveAs(out_png)
    #     print(f"[SAVE] {out_png}")

    # # --- Calibration: peak centers vs electron number (ADU vs e-) ---
    # c_cent = ROOT.TCanvas("c_cent", "Centers vs n_e", 900, 700)
    # ROOT.gStyle.SetOptFit(1111)

    # x  = x_num_electron
    # y  = centers           # ADU
    # ey = e_centers         # ADU
    # ex = [0.0] * len(x)

    # g_cent = ROOT.TGraphErrors(
    #     len(y),
    #     array('d', [float(xx) for xx in x]),
    #     array('d', [float(yy) for yy in y]),
    #     array('d', ex),
    #     array('d', [float(e) for e in ey]),
    # )
    # g_cent.SetMarkerStyle(21)
    # g_cent.SetMarkerSize(0.8)
    # g_cent.Draw("AP E1")
    # c_cent.Update()

    # hframe = g_cent.GetHistogram()
    # if hframe:
    #     hframe.GetXaxis().SetTitle("n_{e^{-}}")
    #     hframe.GetYaxis().SetTitle("#mu_{n_{e^{-}}} [ADU]")  # FIX unit label

    # # Fit over the actually-found range, with generous slope freedom
    # if len(x) >= 2:
    #     xmin_fit, xmax_fit = min(x), max(x)
    #     pol1 = ROOT.TF1("fitfunc_cent", "pol1", xmin_fit, xmax_fit)
    #     pol1.SetLineColor(ROOT.kRed + 1)
    #     pol1.SetLineStyle(9)
    #     pol1.SetParameters(0.0, max(1.0, float(gain)))     # seed intercept=0, slope≈gain
    #     # either remove limits or set them wide:
    #     # pol1.SetParLimits(1, 0, 100*max(1.0, float(gain)))
    #     g_cent.Fit(pol1, "EMR0")
    #     pol1.Draw("SAME")
    # c_cent.Update()

    # Move stats box if present
    try:
        stat = c_cent.GetPrimitive("stats")
        if stat:
            stat.SetX1NDC(0.18); stat.SetY1NDC(0.60)
            stat.SetX2NDC(0.50); stat.SetY2NDC(0.80)
            c_cent.Modified(); c_cent.Update()
    except Exception:
        pass

    if display:
        c_cent.Update()
    else:
        out_png = os.path.join(output, f"{base_tag}_calibration_ne.pdf")
        c_cent.SaveAs(out_png)
        print(f"[SAVE] {out_png}")

    # --- Non-linearity residual: (Reconstructed charge - Peak n.) vs Reconstructed charge (e-) ---
    reco_e = centers_e[:]
    resid_e = [reco_e[i] - x[i] for i in range(len(x))]
    e_resid = e_centers_e[:]
    c_res = ROOT.TCanvas("c_res", "Non-linearity residual", 900, 700)
    g_res = ROOT.TGraphErrors(
        len(x),
        array('d', [float(r) for r in reco_e]),
        array('d', [float(v) for v in resid_e]),
        array('d', [0.0] * len(x)),
        array('d', [float(e) for e in e_resid]),
    )
    g_res.SetMarkerStyle(21)
    g_res.SetMarkerSize(0.8)
    g_res.Draw("AP E1")
    hres = g_res.GetHistogram()
    if hres:
        hres.GetXaxis().SetTitle("Reconstructed pixel charge (e^{-})")
        hres.GetYaxis().SetTitle("(Reconstructed charge - Peak n.) (e^{-})")
        hres.SetTitle("Non-linearity: residual vs reconstructed charge")
    # Parabola fit up to residual_fit_max; optionally extend plot and curve to extrapolate_max (skip if --no-residual-fit)
    x_fit_max = min(residual_fit_max, max(reco_e)) if reco_e else residual_fit_max
    # Always make room on the x-axis to see the extrapolated line up to ~2000 e-
    x_plot_max = max(max(reco_e) if reco_e else 0, 2000.0)
    if extrapolate_max > 0:
        x_plot_max = max(x_plot_max, extrapolate_max)
    def _make_residual_fit_tf1(name, model, xmin, xmax):
        # Use t=x/1000 for stability where requested
        if model == "pol2":
            f = ROOT.TF1(name, "pol2", xmin, xmax)
        elif model == "pol3":
            f = ROOT.TF1(name, "pol3", xmin, xmax)
        elif model == "scaledpol2":
            f = ROOT.TF1(name, "[0] + [1]*(x/1000.0) + [2]*(x/1000.0)*(x/1000.0)", xmin, xmax)
            f.SetParNames("a", "b", "c")
        elif model == "scaledpol3":
            f = ROOT.TF1(name, "[0] + [1]*(x/1000.0) + [2]*(x/1000.0)*(x/1000.0) + [3]*(x/1000.0)*(x/1000.0)*(x/1000.0)", xmin, xmax)
            f.SetParNames("a", "b", "c", "d")
        elif model == "rational":
            # (a + b t + c t^2) / (1 + d t), with t=x/1000
            f = ROOT.TF1(name, "([0] + [1]*(x/1000.0) + [2]*(x/1000.0)*(x/1000.0)) / (1.0 + [3]*(x/1000.0))", xmin, xmax)
            f.SetParNames("a", "b", "c", "d")
            # mild defaults to help the minimizer
            f.SetParameters(0.0, 0.0, 0.0, 0.0)
        else:
            raise RuntimeError(f"Unknown residual fit model: {model}")
        return f

    f_parab = None
    if not no_residual_fit and x_fit_max > (min(reco_e) if reco_e else 0):
        x_fit_min = min(reco_e) if reco_e else 0
        f_parab = _make_residual_fit_tf1("f_res_fit", residual_fit_model, x_fit_min, x_plot_max)
        f_parab.SetLineColor(ROOT.kRed + 1)
        g_res.Fit(f_parab, "EMR0", "", x_fit_min, x_fit_max)
        f_parab.SetRange(x_fit_min, x_plot_max)
        f_parab.Draw("SAME")
    if hres:
        # For TGraph, extend the x-axis limits so the extrapolated parabola segment is visible.
        # Use graph axis limits instead of the histogram frame limits, which were built from the data range only.
        x_min_axis = min(reco_e) if reco_e else 0.0
        x_max_axis = x_plot_max
        # Optional user-specified x-range overrides the default
        if residual_xrange and len(residual_xrange) == 2:
            rx_min, rx_max = residual_xrange
            if rx_max > rx_min:
                x_min_axis, x_max_axis = rx_min, rx_max
        g_res.GetXaxis().SetLimits(x_min_axis, x_max_axis)
        # Optional user-specified y-range
        if residual_yrange and len(residual_yrange) == 2:
            ry_min, ry_max = residual_yrange
            if ry_max > ry_min:
                hres.GetYaxis().SetRangeUser(ry_min, ry_max)
        c_res.Modified()
        c_res.Update()
        # If we extended the curve to extrapolate_max, print its value on the residual axis
        if extrapolate_max > 0 and f_parab:
            r_ext = f_parab.Eval(extrapolate_max)
            print(f"[Residual extrapolation] At {extrapolate_max:.0f} e-: residual = {r_ext:.4f} e-")
    # Extrapolated residual at Fe-55 X-ray line (electron-equivalent) energies
    if f_parab:
        e_Fe55_Ka = keV_to_electrons(FE55_KA_KEV)
        e_Fe55_Kb = keV_to_electrons(FE55_KB_KEV)
        r_Ka = f_parab.Eval(e_Fe55_Ka)
        r_Kb = f_parab.Eval(e_Fe55_Kb)
        print("")
        print("[Residual extrapolation at Fe-55 lines]")
        print(f"  Fe-55 K_alpha ({FE55_KA_KEV} keV = {e_Fe55_Ka:.2f} e-): residual = {r_Ka:.4f} e-")
        print(f"  Fe-55 K_beta  ({FE55_KB_KEV} keV = {e_Fe55_Kb:.2f} e-): residual = {r_Kb:.4f} e-")
        print("")
    c_res.Modified()
    c_res.Update()
    if display:
        c_res.Update()
    else:
        out_res = os.path.join(output, f"{base_tag}_nonlinearity_residual.pdf")
        c_res.SaveAs(out_res)
        print(f"[SAVE] {out_res}")

    # --- Normalized calibration: (μ - μ_0) / n_e  vs  n_e  (no fit params) ---
    c_norm = ROOT.TCanvas("c_norm", "Normalized calibration", 900, 700)

    # Find μ_0 from the data itself
    # Prefer exact n_e == 0; otherwise use the smallest available n_e as an approximation.
    if 0 in x:
        idx0 = x.index(0)
        mu0  = y[idx0]          # ADU at the 0-electron peak
    else:
        # fallback: use the smallest-n point
        imin = int(np.argmin(x))
        mu0  = y[imin]
        print(f"[norm] No n_e=0 peak in data; subtracting μ at n_e={x[imin]} instead (μ0≈{mu0:.3f} ADU).")

    x_norm, y_norm, ex_norm, ey_norm = [], [], [], []
    for xx, yy, ey_val in zip(x, y, ey):   # y and ey are centers/e_centers in ADU
        if xx > 0:  # avoid divide-by-zero
            x_norm.append(xx)
            y_norm.append( (yy - mu0) / xx )     # ADU/e-
            ex_norm.append(0.0)
            ey_norm.append( ey_val / xx )        # subtracting a constant doesn't change σ

    g_norm = ROOT.TGraphErrors(
        len(x_norm),
        array('d', [float(xx) for xx in x_norm]),
        array('d', [float(yy) for yy in y_norm]),
        array('d', ex_norm),
        array('d', [float(e) for e in ey_norm]),
    )

    g_norm.SetMarkerStyle(21)
    g_norm.SetMarkerSize(0.8)
    g_norm.Draw("AP E1")

    hnorm = g_norm.GetHistogram()
    if hnorm:
        hnorm.GetXaxis().SetTitle("n_{e^{-}}")
        hnorm.GetYaxis().SetTitle("(#mu - #mu_{0}) / n_{e^{-}}  [ADU/e^{-}]")
        hnorm.GetYaxis().SetTitleOffset(1.5)
        hnorm.SetTitle("Normalized calibration (data-derived #mu_{0})")

    c_norm.Modified(); c_norm.Update()

    if display:
        c_norm.Update()
    else:
        out_png = os.path.join(output, f"{base_tag}_calibration_norm.pdf")
        c_norm.SaveAs(out_png)
        print(f"[SAVE] {out_png}")

    # Fit (μ - μ0)/n = b + a/n to get the asymptotic gain b
    if x_norm:
        x_min_fit = max(5.0, min(x_norm))   # avoid 1/x blow-up
        x_max_fit = max(x_norm)
        f_bn = ROOT.TF1("f_bn", "[0] + [1]/x", x_min_fit, x_max_fit)
        f_bn.SetParNames("b_asym (ADU/e-)", "a")
        # seed b near the global slope from μ vs n
        f_bn.SetParameter(0, pol1.GetParameter(1) if 'pol1' in globals() else y_norm[-1])
        f_bn.SetParameter(1, 0.0)
        g_norm.Fit(f_bn, "QMR0")
        f_bn.Draw("SAME")
        b_asym = f_bn.GetParameter(0); eb_asym = f_bn.GetParError(0)
        print(f"[normalized] asymptotic gain b = {b_asym:.6f} ± {eb_asym:.6f} ADU/e-")


    # # --- Normalized calibration: μ/n_e vs n_e ---
    # c_norm = ROOT.TCanvas("c_norm", "Normalized calibration", 900, 700)

    # x_norm, y_norm, ex_norm, ey_norm = [], [], [], []
    # for xx, yy, ey_val in zip(x, centers, e_centers):  # centers, e_centers in ADU
    #     if xx > 0:
    #         x_norm.append(xx)
    #         y_norm.append(yy / xx)          # ADU per electron
    #         ex_norm.append(0.0)
    #         ey_norm.append(ey_val / xx)     # propagate error

    # g_norm = ROOT.TGraphErrors(
    #     len(x_norm),
    #     array('d', [float(xx) for xx in x_norm]),
    #     array('d', [float(yy) for yy in y_norm]),
    #     array('d', ex_norm),
    #     array('d', [float(e) for e in ey_norm]),
    # )

    # g_norm.SetMarkerStyle(21)
    # g_norm.SetMarkerSize(0.8)
    # g_norm.Draw("AP E1")

    # hnorm = g_norm.GetHistogram()
    # if hnorm:
    #     hnorm.GetXaxis().SetTitle("n_{e^{-}}")
    #     hnorm.GetYaxis().SetTitle("#mu / n_{e^{-}}  [ADU/e^{-}]")
    #     hnorm.GetYaxis().SetTitleOffset(1.3)
    #     hnorm.SetTitle("Normalized calibration")

    # c_norm.Modified()
    # c_norm.Update()

    # if display:
    #     c_norm.Update()
    # else:
    #     out_png = os.path.join(output, f"{base_tag}_calibration_norm.pdf")
    #     c_norm.SaveAs(out_png)
    #     print(f"[SAVE] {out_png}")

    # --- Spacing between consecutive peaks: Δμ vs n_e ---
    if len(x) >= 2:
        # Build arrays
        x_ne, dmu, edmu = [], [], []
        for i in range(1, len(x)):
            x_ne.append(x[i-1])                          # left peak electron number
            d  = centers[i] - centers[i-1]               # ADU spacing
            ed = (e_centers[i]**2 + e_centers[i-1]**2)**0.5  # error in spacing
            dmu.append(d); edmu.append(ed)

        # Graph
        c_dmu = ROOT.TCanvas("c_dmu", "Consecutive peak spacing", 900, 700)
        g_dmu = ROOT.TGraphErrors(
            len(x_ne),
            array('d', [float(xx) for xx in x_ne]),
            array('d', [float(yy) for yy in dmu]),
            array('d', [0.0]*len(x_ne)),
            array('d', [float(e) for e in edmu]),
        )
        g_dmu.SetMarkerStyle(21)
        g_dmu.SetMarkerSize(0.8)
        g_dmu.Draw("AP E1")

        hdmu = g_dmu.GetHistogram()
        if hdmu:
            hdmu.GetXaxis().SetTitle("n_{e^{-}} (left peak index)")
            hdmu.GetYaxis().SetTitle("#Delta#mu_{i} = #mu_{i+1}-#mu_{i}  [ADU/e^{-}]")
            hdmu.SetTitle("Adjacent peak spacing (anchored to n_{e})")

        # Fit to a constant (expected ~ gain)
        f_const = ROOT.TF1("fitfunc_dmu_const", "pol0", x_ne[0], x_ne[-1])
        f_const.SetLineColor(ROOT.kBlue + 1)
        g_dmu.Fit(f_const, "EMR0")
        f_const.Draw("SAME")

        # Reference line at slope from main fit (b = ADU/e-)
        b_gain = pol1.GetParameter(1)
        ref = ROOT.TF1("gain_ref_line", f"{b_gain}", x_ne[0], x_ne[-1])
        ref.SetLineColor(ROOT.kRed + 1)
        ref.SetLineStyle(9)
        ref.Draw("SAME")

        # Print comparison
        dmu_const  = f_const.GetParameter(0)
        dmu_conste = f_const.GetParError(0)
        print(f"[Δμ] constant fit = {dmu_const:.6f} ± {dmu_conste:.6f} ADU "
            f"(expected gain ~ {b_gain:.6f} ADU/e-)")

        c_dmu.Modified(); c_dmu.Update()
        if display:
            c_dmu.Update()
        else:
            out_png = os.path.join(output, f"{base_tag}_delta_mu_leftpeak.pdf")
            c_dmu.SaveAs(out_png)
            print(f"[SAVE] {out_png}")

    #     # --- Spacing between consecutive peaks: Δμ vs n_e (should be ~ gain) ---
    # if len(x) >= 2:
    #     # Build arrays
    #     x_mid, dmu, edmu = [], [], []
    #     for i in range(1, len(x)):
    #         x_mid.append(0.5 * (x[i] + x[i-1]))                   # electron midpoint
    #         d  = centers[i] - centers[i-1]                         # ADU/e
    #         ed = (e_centers[i]**2 + e_centers[i-1]**2)**0.5        # ADU/e
    #         dmu.append(d); edmu.append(ed)

    #     # Graph
    #     c_dmu = ROOT.TCanvas("c_dmu", "Consecutive peak spacing", 900, 700)
    #     g_dmu = ROOT.TGraphErrors(
    #         len(x_mid),
    #         array('d', [float(xx) for xx in x_mid]),
    #         array('d', [float(yy) for yy in dmu]),
    #         array('d', [0.0]*len(x_mid)),
    #         array('d', [float(e) for e in edmu]),
    #     )
    #     g_dmu.SetMarkerStyle(21)
    #     g_dmu.SetMarkerSize(0.8)
    #     g_dmu.Draw("AP E1")

    #     hdmu = g_dmu.GetHistogram()
    #     if hdmu:
    #         hdmu.GetXaxis().SetTitle("n_{e^{-}} (midpoint between peaks)")
    #         hdmu.GetYaxis().SetTitle("#Delta#mu_{i} = #mu_{i}-#mu_{i-1}  [ADU/e^{-}]")
    #         hdmu.SetTitle("Adjacent peak spacing")

    #     # Fit to a constant (expected ~ gain)
    #     f_const = ROOT.TF1("fitfunc_dmu_const", "pol0", x_mid[0], x_mid[-1])
    #     f_const.SetLineColor(ROOT.kBlue + 1)
    #     g_dmu.Fit(f_const, "EMR0")
    #     f_const.Draw("SAME")

    #     # Reference line at slope from main fit (b = ADU/e-)
    #     b_gain = pol1.GetParameter(1)
    #     ref = ROOT.TF1("gain_ref_line", f"{b_gain}", x_mid[0], max(x_mid[-1], x[-1]))
    #     ref.SetLineColor(ROOT.kRed + 1)
    #     ref.SetLineStyle(9)
    #     ref.Draw("SAME")

    #     # Print comparison
    #     dmu_const  = f_const.GetParameter(0)
    #     dmu_conste = f_const.GetParError(0)
    #     print(f"[Δμ] constant fit = {dmu_const:.6f} ± {dmu_conste:.6f} ADU (expected gain ~ {b_gain:.6f} ADU/e-)")

    #     c_dmu.Modified(); c_dmu.Update()
    #     if display:
    #         c_dmu.Update()
    #     else:
    #         out_png = os.path.join(output, f"{base_tag}_delta_mu.pdf")
    #         c_dmu.SaveAs(out_png)
    #         print(f"[SAVE] {out_png}")

    #             # --- Residuals of spacing relative to a reference gain, WITH error bars ---
    #     # Choose reference: "slope" uses b from main fit; "spacing_fit" uses the blue pol0 on Δμ
    #     residual_ref = "slope"  # change to "spacing_fit" if you prefer that reference

    #     if residual_ref == "spacing_fit":
    #         G_ref     = f_const.GetParameter(0)      # ADU/e- (blue line)
    #         G_ref_err = f_const.GetParError(0)
    #         ref_label = "Δμ constant fit"
    #     else:
    #         G_ref     = pol1.GetParameter(1)         # ADU/e- (main slope)
    #         G_ref_err = pol1.GetParError(1)
    #         ref_label = "slope b from μ vs n_e"

    #     # Build residuals (%). Error bars from edmu; (optionally include G_ref_err).
    #     res_y, res_ey = [], []
    #     for d, ed in zip(dmu, edmu):
    #         # residual = (Δμ - G_ref)/G_ref * 100
    #         r = 100.0 * (d - G_ref) / G_ref
    #         # σ_resid ≈ (σ_Δμ / G_ref) * 100  (ignore σ_G_ref to keep it simple and local)
    #         er = 100.0 * (ed / max(1e-12, G_ref))
    #         res_y.append(r); res_ey.append(er)

    #     c_dmu_res = ROOT.TCanvas("c_dmu_res", "Δμ residuals", 900, 600)
    #     g_dmu_res = ROOT.TGraphErrors(
    #         len(x_mid),
    #         array('d', [float(xx) for xx in x_mid]),
    #         array('d', [float(rr) for rr in res_y]),
    #         array('d', [0.0]*len(x_mid)),
    #         array('d', [float(e) for e in res_ey]),
    #     )
    #     g_dmu_res.SetMarkerStyle(20)
    #     g_dmu_res.Draw("AP E1")
    #     g_dmu_res.SetTitle(f"Spacing residuals vs {ref_label}; n_{{e^-}} (midpoint); (#Delta#mu - G_{{ref}})/G_{{ref}} [%]")

    #     zero = ROOT.TF1("zero", "0.0", x_mid[0], x_mid[-1])
    #     zero.SetLineStyle(7); zero.SetLineColor(ROOT.kGray+2); zero.Draw("SAME")

    #     # Print a quick summary
    #     print(f"[Δμ residuals] Reference gain ({ref_label}): {G_ref:.6f} ± {G_ref_err:.6f} ADU/e-")
    #     print(f"[Δμ residuals] Mean residual = {np.mean(res_y):.3f}%   RMS = {np.std(res_y, ddof=1):.3f}%")

    #     c_dmu_res.Modified(); c_dmu_res.Update()
    #     if display:
    #         c_dmu_res.Update()
    #     else:
    #         out_png = os.path.join(output, f"{base_tag}_delta_mu_residuals.png")
    #         c_dmu_res.SaveAs(out_png)
    #         print(f"[SAVE] {out_png}")
    # else:
    #     print("[Δμ] Not enough peaks to compute consecutive spacings.")

    # --- Absolute adjacent spacing: Δμ (ADU) vs n_e(midpoint) ---
    if len(x) >= 2:
        x_mid, dmu, edmu = [], [], []
        for i in range(1, len(x)):
            x_mid.append(0.5 * (x[i] + x[i-1]))
            dm  = centers[i] - centers[i-1]                   # ADU
            edm = (e_centers[i]**2 + e_centers[i-1]**2)**0.5  # ADU
            dmu.append(dm); edmu.append(edm)

        # 1) Raw absolute spacing (ADU)
        c_abs = ROOT.TCanvas("c_abs", "Absolute adjacent spacing", 900, 700)
        g_abs = ROOT.TGraphErrors(
            len(x_mid),
            array('d', [float(xx) for xx in x_mid]),
            array('d', [float(yy) for yy in dmu]),
            array('d', [0.0]*len(x_mid)),
            array('d', [float(e) for e in edmu]),
        )
        g_abs.SetMarkerStyle(21); g_abs.SetMarkerSize(0.8)
        g_abs.Draw("AP E1")
        h_abs = g_abs.GetHistogram()
        if h_abs:
            h_abs.SetTitle("Adjacent peak spacing; n_{e^{-}} (midpoint); #Delta#mu_{i} [ADU]")

        # Optional: constant fit for the spacing itself
        f_abs_const = ROOT.TF1("f_abs_const", "pol0", x_mid[0], x_mid[-1])
        f_abs_const.SetLineColor(ROOT.kBlue+1)
        g_abs.Fit(f_abs_const, "QMR0")
        f_abs_const.Draw("SAME")

        # Reference line at global slope b (from μ vs n_e fit)
        b_gain = pol1.GetParameter(1)
        ref_line = ROOT.TF1("ref_gain_line", f"{b_gain}", x_mid[0], x_mid[-1])
        ref_line.SetLineColor(ROOT.kRed+1); ref_line.SetLineStyle(7)
        ref_line.Draw("SAME")

        c_abs.Modified(); c_abs.Update()

        # 2) Absolute deviation from reference gain: |Δμ - G_ref| (ADU)
        # Choose reference: "slope" (global b) or "spacing_fit" (blue constant above)
        use_ref = "slope"  # change to "spacing_fit" if you prefer that
        if use_ref == "spacing_fit":
            G_ref     = f_abs_const.GetParameter(0)
            G_ref_err = f_abs_const.GetParError(0)
            ref_name  = "Δμ constant fit"
        else:
            G_ref     = b_gain
            G_ref_err = pol1.GetParError(1)
            ref_name  = "global slope b"

        abs_dev  = [abs(dm - G_ref) for dm in dmu]
        e_absdev = edmu[:]  # same uncertainties in ADU

        c_dev = ROOT.TCanvas("c_absdev", "Absolute deviation |Δμ - G_ref|", 900, 700)
        g_dev = ROOT.TGraphErrors(
            len(x_mid),
            array('d', [float(xx) for xx in x_mid]),
            array('d', [float(yy) for yy in abs_dev]),
            array('d', [0.0]*len(x_mid)),
            array('d', [float(e) for e in e_absdev]),
        )
        g_dev.SetMarkerStyle(20); g_dev.SetMarkerSize(0.8)
        g_dev.Draw("AP E1")
        h_dev = g_dev.GetHistogram()
        if h_dev:
            h_dev.SetTitle(f"Absolute deviation from {ref_name}; n_{'{e^-}'} (midpoint); |#Delta#mu - G_{{ref}}| [ADU]")

        # Optional: fit a constant to |Δμ - G_ref| to summarize typical absolute deviation
        f_absdev_const = ROOT.TF1("f_absdev_const", "pol0", x_mid[0], x_mid[-1])
        f_absdev_const.SetLineColor(ROOT.kMagenta+2)
        g_dev.Fit(f_absdev_const, "QMR0")
        f_absdev_const.Draw("SAME")

        print(f"[abs Δμ] mean spacing (blue) = {f_abs_const.GetParameter(0):.6f} ± {f_abs_const.GetParError(0):.6f} ADU")
        print(f"[abs dev] reference gain = {G_ref:.6f} ± {G_ref_err:.6f} ADU/e-")
        print(f"[abs dev] mean |Δμ - G_ref| = {f_absdev_const.GetParameter(0):.6f} ± {f_absdev_const.GetParError(0):.6f} ADU")

        c_dev.Modified(); c_dev.Update()
    else:
        print("[abs Δμ] Not enough peaks to compute adjacent spacings.")




    # # --- Normalized calibration: (ADU / n_e) vs n_e ---
    # c_norm = ROOT.TCanvas("c_norm", "Normalized calibration", 900, 700)
    # ROOT.gStyle.SetOptFit(1111)

    # # Normalize by electron number
    # x_norm = []
    # y_norm = []
    # ex_norm = []
    # ey_norm = []
    # for i, (xx, yy, ey_val) in enumerate(zip(x, centers, e_centers)):
    #     if xx > 0:  # avoid division by zero at n_e=0
    #         x_norm.append(xx)
    #         y_norm.append(yy / xx)
    #         ex_norm.append(0.0)
    #         ey_norm.append(ey_val / xx)

    # g_norm = ROOT.TGraphErrors(
    #     len(x_norm),
    #     array('d', [float(xx) for xx in x_norm]),
    #     array('d', [float(yy) for yy in y_norm]),
    #     array('d', ex_norm),
    #     array('d', [float(e) for e in ey_norm]),
    # )
    # g_norm.SetMarkerStyle(21)
    # g_norm.SetMarkerSize(0.8)
    # g_norm.Draw("AP E1")

    # hnorm = g_norm.GetHistogram()
    # if hnorm:
    #     hnorm.GetXaxis().SetTitle("n_{e^{-}}")
    #     hnorm.GetYaxis().SetTitle("#mu / n_{e^{-}} [ADU / e^{-}]")
    #     hnorm.GetYaxis().SetTitleOffset(1.3)
    #     hnorm.SetTitle("Normalized calibration vs electron peak number")


    # # Fit to a constant (should be flat if perfectly linear)
    # # pol0 = ROOT.TF1("fitfunc_norm", "pol0", min(x_norm), max(x_norm))
    # # pol0.SetLineColor(ROOT.kRed + 1)
    # # g_norm.Fit(pol0, "EMR0")
    # # pol0.Draw("SAME")
    #     # --- Fit normalized data with b + a/x (not a constant) ---
    # # choose a safe lower bound to avoid 1/x blow-up; e.g., start at n_e >= 5
    # x_min_fit = max(5.0, min(x_norm) if x_norm else ne_ini)
    # x_max_fit = max(x_norm) if x_norm else ne_max

    # f_norm = ROOT.TF1("fitfunc_norm_ax", "[0] + [1]/x", x_min_fit, x_max_fit)
    # # Good initial values: [0] ~ slope b from the main (pol1) fit; [1] ~ intercept a
    # f_norm.SetParameter(0, pol1.GetParameter(1))  # ~ b (ADU/e-)
    # f_norm.SetParameter(1, pol1.GetParameter(0))  # ~ a (ADU)

    # f_norm.SetParNames("b (ADU/e-)", "a (ADU)")
    # f_norm.SetLineColor(ROOT.kRed + 1)
    # f_norm.SetLineStyle(1)

    # g_norm.Fit(f_norm, "EMR0")
    # f_norm.Draw("SAME")

    # # Report results in a useful way
    # b_asym    = f_norm.GetParameter(0)
    # eb_asym   = f_norm.GetParError(0)
    # a_offset  = f_norm.GetParameter(1)
    # ea_offset = f_norm.GetParError(1)
    # print(f"[norm fit] asymptotic gain b = {b_asym:.6f} ± {eb_asym:.6f} ADU/e-")
    # print(f"[norm fit] offset a = {a_offset:.3f} ± {ea_offset:.3f} ADU (from normalized fit)")


    # c_norm.Modified()
    # c_norm.Update()

    # if display:
    #     c_norm.Update()
    # else:
    #     out_png = os.path.join(output, f"{base_tag}_calibration_norm.pdf")
    #     c_norm.SaveAs(out_png)
    #     print(f"[SAVE] {out_png}")

    

    # # ---------- Extrapolate normalized curve to Fe-55 and annotate ----------
    # e_per_electron = 3.74  # eV / e-
    # E_targets_eV = {"K#alpha": 5895.0, "K#beta": 6490.0}
    # ne_targets   = {name: (E / e_per_electron) for name, E in E_targets_eV.items()}

    # # Use your linear calibration (pol1) to predict μ(ADU), then normalize by n_e
    # mu_pred_ADU   = {name: pol1.Eval(ne)        for name, ne in ne_targets.items()}
    # mu_pred_norm  = {name: mu_pred_ADU[name]/ne for name, ne in ne_targets.items()}  # ADU/e-

    # # Extend axes robustly for TGraph
    # xmax_needed = 1.05 * max(ne_targets.values())
    # current_ymin = g_norm.GetMinimum()
    # current_ymax = g_norm.GetMaximum()
    # ymax_needed  = 1.10 * max(current_ymax, max(mu_pred_norm.values()))
    # ymin_needed  = min(current_ymin, min(y_norm)) if y_norm else current_ymin

    # g_norm.GetXaxis().SetLimits(ne_ini, xmax_needed)
    # g_norm.SetMaximum(ymax_needed)
    # g_norm.SetMinimum(ymin_needed)
    # c_norm.Modified(); c_norm.Update()

    # # Draw the constant fit solid over data range, dashed to x-ray range
    # if x_norm:
    #     x_lo_data = min(x_norm)
    #     x_hi_data = max(x_norm)
    # else:
    #     x_lo_data = ne_ini
    #     x_hi_data = ne_max

    # # solid over data
    # pol0.SetRange(x_lo_data, x_hi_data)
    # pol0.SetLineStyle(1)
    # pol0.Draw("SAME")

    # # dashed extension
    # pol0_ex = ROOT.TF1("fitfunc_norm_ex", "pol0", x_hi_data, xmax_needed)
    # pol0_ex.SetParameter(0, pol0.GetParameter(0))
    # pol0_ex.SetLineColor(pol0.GetLineColor())
    # pol0_ex.SetLineStyle(2)
    # pol0_ex.Draw("SAME")

    # # Guides/markers/labels at Fe-55 points (in ADU/e-)
    # col_map = {"K#alpha": ROOT.kBlue + 1, "K#beta": ROOT.kGreen + 2}
    # for name, ne_val in ne_targets.items():
    #     y_val = mu_pred_norm[name]   # ADU/e-
    #     col   = col_map.get(name, ROOT.kBlue + 1)

    #     v = ROOT.TLine(ne_val, ymin_needed, ne_val, y_val)
    #     v.SetLineStyle(7); v.SetLineColor(col); v.Draw()

    #     h = ROOT.TLine(ne_ini, y_val, ne_val, y_val)
    #     h.SetLineStyle(7); h.SetLineColor(col); h.Draw()

    #     m = ROOT.TMarker(ne_val, y_val, 20)
    #     m.SetMarkerColor(col); m.SetMarkerSize(1.0); m.Draw()

    #     txt = ROOT.TLatex()
    #     txt.SetTextSize(0.030)
    #     txt.SetTextColor(col)
    #     txt.DrawLatex(ne_val * 0.70, y_val * 1.02,
    #                   f"{name}: n_{{e}}={ne_val:.1f},  #mu/n_{{e}}={y_val:.3f} ADU/e^-")

    # c_norm.Modified()
    # c_norm.Update()

    # if display:
    #     c_norm.Update()
    # else:
    #     out_png = os.path.join(output, f"{base_tag}_calibration_norm.pdf")
    #     c_norm.SaveAs(out_png)
    #     print(f"[SAVE] {out_png}")




    # --- Resolution: sigma vs electron number (sigma in electrons) ---
    c_sig = ROOT.TCanvas("c_sig", "Sigma vs n_e", 900, 700)

    y_sig = sigma
    ey_sig = e_sigma

    g_sig = ROOT.TGraphErrors(
        len(y_sig),
        array('d', [float(xx) for xx in x]),
        array('d', [float(yy) for yy in y_sig]),
        array('d', ex),
        array('d', [float(e) for e in ey_sig]),
    )
    g_sig.SetMarkerStyle(21)
    g_sig.SetMarkerSize(0.8)
    g_sig.Draw("AP E1")
    c_sig.Update()
    hsig = g_sig.GetHistogram()
    if hsig:
        hsig.GetXaxis().SetTitle("n_{e^{-}}")
        hsig.GetYaxis().SetTitle("#sigma_{n_{e^{-}}} [e^{-}]")

    # Simple linear fit (descriptive); consider sigma^2 model for physics
    pol1s = ROOT.TF1("fitfunc_sig", "pol1", ne_ini, ne_max)
    pol1s.SetLineColor(ROOT.kRed + 1)
    pol1s.SetLineStyle(9)
    pol1s.SetParameter(0, 0.0)
    pol1s.SetParLimits(1, 0.0, 10.0)  # slope in e-/e-
    pol1s.SetParameter(1, 0.1)

    g_sig.Fit(pol1s, "EMR0")
    pol1s.Draw("SAME")
    c_sig.Update()

    try:
        stat = c_sig.GetPrimitive("stats")
        if stat:
            stat.SetX1NDC(0.18); stat.SetY1NDC(0.60)
            stat.SetX2NDC(0.50); stat.SetY2NDC(0.80)
            c_sig.Modified(); c_sig.Update()
    except Exception:
        pass

    if display:
        c_sig.Update()
        input("Press <Enter> to close…")
    else:
        out_png = os.path.join(output, f"{base_tag}_sigma_ne.pdf")
        c_sig.SaveAs(out_png)
        print(f"[SAVE] {out_png}")


if __name__ == '__main__':
    main()
