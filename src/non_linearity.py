#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit Gaussian peaks to CCD charge spectra (from FITS images or ROOT clusters).

Usage examples:
  # From FITS (ADU -> electrons via --gain if --inadu is set)
  python fit_peaks.py -f '/path/to/*.fits' --inadu true -k 3.74 --ne 0 50 --de 0.05 --ext 1 -o out

  # From ROOT tree (expects TTree "clustersRec" with pixels_E in keV)
  python fit_peaks.py -r clusters.root -k 3.74 --ne 0 50 --de 0.05 -o out
"""

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
else:
    import ROOT  # type: ignore
    from array import array  # type: ignore
    import numpy as np  # type: ignore
    from glob import glob  # type: ignore
    from astropy.io import fits  # type: ignore


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
@click.option("-r", "--root", type=str,
              help="Input ROOT file with reconstructed clusters (TTree 'clustersRec').")
@click.option("-k", "--gain", type=float, required=True,
              help="Gain (ADU per e-). Used for calibration and plotting.")
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
@click.option("--ext", type=int, default=1, help="FITS extension for data.")
@click.option("-o", "--output", type=str, default='.', help="Output directory.")
def main(infile, inadu, root, gain, ne, de, display, cols, rows, erange, ext, output):
    # Basic checks
    if (infile is None) and (root is None):
        raise click.UsageError("Provide either --infile (FITS) or --root (ROOT).")
    if (infile is not None) and (root is not None):
        raise click.UsageError("Provide only one of --infile or --root, not both.")

    os.makedirs(output, exist_ok=True)

    # Batch mode handling
    ROOT.gROOT.SetBatch(not display)

    ne_ini, ne_max = ne

    # Build histogram (either from ROOT TTree draw or manual fill from FITS)
    hist_name = "hist_e"
    hist_e = None
    base_tag = "spectra"

    if infile is None:
        # From ROOT clusters: pixels_E in keV -> electrons via (keV*1000)/3.74
        f = ROOT.TFile.Open(root)
        if not f or f.IsZombie():
            raise RuntimeError(f"Failed to open ROOT file: {root}")

        clusters = f.Get("clustersRec")
        if clusters is None:
            raise RuntimeError("Could not find TTree 'clustersRec' in ROOT file.")

        nbins = int((ne_max - ne_ini) / de)
        draw_expr = f"pixels_E*1000/3.74>>{hist_name}({nbins},{ne_ini},{ne_max})"
        _ = clusters.Draw(draw_expr, "", "GOFF norm")
        hist_e = ROOT.gDirectory.Get(hist_name)
        base_tag = os.path.splitext(os.path.basename(root))[0]
    else:
        # From FITS images: convert to electrons if --inadu
        nbins = int((ne_max - ne_ini) / de)
        hist_e = ROOT.TH1D(hist_name, hist_name, nbins, ne_ini - 3, ne_max + 3)
        hist_e.Sumw2()
        get_histe(infile, hist_e, inadu, gain, cols=cols, rows=rows, ext=ext)
        globbed = glob(infile)
        base_tag = os.path.splitext(os.path.basename(globbed[0]))[0] if globbed else "fits_input"

    if not hist_e:
        raise RuntimeError("Failed to create/fill histogram.")

    # Style
    hist_e.GetXaxis().SetTitle("n_{e^{-}}")
    hist_e.GetYaxis().SetTitle("Counts")
    hist_e.SetDirectory(0)

    # --- Fit individual Gaussian peaks window-by-window ---
    gauss = []
    centers, e_centers = [], []   # in ADU
    sigma, e_sigma = [], []       # in electrons
    x_num_electron = []

    mu_running = ne_ini - 1
    fit_opts = "LEMRQ0"

    for i in range(ne_ini, ne_max + 1):
        x_min = (mu_running + 1) - erange[0]
        x_max = (mu_running + 1) + erange[1]

        fgi = ROOT.TF1(f"ne_{i}", "gaus", x_min, x_max)
        fgi.SetParLimits(1, x_min, x_max)   # mean
        fgi.SetParLimits(2, 0.01, 1.0)      # sigma in electrons
        fgi.SetLineColor(ROOT.kRed + 1)

        hist_e.Fit(fgi, fit_opts)
        # Retry once if the fit latched onto a wrong region
        if (fgi.GetParameter(1) - mu_running) > 2:
            fgi.SetParLimits(1, x_min, x_max)
            hist_e.Fit(fgi, fit_opts)

        mu_running = fgi.GetParameter(1)

        # Store: centers in ADU, sigma in electrons
        centers.append(gain * fgi.GetParameter(1))
        e_centers.append(gain * fgi.GetParError(1))
        sigma.append(fgi.GetParameter(2))
        e_sigma.append(fgi.GetParError(2))
        x_num_electron.append(i)
        gauss.append(fgi)

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

    # --- Calibration: peak centers vs electron number (ADU vs e-) ---
    c_cent = ROOT.TCanvas("c_cent", "Centers vs n_e", 900, 700)
    ROOT.gStyle.SetOptFit(1111)

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

    # ---------- Extend calibration to Fe-55 and annotate predicted ADU ----------
    e_per_electron = 3.74  # eV / e-
    E_targets_eV = {"K#alpha": 5895.0, "K#beta": 6490.0}
    ne_targets = {name: (E / e_per_electron) for name, E in E_targets_eV.items()}
    mu_pred_ADU = {name: pol1.Eval(ne) for name, ne in ne_targets.items()}

    xmax_needed = 1.05 * max(ne_targets.values())
    hframe = g_cent.GetHistogram()
    if hframe:
        hframe.GetXaxis().SetRangeUser(ne_ini, max(xmax_needed, hframe.GetXaxis().GetXmax()))
        ymax_needed = 1.10 * max(hframe.GetMaximum(), max(mu_pred_ADU.values()))
        hframe.SetMaximum(ymax_needed)

    # Draw the fitted line out to xmax_needed
    pol1.SetRange(ne_ini, xmax_needed)
    pol1.Draw("SAME")

    col_map = {"K#alpha": ROOT.kBlue + 1, "K#beta": ROOT.kGreen + 2}
    for name, ne_val in ne_targets.items():
        mu_val = mu_pred_ADU[name]
        col = col_map.get(name, ROOT.kBlue + 1)

        v = ROOT.TLine(ne_val, 0, ne_val, mu_val);  v.SetLineStyle(7); v.SetLineColor(col); v.Draw()
        h = ROOT.TLine(ne_ini, mu_val, ne_val, mu_val); h.SetLineStyle(7); h.SetLineColor(col); h.Draw()
        m = ROOT.TMarker(ne_val, mu_val, 20); m.SetMarkerColor(col); m.SetMarkerSize(1.0); m.Draw()

        txt = ROOT.TLatex()
        txt.SetTextSize(0.030)
        txt.DrawLatex(ne_val * 0.70, mu_val * 1.02,
                      f"{name}: n_{{e}}={ne_val:.1f},  #mu={mu_val:.1f} ADU")

    c_cent.Modified(); c_cent.Update()

    # ---------- Fe-55 guides, labels, and axis extension (robust for TGraph) ----------
    e_per_electron = 3.74  # eV / e-
    E_targets_eV = {"K#alpha": 5895.0, "K#beta": 6490.0}

    # n_e for each target and predicted ADU from the linear fit
    ne_targets  = {name: (E / e_per_electron) for name, E in E_targets_eV.items()}
    mu_pred_ADU = {name: pol1.Eval(ne)        for name, ne in ne_targets.items()}

    # Determine needed axis limits
    xmax_needed = 1.05 * max(ne_targets.values())
    # Keep current ymin; expand ymax to include predicted points
    # (If you prefer to keep the lower part tight, just don't set SetMinimum.)
    current_ymin = g_cent.GetMinimum()
    current_ymax = g_cent.GetMaximum()
    ymax_needed  = 1.10 * max(current_ymax, max(mu_pred_ADU.values()))

    # >>> IMPORTANT: use TGraph axis setters, not the frame's SetRangeUser <<<
    g_cent.GetXaxis().SetLimits(ne_ini, xmax_needed)  # extend X
    g_cent.SetMaximum(ymax_needed)                    # extend Y up
    # g_cent.SetMinimum(current_ymin)                # (optional) keep same bottom
    c_cent.Modified(); c_cent.Update()

    # Draw solid fit over the measured range [ne_ini, ne_max]
    pol1.SetRange(ne_ini, ne_max)
    pol1.SetLineStyle(1)   # solid
    pol1.Draw("SAME")

    # Draw dashed fit for the extrapolated segment (ne_max -> xmax_needed)
    pol1_ex = ROOT.TF1("fitfunc_cent_ex", "pol1", ne_max, xmax_needed)
    pol1_ex.SetParameters(pol1.GetParameter(0), pol1.GetParameter(1))
    pol1_ex.SetLineColor(pol1.GetLineColor())
    pol1_ex.SetLineStyle(2)  # dashed
    pol1_ex.Draw("SAME")

    # Guides, markers, labels
    col_map = {"K#alpha": ROOT.kBlue + 1, "K#beta": ROOT.kGreen + 2}
    for name, ne_val in ne_targets.items():
        mu_val = mu_pred_ADU[name]
        col = col_map.get(name, ROOT.kBlue + 1)

        v = ROOT.TLine(ne_val, current_ymin if current_ymin is not None else 0.0, ne_val, mu_val)
        v.SetLineStyle(7); v.SetLineColor(col); v.Draw()

        h = ROOT.TLine(ne_ini, mu_val, ne_val, mu_val)
        h.SetLineStyle(7); h.SetLineColor(col); h.Draw()

        m = ROOT.TMarker(ne_val, mu_val, 20)
        m.SetMarkerColor(col); m.SetMarkerSize(1.0); m.Draw()

        txt = ROOT.TLatex()
        txt.SetTextSize(0.030)
        txt.SetTextColor(col)
        txt.DrawLatex(ne_val * 0.70, mu_val * 1.02,
                      f"{name}: n_{{e}}={ne_val:.1f},  #mu={mu_val:.1f} ADU/e-")

    c_cent.Modified(); c_cent.Update()



    # Extract linear fit parameters and print with correct unit conversions
    a = pol1.GetParameter(0)   # ADU
    b = pol1.GetParameter(1)   # ADU/e-
    print(f"[poly1] Calibration: offset a = {a:.3f} ADU, gain b = {b:.6f} ADU/e-")
    adu_per_keV = b * (1000.0 / e_per_electron)    # ADU/keV
    eV_per_ADU  = e_per_electron / b               # eV/ADU
    print(f"[poly1] => {adu_per_keV:.3f} ADU/keV, {eV_per_ADU:.3f} eV/ADU")

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
    ROOT.gStyle.SetOptFit(1111)

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
