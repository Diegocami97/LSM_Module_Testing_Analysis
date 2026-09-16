#!/usr/bin/env python3
# Diego Venegas-Vargas
# DAMIC-M Collaboration
# Builds one FITS image per skip from an L1 CSV via CDS (ramp_1 - ramp_0).
"""Build one FITS image per skip from an L1 CSV: each image is 480×1280, pixel = ramp_1 - ramp_0 (CDS).

Reads the CSV in chunks. For each skip, fills pedestal (ramp=0) and signal (ramp=1)
within 480 rows × 1280 columns (pixel index rolls over per image), then outputs
one image per skip with value = signal - pedestal. Saves as .fz. First 10 skips by default.

Usage:
  python csv_l1_to_fits_per_skip.py input_L1.csv [-o outdir] [--header reference.fits] [--max 10]
"""

import argparse
import os
import shutil
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits


# Image size per skip: 480 rows x 1280 columns = 614400 pixels per frame.
# For each skip we build one frame (614400 pixels). Line and pix reset per skip:
# within each (skip, ramp), line runs 0..479 and pix runs 0..614399 (unique, no repeat).
# Column in image = pix % NCOL; row = line.
NROW = 480
NCOL = 1280


def get_dtype_for_sum():
    """L1 'sum' is int32; CDS (ramp1 - ramp0) can be int32."""
    return np.int32


def report_dimensions(csv_path, chunksize=int(5e6)):
    """Scan CSV once: report max line/pix, distinct counts, and whether (line, pix) is unique per (skip, ramp)."""
    max_line = -1
    max_pix = -1
    lines_seen = set()
    pix_seen = set()
    skip_ramps = set()
    first_skip, first_ramp = None, None
    rows_first = 0
    lp_first = set()

    for chunk in pd.read_csv(csv_path, chunksize=chunksize, dtype={"pix": np.uint32, "line": np.uint32, "skip": np.uint16, "ramp": np.uint8, "idx": np.uint8, "n": np.uint32, "sum": np.int32}):
        if len(chunk) == 0:
            continue
        max_line = max(max_line, int(chunk["line"].max()))
        max_pix = max(max_pix, int(chunk["pix"].max()))
        lines_seen.update(chunk["line"].unique())
        pix_seen.update(chunk["pix"].unique())
        for _, row in chunk[["skip", "ramp"]].drop_duplicates().iterrows():
            skip_ramps.add((int(row["skip"]), int(row["ramp"])))

        if first_skip is None and len(skip_ramps) > 0:
            first_skip, first_ramp = sorted(skip_ramps)[0]
        if first_skip is not None:
            m = (chunk["skip"] == first_skip) & (chunk["ramp"] == first_ramp)
            if m.any():
                sub = chunk.loc[m]
                rows_first += len(sub)
                for _, r in sub.iterrows():
                    lp_first.add((int(r["line"]), int(r["pix"])))

    nrow = max_line + 1
    ncol = max_pix + 1
    print("From CSV (one chunked pass):")
    print("  max(line)       =", max_line)
    print("  max(pix)        =", max_pix)
    print("  nrow (max+1)    =", nrow)
    print("  ncol (max+1)    =", ncol)
    print("  n_distinct_line =", len(lines_seen))
    print("  n_distinct_pix  =", len(pix_seen))
    print("  n (skip, ramp)  =", len(skip_ramps))
    if first_skip is not None and rows_first > 0:
        n_unique_lp = len(lp_first)
        print("  For first (skip, ramp)=({}, {}): rows={}, unique (line,pix)={}".format(first_skip, first_ramp, rows_first, n_unique_lp))
        if n_unique_lp == rows_first:
            print("  -> (line, pix) is UNIQUE per (skip, ramp); pix does not repeat within a frame.")
        else:
            print("  -> (line, pix) REPEATS (duplicate keys); last write wins in image.")
    print("Use these to set image shape (e.g. nrow x ncol) for building FITS.")
    return nrow, ncol, sorted(skip_ramps)


def first_pass_skip_list(csv_path, chunksize=int(5e6)):
    """Scan CSV to get sorted list of unique skip indices."""
    skips = set()
    for chunk in pd.read_csv(csv_path, chunksize=chunksize, dtype={"pix": np.uint32, "line": np.uint32, "skip": np.uint16, "ramp": np.uint8, "idx": np.uint8, "n": np.uint32, "sum": np.int32}):
        skips.update(chunk["skip"].unique().astype(int).tolist())
    return sorted(skips)


def build_one_cds_image(csv_path, skip, nrow=NROW, ncol=NCOL, chunksize=int(5e6)):
    """Build one 480×1280 CDS image for this skip.
    Raw ADC sums (sig, ped) are typically negative; decoder uses cds_avg = -1.0 * cds_sum/n
    so that more signal → more positive. We output image = -(sig - ped) = ped - sig."""
    ped = np.zeros((nrow, ncol), dtype=np.int64)
    sig = np.zeros((nrow, ncol), dtype=np.int64)
    for chunk in pd.read_csv(csv_path, chunksize=chunksize, dtype={"pix": np.uint32, "line": np.uint32, "skip": np.uint16, "ramp": np.uint8, "idx": np.uint8, "n": np.uint32, "sum": np.int32}):
        mask_skip = chunk["skip"] == skip
        if not mask_skip.any():
            continue
        in_bounds = (chunk["line"] < nrow) & (chunk["pix"] < nrow * ncol)
        for ramp_val, arr in [(0, ped), (1, sig)]:
            m = mask_skip & (chunk["ramp"] == ramp_val) & in_bounds
            if not m.any():
                continue
            sub = chunk.loc[m]
            line = sub["line"].values
            pix = sub["pix"].values
            col = np.mod(pix, ncol)
            sum_ = sub["sum"].values.astype(np.int64)
            arr[line, col] = sum_
    cds = ped - sig
    # cds = sig - ped
    return cds.astype(get_dtype_for_sum())


def load_header_from_fits(header_path):
    """Copy header from a reference FITS file (pr imary HDU)."""
    with fits.open(header_path) as hdul:
        return hdul[0].header.copy()


def write_one_fits(image, skip, out_dir, base_name, header=None, ndcm=None, compress_fz=True):
    """Write one CDS image (ramp1 - ramp0) per skip to a FITS file and optionally compress to .fz."""
    nrow, ncol = image.shape
    name = f"{base_name}_skip{skip:04d}.fits"
    path = os.path.join(out_dir, name)
    hdu = fits.PrimaryHDU(data=image)
    if header is not None:
        for k, v in header.items():
            if k in ("SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"):
                continue
            try:
                hdu.header[k] = v
            except Exception:
                pass
    hdu.header["NAXIS1"] = ncol
    hdu.header["NAXIS2"] = nrow
    hdu.header["SKIP"] = skip
    if header is None:
        hdu.header["BUNIT"] = ("CDS (ped-sig)", "pedestal minus signal sum")
    if ndcm is not None:
        hdu.header["NDCM"] = ndcm
    with warnings.catch_warnings():
        try:
            warnings.simplefilter("ignore", fits.verify.VerifyWarning)
        except AttributeError:
            pass
        hdu.writeto(path, overwrite=True, output_verify="silentfix")
    if compress_fz:
        fpack_cmd = shutil.which("fpack")
        if not fpack_cmd:
            for cand in ("/opt/homebrew/bin/fpack", "/usr/local/bin/fpack"):
                if os.path.isfile(cand) and os.access(cand, os.X_OK):
                    fpack_cmd = cand
                    break
        try:
            if not fpack_cmd:
                raise FileNotFoundError("fpack not found")
            subprocess.run([fpack_cmd, "-D", "-Y", path], check=True, capture_output=True, timeout=60)
            path_fits_fz = path + ".fz"  # fpack creates file.fits.fz
            path_fz = os.path.splitext(path)[0] + ".fz"  # we want file.fz
            os.rename(path_fits_fz, path_fz)
            os.remove(path)
            print("  Written", path_fz)
        except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print("  Written (uncompressed)", path)
            if isinstance(e, FileNotFoundError):
                print("    (fpack not found; install with: brew install cfitsio - or use --no-fz)")
    else:
        print("  Written", path)


def main():
    parser = argparse.ArgumentParser(description="Build one CDS image (480×1280, ramp1-ramp0) per skip from L1 CSV; save as .fz")
    parser.add_argument("csv", help="Path to L1 CSV (e.g. from run_decoder_l1.py)")
    parser.add_argument("-o", "--outdir", default=None, help="Output directory (default: same as CSV)")
    parser.add_argument("--header", default=None, help="Reference FITS file to copy header from (optional)")
    parser.add_argument("--ndcm", type=int, default=None, help="NDCM value for header (overrides reference; set per run)")
    parser.add_argument("--no-fz", action="store_true", help="Do not compress to .fz (keep .fits only)")
    parser.add_argument("--max", type=int, default=10, help="Max number of skips to output (default 10); ignored if --skips is set")
    parser.add_argument("--skips", type=int, nargs="+", default=None, metavar="N", help="Output only these skip numbers (e.g. --skips 10 100 200). Overrides --max.")
    parser.add_argument("--print-image", action="store_true", help="Print the 2D CDS image array to the screen (use with --max 1 for one image)")
    parser.add_argument("--flip-cds-after-first", action="store_true", help="Use opposite CDS sign for skip > 0 (sig - ped instead of ped - sig) so later skips are positive like the first")
    parser.add_argument("--report-dims", action="store_true", help="Only report max line/pix and distinct counts from CSV, then exit")
    parser.add_argument("--chunksize", type=int, default=int(5e6), help="CSV read chunk size (default 5e6)")
    args = parser.parse_args()

    csv_path = os.path.abspath(args.csv)
    if not os.path.isfile(csv_path):
        print("Error: CSV file not found:", csv_path, file=sys.stderr)
        sys.exit(1)

    out_dir = args.outdir
    if out_dir is None:
        out_dir = os.path.dirname(csv_path)
    os.makedirs(out_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(csv_path))[0]
    if base_name.endswith("_L1"):
        base_name = base_name[:-3]

    if args.report_dims:
        report_dimensions(csv_path, chunksize=args.chunksize)
        return

    header = None
    if args.header:
        if os.path.isfile(args.header):
            header = load_header_from_fits(args.header)
            print("Using header from", args.header)
        else:
            print("Warning: header file not found, using minimal header", file=sys.stderr)

    if args.skips is not None:
        skip_list = sorted(set(args.skips))
        print(f"  Image shape: {NROW} x {NCOL}. User-specified skips: {skip_list} ({len(skip_list)} images)")
    else:
        print("First pass: scanning CSV for skip list ...")
        skip_list = first_pass_skip_list(csv_path, chunksize=args.chunksize)
        skip_list = skip_list[: args.max]
        print(f"  Image shape: {NROW} x {NCOL} (CDS = ramp_1 - ramp_0), outputting first {len(skip_list)} skips (max={args.max})")

    if args.ndcm is not None:
        print("Using NDCM =", args.ndcm)
    if not skip_list:
        print("No skip indices found in CSV. Nothing to do.")
        return

    print("Building and writing ...")
    for i, skip in enumerate(skip_list):
        image = build_one_cds_image(csv_path, skip, nrow=NROW, ncol=NCOL, chunksize=args.chunksize)
        if args.print_image:
            n_filled = np.count_nonzero(image)
            n_total = image.size
            print(f"\n--- 2D CDS image (skip={skip}), shape {image.shape} ---")
            if n_filled == 0:
                print("  WARNING: image all zeros — check that CSV has both ramp=0 and ramp=1 for this skip.")
            with np.printoptions(threshold=np.inf, linewidth=200):
                print(image)
            print("--- end image ---\n")
        write_one_fits(image, skip, out_dir, base_name, header=header, ndcm=args.ndcm, compress_fz=not args.no_fz)
        del image
        if (i + 1) % 5 == 0 or (i + 1) == len(skip_list):
            print(f"  Progress: {i + 1}/{len(skip_list)}")
    print("Done.")


if __name__ == "__main__":
    main()
