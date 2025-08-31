import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from skimage import filters, measure, morphology
from skimage.transform import probabilistic_hough_line
from pathlib import Path
import argparse
import os

plt.rcParams.update({'font.size': 18})

overscan_start = 309
prescan = 8

# def subtract_pedestal(data):
#     row_medians = np.median(data[:, overscan_start:], axis=1)
#     pedestal_subtracted = data - row_medians[:, np.newaxis]
#     print("Noise [ADU]:", round(np.std(pedestal_subtracted[:, overscan_start:]), 1))
#     return pedestal_subtracted

def subtract_pedestal(data, clip_sigma=5.0, max_iter=2):
    """
    Robust per-row pedestal subtraction using X-overscan with MAD clipping.
    Avoids 2D->1D collapse by masking with NaNs (keeps axes intact).
    """
    # Overscan region in X
    os = data[:, overscan_start:]  # shape (rows, overscan_cols)

    # --- initial row medians from overscan ---
    row_med = np.nanmedian(os, axis=1)  # (rows,)

    for _ in range(max_iter):
        # residuals to current pedestal
        res = os - row_med[:, None]  # (rows, overscan_cols)
        # per-row robust sigma via MAD
        mad = 1.4826 * np.nanmedian(np.abs(res), axis=1)  # (rows,)
        mad = np.maximum(mad, 1e-8)  # avoid div by zero
        # mask outliers (keep within clip_sigma)
        keep = np.abs(res) <= (clip_sigma * mad[:, None])  # 2D bool
        # recompute per-row medians on kept pixels using NaN masking
        os_kept = np.where(keep, os, np.nan)
        new_row_med = np.nanmedian(os_kept, axis=1)
        # if a row loses all kept pixels, fall back to previous value
        row_med = np.where(np.isfinite(new_row_med), new_row_med, row_med)

    # subtract pedestal
    pedestal_subtracted = data - row_med[:, None]

    # --- compute noise from overscan after subtraction ---
    os_sub = pedestal_subtracted[:, overscan_start:]
    res2 = os_sub - np.nanmedian(os_sub, axis=1)[:, None]
    mad2 = 1.4826 * np.nanmedian(np.abs(res2), axis=1)
    mad2 = np.maximum(mad2, 1e-8)
    keep2 = np.abs(res2) <= (clip_sigma * mad2[:, None])
    # std over kept pixels (NaN-mask to keep 2D)
    noise = float(np.nanstd(np.where(keep2, res2, np.nan)))
    print("Noise [ADU] :", round(noise, 2))
    return pedestal_subtracted



def detect_diagonal_lines(img, threshold_abs=50, min_len=150):
    binary = img > threshold_abs
    binary = morphology.remove_small_objects(binary, min_size=20)
    edges = filters.sobel(binary)
    lines = probabilistic_hough_line(edges, threshold=10, line_length=min_len, line_gap=5)

    filtered_lines = []
    print("Detected line angles:")
    for p0, p1 in lines:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        angle = 90.0 if dx == 0 else np.degrees(np.arctan2(dy, dx))
        print(f"  angle = {angle:.1f}°")
        if 120 < abs(angle) < 135:  # keep diagonals
            filtered_lines.append((p0, p1))

    print(f"  → {len(filtered_lines)} kept as diagonals\n")
    return filtered_lines

def main():
    # ---------------- Argument parsing ----------------
    parser = argparse.ArgumentParser(description="CCD plot generator")
    parser.add_argument("--files", required=True, help="Input FITS/FZ file")
    parser.add_argument("--output_dir", required=True, help="Path for output directory.")
    parser.add_argument("--module", required=True, help="CCD Module name for output file")
    args = parser.parse_args()

    fits_file = args.files
    module_name = args.module

    hdul = fits.open(fits_file)
    images = []

    # Load first 4 extensions
    for ext in range(1, 5):
        data = np.array(hdul[ext].data)
        image = subtract_pedestal(data)
        images.append(image)
    hdul.close()

    # Find common vmin and vmax
    all_pixels = np.concatenate([img[:, prescan:overscan_start].flatten() for img in images])
    vmin = np.percentile(all_pixels, 0.5)
    vmax = np.percentile(all_pixels, 99.5)

    # Create single figure
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))
    axs = axs.flatten()
    reordered = [images[0], images[1], images[2], images[3]]
    titles = [1, 2, 3, 4]

    for i, (ax, img) in enumerate(zip(axs, reordered)):
        im = ax.imshow(img, cmap='cividis', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(module_name + f" - Channel {titles[i]-1}", fontsize=14, fontweight="bold")
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("Column Number", fontsize=12)
        ax.set_ylabel("Row Number", fontsize=12)

    # Reduce vertical space
    plt.subplots_adjust(hspace=-0.6, wspace=0.1, left=0.05, right=0.9, top=0.971, bottom=0.05)

    # Shared colorbar
    cbar_ax = fig.add_axes([0.92, 0.30, 0.02, 0.42])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Charge [ADU]', fontsize=18)
    cbar.ax.tick_params(labelsize=12)

    # Save with module name
    output_dir = Path(args.output_dir)   # <-- make it a Path
    out_png = output_dir / f"{module_name}.png"
    out_pdf = output_dir / f"{module_name}.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    # plt.show()

if __name__ == "__main__":
    main()
