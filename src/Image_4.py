import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
from pathlib import Path
from scipy import stats
import sys
# Configuration constants
overscan_start = 6152
y_overscan_start = 1536
prescan = 8
running_window = 40
col_threshold = 6
cti_window = 100

def robust_noise_from_overscan(sub, overscan_start, clip_sigma=5.0, max_iter=2):
    """
    Robust RMS from the X-overscan of a pedestal-subtracted image.
    Centers per row, then iteratively MAD-clips outliers (tracks, spikes).
    """
    os_sub = sub[:, overscan_start:]                           # overscan region
    res = os_sub - np.nanmedian(os_sub, axis=1)[:, None]       # row-centered

    mask = np.isfinite(res)
    for _ in range(max_iter):
        mad = 1.4826 * np.nanmedian(np.abs(np.where(mask, res, np.nan)))
        mad = max(float(mad), 1e-8)
        new_mask = np.abs(res) <= (clip_sigma * mad)
        if new_mask.sum() == mask.sum():
            break
        mask = new_mask

    return float(np.nanstd(np.where(mask, res, np.nan)))

def subtract_pedestal(data, clip_sigma=5.0, max_iter=2):
    """
    Robust per-row pedestal using X-overscan with iterative MAD clipping,
    plus robust noise estimate (printed) from the overscan after subtraction.
    """
    os = data[:, overscan_start:]
    row_med = np.nanmedian(os, axis=1)

    for _ in range(max_iter):
        res = os - row_med[:, None]
        # use a global MAD for stability against sparse tracks
        mad = 1.4826 * np.nanmedian(np.abs(res))
        mad = max(float(mad), 1e-8)
        keep = np.abs(res) <= (clip_sigma * mad)
        os_kept = np.where(keep, os, np.nan)
        new_row_med = np.nanmedian(os_kept, axis=1)
        row_med = np.where(np.isfinite(new_row_med), new_row_med, row_med)

    subtracted = data - row_med[:, None]

    noise = robust_noise_from_overscan(subtracted, overscan_start,
                                       clip_sigma=clip_sigma, max_iter=max_iter)
    print("Noise [ADU]:", round(noise, 2))
    return subtracted

# def subtract_pedestal(data):
#     row_medians = np.median(data[:, overscan_start:], axis=1)
#     subtracted = data - row_medians[:, np.newaxis]
#     noise_os = np.std(subtracted[:, overscan_start:])
#     mean_os = np.mean(subtracted[:, overscan_start:])
#     threshold = mean_os + 5 * noise_os
#     valid_pixels = subtracted[:, overscan_start:][subtracted[:, overscan_start:] < threshold]
#     print("Noise [ADU] : ", round(np.std(valid_pixels), 1))
#     return subtracted

def compute_median_images(fz_files):
    all_images = [[] for _ in range(4)]
    all_noises = [[] for _ in range(4)]

    for f in fz_files:
        with fits.open(f) as file:
            for channel in range(1, 5):
                raw = np.array(file[channel].data)
                print(f"Loaded {Path(f).name} channel {channel-1}")
                sub = subtract_pedestal(raw)
                all_images[channel - 1].append(sub)
                # noise_os = np.std(sub[:, overscan_start:])
                # mean_os = np.mean(sub[:, overscan_start:])
                # threshold = mean_os + 5 * noise_os
                # valid = sub[:, overscan_start:][sub[:, overscan_start:] < threshold]
                # all_noises[channel - 1].append(np.std(valid))
                all_noises[channel - 1].append(robust_noise_from_overscan(sub, overscan_start))


    median_images = [np.median(np.stack(imgs), axis=0) for imgs in all_images]
    median_noises = [round(np.median(noises), 1) for noises in all_noises]
    return median_images, median_noises

def running_median_absolute_deviation(arr, window_size):
    medians, mads = [], []
    arr = np.where(np.ma.getmask(arr), np.nan, arr)
    for i in range(len(arr)):
        start = max(0, i - window_size // 2)
        end = min(len(arr), i + window_size // 2 + 1)
        medians.append(np.nanmedian(arr[start:end]))
        mads.append(stats.median_abs_deviation(arr[start:end], nan_policy='omit'))
    return np.ma.masked_array(medians), np.ma.masked_array(mads)

def find_col_traps(image, threshold=col_threshold):
    masked = np.ma.masked_array(image)
    masked[:, :prescan] = np.ma.masked
    masked[:, overscan_start - prescan:] = np.ma.masked
    col_medians = np.ma.median(masked, axis=0)
    mask = np.ma.getmaskarray(col_medians)
    running_medians, running_mads = running_median_absolute_deviation(col_medians, window_size=running_window)
    high = running_medians + threshold * running_mads
    low = running_medians - threshold * running_mads
    traps = ~mask & ((col_medians > high) | (col_medians < low))
    return high, low, traps, np.where(traps)[0]

def calculate_cti_x(col_medians):
    pixels = col_medians[overscan_start - cti_window:overscan_start]
    overscan_val = col_medians[overscan_start]
    overscan = col_medians[overscan_start + cti_window:]
    px_mean = np.mean(pixels)
    px_std = np.std(pixels)
    os_mean = np.mean(overscan)
    os_std = np.std(overscan)
    val = px_mean - os_mean
    unc_val = np.sqrt(px_std**2 + os_std**2)
    os_step = overscan_val - os_mean
    os_unc = np.sqrt(2) * os_std
    cti = os_step / val
    cti_unc = cti * np.sqrt((unc_val / val)**2 + (os_unc / os_step)**2)
    return cti * 100, cti_unc * 100

def calculate_cti_y(row_medians):
    pixels = row_medians[y_overscan_start - cti_window:y_overscan_start]
    overscan_val = row_medians[y_overscan_start + 1]
    overscan = row_medians[y_overscan_start + 1:]
    px_mean = np.mean(pixels)
    px_std = np.std(pixels)
    os_mean = np.mean(overscan)
    os_std = np.std(overscan)
    val = px_mean - os_mean
    unc_val = np.sqrt(px_std**2 + os_std**2)
    os_step = overscan_val - os_mean
    os_unc = np.sqrt(2) * os_std
    cti = os_step / val
    cti_unc = cti * np.sqrt((unc_val / val)**2 + (os_unc / os_step)**2)
    return cti * 100, cti_unc * 100

def compute_display_vrange(median_images, prescan, overscan_start,
                           pct_low=0.5, pct_high=99.5):
    """
    Return (vmin, vmax) as robust percentiles over the data region of the
    four composite images (columns [prescan:overscan_start]).
    Floors: vmin >= 0, vmax >= 10.
    """
    data_pixels = []
    for img in median_images:
        roi = img[:, prescan:overscan_start]   # data region only
        data_pixels.append(roi.ravel())

    if not data_pixels:
        return 0.0, 10.0

    all_pixels = np.concatenate(data_pixels)
    # Keep only finite values
    all_pixels = all_pixels[np.isfinite(all_pixels)]
    if all_pixels.size == 0:
        return 0.0, 10.0

    vmin = np.nanpercentile(all_pixels, pct_low)
    vmax = np.nanpercentile(all_pixels, pct_high)

    # Fallbacks in degenerate cases
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin = float(np.nanmin(all_pixels))
        vmax = float(np.nanmax(all_pixels))

    # Clamp floors
    vmin = max(0.0, float(vmin))
    vmax = max(10.0, float(vmax))

    # Ensure non-degenerate range
    if vmax <= vmin:
        vmax = max(vmin + 1.0, 10.0)

    return float(vmin), float(vmax)


def _finite_vals(x):
    # Extract finite, unmasked values from ndarray or masked array
    if np.ma.isMaskedArray(x):
        return x.compressed()
    v = np.asarray(x).ravel()
    return v[np.isfinite(v)]

# def compute_colpanel_yrange(median_images, prescan, overscan_start,
#                             pct_low=0.5, pct_high=99.5, pad_frac=0.05):
#     """
#     Compute common (ymin, ymax) for the 'Column Median' panel by scanning
#     all 4 composite images. Uses robust percentiles over:
#       - column medians in data region (excludes prescan & overscan)
#       - running thresholds (Median ± N*MAD) so curves aren't clipped.

#     Final bounds are clamped so ymin >= 0 and ymax >= 10.
#     """
#     vals = []
#     for img in median_images:
#         # Same masking as find_col_traps to stay consistent
#         masked = np.ma.masked_array(img)
#         masked[:, :prescan] = np.ma.masked
#         masked[:, overscan_start - prescan:] = np.ma.masked

#         col_medians = np.ma.median(masked, axis=0)
#         # Reuse your existing logic to get high/low curves
#         high, low, _, _ = find_col_traps(img)

#         vals.append(_finite_vals(col_medians))
#         vals.append(_finite_vals(high))
#         vals.append(_finite_vals(low))

#     all_vals = np.concatenate([v for v in vals if v.size > 0]) if vals else np.array([0.0, 1.0])

#     vmin = np.nanpercentile(all_vals, pct_low)
#     vmax = np.nanpercentile(all_vals, pct_high)
#     if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
#         vmin, vmax = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))

#     # Add a small padding so curves don’t touch the frame
#     span = (vmax - vmin) if np.isfinite(vmax - vmin) and (vmax - vmin) > 0 else 1.0
#     pad = pad_frac * span
#     y0 = vmin - pad
#     y1 = vmax + pad

#     # Clamp floors: ymin >= 0, ymax >= 10
#     y0 = max(0.0, y0)
#     y1 = max(10.0, y1)

#     # Ensure non-degenerate range
#     if y1 <= y0:
#         y1 = y0 + max(1.0, 10.0)

#     return float(y0), float(y1)

def compute_colpanel_yrange(median_images, prescan, overscan_start,
                            pct_low=0.5, pct_high=99.5, pad_frac=0.05):
    """
    Compute common (ymin, ymax) for the 'Column Median' panel by scanning
    all 4 composite images. Uses robust percentiles over:
      - column medians in data region (excludes prescan & overscan)
      - running thresholds (Median ± N*MAD) so curves aren't clipped.
    """
    vals = []
    for img in median_images:
        # Same masking as find_col_traps to stay consistent
        masked = np.ma.masked_array(img)
        masked[:, :prescan] = np.ma.masked
        masked[:, overscan_start - prescan:] = np.ma.masked

        col_medians = np.ma.median(masked, axis=0)
        # Reuse your existing logic to get high/low curves
        high, low, _, _ = find_col_traps(img)

        vals.append(_finite_vals(col_medians))
        vals.append(_finite_vals(high))
        vals.append(_finite_vals(low))

    all_vals = np.concatenate([v for v in vals if v.size > 0]) if vals else np.array([0.0, 1.0])

    vmin = np.nanpercentile(all_vals, pct_low)
    vmax = np.nanpercentile(all_vals, pct_high)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin, vmax = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))

    # Add a small padding so curves don’t touch the frame
    pad = pad_frac * (vmax - vmin if np.isfinite(vmax - vmin) and (vmax - vmin) > 0 else 1.0)
    return float(vmin - pad), float(vmax + pad)




def save_image_detailed(images, ccd_image_name, log_file,vmin, vmax, colpanel_ylim, output_dir = '.'):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    axs = axs.flatten()


    for i, (image, channel, noise) in enumerate(images):
        ax = axs[i]
        # im = ax.imshow(image, cmap='cividis', origin='lower', vmin=-1, vmax=100, aspect='auto')
        im = ax.imshow(image, cmap='cividis', origin='lower', vmin=vmin, vmax=vmax, aspect='auto')
        ax.set_title(f'{ccd_image_name} – Channel {channel - 1}', fontsize=14, fontweight="bold")
        ax.set_xlabel("Column Number", fontsize=12)
        ax.set_ylabel("Row Number", fontsize=12)

        divider = make_axes_locatable(ax)
        axbottom = divider.append_axes("bottom", size="30%", pad=0.8, sharex=ax)
        axright = divider.append_axes("right", size="15%", pad=0.8, sharey=ax)
        cax = divider.append_axes("right", size="5%", pad=0.3)

        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label("Signal [ADU]", fontsize=12)

        col_medians = np.ma.median(image, axis=0)
        row_medians = np.ma.median(image, axis=1)
        col_medians[:prescan] = np.nan

        condition_high, condition_low, _, trap_locations = find_col_traps(image)

        # axbottom.step(np.arange(len(col_medians)), col_medians, where='mid', label='Median')
        # axbottom.step(np.arange(len(condition_high)), condition_high, color='red', linestyle='dashed', label='Median + 6 MAD')
        # axbottom.step(np.arange(len(condition_low)), condition_low, color='red', linestyle='dotted', label='Median - 6 MAD')
        # axbottom.set_ylabel("ADU", fontsize=10)
        # axbottom.set_xlabel("Column Number", fontsize=10)
        # axbottom.set_title("Column Median", fontsize=10)

        # axbottom.legend(loc='lower left', bbox_to_anchor=(1.0, 0.5), borderaxespad=0., fontsize=8, frameon=True)

        # --- bottom panel (shared y-limits across channels) ---
        x_idx = np.arange(len(col_medians))
        axbottom.step(x_idx, col_medians, where='mid', label='Median')
        axbottom.step(np.arange(len(condition_high)), condition_high, color='red',
                      linestyle='dashed', label='Median + 6 MAD')
        axbottom.step(np.arange(len(condition_low)),  condition_low,  color='red',
                      linestyle='dotted', label='Median - 6 MAD')
        axbottom.set_ylabel("ADU", fontsize=10)
        axbottom.set_xlabel("Column Number", fontsize=10)
        axbottom.set_title("Column Median", fontsize=10)

        # >>> apply common y-limits here <<<
        axbottom.set_ylim(colpanel_ylim)

        axbottom.legend(loc='lower left', bbox_to_anchor=(1.0, 0.5),
                        borderaxespad=0., fontsize=8, frameon=True)


        axright.step(row_medians, np.arange(image.shape[0]), where='mid')
        axright.set_xlabel("ADU", fontsize=10)
        axright.set_title("Row Median", fontsize=10)

        cti_x, cti_x_err = calculate_cti_x(col_medians)
        cti_y, cti_y_err = calculate_cti_y(row_medians)

        log_file.write("-" * 72 + "\n")
        log_file.write(f"---------------------   {ccd_image_name} - Channel : {i}   -------------------------\n")
        log_file.write(f"Median Noise : {noise} ADU\n")
        log_file.write(f"Column Defects Location: {trap_locations.tolist()}\n")
        log_file.write(f"Number of Column Defects: {len(trap_locations)}\n")
        log_file.write(f"CTI X [%]: {cti_x:.2f} ± {cti_x_err:.2f}\n")
        log_file.write(f"CTI Y [%]: {cti_y:.2f} ± {cti_y_err:.2f}\n")

    fig.suptitle(f'Column Defect Maps – {ccd_image_name}', fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / f"{ccd_image_name}.png", dpi=300)
    plt.close()

def main():
    # argparse interface
    parser = argparse.ArgumentParser(description="Compose CCD median images from multiple FZ files")
    parser.add_argument("--files", nargs="+", required=True, help="Input .fz files (space-separated)")
    parser.add_argument("--output_dir", required=True, help="Path for output directory.")
    parser.add_argument("--module", required=True, help="Module name (used for outputs: <module>.png and <module>.log)")
    args = parser.parse_args()

    fz_files = [str(Path(f)) for f in args.files if Path(f).suffix == ".fz"]

    if not fz_files:
        print("[ERROR] No .fz files provided or matched.")
        sys.exit(1)

    # check existence
    missing = [f for f in fz_files if not Path(f).exists()]
    if missing:
        for f in missing:
            print(f"[ERROR] File not found: {f}")
        sys.exit(1)

    ccd_image_name = args.module

    output_dir = Path(args.output_dir)
    median_images, median_noises = compute_median_images(fz_files)
    images = [(median_images[i], i + 1, median_noises[i]) for i in range(4)]
    # NEW: dynamic display range over the 4 composites
    vmin, vmax = compute_display_vrange(median_images, prescan, overscan_start,
                                        pct_low=0.5, pct_high=99.5)
    
    # NEW: common y-limits for the bottom 'Column Median' panels
    ymin_cm, ymax_cm = compute_colpanel_yrange(median_images, prescan, overscan_start,
                                           pct_low=0.5, pct_high=99.5)


    log_file_path = Path(output_dir / f"{ccd_image_name}.log")
    with log_file_path.open('w') as log_file:
        save_image_detailed(images, ccd_image_name, log_file, vmin, vmax, (ymin_cm, ymax_cm), output_dir=output_dir)

   

    print(f"[DONE] Saved: {ccd_image_name}.png and {ccd_image_name}.log")

if __name__ == "__main__":
    main()
