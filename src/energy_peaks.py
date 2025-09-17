# use: python3 src/energy_peaks.py --id 19 --csvpath Data/MODULE_UNDERGROUND2.csv
# Plot energy peaks over the modules
# Cinyu Zhu, JHU, Sept 2025

# todo: move the script to ccdqc server and make it automatic

# ## contents
# - prepare a function to convert a string of fitted result to [center, error], and also for the case of one number/NULL
# - prepare a numpy array[module, ccd, image4-7, peak, content/error] (shape = 28 * 4 * 4 * 2 * 2)
# - prepare a function to query csv given a id + column name 
# - query the sql db according to the column name assembled from the above, given specific indiced
# - fill the numpy array with numbers
# - plot

### csv is obtained from exporting csv from phpMyAdmin (direct connection to the database is not allowed)
# ````
# SELECT id, Name, Image4_Low_Peak1_A, Image4_Low_Peak1_B, Image4_Low_Peak1_C, Image4_Low_Peak1_D, Image4_Low_Peak2_A, Image4_Low_Peak2_B, Image4_Low_Peak2_C, Image4_Low_Peak2_D, Image5_Low_Peak1_A, Image5_Low_Peak1_B, Image5_Low_Peak1_C, Image5_Low_Peak1_D, Image5_Low_Peak2_A, Image5_Low_Peak2_B, Image5_Low_Peak2_C, Image5_Low_Peak2_D, Image6_Low_Peak1_A, Image6_Low_Peak1_B, Image6_Low_Peak1_C, Image6_Low_Peak1_D, Image6_Low_Peak2_A, Image6_Low_Peak2_B, Image6_Low_Peak2_C, Image6_Low_Peak2_D, Image7_Low_Peak1_A, Image7_Low_Peak1_B, Image7_Low_Peak1_C, Image7_Low_Peak1_D, Image7_Low_Peak2_A, Image7_Low_Peak2_B, Image7_Low_Peak2_C, Image7_Low_Peak2_D
# FROM MODULE_UNDERGROUND2
# ````

import numpy as np
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_value_error(s: str) -> np.ndarray:
    """Parse 'x +/- y', single number, or empty → [x, y] with -1 fallbacks."""
    if s is None:
        return np.array([0, 0])
    if isinstance(s, float) and np.isnan(s):
        return np.array([0, 0])
    s = str(s).strip()
    if s == "":
        return np.array([0, 0])
    if "+/-" in s:
        try:
            left, right = s.split("+/-", 1)
            return np.array([float(left.strip()), float(right.strip())], dtype=float)
        except Exception:
            print(f"Error parsing string: {s}")
            return np.array([0, 0])
    try:
        return np.array([float(s), 0], dtype=float)
    except Exception:
        print(f"Error parsing string: {s}")
        return np.array([0, 0])
    
def build_array_from_patterned_csv(
    csv_path: str,
    images=(4,5,6,7),
    peaks=(1,2),
    ccd_letters=("A","B","C","D"),
    fill_value: float = 0,
) -> np.ndarray:
    """
    Build array with shape (28, 4, 4, 2, 2) ordered as:
    [module(0..27), ccd(0..3), image_idx(0..3 for 4..7), peak(0..1 for 1..2), content/error].
    CSV must have columns named like: 'Image{4..7}_Low_Peak{1|2}_{A|B|C|D}'.
    Also requires an 'id' column giving module numbers 1..28.
    """
    df = pd.read_csv(csv_path)

    if "id" not in df.columns:
        raise ValueError("CSV must contain a column named 'id' with values 1..28")

    # output array
    arr = np.full((28, 4, 4, 2, 2), fill_value, dtype=float)

    # precompute index maps
    img_to_axis = {img: i for i, img in enumerate(images)}
    ccd_to_axis = {letter: i for i, letter in enumerate(ccd_letters)}
    peak_to_axis = {p: i for i, p in enumerate(peaks)}

    # verify expected columns exist; warn if not
    expected_cols = []
    for img in images:
        for p in peaks:
            for c in ccd_letters:
                expected_cols.append(f"Image{img}_Low_Peak{p}_{c}")
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        print(f"Warning: {len(missing)} expected columns not found; examples: {missing[:5]}")

    # fill array
    for _, row in df.iterrows():
        try:
            mval = int(row["id"])
            if not (1 <= mval <= 28):
                continue
            m_idx = mval - 1
        except Exception:
            continue

        for img in images:
            i_idx = img_to_axis[img]
            for p in peaks:
                p_idx = peak_to_axis[p]
                for c in ccd_letters:
                    c_idx = ccd_to_axis[c]
                    col = f"Image{img}_Low_Peak{p}_{c}"
                    if col not in df.columns:
                        continue
                    parsed = parse_value_error(row[col])
                    arr[m_idx, c_idx, i_idx, p_idx, 0] = parsed[0]
                    arr[m_idx, c_idx, i_idx, p_idx, 1] = parsed[1]

    return arr

# ---------- plotting utility ----------
def plot_module_energy_peaks(
    arr: np.ndarray,
    module_number: int = 20,
    images=(4,5,6),
    peaks=(1,2),
):
    """
    Plot energy peaks with error bars for a given module.
    X-axis: ext1..ext4 (CCD A..D).
    Marker by image (4:o, 5:s, 6:^, 7:v),
    color by peak (1:red, 2:blue).
    Each image is horizontally offset for clarity.
    """
    if not (1 <= module_number <= 28):
        raise ValueError("module_number must be in 1..28")
    m_idx = module_number - 1

    # x positions for CCDs
    x = np.arange(4)  # 0..3
    x_labels = [f"ext{i}" for i in range(1,5)]

    # style maps
    marker_map = {4: "o", 5: "s", 6: "^", 7: "v"}
    color_map = {1: "red", 2: "blue"}

    # horizontal offsets for images (spread evenly around 0)
    offsets = np.linspace(-0.15, 0.15, num=len(images))
    offset_map = {img: offsets[i] for i, img in enumerate(images)}

    plt.figure(figsize=(9, 6))
    for pk in peaks:
        if pk == 1:
            plt.axhline(y=5.9, color=color_map[pk], linestyle='--', label="_nolegend_")
        elif pk == 2:
            plt.axhline(y=6.49, color=color_map[pk], linestyle='--', label="_nolegend_")
        for img in images:
            i_idx = images.index(img)
            dx = offset_map[img]
            p_idx = peaks.index(pk)

            y = arr[m_idx, :, i_idx, p_idx, 0]
            yerr = arr[m_idx, :, i_idx, p_idx, 1]

            mask = y != -1
            if not np.any(mask):
                continue

            plt.errorbar(
                x[mask] + dx,     # apply offset
                y[mask],
                yerr=yerr[mask],
                fmt=marker_map[img],
                linestyle="none",
                color=color_map[pk],
                label=f"Img{img} Peak{pk}",
                capsize=2,
                linewidth=1,
                markersize=5,
            )

    plt.xticks(x, x_labels)
    plt.xlabel("CCD (extensions)")
    plt.ylabel("Energy peak center [KeV]")
    plt.ylim(5.4, 6.6)
    plt.title(f"Module {module_number}: Energy peaks with error bars")

    # deduplicate legend
    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    plt.legend(uniq.values(), uniq.keys(), ncol=2, frameon=True,title="5.9 (KeV)    6.49 (KeV)", fontsize=8)

def main():
    # ---------------- Argument parsing ----------------
    parser = argparse.ArgumentParser(description="plot Fe55 evergy peaks")
    parser.add_argument("--id", type=int, required=True, help="Module number to plot")
    parser.add_argument("--csvpath", required=True, help="path to the csv file of sql export")
    args = parser.parse_args()

    csv_path=args.csvpath
    arr = build_array_from_patterned_csv(csv_path=csv_path, images=(4,5,6))
    plot_module_energy_peaks(arr, module_number=args.id)
    plt.savefig(f"Module{args.id}_Fe55Energy.png", dpi=300)



if __name__ == "__main__":
    main()
