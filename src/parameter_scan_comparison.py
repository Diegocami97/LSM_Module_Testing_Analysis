import numpy as np
import matplotlib.pyplot as plt

# =========================
# USER INPUT (manual values)
# =========================
modules = ["DM01", "DM02"]
exts    = ["Ext1", "Ext2", "Ext3", "Ext4"]

scan_ids = ["S1", "S2", "S3", "S4", "S5"]  # 5 scan steps

# True line energies (keV)
E_KA = 5.898  # Mn Kα
E_KB = 6.490  # Mn Kβ

# -------------------------------------------------------------------
# FRONT EVENTS: manual values (replace placeholders with your numbers)
# -------------------------------------------------------------------
mu_KA_front = {
    "S1": {"DM01": {"Ext1": 5.22, "Ext2": 5.727, "Ext3": 5.82, "Ext4": 5.748},
           "DM02": {"Ext1": 5.870, "Ext2": 5.807, "Ext3": 5.862, "Ext4": 5.791},
           "DM03": {"Ext1": 5.82, "Ext2": 5.83, "Ext3": 5.84, "Ext4": 5.82}},
    "S2": {"DM01": {"Ext1": 5.22, "Ext2": 5.848, "Ext3": 5.865, "Ext4": 5.805},
           "DM02": {"Ext1": 5.841, "Ext2": 5.852, "Ext3": 5.844, "Ext4": 5.837},
           "DM03": {"Ext1": 5.86, "Ext2": 5.86, "Ext3": 5.87, "Ext4": 5.86}},
    "S3": {"DM01": {"Ext1": 5.22, "Ext2": 5.946, "Ext3": 5.839, "Ext4": 5.853},
           "DM02": {"Ext1": 5.870, "Ext2": 5.852, "Ext3": 5.863, "Ext4": 5.843},
           "DM03": {"Ext1": 5.88, "Ext2": 5.90, "Ext3": 5.89, "Ext4": 5.89}},
    "S4": {"DM01": {"Ext1": 5.22, "Ext2": 5.973, "Ext3": 5.912, "Ext4": 5.874},
           "DM02": {"Ext1": 5.882, "Ext2": 5.867, "Ext3": 5.894, "Ext4": 5.863},
           "DM03": {"Ext1": 5.880, "Ext2": 5.92, "Ext3": 5.94, "Ext4": 5.94}},
    "S5": {"DM01": {"Ext1": 5.371, "Ext2": 6.040, "Ext3": 5.752, "Ext4": 5.945},
           "DM02": {"Ext1": 5.880, "Ext2": 6.148, "Ext3": 5.883, "Ext4": 6.292},
           "DM03": {"Ext1": 5.95, "Ext2": 5.97, "Ext3": 5.98, "Ext4": 5.97}},
}

mu_KB_front = {
    "S1": {"DM01": {"Ext1": 4.187, "Ext2": 6.297, "Ext3": 6.405, "Ext4": 6.322 },
           "DM02": {"Ext1": 6.440, "Ext2": 6.386, "Ext3": 6.432, "Ext4": 6.372},
           "DM03": {"Ext1": 6.41, "Ext2": 6.41, "Ext3": 6.44, "Ext4": 6.42}},
    "S2": {"DM01": {"Ext1": 4.187, "Ext2": 6.411, "Ext3": 6.453, "Ext4": 6.404},
           "DM02": {"Ext1": 6.400, "Ext2": 6.421, "Ext3": 6.430, "Ext4": 6.409},
           "DM03": {"Ext1": 6.45, "Ext2": 6.46, "Ext3": 6.47, "Ext4": 6.45}},
    "S3": {"DM01": {"Ext1": 4.187, "Ext2": 6.522, "Ext3": 6.423, "Ext4": 6.448},
           "DM02": {"Ext1": 6.451, "Ext2": 6.428, "Ext3": 6.455, "Ext4": 6.423},
           "DM03": {"Ext1": 6.48, "Ext2": 6.49, "Ext3": 6.49, "Ext4": 6.48}},
    "S4": {"DM01": {"Ext1": 6.50, "Ext2": 6.554, "Ext3": 6.507, "Ext4": 6.450},
           "DM02": {"Ext1": 6.463, "Ext2": 6.445, "Ext3": 6.479, "Ext4": 6.418},
           "DM03": {"Ext1": 6.52, "Ext2": 6.52, "Ext3": 6.54, "Ext4": 6.52}},
    "S5": {"DM01": {"Ext1": 5.893, "Ext2": 6.617, "Ext3": 6.315, "Ext4": 6.530},
           "DM02": {"Ext1": 6.467, "Ext2": 6.766, "Ext3": 6.482, "Ext4": 6.947},
           "DM03": {"Ext1": 6.56, "Ext2": 6.56, "Ext3": 6.57, "Ext4": 6.57}},
}

# ------------------------------------------------------------------
# BACK EVENTS: manual values (placeholders — replace with your numbers)
# ------------------------------------------------------------------
# Tip: start by copying the front dicts and then editing values.
mu_KA_back = {
    "S1": {"DM01": {"Ext1": 5.367, "Ext2": 5.944, "Ext3": 5.960, "Ext4": 5.964},
           "DM02": {"Ext1": 6.016, "Ext2": 5.975, "Ext3": 6.015, "Ext4": 5.993},
           "DM03": {"Ext1": 5.90, "Ext2": 5.91, "Ext3": 5.90, "Ext4": 5.90}},
    "S2": {"DM01": {"Ext1": 5.367, "Ext2": 6.012, "Ext3": 5.988, "Ext4": 5.977},
           "DM02": {"Ext1": 5.996, "Ext2": 5.998, "Ext3": 5.974, "Ext4": 5.997},
           "DM03": {"Ext1": 5.92, "Ext2": 5.92, "Ext3": 5.93, "Ext4": 5.92}},
    "S3": {"DM01": {"Ext1": 5.367, "Ext2": 6.070, "Ext3": 5.957, "Ext4": 5.993},
           "DM02": {"Ext1": 5.999, "Ext2": 5.988, "Ext3": 5.999, "Ext4": 6.000},
           "DM03": {"Ext1": 5.94, "Ext2": 5.95, "Ext3": 5.95, "Ext4": 5.94}},
    "S4": {"DM01": {"Ext1": 5.92, "Ext2": 6.084, "Ext3": 6.022, "Ext4": 6.007},
           "DM02": {"Ext1": 6.010, "Ext2": 6.012, "Ext3": 6.031, "Ext4": 5.982},
           "DM03": {"Ext1": 5.96, "Ext2": 5.97, "Ext3": 5.98, "Ext4": 5.96}},
    "S5": {"DM01": {"Ext1": 4.877, "Ext2": 6.314, "Ext3": 5.911, "Ext4": 6.087},
           "DM02": {"Ext1": 5.992, "Ext2": 6.063, "Ext3": 6.028, "Ext4": 6.567},
           "DM03": {"Ext1": 5.98, "Ext2": 5.99, "Ext3": 6.00, "Ext4": 5.99}},
}

mu_KB_back = {
    "S1": {"DM01": {"Ext1": 5.367, "Ext2": 6.577, "Ext3": 6.565, "Ext4": 6.532},
           "DM02": {"Ext1": 6.567, "Ext2": 6.519, "Ext3": 6.636, "Ext4": 6.591},
           "DM03": {"Ext1": 6.52, "Ext2": 6.52, "Ext3": 6.53, "Ext4": 6.52}},
    "S2": {"DM01": {"Ext1": 5.367, "Ext2": 6.636, "Ext3": 6.611, "Ext4": 6.548},
           "DM02": {"Ext1": 6.570, "Ext2": 6.563, "Ext3": 6.555, "Ext4": 6.582},
           "DM03": {"Ext1": 6.53, "Ext2": 6.54, "Ext3": 6.54, "Ext4": 6.53}},
    "S3": {"DM01": {"Ext1": 5.367, "Ext2": 6.663, "Ext3": 6.526, "Ext4": 6.573},
           "DM02": {"Ext1": 6.596, "Ext2": 6.477, "Ext3": 6.561, "Ext4": 6.563},
           "DM03": {"Ext1": 6.55, "Ext2": 6.56, "Ext3": 6.56, "Ext4": 6.55}},
    "S4": {"DM01": {"Ext1": 5.367, "Ext2": 6.674, "Ext3": 6.669, "Ext4": 6.582},
           "DM02": {"Ext1": 6.540, "Ext2": 6.570, "Ext3": 6.611, "Ext4": 6.606},
           "DM03": {"Ext1": 6.57, "Ext2": 6.58, "Ext3": 6.58, "Ext4": 6.57}},
    "S5": {"DM01": {"Ext1": 5.385, "Ext2": 6.679, "Ext3": 6.6422, "Ext4": 6.670},
           "DM02": {"Ext1": 6.613, "Ext2": 6.521, "Ext3": 6.568, "Ext4": 6.783},
           "DM03": {"Ext1": 6.59, "Ext2": 6.60, "Ext3": 6.60, "Ext4": 6.59}},
}

# =========================
# PLOTTING HELPERS
# =========================
def make_x_positions(modules, exts, gap=1):
    x_positions = {}
    x_all, xtick_labels = [], []
    for mi, m in enumerate(modules):
        base = mi * (len(exts) + gap)
        for ei, e in enumerate(exts):
            x = base + ei
            x_positions[(m, e)] = x
            x_all.append(x)
            xtick_labels.append(f"{m}\n{e}")
    return x_positions, x_all, xtick_labels

def plot_scan_summary(mu_KA, mu_KB, title, modules, exts, scan_ids, E_KA, E_KB, gap=1):
    x_positions, x_all, xtick_labels = make_x_positions(modules, exts, gap=gap)

    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]  # supports >5 scans if needed
    if len(scan_ids) > len(markers):
        raise ValueError(f"Too many scan steps ({len(scan_ids)}). Add more marker styles.")
    jitter = np.linspace(-0.15, 0.15, len(scan_ids))

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11, 6.8), sharex=True, constrained_layout=True
    )

    # TOP: Kβ
    for si, sid in enumerate(scan_ids):
        xs, ys = [], []
        for m in modules:
            for e in exts:
                xs.append(x_positions[(m, e)] + jitter[si])
                ys.append(mu_KB[sid][m][e])
        ax_top.scatter(xs, ys, marker=markers[si], label=sid)

    ax_top.axhline(E_KB, linestyle="--")
    ax_top.set_ylabel(r"$\mu$ [keV]")
    ax_top.text(0.01, 0.92, "Kβ", transform=ax_top.transAxes, fontsize=12, va="top")

    # BOTTOM: Kα (requested on bottom)
    for si, sid in enumerate(scan_ids):
        xs, ys = [], []
        for m in modules:
            for e in exts:
                xs.append(x_positions[(m, e)] + jitter[si])
                ys.append(mu_KA[sid][m][e])
        ax_bot.scatter(xs, ys, marker=markers[si], label=sid)

    ax_bot.axhline(E_KA, linestyle="--")
    ax_bot.set_ylabel(r"$\mu$ [keV]")
    ax_bot.text(0.01, 0.92, "Kα", transform=ax_bot.transAxes, fontsize=12, va="top")

    ax_bot.set_xticks(x_all)
    ax_bot.set_xticklabels(xtick_labels)

    # module separators
    for mi in range(len(modules) - 1):
        sep_x = (mi + 1) * (len(exts) + gap) - gap / 2 - 0.5
        ax_top.axvline(sep_x, linestyle=":", linewidth=1)
        ax_bot.axvline(sep_x, linestyle=":", linewidth=1)

    ax_top.legend(
        ncol=min(len(scan_ids), 8),
        loc="upper center",
        bbox_to_anchor=(0.5, 1.25),
        frameon=True,
    )

    fig.suptitle(title, y=1.02)
    plt.show()

# =========================
# MAKE TWO FIGURES: FRONT + BACK
# =========================
plot_scan_summary(
    mu_KA_front, mu_KB_front,
    title="Front events — peak mean vs (Module, Extension)",
    modules=modules, exts=exts, scan_ids=scan_ids,
    E_KA=E_KA, E_KB=E_KB, gap=1
)

plot_scan_summary(
    mu_KA_back, mu_KB_back,
    title="Back events — peak mean vs (Module, Extension)",
    modules=modules, exts=exts, scan_ids=scan_ids,
    E_KA=E_KA, E_KB=E_KB, gap=1
)
