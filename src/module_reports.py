#!/usr/bin/env python3
# Diego Venegas-Vargas
# DAMIC-M Collaboration
# Reads module data from CSV and generates module-wise plots for CCD metrics.
"""
Module Reports Script
Reads module data from CSV and generates module-wise plots for CCD metrics.
"""

import argparse
import os
import sys
import re
import json
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def ensure_interactive_backend():
    """Switch off Agg so plt.show() can open a GUI window (pan/zoom/axis limits)."""
    current = matplotlib.get_backend()
    if current and current.lower() != 'agg':
        return current
    for name in ('MacOSX', 'QtAgg', 'Qt5Agg', 'TkAgg'):
        try:
            plt.switch_backend(name)
            return name
        except Exception:
            continue
    print("Warning: no GUI matplotlib backend found; --show will not open a window.",
          file=sys.stderr)
    return matplotlib.get_backend()

DEFAULT_MODULES_EXCLUDED_FROM_PLOT = []
DEFAULT_EXCLUDE_FROM_MEAN = ['DM-08', 'DM-12']
DEFAULT_OUTLIER_SIGMA = 3.0


def outlier_modules_mad(df_pivot, n_sigma=DEFAULT_OUTLIER_SIGMA):
    """
    Return module IDs whose per-CCD mean is an outlier relative to other modules.

    Uses median and MAD of the per-module means: a module is an outlier if
    |mean - median| > n_sigma * 1.4826 * MAD. Those modules stay on the plot
    but can be dropped from the ensemble Mean ± 1σ band.
    """
    mod_means = df_pivot.mean(axis=1).dropna()
    if len(mod_means) < 3:
        return []
    values = mod_means.values.astype(float)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    sigma = 1.4826 * mad
    if sigma <= 0:
        return []
    z = np.abs(values - median) / sigma
    return [str(mid) for mid, zi in zip(mod_means.index, z) if zi > n_sigma]


def parse_value_with_uncertainty(value):
    """
    Parse a value that may contain uncertainty information.
    
    Handles formats like:
    - "5.73442 +/- 0.00165312" -> (5.73442, 0.00165312)
    - "1.23 ± 0.05" -> (1.23, 0.05)
    - "1.23(5)" -> (1.23, 0.05)
    - "1.23" or 1.23 -> (1.23, None)
    
    Returns:
        tuple: (value, uncertainty) where uncertainty is None if not present
    """
    if pd.isna(value):
        return (np.nan, None)
    
    # If already numeric, return as-is
    if isinstance(value, (int, float)):
        return (float(value), None)
    
    # Convert to string
    val_str = str(value).strip()
    
    # Try to parse as float first
    try:
        return (float(val_str), None)
    except ValueError:
        pass
    
    # Try patterns with uncertainty
    patterns = [
        r'([+-]?\d+\.?\d*)\s*\+\/-\s*([+-]?\d+\.?\d*)',  # "5.73442 +/- 0.00165312"
        r'([+-]?\d+\.?\d*)\s*±\s*([+-]?\d+\.?\d*)',        # "1.23 ± 0.05"
        r'([+-]?\d+\.?\d*)\((\d+)\)',                       # "1.23(5)" -> uncertainty in last digit
    ]
    
    for pattern in patterns:
        match = re.match(pattern, val_str)
        if match:
            val = float(match.group(1))
            err_str = match.group(2)
            
            # For format "1.23(5)", convert to absolute uncertainty
            if pattern == patterns[2]:
                # "1.23(5)" means uncertainty of 5 in the last decimal place
                # So if value is 1.23, uncertainty is 0.05
                decimal_places = len(str(val).split('.')[-1]) if '.' in str(val) else 0
                uncertainty = float(err_str) * (10 ** (-decimal_places))
            else:
                uncertainty = float(err_str)
            
            return (val, uncertainty)
    
    # If no pattern matches, try to extract just the number
    num_match = re.search(r'([+-]?\d+\.?\d*)', val_str)
    if num_match:
        return (float(num_match.group(1)), None)
    
    # If all else fails, return NaN
    return (np.nan, None)


def get_gain_metric_for_res(res_metric_name):
    """
    Return the corresponding Gain metric name for a Res metric, or None.
    E.g. Image31_Low_Res -> Image31_Low_Gain.
    """
    if not res_metric_name.endswith('_Res'):
        return None
    return res_metric_name[:-4] + '_Gain'  # replace _Res with _Gain


def get_gain_metric_for_dark_current(dark_metric_name):
    """
    Return the corresponding Gain metric name for a Dark_Current metric, or None.
    E.g. Image31_Low_Dark_Current -> Image31_Low_Gain.
    """
    if 'Ratio' in dark_metric_name:
        return None
    if not dark_metric_name.endswith('_Dark_Current'):
        return None
    return dark_metric_name[:-13] + '_Gain'  # replace _Dark_Current with _Gain (len("_Dark_Current")==13)


# Time/unit factor to e-/pix/day. Serial-binning (NSBIN) is applied separately:
# Low-T Image 3 uses NSBIN=10; High-T Image 3 is 1x1 (NSBIN=1).
DARK_CURRENT_TIME_FACTOR = 232.9
LOW_DC_SERIAL_BINNING = 10
HIGH_DC_SERIAL_BINNING = 1


def serial_binning_for_dc_metric(metric_name):
    """NSBIN used when converting a Dark_Current metric. None if not a DC metric."""
    if 'Ratio' in metric_name or not metric_name.endswith('_Dark_Current'):
        return None
    if '_High_Dark_Current' in metric_name:
        return HIGH_DC_SERIAL_BINNING
    if '_Low_Dark_Current' in metric_name:
        return LOW_DC_SERIAL_BINNING
    return LOW_DC_SERIAL_BINNING

IMAGE31_LOW_DC_METRIC = 'Image31_Low_Dark_Current'
IMAGE31_HIGH_DC_METRIC = 'Image31_High_Dark_Current'
IMAGE31_DC_RATIO_METRIC = 'Image31_Ratio_Low_High_Dark_Current'
IMAGE3_LOW_DC_METRIC = 'Image3_Low_Dark_Current'
IMAGE3_HIGH_DC_METRIC = 'Image3_High_Dark_Current'
IMAGE3_DC_RATIO_METRIC = 'Image3_Ratio_Low_High_Dark_Current'

DC_RATIO_SPECS = [
    (IMAGE31_LOW_DC_METRIC, IMAGE31_HIGH_DC_METRIC, IMAGE31_DC_RATIO_METRIC),
    (IMAGE3_LOW_DC_METRIC, IMAGE3_HIGH_DC_METRIC, IMAGE3_DC_RATIO_METRIC),
]


def apply_dark_current_conversion(df_long, df_metric, metric_name):
    """
    Convert Dark_Current to e-/pix/day: value * (1/gain) * (1/NSBIN) * 232.9.

    NSBIN is 10 for Low (10-column serial binning) and 1 for High (1x1).
    """
    gain_metric = get_gain_metric_for_dark_current(metric_name)
    if gain_metric is None:
        return df_metric
    nsbin = serial_binning_for_dc_metric(metric_name)
    if nsbin is None:
        nsbin = 1

    df_metric = df_metric.copy()
    df_metric['module_id'] = df_metric['module_id'].astype(str).str.strip()
    df_metric['ccd_id'] = df_metric['ccd_id'].astype(str).str.strip()
    df_gain = df_long.loc[df_long['metric_name'] == gain_metric, ['module_id', 'ccd_id', 'value']].copy()
    if df_gain.empty:
        print(f"    Skipping dark current conversion: no rows for gain metric '{gain_metric}'")
        return df_metric

    df_gain = df_gain.drop_duplicates(subset=['module_id', 'ccd_id'], keep='first')
    df_gain['module_id'] = df_gain['module_id'].astype(str).str.strip()
    df_gain['ccd_id'] = df_gain['ccd_id'].astype(str).str.strip()
    df_gain = df_gain.rename(columns={'value': 'gain'})
    df_gain['gain'] = pd.to_numeric(df_gain['gain'], errors='coerce')
    df_metric = df_metric.merge(df_gain[['module_id', 'ccd_id', 'gain']], on=['module_id', 'ccd_id'], how='left')
    gain_vals = pd.to_numeric(df_metric['gain'], errors='coerce')
    valid = gain_vals.notna() & (gain_vals != 0)
    scale = np.where(valid, DARK_CURRENT_TIME_FACTOR / (gain_vals.values * nsbin), np.nan)
    vals = pd.to_numeric(df_metric['value'], errors='coerce')
    df_metric['value'] = vals.values * scale
    if 'uncertainty' in df_metric.columns:
        unc = pd.to_numeric(df_metric['uncertainty'], errors='coerce')
        df_metric['uncertainty'] = unc.values * scale
    df_metric = df_metric.drop(columns=['gain'])
    n_applied = int(valid.sum())
    print(f"    Applied dark current conversion: {metric_name} * (1/{gain_metric})*(1/{nsbin})*{DARK_CURRENT_TIME_FACTOR}  ({n_applied} values)")
    return df_metric


def add_derived_dc_ratio_metrics(df_long):
    """
    Build Low/High dark-current ratio metrics from CSV pairs.

    Supports Image31_* (underground) and Image3_* (surface). Both Low and High
    are converted first, then ratio = Low / High per (module, CCD).
    """
    available = set(df_long['metric_name'].unique())
    derived_names = []
    keys = ['module_id', 'ccd_id']

    for low_name, high_name, ratio_name in DC_RATIO_SPECS:
        if low_name not in available or high_name not in available:
            continue

        df_low = apply_dark_current_conversion(df_long, df_long[df_long['metric_name'] == low_name], low_name)
        df_high = apply_dark_current_conversion(df_long, df_long[df_long['metric_name'] == high_name], high_name)

        low_cols = keys + ['value'] + (['uncertainty'] if 'uncertainty' in df_low.columns else [])
        high_cols = keys + ['value'] + (['uncertainty'] if 'uncertainty' in df_high.columns else [])
        left = df_low[low_cols].rename(columns={'value': 'low', 'uncertainty': 'low_unc'})
        right = df_high[high_cols].rename(columns={'value': 'high', 'uncertainty': 'high_unc'})
        left = left.drop_duplicates(subset=keys, keep='first')
        right = right.drop_duplicates(subset=keys, keep='first')
        merged = left.merge(right, on=keys, how='inner')

        low_vals = pd.to_numeric(merged['low'], errors='coerce')
        high_vals = pd.to_numeric(merged['high'], errors='coerce')
        valid = low_vals.notna() & high_vals.notna() & (high_vals != 0)
        merged['value'] = np.where(valid, low_vals / high_vals, np.nan)

        merged['uncertainty'] = np.nan
        if 'low_unc' in merged.columns and 'high_unc' in merged.columns:
            low_unc = pd.to_numeric(merged['low_unc'], errors='coerce')
            high_unc = pd.to_numeric(merged['high_unc'], errors='coerce')
            rel = valid & low_vals.notna() & (low_vals != 0) & low_unc.notna() & high_unc.notna()
            merged.loc[rel, 'uncertainty'] = np.abs(merged.loc[rel, 'value']) * np.sqrt(
                (low_unc[rel] / low_vals[rel]) ** 2 + (high_unc[rel] / high_vals[rel]) ** 2
            )

        merged['metric_name'] = ratio_name
        derived = merged[keys + ['metric_name', 'value', 'uncertainty']].copy()
        n_ok = int(valid.sum())
        print(f"    Derived {ratio_name}: Low/High after DC conversion ({n_ok} values)")
        df_long = pd.concat([df_long, derived], ignore_index=True)
        derived_names.append(ratio_name)

    if not derived_names:
        print("    Skipping DC ratio: no Image31 or Image3 Low/High Dark_Current pair found")
    return df_long, derived_names


def dc_plot_metrics_for_csv(available_metrics):
    """Pick Image31 (underground) or Image3 (surface) DC plot metrics."""
    available = set(available_metrics)
    if IMAGE31_LOW_DC_METRIC in available and IMAGE31_HIGH_DC_METRIC in available:
        return [IMAGE31_LOW_DC_METRIC, IMAGE31_HIGH_DC_METRIC, IMAGE31_DC_RATIO_METRIC]
    if IMAGE3_LOW_DC_METRIC in available and IMAGE3_HIGH_DC_METRIC in available:
        return [IMAGE3_LOW_DC_METRIC, IMAGE3_HIGH_DC_METRIC, IMAGE3_DC_RATIO_METRIC]
    return None


def get_metric_config(metric_name, config_file=None):
    """
    Get plot configuration for a specific metric.
    
    Configuration options:
    - ylabel: Y-axis label (with units if needed)
    - ylabel_unit: Unit to append to metric name
    - value_format: Format string for values (e.g., "{:.2f}", "{:.3e}")
    - log_scale: Use log scale for y-axis (True/False)
    - show_error_bars: Show error bars if uncertainties are present (True/False)
    - decimal_places: Number of decimal places for value labels
    
    Args:
        metric_name: Name of the metric
        config_file: Optional path to JSON config file
    
    Returns:
        dict: Configuration dictionary with default values if not specified
    """
    # Default configuration
    default_config = {
        'ylabel': None,  # If None, use metric name
        'ylabel_unit': None,  # Unit to append (e.g., "eV", "e-/ADU", "K")
        'value_format': None,  # If None, auto-format based on value
        'log_scale': False,
        'show_error_bars': True,
        'decimal_places': None,  # If None, auto-determine
        'plot_title': None,  # If None, use default title format
        'decay_energy': None,  # If set (e.g. 5.9), draw horizontal line at this energy (keV)
        'decay_line_label': None,  # Label for decay line (e.g. "Fe-55 K#alpha")
    }
    
    # Load custom config from file if provided
    custom_configs = {}
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                custom_configs = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load config file {config_file}: {e}")
    
    # Check for exact metric match first
    if metric_name in custom_configs:
        config = {**default_config, **custom_configs[metric_name]}
    else:
        # Check for pattern matches (metric names containing pattern)
        config = default_config.copy()
        for pattern, pattern_config in custom_configs.items():
            if pattern.lower() in metric_name.lower():
                config.update(pattern_config)
                break
    
    # Apply metric-specific defaults based on metric name patterns
    if metric_name in (IMAGE31_LOW_DC_METRIC, IMAGE3_LOW_DC_METRIC):
        config.setdefault('ylabel', 'Dark Current [e-/pix/day]')
        config.setdefault('decimal_places', 4)
        config.setdefault('value_format', '{:.4e}')
        config['mean_format'] = '{:.4f}'
        config.setdefault('log_scale', True)
        config.setdefault('plot_title', 'Dark Current from 500-Skip SR Image at 135K')
        config.setdefault('y_min', 1e-2)
        config.setdefault('y_max', 3e1)
        config.setdefault('dm_modules_only', True)
    elif metric_name in (IMAGE31_HIGH_DC_METRIC, IMAGE3_HIGH_DC_METRIC):
        config.setdefault('ylabel', 'Dark Current [e-/pix/day]')
        config.setdefault('decimal_places', 4)
        config.setdefault('value_format', '{:.4e}')
        config['mean_format'] = '{:.4f}'
        config.setdefault('log_scale', True)
        config.setdefault('plot_title', 'Dark Current from 500-Skip SR Image at 180K')
        config.setdefault('y_min', 1e1)
        config.setdefault('y_max', 1e3)
        config.setdefault('dm_modules_only', True)
    elif 'Ratio_Low_High_Dark_Current' in metric_name:
        config.setdefault('ylabel', 'Low / High Dark Current')
        config.setdefault('decimal_places', 4)
        config.setdefault('value_format', '{:.4e}')
        config['mean_format'] = '{:.4f}'
        config.setdefault('log_scale', True)
        config.setdefault('line_at_1', False)
        config.setdefault('plot_title', 'Dark Current Ratio (Low / High) from 500-Skip SR Image')
        config.setdefault('y_min', 1e-4)
        config.setdefault('y_max', 1e-2)
        config.setdefault('dm_modules_only', True)
    elif 'Res' in metric_name:
        # Resolution metrics (typically in eV)
        config.setdefault('ylabel_unit', 'eV')
        config.setdefault('decimal_places', 3)
        config.setdefault('value_format', '{:.3f}')
    elif 'Gain' in metric_name:
        # Gain metrics (typically in e-/ADU)
        config.setdefault('ylabel_unit', 'e-/ADU')
        config.setdefault('decimal_places', 3)
        config.setdefault('value_format', '{:.3f}')
    elif 'Dark_Current' in metric_name:
        # Dark current after (1/gain)*(1/NSBIN)*232.9  (e-/pix/day)
        config.setdefault('ylabel_unit', 'e-/pix/day')
        config.setdefault('decimal_places', 4)
        config.setdefault('value_format', '{:.4e}')
        config.setdefault('log_scale', True)
    elif 'Temp' in metric_name:
        # Temperature (typically in K or °C)
        config.setdefault('ylabel_unit', 'K')
        config.setdefault('decimal_places', 1)
        config.setdefault('value_format', '{:.1f}')
    elif 'Peak' in metric_name:
        # Peak values (may have uncertainties)
        config.setdefault('decimal_places', 3)
        config.setdefault('value_format', '{:.3f}')
    elif 'CTI' in metric_name:
        # CTI metrics (typically fractions or small numbers)
        config.setdefault('decimal_places', 4)
        config.setdefault('value_format', '{:.4f}')
    elif 'Noise' in metric_name:
        # Noise metrics
        config.setdefault('ylabel_unit', 'ADU')
        config.setdefault('decimal_places', 2)
        config.setdefault('value_format', '{:.2f}')
    
    # Build ylabel if unit is specified
    if config['ylabel'] is None and config['ylabel_unit']:
        config['ylabel'] = f"{metric_name} [{config['ylabel_unit']}]"
    elif config['ylabel'] is None:
        config['ylabel'] = metric_name
    
    return config


def read_and_process_data(filepath):
    """
    Read CSV and reshape from wide format to long format.
    
    Returns:
        df_long: DataFrame with columns [module_id, ccd_id, metric_name, value]
        metrics: list of unique metric base names
    """
    print(f"Reading data from: {filepath}")
    df = pd.read_csv(filepath)
    print(f"Loaded {len(df)} rows, {len(df.columns)} columns")
    
    # Identify module ID column (usually 'Name')
    module_col = None
    for col in ['Name', 'name', 'Module', 'module', 'Module_ID', 'module_id']:
        if col in df.columns:
            module_col = col
            break
    
    if module_col is None:
        raise ValueError("Could not find module ID column. Expected: Name, name, Module, module, Module_ID, or module_id")
    
    # Identify CCD-specific columns (ending in _A, _B, _C, _D)
    # A "metric" is the base name of a measurement that exists for each CCD.
    # Example: Column "Image31_High_Res_A" -> metric = "Image31_High_Res"
    #          This metric has 4 columns: _A, _B, _C, _D (one per CCD)
    ccd_suffixes = ['_A', '_B', '_C', '_D']
    ccd_columns = []
    metrics = set()
    
    for col in df.columns:
        for suffix in ccd_suffixes:
            if col.endswith(suffix):
                base_name = col[:-len(suffix)]  # Extract metric name (e.g., "Image31_High_Res")
                ccd_columns.append(col)          # Keep full column name (e.g., "Image31_High_Res_A")
                metrics.add(base_name)           # Store unique metric name
                break
    
    metrics = sorted(list(metrics))  # List of all available metrics (e.g., ["Image31_High_Res", "Image31_Low_Gain", ...])
    print(f"Found {len(set(df[module_col]))} modules: {sorted(set(df[module_col]))[:10]}...")
    print(f"Found {len(metrics)} metrics with CCD suffixes (_A, _B, _C, _D)")
    
    # Reshape to long format
    id_vars = [module_col]
    value_vars = ccd_columns
    
    df_long = pd.melt(
        df,
        id_vars=id_vars,
        value_vars=value_vars,
        var_name='metric_ccd',
        value_name='value'
    )
    
    # Extract metric name and CCD ID
    df_long['ccd_id'] = df_long['metric_ccd'].str[-1]  # Last character (_A -> A)
    df_long['metric_name'] = df_long['metric_ccd'].str[:-2]  # Everything except _A
    
    # Rename module column
    df_long = df_long.rename(columns={module_col: 'module_id'})
    
    # Select relevant columns
    df_long = df_long[['module_id', 'ccd_id', 'metric_name', 'value']].copy()
    
    # Parse values that may contain uncertainties
    # Store both the value and uncertainty separately
    parsed = df_long['value'].apply(parse_value_with_uncertainty)
    df_long['value'] = [v[0] for v in parsed]
    df_long['uncertainty'] = [v[1] for v in parsed]

    df_long, derived_metrics = add_derived_dc_ratio_metrics(df_long)
    if derived_metrics:
        metrics = sorted(list(metrics) + derived_metrics)
    
    return df_long, metrics


def plot_metrics(df_long, metrics, output_dir, pattern=None, config_file=None,
                 exclude_from_mean=None, exclude_from_plot=None, show_module_averages=False,
                 auto_exclude_from_mean=False, outlier_sigma=DEFAULT_OUTLIER_SIGMA,
                 show=False, log_scale=None):
    """
    Generate plots for specified metrics.
    
    Args:
        df_long: Long-format DataFrame with columns [module_id, ccd_id, metric_name, value, uncertainty]
        metrics: List of metric names to plot
        output_dir: Output directory for plots
        pattern: Optional pattern to filter metrics (substring match, case-insensitive)
                 Example: pattern="Image31" matches "Image31_High_Res", "Image31_Low_Gain", etc.
        config_file: Optional path to JSON file with per-metric plot configurations
        exclude_from_mean: List of module IDs to exclude from mean ± 1σ calculation (still plotted)
        exclude_from_plot: List of module IDs to exclude entirely from plots
        show_module_averages: If True, print per-module average (mean across CCDs) to the screen
        auto_exclude_from_mean: If True, also drop MAD-outlier modules from the mean band
        outlier_sigma: MAD threshold used when auto_exclude_from_mean is True
        show: If True, open each figure in an interactive matplotlib window
        log_scale: If True, force a logarithmic y-axis
    """
    if exclude_from_mean is None:
        exclude_from_mean = []
    if exclude_from_plot is None:
        exclude_from_plot = []

    # Filter metrics if pattern provided
    # A "pattern" is a substring to match against metric names (case-insensitive)
    # Example: pattern="Image31" will match all metrics containing "Image31"
    #          pattern="Res" will match "Image31_High_Res", "Image31_Low_Res", "Image32_High_Res", etc.
    if pattern:
        metrics = [m for m in metrics if pattern.lower() in m.lower()]
    
    if not metrics:
        print(f"No metrics found matching pattern: {pattern}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\nGenerating plots for {len(metrics)} metric(s)...")
    
    for metric in metrics:
        print(f"  Processing: {metric}")
        
        # Get configuration for this metric
        config = get_metric_config(metric, config_file)
        if log_scale:
            config['log_scale'] = True
        
        # Filter data for this specific metric
        # A "metric" here is an exact metric name (e.g., "Image31_High_Res")
        # This extracts all rows for this metric across all modules and CCDs
        df_metric = df_long[df_long['metric_name'] == metric].copy()
        
        if df_metric.empty:
            print(f"    Warning: No data found for {metric}")
            continue
        
        # If this is a Res metric, divide by the corresponding Gain metric when available
        gain_metric = get_gain_metric_for_res(metric)
        if gain_metric is not None:
            gain_in_long = (df_long['metric_name'] == gain_metric).any()
            if gain_in_long:
                df_gain = df_long.loc[df_long['metric_name'] == gain_metric, ['module_id', 'ccd_id', 'value']].copy()
                df_gain = df_gain.rename(columns={'value': 'gain'})
                df_metric = df_metric.merge(df_gain, on=['module_id', 'ccd_id'], how='left')
                # value = Res / Gain (avoid div by zero)
                df_metric['value'] = np.where(
                    (df_metric['gain'].notna()) & (df_metric['gain'] != 0),
                    df_metric['value'] / df_metric['gain'],
                    np.nan
                )
                if 'uncertainty' in df_metric.columns:
                    df_metric['uncertainty'] = np.where(
                        (df_metric['gain'].notna()) & (df_metric['gain'] != 0),
                        df_metric['uncertainty'] / df_metric['gain'],
                        np.nan
                    )
                df_metric = df_metric.drop(columns=['gain'])
                print(f"    Applied gain correction: {metric} / {gain_metric}")
        
        # For *_Dark_Current metrics: multiply by (1/gain)*(1/NSBIN)*232.9
        # Low uses NSBIN=10; High is 1x1 (NSBIN=1). Skip derived ratio metrics.
        if 'Ratio' not in metric:
            df_metric = apply_dark_current_conversion(df_long, df_metric, metric)

        # Exact zeros in the QC table are unfilled / failed entries, not a measured rate.
        # Drop them so they do not appear on the log axis or pull down the mean.
        if 'Dark_Current' in metric:
            zero = pd.to_numeric(df_metric['value'], errors='coerce') == 0
            if zero.any():
                bad = df_metric.loc[zero, ['module_id', 'ccd_id']]
                pairs = [f"{m}-{c}" for m, c in zip(bad['module_id'], bad['ccd_id'])]
                print(f"    Treating value=0 as missing: {', '.join(pairs)}")
                df_metric.loc[zero, 'value'] = np.nan
                if 'uncertainty' in df_metric.columns:
                    df_metric.loc[zero, 'uncertainty'] = np.nan
         
        # Pivot to wide format for plotting (values)
        df_pivot = df_metric.pivot_table(
            index='module_id',
            columns='ccd_id',
            values='value',
            aggfunc='first'
        )
        
        # Pivot uncertainties if available
        df_uncertainty = None
        if 'uncertainty' in df_metric.columns:
            df_uncertainty = df_metric.pivot_table(
                index='module_id',
                columns='ccd_id',
                values='uncertainty',
                aggfunc='first'
            )
            # Check if we have any uncertainties
            has_uncertainties = df_uncertainty.notna().any().any() and config['show_error_bars']
        else:
            has_uncertainties = False
        
        # Sort modules naturally (DM-01, DM-02, ...)
        module_order = sorted(df_pivot.index, key=lambda x: int(str(x).split('-')[1]) if '-' in str(x) else 0)
        df_pivot = df_pivot.reindex(module_order)
        # Exclude specified modules from plot
        if exclude_from_plot:
            df_pivot = df_pivot.drop(index=exclude_from_plot, errors='ignore')
        if config.get('dm_modules_only'):
            keep = [m for m in df_pivot.index if str(m).startswith('DM-')]
            dropped = [m for m in df_pivot.index if m not in keep]
            if dropped:
                print(f"    Excluding non-DM modules from plot: {', '.join(map(str, dropped))}")
            df_pivot = df_pivot.reindex(keep)

        # Align uncertainty pivot with value pivot
        if df_uncertainty is not None:
            df_uncertainty = df_uncertainty.reindex(df_pivot.index)
        
        # Create plot with reasonable size limits (larger aspect: bigger figure, bars, text)
        scale = 1.4  # overall scale factor for figure and elements
        fig_width = min(max(16, len(df_pivot) * 0.7), 60) * scale
        fig_height = 8 * scale
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Get CCDs present in data
        ccds = sorted([c for c in ['A', 'B', 'C', 'D'] if c in df_pivot.columns])
        n_ccds = len(ccds)
        
        if n_ccds == 0:
            print(f"    Warning: No valid CCD data found for {metric}")
            continue
        
        if len(df_pivot) == 0:
            print(f"    Warning: No valid module data found for {metric}")
            continue

        # Optionally print per-module average (mean across CCDs) to the screen
        if show_module_averages:
            mod_means = df_pivot.mean(axis=1)
            print(f"    Per-module average ({metric}):")
            for mid in mod_means.index:
                val = mod_means[mid]
                if pd.isna(val):
                    print(f"      {mid}: (no data)")
                else:
                    print(f"      {mid}: {val}")
            print()

        width = 1 / n_ccds  # wider bars
        # Slightly space out modules so tick marks and bar clusters have a bit more separation
        module_spacing = 1.5
        x = np.arange(len(df_pivot)) * module_spacing

        # Plot bars for each CCD
        # colors = plt.cm.tab(np.linspace(0, 1, 4))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        for i, ccd in enumerate(ccds):
            values = df_pivot[ccd].values
            uncertainties = None
            if has_uncertainties and df_uncertainty is not None and ccd in df_uncertainty.columns:
                uncertainties = df_uncertainty[ccd].values
            
            # Offset the bars to the center of the plot
            # This centers the bars around the x-axis
            # i is the index of the CCD
            # n_ccds is the number of CCDs
            # width is the width of the bars
            # The offset is the distance between the bars
            # The offset is calculated as the distance between the center of the plot and the center of the bar
            # I want the distance between modules to be larger 
            offset = (i - n_ccds/2 + 0.5) * width
            
            # Plot with or without error bars
            if has_uncertainties and uncertainties is not None:
                bars = ax.bar(x + offset, values, width, yerr=uncertainties,
                             label=f'CCD {ccd}', color=colors[i], capsize=5, error_kw={'elinewidth': 2, 'capthick': 1.5})
            else:
                bars = ax.bar(x + offset, values, width, label=f'CCD {ccd}', color=colors[i])
            
            # Add value labels on bars
            for j, (bar, val) in enumerate(zip(bars, values)):
                if not pd.isna(val):
                    height = bar.get_height()
                    # Format value based on config
                    if config['value_format']:
                        label_text = config['value_format'].format(val)
                    elif config['decimal_places'] is not None:
                        label_text = f'{val:.{config["decimal_places"]}f}'
                    else:
                        # Auto-format
                        if abs(val) < 0.01 or abs(val) > 1000:
                            label_text = f'{val:.2e}'
                        else:
                            label_text = f'{val:.2f}'
                    
                    # Place text label with clipping to prevent figure size issues
                    # For log scale, use a small offset to avoid issues
                    if config['log_scale']:
                        # For log scale, place label slightly above bar
                        y_pos = height * 1.1 if height > 0 else height * 0.9
                    else:
                        y_pos = height
                    
                    # ax.text(bar.get_x() + bar.get_width()/2., y_pos,
                    #        label_text, ha='center', va='bottom', fontsize=20 
                    #        rotation=90, clip_on=True)
        
        # Calculate mean and standard deviation (optionally excluding some modules from mean only)
        mean_exclude = list(exclude_from_mean)
        if auto_exclude_from_mean:
            auto_outliers = outlier_modules_mad(df_pivot, n_sigma=outlier_sigma)
            extra = [m for m in auto_outliers if m not in mean_exclude]
            if extra:
                print(f"    Auto-excluded from mean (MAD > {outlier_sigma}σ): {', '.join(extra)}")
            mean_exclude = list(dict.fromkeys(mean_exclude + auto_outliers))
        df_for_mean = df_pivot.drop(index=mean_exclude, errors='ignore')
        all_values = df_for_mean.values.flatten()
        all_values = all_values[~pd.isna(all_values)]  # Remove NaN values
        
        if len(all_values) > 0:
            mean_val = np.mean(all_values)
            std_val = np.std(all_values)
            
            # Format mean and std for legend (mean_format wins so the legend can be
            # more precise than bar labels / scientific notation).
            mean_fmt = config.get('mean_format') or config.get('value_format')
            if mean_fmt:
                mean_str = mean_fmt.format(mean_val)
                std_str = mean_fmt.format(std_val)
            elif config['decimal_places'] is not None:
                mean_str = f'{mean_val:.{config["decimal_places"]}f}'
                std_str = f'{std_val:.{config["decimal_places"]}f}'
            else:
                # Auto-format
                if abs(mean_val) < 0.01 or abs(mean_val) > 1000:
                    mean_str = f'{mean_val:.4e}'
                    std_str = f'{std_val:.4e}'
                else:
                    mean_str = f'{mean_val:.4f}'
                    std_str = f'{std_val:.4f}'
            print(f"    Mean = {mean_val:.6g}  ± 1σ = {std_val:.6g}  (legend: {mean_str})")
            
            # Add grey band for ±1 sigma (behind bars, zorder=0)
            # For log scale, ensure values are positive
            if config['log_scale']:
                band_low = max(mean_val - std_val, all_values.min() * 0.1) if mean_val - std_val > 0 else mean_val * 0.5
                band_high = mean_val + std_val
            else:
                band_low = mean_val - std_val
                band_high = mean_val + std_val
            
            ax.axhspan(band_low, band_high, 
                      alpha=0.2, color='grey', zorder=0,
                      label=f'Mean ± 1σ')
            
            # Add dashed black line for mean
            ax.axhline(y=mean_val, color='grey', linestyle=':', linewidth=2.5, 
                      zorder=1, label=f'Mean: {mean_str}')

            if config.get('decay_energy') is not None:
                label = config.get('decay_line_label') or 'Decay'
                ax.axhline(
                    y=config['decay_energy'],
                    color='magenta',
                    linestyle=':',
                    linewidth=2.5,
                    zorder=1,
                    label=f"{label} Energy: {config['decay_energy']} keV"
                )

            if config.get('line_at_1'):
                ax.axhline(y=1, color='magenta', linestyle=':', linewidth=2.5,
                          zorder=1, label='Ratio = 1')

        ax.set_xlabel('Module', fontsize=30)
        ax.set_ylabel(config['ylabel'], fontsize=30)
        # Use custom title if provided, otherwise use default
        if config['plot_title']:
            plot_title = config['plot_title']
        else:
            plot_title = f'{metric} by Module and CCD'
        ax.set_title(plot_title, fontsize=32, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(df_pivot.index, rotation=45, ha='right', fontsize=30)
        ax.tick_params(axis='y', labelsize=30)
        handles, labels = ax.get_legend_handles_labels()
        mean_idx = [i for i, lab in enumerate(labels) if lab.startswith('Mean')]
        other_idx = [i for i in range(len(labels)) if i not in mean_idx]
        order = mean_idx + other_idx
        ax.legend([handles[i] for i in order], [labels[i] for i in order],
                  loc='upper left', fontsize=26, frameon=True, bbox_to_anchor=(1, 1),
                  handlelength=2, handleheight=1.5, labelspacing=0.8, borderpad=1)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Apply log scale if configured
        if config['log_scale']:
            ax.set_yscale('log')
            # Set reasonable y-axis limits for log scale to prevent huge figures
            if config.get('y_min') and config.get('y_max') is not None:
                ax.set_ylim(config['y_min'], config['y_max'])
            else:
                y_min = df_pivot.min().min()
                y_max = df_pivot.max().max()
                if y_min > 0 and y_max > 0:
                    # Add small margins for log scale
                    ax.set_ylim(y_min * 0.5, y_max * 2.0)
        else:
            # For linear scale, set reasonable limits with small margin
            # y_min = df_pivot.min().min()
            # y_max = df_pivot.max().max()
            # if not (pd.isna(y_min) or pd.isna(y_max)):
            #     margin = (y_max - y_min) * 0.1 if y_max != y_min else abs(y_max) * 0.1
                # ax.set_ylim(max(0, y_min - margin), y_max + margin)
            if config.get('y_min') and config.get('y_max') is not None:
                ax.set_ylim(config['y_min'], config['y_max'])
            elif config.get('y_min') is not None:
                ax.set_ylim(config['y_min'], band_high + band_high * 0.2)
            elif config.get('y_max') is not None:
                ax.set_ylim(0, config['y_max'])
            else:
                ax.set_ylim(0, band_high + band_high * 0.2)

        
        # Try to use constrained_layout, fall back to tight_layout if not available
        try:
            fig.set_constrained_layout(True)
        except AttributeError:
            # Older matplotlib versions don't have constrained_layout
            plt.tight_layout()
        
        # Save plot with reasonable DPI and size constraints
        filename = metric.lower().replace(' ', '_').replace('/', '_')
        filepath = os.path.join(output_dir, f'{filename}_all_modules.png')
        filepath = os.path.join(output_dir, f'{filename}_all_modules.pdf')
        
        # Limit DPI to prevent huge file sizes, but ensure minimum quality
        save_dpi = 150

        saved = {'done': False}

        def _save_current_view(_event=None):
            if saved['done']:
                return
            saved['done'] = True
            fig.savefig(filepath, dpi=save_dpi, bbox_inches='tight', facecolor='white',
                        pad_inches=0.2)

        if show:
            print(f"    Interactive view: {metric}")
            print("    Use the toolbar to pan/zoom. Close the window to save that view and continue.")
            fig.canvas.mpl_connect('close_event', _save_current_view)
            plt.show()
        else:
            _save_current_view()
        print(f"    Saved: {filepath}")
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate module-wise plots for CCD metrics from CSV data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('csv_file', help='Path to CSV file with module data')
    parser.add_argument('--output-dir', '-o', default='plots',
                       help='Output directory for plots (default: plots)')
    parser.add_argument('--metric', '-m', action='append', dest='metrics',
                       help='Metric name to plot (can be used multiple times)')
    parser.add_argument('--pattern', '-p', 
                       help='Pattern to match metric names (e.g., "Image31" or "Res")')
    parser.add_argument('--list-metrics', action='store_true',
                       help='List all available metrics and exit')
    parser.add_argument('--config', '-c',
                       help='Path to JSON file with per-metric plot configurations')
    parser.add_argument('--exclude-from-mean',
                       help='Comma-separated module IDs to exclude from mean ± 1σ calculation (e.g. PD-06,PD-07,PD-08)')
    parser.add_argument('--exclude-from-plot',
                       help='Comma-separated module IDs to exclude entirely from plots (default: none)')
    parser.add_argument('--auto-exclude-from-mean', action='store_true',
                       help='Exclude MAD-outlier modules from mean ± 1σ (still plotted)')
    parser.add_argument('--no-auto-exclude-from-mean', action='store_true',
                       help='Do not auto-exclude outlier modules from the mean')
    parser.add_argument('--outlier-sigma', type=float, default=DEFAULT_OUTLIER_SIGMA,
                       help=f'MAD sigma threshold for --auto-exclude-from-mean (default: {DEFAULT_OUTLIER_SIGMA})')
    parser.add_argument('--log-scale', action='store_true',
                       help='Use a logarithmic y-axis (already the default for Image31 DC plots)')
    parser.add_argument('--show', action='store_true',
                       help='Open each plot in an interactive matplotlib window (pan/zoom/axis limits). Close the window to save and continue')
    parser.add_argument('--show-module-averages', action='store_true',
                       help='Print the average value of the metric per module (mean across CCDs) to the screen')
    parser.add_argument('--image31-dc-plots', action='store_true',
                       help='Generate low/high dark-current and Low/High DC ratio plots (Image31 underground or Image3 surface)')

    args = parser.parse_args()
    
    # Read and process data
    try:
        df_long, metrics = read_and_process_data(args.csv_file)
    except Exception as e:
        print(f"Error reading data: {e}", file=sys.stderr)
        sys.exit(1)
    
    # List metrics if requested
    if args.list_metrics:
        print(f"\nAvailable metrics ({len(metrics)} total):")
        print("=" * 80)
        for i, metric in enumerate(metrics, 1):
            print(f"{i:3d}. {metric}")
        return
    
    # Determine which metrics to plot based on user input
    # --metric: User specifies exact metric name(s) (e.g., --metric Image31_High_Res)
    # --pattern: User specifies substring to match (e.g., --pattern Image31 matches all Image31 metrics)
    # --image31-dc-plots: Image31 low DC, high DC, and derived Low/High DC ratio
    metrics_to_plot = []
    if getattr(args, 'image31_dc_plots', False):
        chosen = dc_plot_metrics_for_csv(metrics)
        if not chosen:
            print("Error: --image31-dc-plots needs Image31_* or Image3_* Low/High Dark_Current columns",
                  file=sys.stderr)
            sys.exit(1)
        print(f"DC plots using metrics: {', '.join(chosen)}")
        metrics_to_plot.extend(chosen)
    if args.metrics:
        # Specific metrics: exact names provided by user
        # Example: --metric Image31_High_Res --metric Image31_Low_Res
        metrics_to_plot.extend(args.metrics)
    elif args.pattern and not metrics_to_plot:
        # Pattern matching: find all metrics containing the pattern (case-insensitive)
        # Example: --pattern Image31 matches "Image31_High_Res", "Image31_Low_Gain", etc.
        #          --pattern Res matches all metrics with "Res" in the name
        metrics_to_plot = [m for m in metrics if args.pattern.lower() in m.lower()]
        if not metrics_to_plot:
            print(f"No metrics found matching pattern: {args.pattern}")
            return
    if not metrics_to_plot:
        print("Error: Must specify --metric, --pattern, --image31-dc-plots, or --list-metrics")
        parser.print_help()
        sys.exit(1)
    
    # Parse exclude-from-mean list
    exclude_from_mean = []
    if getattr(args, 'image31_dc_plots', False):
        exclude_from_mean.extend(DEFAULT_EXCLUDE_FROM_MEAN)
    if getattr(args, 'exclude_from_mean', None):
        exclude_from_mean.extend(
            [m.strip() for m in args.exclude_from_mean.split(',') if m.strip()]
        )
    exclude_from_mean = list(dict.fromkeys(exclude_from_mean))
    if exclude_from_mean:
        print(f"Excluding from mean ± 1σ (still plotted): {', '.join(exclude_from_mean)}")
    exclude_from_plot = []
    if getattr(args, 'exclude_from_plot', None):
        exclude_from_plot = [m.strip() for m in args.exclude_from_plot.split(',') if m.strip()]

    auto_exclude_from_mean = bool(getattr(args, 'auto_exclude_from_mean', False))
    if getattr(args, 'no_auto_exclude_from_mean', False):
        auto_exclude_from_mean = False

    show = bool(getattr(args, 'show', False))
    if show:
        backend = ensure_interactive_backend()
        print(f"Interactive matplotlib backend: {backend}")

    # Generate plots
    try:
        plot_metrics(df_long, metrics_to_plot, args.output_dir, args.pattern, args.config,
                     exclude_from_mean=exclude_from_mean,
                     exclude_from_plot=exclude_from_plot,
                     show_module_averages=getattr(args, 'show_module_averages', False),
                     auto_exclude_from_mean=auto_exclude_from_mean,
                     outlier_sigma=getattr(args, 'outlier_sigma', DEFAULT_OUTLIER_SIGMA),
                     show=show,
                     log_scale=getattr(args, 'log_scale', False))
        print("\nDone!")
    except Exception as e:
        print(f"Error generating plots: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
