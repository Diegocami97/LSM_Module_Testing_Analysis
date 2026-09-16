#!/usr/bin/env python3
# Diego Venegas-Vargas
# DAMIC-M Collaboration
# Loads an ACM image with acmpy Decoder and writes its Level-1 dataframe to CSV.
"""Load ACM image with acmpy Decoder and expose L1 dataframe (df1_view).

CSV output columns (Level-1 decoded data, one row per pixel/sample):
  pix   : pixel index (column position in sensor)
  line  : line index (row position in sensor)
  skip  : sample index within pedestal (ramp=0) or signal (ramp=1) phase
  ramp  : 0 = pedestal (RD), 1 = signal (RU)
  idx   : cyclic index 0–127 from firmware
  n     : number of ADC samples in the integration window
  sum   : sum of ADC values (raw; CDS = signal_sum - pedestal_sum per pixel)
"""

import sys
import os

from acmpy import Decoder

# Set your image path here, or pass as first command-line argument
img_path = sys.argv[1] if len(sys.argv) > 1 else "/path/to/img"

data = Decoder(img_path)
L1_df = data.df1_view

# Save full dataframe to CSV (same location as input, with _L1.csv suffix)
base, _ = os.path.splitext(img_path)
out_csv = base + "_L1.csv"
print("L1_df shape:", L1_df.shape)
print("CSV columns:", list(L1_df.columns))
print("Writing to", out_csv, "...")
L1_df.to_csv(out_csv, index=False)
print("Done.")
