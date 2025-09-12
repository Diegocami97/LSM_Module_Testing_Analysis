# usecase: wildcard pattern search, take the first 4 matches
# python plot_pdfs.py --files "trace_*.pdf" outputname
#
# compatible to the old version: take the first 4 files in a directory ()
# Example: python plot_pdfs.py --directory ./pdf_outputs outputname
import pdfplumber
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
import argparse
from natsort import natsorted
import glob
import os
import re


def plot_pdf(pdf_paths, ccd_name, output_dir):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    axs = axs.flatten()
    
    for i, pdf_path in enumerate(pdf_paths):
        if i >= 4:  # We only have 4 subplots, so limit to 4 files
            break
        ax = axs[i]
        with pdfplumber.open(pdf_path) as pdf:
            page = pdf.pages[0]
        im = page.to_image(resolution=600)
        image_np = np.array(im.original)
        ax.imshow(image_np, cmap='gray')
        ax.axis('off')
        ax.set_title(pdf_path.split('/')[-1],fontweight='bold',fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir + '/' + f'{ccd_name}.png', dpi=600)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Process and plot PDF files.")
    parser.add_argument('--files', nargs=4, metavar='PDF', help="Specify 4 individual PDF files.")
    parser.add_argument('--input_dir', type=str, help="Directory containing PDF files.")
    parser.add_argument('--output_dir', type=str, help="Path for output directory.")
    parser.add_argument('--module', type=str, help="Name of the CCD module for output image.")
    
    args = parser.parse_args()

    pdf_paths = []
    if args.files:
        # Ensure exactly 4 files are provided
        if len(args.files) != 4:
            sys.exit("You must provide exactly 4 PDF files.")
        pdf_paths = [Path(pdf) for pdf in args.files]
        for pdf_path in pdf_paths:
            if not pdf_path.exists():
                sys.exit(f"File {pdf_path} does not exist.")
    elif args.input_dir:
        # Get all PDFs from the directory, limit to 4
        pdf_dir = Path(args.input_dir)
        if not pdf_dir.exists() or not pdf_dir.is_dir():
            sys.exit(f"Directory {pdf_dir} does not exist or is not a directory.")
        pdf_paths = []
        for filename in natsorted(os.listdir(pdf_dir)):
            f = os.path.join(pdf_dir, filename)
            if os.path.isfile(f) and filename.lower().endswith('.pdf'):
                pdf_paths.append(f)    
        #pdf_paths = list(pdf_dir.glob("*.pdf"))[:4]
        if len(pdf_paths) < 4:
            sys.exit("The directory must contain at least 4 PDF files.")
    else:
        sys.exit("You must provide either 4 PDF files or a directory with PDFs.")

    # Plot the PDFs
    plot_pdf(pdf_paths, args.module,args.output_dir)

if __name__ == "__main__":
    main()