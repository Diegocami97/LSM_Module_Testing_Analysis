#!/usr/bin/env python3
import sys, glob, argparse
from PIL import Image

def main():
    parser = argparse.ArgumentParser(description="Make a 2x2 PNG collage from 4 images.")
    parser.add_argument("pattern", help="Glob pattern for input PNGs, e.g. 'trace_ch3_109_*.png'")
    parser.add_argument("--outname", required=True, help="naming output file")
    args = parser.parse_args()

    files = sorted(glob.glob(args.pattern))
    print("Input filenames:", files)
    if len(files) < 4:
        raise SystemExit(f"Need at least 4 files, found {len(files)} for pattern: {args.pattern}")

    imgs = [Image.open(f).convert("RGBA") for f in files[:4]]
    w, h = imgs[0].size

    for f, im in zip(files[:4], imgs):
        if im.size != (w, h):
            raise SystemExit(f"Image {f} size {im.size} != {w}x{h}")

    out = Image.new("RGBA", (2*w, 2*h))
    out.paste(imgs[0], (0,   0))
    out.paste(imgs[1], (w,   0))
    out.paste(imgs[2], (0,   h))
    out.paste(imgs[3], (w,   h))

    out_name = f"trace_{args.outname}.png"
    out.convert("RGB").save(out_name, "PNG")
    print(f"output filenames: {out_name}")

if __name__ == "__main__":
    main()
