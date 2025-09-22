#!/usr/bin/env python3
import os
import argparse
from ase.io import write
from ase.io.vasp import read_vasp_out


def outcar_to_xyz(outcar_path, xyz_out):
    """Convert VASP OUTCAR trajectory into XYZ for MLIP training."""
    try:
        images = read_vasp_out(outcar_path, index=":")
    except Exception as e:
        raise RuntimeError(f"Failed to read {outcar_path}: {e}")
    
    if not images:
        raise RuntimeError(f"No frames found in {outcar_path}")
    
    write(xyz_out, images, format="extxyz")
    print(f"{os.path.basename(outcar_path)} → {xyz_out} ({len(images)} frames)")


def find_outcars(root_dir):
    """Recursively find all OUTCAR files under root_dir."""
    outcar_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname == "OUTCAR":
                outcar_files.append(os.path.join(dirpath, fname))
    return outcar_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Recursively convert all VASP OUTCAR files to XYZ"
    )
    parser.add_argument("rootdir", help="Root directory to search for OUTCAR files")
    parser.add_argument(
        "-o", "--outdir",
        default="xyz_outputs",
        help="Directory where XYZ files will be saved (default: xyz_outputs)"
    )

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    outcar_files = find_outcars(args.rootdir)
    if not outcar_files:
        print("No OUTCAR files found.")
        exit(0)

    for outcar in outcar_files:
        rel_path = os.path.relpath(os.path.dirname(outcar), args.rootdir)
        parts = rel_path.split(os.sep)

        # avoid duplication if top and last folder are same
        # take top-level mineral folder + last subfolder only, avoid repeating
        if len(parts) >= 2:
            if parts[-1] in parts[0]:
                filename = f"{parts[0]}.xyz"
            else:
                filename = f"{parts[0]}_{parts[-1]}.xyz"
        else:
            filename = f"{parts[0]}.xyz"


        xyz_file = os.path.join(args.outdir, filename)

        try:
            outcar_to_xyz(outcar, xyz_file)
        except Exception as e:
            print(f"Skipping {outcar}: {e}")
