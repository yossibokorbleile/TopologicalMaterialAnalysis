##
# @internal
# @file max_flow_cli.py
# @cli script for computing circular max flow.
# @version 1.3.1
# @date May 2026
# @author: Yossi Bokor Bleile
# @author: Matteo Pegoraro
# @author: Rasmus Christensen
# @copyright 2026 GPL


"""
Max Flow Analysis — config-driven CLI
======================================
Reads all parameters from a YAML config file so the same script runs
unchanged across different systems, backbones, and trajectories.

Usage
-----
    python max_flow_cli.py --config nasio.yaml
    python max_flow_cli.py --config lips.yaml --verbose   # override verbosity
    python max_flow_cli.py --help

YAML keys
---------
Required
    trajectory_path   str        Path to the LAMMPS dump file.
    atom_type_map     dict       LAMMPS int type → element symbol.
    backbone_atoms    list[str]  Elements forming the rigid backbone.
    flow_atoms        list[str]  Mobile ion elements.

Optional — trajectory slicing (applied before analysis stride)
    frame_start  int   First frame index to load (default 0).
    frame_end    int   Last frame index, exclusive (default: all frames).
    frame_step   int   Read every Nth frame from the slice (default 1).

Optional — calculation
    grid_spacing      float         Reeb-graph resolution in Å (default 0.3).
    fat               float         Backbone construction factor (default 1.0).
    reeb_stride       int           Reeb-graph sub-sampling stride (default 2).
    stride            int           Analysis frame sub-sampling (default 20).
    periodic          bool          Use PBC (default true).
    common_backbone   bool          Same radii per atom type (default true).
    apply_frechet_mean bool         Use Fréchet mean for backbone, instead of Arithmetic mean. Computational cost is higher. (default false).
    print_mean_stucture bool         Print the mean structure to a file. (default false)
    repeat            [int,int,int] Supercell repeat (default [1,1,1]).

Optional — output
    output_csv          str   CSV file path (default "flow_results.csv").
    verbose             bool  Print progress (default false).
    save_reeb_backbone  bool  Serialise Reeb backbone to disk (default false).
"""

import argparse
import sys
from dataclasses import dataclass, field, fields
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from ase.io import read

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from max_flow import compute_backbone, compute_max_flow


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
@dataclass
class FlowConfig:
    """All parameters for one max-flow run, loaded from YAML."""

    # Required
    trajectory_path: str = ""
    file_format: str = "lammps-dump-text"
    atom_type_map: dict = field(default_factory=dict)
    backbone_atoms: list = field(default_factory=lambda: ["Si", "O"])
    flow_atoms: list = field(default_factory=lambda: ["Na"])

    # Trajectory slicing
    frame_start: int = 0
    frame_end: int | None = None  # None → load to end
    frame_step: int = 1

    # Calculation
    grid_spacing: float = 0.3
    fat: float = 1.0
    reeb_stride: int = 2
    stride: int = 20
    periodic: bool = True
    common_backbone: bool = True
    apply_frechet_mean: bool = False
    print_mean_stucture: bool = False
    repeat: tuple = (1, 1, 1)

    # Output
    output_csv: str = "flow_results.csv"
    verbose: bool = False
    save_reeb_backbone: bool = False
    radii_figures: bool = False


def load_config(yaml_path: str, cli_overrides: dict) -> FlowConfig:
    """Parse YAML and apply any CLI flag overrides."""
    with open(yaml_path) as fh:
        raw = yaml.safe_load(fh)

    # atom_type_map keys must be ints (YAML may parse them as strings)
    raw["atom_type_map"] = {int(k): v for k, v in raw["atom_type_map"].items()}

    # repeat as tuple
    if "repeat" in raw:
        raw["repeat"] = tuple(raw["repeat"])

    valid_keys = {f.name for f in fields(FlowConfig)}
    cfg = FlowConfig(**{k: raw[k] for k in raw if k in valid_keys})

    # CLI flags win over YAML
    for key, val in cli_overrides.items():
        if val is not None:
            setattr(cfg, key, val)

    return cfg


# ---------------------------------------------------------------------------
# Trajectory loading
# ---------------------------------------------------------------------------


def remap_atom_types(frames: list, index_to_symbol: dict) -> list:
    """Replace numeric atom types with chemical symbols in every frame.

    ASE returns np.int64 from get_atomic_numbers(); cast to plain int so
    the lookup works regardless of how the map keys were constructed.
    """
    for atoms in frames:
        try:
            symbols = [index_to_symbol[int(i)] for i in atoms.get_atomic_numbers()]
        except KeyError as exc:
            known = sorted(index_to_symbol.keys())
            raise KeyError(
                f"Atom type {exc} not found in atom_type_map. "
                f"Mapped types: {known}. Check your YAML atom_type_map."
            ) from exc
        atoms.set_chemical_symbols(symbols)
    return frames


def load_trajectory(cfg: FlowConfig) -> list:
    """Read a LAMMPS dump, apply frame slice, wrap, and fix atom labels."""
    # Build ASE index string for the requested slice
    start = cfg.frame_start
    stop = cfg.frame_end if cfg.frame_end is not None else ""
    step = cfg.frame_step
    index = f"{start}:{stop}:{step}"

    if cfg.verbose:
        label = f"frames {start}:{stop or 'end'}" + (
            f" step {step}" if step != 1 else ""
        )
        print(f"Loading trajectory: {cfg.trajectory_path}  ({label})")

    frames = read(cfg.trajectory_path, format=cfg.file_format, index=index)

    if not isinstance(frames, list):
        frames = [frames]

    [atom.wrap() for atom in frames]
    remap_atom_types(frames, cfg.atom_type_map)

    if cfg.verbose:
        print(f"  Loaded {len(frames)} frames.")

    return frames


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------


def run_flow_analysis(cfg: FlowConfig) -> dict:
    """Run the full max-flow analysis and return a results dict."""
    frames = load_trajectory(cfg)

    present_types = set(frames[0].get_chemical_symbols())
    backbone_present = [el for el in cfg.backbone_atoms if el in present_types]
    flow_present = [el for el in cfg.flow_atoms if el in present_types]

    if cfg.verbose:
        print(f"  Backbone elements found: {backbone_present}")
        print(f"  Flow elements found:     {flow_present}")

    mean_volume = np.mean([f.get_volume() for f in frames])

    backbone_mean, backbone_radii, cell, m, rdf_figures = compute_backbone(
        frames,
        backbone_present,
        flow_present,
        cfg.reeb_stride,
        common_backbone=cfg.common_backbone,
        repeat=cfg.repeat,
        apply_frechet_mean=cfg.apply_frechet_mean,
        print_mean_stucture=cfg.print_mean_stucture,
        radii_figures=cfg.radii_figures,
    )

    flow = compute_max_flow(
        backbone_mean,
        backbone_radii,
        cell,
        m,
        cfg.grid_spacing,
        mean_volume,
        cfg.fat,
        cfg.reeb_stride,
        cfg.stride,
        cfg.periodic,
        cfg.save_reeb_backbone,
        verbose=cfg.verbose,
    )

    slice_area = cell[0] * cell[1]

    output = {
        "flow": flow,
        "slice_area": slice_area,
        "norm_flow": flow / slice_area,
    }

    if cfg.radii_figures:
        output["rdf_figures"] = rdf_figures
    
    return output


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description="Config-driven max-flow analysis for LAMMPS trajectories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--config",
        "-c",
        required=True,
        metavar="FILE",
        help="Path to a YAML config file.",
    )
    # Optional CLI overrides (all optional; they win over YAML values)
    p.add_argument(
        "--trajectory", metavar="FILE", help="Override trajectory_path from YAML."
    )
    p.add_argument("--output", metavar="FILE", help="Override output_csv from YAML.")
    p.add_argument("--frame-start", type=int, metavar="N", help="Override frame_start.")
    p.add_argument("--frame-end", type=int, metavar="N", help="Override frame_end.")
    p.add_argument("--frame-step", type=int, metavar="N", help="Override frame_step.")
    p.add_argument("--stride", type=int, metavar="N", help="Override analysis stride.")
    p.add_argument(
        "--verbose",
        action="store_true",
        default=None,
        help="Enable verbose output (overrides YAML).",
    )
    return p.parse_args()


def main():
    args = parse_args()

    cli_overrides = {
        "trajectory_path": args.trajectory,
        "output_csv": args.output,
        "frame_start": args.frame_start,
        "frame_end": args.frame_end,
        "frame_step": args.frame_step,
        "stride": args.stride,
        "verbose": args.verbose if args.verbose else None,
    }

    cfg = load_config(args.config, cli_overrides)

    if cfg.verbose:
        print("=" * 60)
        print(f"Config : {args.config}")
        print(f"System : backbone={cfg.backbone_atoms}  flow={cfg.flow_atoms}")
        print(f"Output : {cfg.output_csv}")
        print("=" * 60)

    results = run_flow_analysis(cfg)

    if cfg.radii_figures:
        for fig in results.get("rdf_figures", []):
            plt.close(fig)

    pd.DataFrame([results]).rename(
        columns={
            "flow": "Flow",
            "slice_area": "Slice Area",
            "norm_flow": "Normalized Flow",
        }
    ).to_csv(cfg.output_csv, index=False)

    print(f"\nResults saved to '{cfg.output_csv}'")
    print(f"  Flow            : {results['flow']:.4f}")
    print(f"  Slice area      : {results['slice_area']:.4f} Å²")
    print(f"  Normalised flow : {results['norm_flow']:.6f}")


if __name__ == "__main__":
    main()
