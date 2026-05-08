##
# @internal
# @file Circular_Max_Flow.py
# @brief Streamlit page for computing circular max flow.
# @version 1.3.1
# @date May 2026
# @author: Yossi Bokor Bleile
# @author: Matteo Pegoraro
# @author: Rasmus Christensen


import streamlit as st
from utils import *

from max_flow_cli import *

st.header("Circular Max Flow")

st.markdown("""
This page allows you to compute the circular max flow of a configuration, as implemented in the [Circular Max Flow](https://github.com/pego91C/ircular_Max_Flow) repository.
""")

st.markdown("""
It creates a static approximation of the backbone, by taking approximate locations of atoms in the backbone, and for each atom type a radius to thicken these atoms by. 
""")

st.checkbox(
    "Load settings from yaml config file (overrides manual settings)",
    key="load_config",
    value=False,
)

if st.session_state.load_config:
    st.text_input("Config file", key="config_file", placeholder="path/to/config.yaml")

st.subheader("Trajectory")

st.text_input("Input file", key="input_file", placeholder="path/to/input.xyz (any format supported by ase.io.read)")
st.text_input("File format (Plain text of given to ase.io.read):", key="file_format",placeholder="Auto") #specify format of the initial strutcure file

st.text_area(
    "Atom type map  (one  `<int>: <symbol>`  pair per line)",
    key="atom_type_map_input",
    placeholder="1: Si\n2: O\n3: Na",
    help="Maps the integer atom types in the LAMMPS dump to chemical symbols.",
)

col_start, col_end, col_step = st.columns(3)
with col_start:
    st.number_input("Frame start", key="frame_start", value=0, min_value=0, step=1)
with col_end:
    st.text_input("Frame end", key="frame_end_input", placeholder="(all frames)")
with col_step:
    st.number_input("Frame step", key="frame_step", value=1, min_value=1, step=1)

st.subheader("Atom selection")

col_bb, col_fl = st.columns(2)

with col_bb:
    st.text_input(
        "Backbone atoms",
        key="backbone_atoms_input",
        placeholder="Si, O",
        help="Comma-separated list of backbone element symbols.",
    )
with col_fl:
    st.text_input(
        "Flow atoms",
        key="flow_atoms_input",
        placeholder="Na",
        help="Comma-separated list of mobile-ion element symbols.",
    )

st.subheader("Calculation settings")

col_gs, col_fat = st.columns(2)
col_gs, col_fat = st.columns(2)
with col_gs:
    st.number_input(
        "Grid spacing (Å)",
        key="grid_spacing",
        value=0.3,
        min_value=0.01,
        step=0.05,
        format="%.3f",
        help="Reeb-graph resolution in Ångströms.",
    )
with col_fat:
    st.number_input(
        "Fat radius",
        key="fat",
        value=1.0,
        min_value=0.0,
        step=0.1,
        format="%.2f",
        help="Backbone construction factor.",
    )
 
col_rs, col_st = st.columns(2)
with col_rs:
    st.number_input(
        "Reeb stride",
        key="reeb_stride",
        value=2,
        min_value=1,
        step=1,
        help="Sub-sampling stride used when building the Reeb graph.",
    )
with col_st:
    st.number_input(
        "Analysis stride",
        key="stride",
        value=20,
        min_value=1,
        step=1,
        help="Sub-sampling stride used over the full trajectory.",
    )
 
col_rep1, col_rep2, col_rep3 = st.columns(3)
with col_rep1:
    st.number_input("Repeat x", key="repeat_x", value=1, min_value=1, step=1)
with col_rep2:
    st.number_input("Repeat y", key="repeat_y", value=1, min_value=1, step=1)
with col_rep3:
    st.number_input("Repeat z", key="repeat_z", value=1, min_value=1, step=1)
 
st.checkbox(
    "Periodic boundary conditions",
    key="periodic",
    value=True,
)
st.checkbox(
    "Common backbone radii (same radius per atom type)",
    key="common_backbone",
    value=True,
)

if st.session_state.common_backbone:
    st.checkbox(
    "Generate RDF figures with cutoff radius estimation",
    key="radii_figures",
    value=False,
)

    
st.checkbox(
    "Apply Fréchet mean for backbone  (slower, more accurate)",
    key="apply_frechet_mean",
    value=False,
)

st.subheader("Output")
 
st.text_input(
    "Output CSV",
    key="output_csv",
    value="flow_results.csv",
    help="Results are always written here after a successful run.",
)
st.checkbox("Verbose output", key="verbose", value=False)
st.checkbox("Save Reeb backbone to disk", key="save_reeb_backbone", value=False)

def _parse_atom_type_map(raw: str) -> dict:
    """Parse 'int: symbol' lines into {int: str}."""
    result = {}
    for line in raw.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        k, _, v = line.partition(":")
        result[int(k.strip())] = v.strip()
    return result

def _parse_atom_list(raw: str) -> list:
    return [s.strip() for s in raw.split(",") if s.strip()]

def _build_config() -> FlowConfig:
    frame_end_raw = st.session_state.frame_end_input.strip()
    frame_end = int(frame_end_raw) if frame_end_raw else None
 
    return FlowConfig(
        trajectory_path=st.session_state.trajectory_path,
        file_format=st.session_state.file_format,
        atom_type_map=_parse_atom_type_map(st.session_state.atom_type_map_input),
        backbone_atoms=_parse_atom_list(st.session_state.backbone_atoms_input),
        flow_atoms=_parse_atom_list(st.session_state.flow_atoms_input),
        frame_start=int(st.session_state.frame_start),
        frame_end=frame_end,
        frame_step=int(st.session_state.frame_step),
        grid_spacing=float(st.session_state.grid_spacing),
        fat=float(st.session_state.fat),
        reeb_stride=int(st.session_state.reeb_stride),
        stride=int(st.session_state.stride),
        periodic=st.session_state.periodic,
        common_backbone=st.session_state.common_backbone,
        apply_frechet_mean=st.session_state.apply_frechet_mean,
        print_mean_stucture=st.session_state.print_mean_stucture,
        repeat=(
            int(st.session_state.repeat_x),
            int(st.session_state.repeat_y),
            int(st.session_state.repeat_z),
        ),
        output_csv=st.session_state.output_csv,
        verbose=st.session_state.verbose,
        save_reeb_backbone=st.session_state.save_reeb_backbone,
        radii_figures=st.session_state.radii_figures,)

def compute_circular_max_flow():
    try:       
        if st.session_state.load_config:
            try:
                cfg = load_config(yaml_path=st.session_state.config_file, cli_overrides={})
                st.success(f"Loaded config from {st.session_state.config_file}")
            except Exception as exc:
                st.error(f"Failed to load config: {exc}")
                return
        else:
            cfg = _build_config()

        st.session_state.last_cfg = cfg
        results = run_flow_analysis(cfg)
 
        # Save to CSV
        df = pd.DataFrame([results]).rename(
            columns={
                "flow": "Flow",
                "slice_area": "Slice Area",
                "norm_flow": "Normalized Flow",
            }
        )
        df.to_csv(cfg.output_csv, index=False)
 
        st.session_state.results = results
        st.session_state.results_csv = df.to_csv(index=False)
        st.session_state.run_error = None
 
    except Exception as exc:
        st.session_state.results = None
        st.session_state.run_error = str(exc)


def load_test_values():
    """Populate inputs with the standard test case."""
    st.session_state.trajectory_path = "../examples/max_flow/md_stats.lammpstrj"
    st.session_state.file_format = "lammps-dump-text"
    st.session_state.atom_type_map_input = "1: Si\n5: O\n8: Na"
    st.session_state.backbone_atoms_input = "Si, O"
    st.session_state.flow_atoms_input = "Na"
    st.session_state.frame_start = 0
    st.session_state.frame_end_input = ""
    st.session_state.frame_step = 1
    st.session_state.grid_spacing = 0.2
    st.session_state.fat = 1.0
    st.session_state.reeb_stride = 2
    st.session_state.stride = 20
    st.session_state.periodic = True
    st.session_state.common_backbone = True
    st.session_state.apply_frechet_mean = True
    st.session_state.print_mean_stucture = False
    st.session_state.repeat_x = 1
    st.session_state.repeat_y = 1
    st.session_state.repeat_z = 1
    st.session_state.radii_figures = True
    st.session_state.output_csv = "./examples/max_flow/flow_example.csv"


col_test, col_run = st.columns([1, 3])
with col_test:
    st.button("Load test values", on_click=load_test_values)
with col_run:
    st.button("Compute circular max flow", on_click=compute_circular_max_flow, type="primary")
    

if "run_error" in st.session_state and st.session_state.run_error:
    st.error(f"Run failed: {st.session_state.run_error}")
    

if "results" in st.session_state and st.session_state.results:
    r = st.session_state.results
    st.success("Computation complete.")
 
    col_f, col_a, col_n = st.columns(3)
    col_f.metric("Flow", f"{r['flow']:.4f}")
    col_a.metric("Slice area (Å^2)", f"{r['slice_area']:.4f}")
    col_n.metric("Normalised flow", f"{r['norm_flow']:.6f}")
 
    cfg = st.session_state.last_cfg
    st.caption(f"Results saved to `{cfg.output_csv}`")
 
    st.download_button(
        label="Download CSV",
        data=st.session_state.results_csv,
        file_name=cfg.output_csv,
        mime="text/csv",
    )
    
    rdf_figures = r.get("rdf_figures", [])
    if rdf_figures:
        st.subheader("RDF onset radii")
        st.markdown(
			"The dashed red line marks the onset of the first RDF peak, "
			"used as the backbone radius for each atom type."
		)
        cols_per_row = min(len(rdf_figures), 3)
        for row_start in range(0, len(rdf_figures), cols_per_row):
            row_figs = rdf_figures[row_start : row_start + cols_per_row]
            cols = st.columns(len(row_figs))
            for col, fig in zip(cols, row_figs):
                col.pyplot(fig)