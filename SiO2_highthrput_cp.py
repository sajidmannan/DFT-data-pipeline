#!/usr/bin/env python3
import subprocess
import shutil
from pathlib import Path

# ================================================================
# 1. Configuration
# ================================================================
# Detect CP2K executable
for exe_name in ["cp2k.psmp", "cp2k.popt", "cp2k"]:
    path = shutil.which(exe_name)
    if path:
        CP2K_CMD = path
        break
else:
    raise FileNotFoundError("No CP2K executable found in PATH — please load CP2K before running.")

print(f"Using CP2K executable: {CP2K_CMD}")

BASE_DIR = Path.cwd()
ORIGINAL_XYZ = BASE_DIR / "SiO2.xyz"
CLEANED_XYZ = BASE_DIR / "SiO2_cleaned.xyz"

FUNCTIONAL = "PBE_D3"
CUTOFF = 600
REL_CUTOFF = 40

# ================================================================
# 2. XYZ Cleaning
# ================================================================
def clean_xyz_file(original_xyz: Path, cleaned_xyz: Path):
    lines = original_xyz.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError("XYZ file too short.")
    cleaned_xyz.write_text("\n".join(lines[2:]) + "\n")

# ================================================================
# 3. CP2K Blocks
# ================================================================
def generate_xc_block(functional: str):
    pot_func = functional.replace("_D3", "")
    xc_content = f"""    &XC
      &XC_FUNCTIONAL {pot_func}
      &END XC_FUNCTIONAL"""
    if "_D3" in functional:
        xc_content += """
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF 16.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL"""
    return xc_content

def generate_cell_block():
    return """    &CELL
      A 3.7528 0.0    0.0
      B 0.0    3.7528 0.0
      C 0.0    0.0    3.7528
      PERIODIC XYZ
    &END CELL
"""

def generate_subsys_block(xyz_file: str, functional: str):
    pot_func = functional.replace("_D3", "")
    return f"""{generate_cell_block()}
    &COORD
      @INCLUDE {xyz_file}
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-{pot_func}-q4
    &END KIND

    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-{pot_func}-q6
    &END KIND
"""

def cell_opt_block():
    """Cell optimization block with stress and force output"""
    return """&MOTION
  &CELL_OPT
    OPTIMIZER BFGS
    MAX_ITER 200
    MAX_DR 0.003
    MAX_FORCE 0.00045
  &END CELL_OPT

  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      FILENAME cell_opt-pos
    &END TRAJECTORY

    &FORCES
      FORMAT XYZ
      FILENAME cell_opt-forces
    &END FORCES

    &CELL
      FILENAME cell_opt-cell
      &EACH
        CELL_OPT 1
      &END EACH
    &END CELL

    &STRESS
      FILENAME cell_opt-stress
      &EACH
        CELL_OPT 1
      &END EACH
    &END STRESS
  &END PRINT
&END MOTION
"""

def nvt_md_block(temp: float, steps: int = 5000, timecon: float = 100):
    """NVT equilibration block with stress and force output"""
    return f"""&MOTION
  &MD
    ENSEMBLE NVT
    STEPS {steps}
    TIMESTEP 0.5
    TEMPERATURE {temp}
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON {timecon}
      &END CSVR
    &END THERMOSTAT
  &END MD

  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      FILENAME nvt-pos
    &END TRAJECTORY

    &VELOCITIES
      FORMAT XYZ
      FILENAME nvt-vel
    &END VELOCITIES

    &FORCES
      FORMAT XYZ
      FILENAME nvt-forces
    &END FORCES

    &CELL
      FILENAME nvt-cell
      &EACH
        MD 1
      &END EACH
    &END CELL

    &STRESS
      FILENAME nvt-stress
      &EACH
        MD 1
      &END EACH
    &END STRESS
  &END PRINT
&END MOTION
"""

def npt_md_block(temp: float, pressure: float, steps: int = 10000, timecon: float = 100):
    """NPT production block with stress and force output"""
    return f"""&MOTION
  &MD
    ENSEMBLE NPT_I
    STEPS {steps}
    TIMESTEP 0.5
    TEMPERATURE {temp}
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON {timecon}
      &END CSVR
    &END THERMOSTAT
    &BAROSTAT
      PRESSURE {pressure}
      TIMECON 500
    &END BAROSTAT
  &END MD

  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      FILENAME npt-pos
    &END TRAJECTORY

    &VELOCITIES
      FORMAT XYZ
      FILENAME npt-vel
    &END VELOCITIES

    &FORCES
      FORMAT XYZ
      FILENAME npt-forces
    &END FORCES

    &CELL
      FILENAME npt-cell
      &EACH
        MD 1
      &END EACH
    &END CELL

    &STRESS
      FILENAME npt-stress
      &EACH
        MD 1
      &END EACH
    &END STRESS
  &END PRINT
&END MOTION
"""

# ================================================================
# 4. Input Writer
# ================================================================
def write_input_file(run_type: str, dir_path: Path, encut: int, relcut: int, motion_section: str = ""):
    xc_block = generate_xc_block(FUNCTIONAL)
    subsys_block = generate_subsys_block(CLEANED_XYZ.name, FUNCTIONAL)
    input_content = f"""&GLOBAL
  PROJECT {run_type}_{FUNCTIONAL}_encut{encut}_relcut{relcut}
  RUN_TYPE {run_type}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
{xc_block}
    &END XC

    &POISSON
      PERIODIC XYZ
      POISSON_SOLVER PERIODIC
    &END POISSON

    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 50
    &END SCF

    &MGRID
      CUTOFF {encut}
      REL_CUTOFF {relcut}
      NGRIDS 4
    &END MGRID
  &END DFT

  &SUBSYS
{subsys_block}
  &END SUBSYS
&END FORCE_EVAL
{motion_section}
"""
    (dir_path / "input.inp").write_text(input_content)

# ================================================================
# 5. Run CP2K
# ================================================================
def run_cp2k(run_type: str, encut: int, relcut: int, motion_section: str = "", extra_tag: str = "", input_xyz: Path = None):
    folder_name = f"{run_type.lower()}_{FUNCTIONAL}_encut{encut}_relcut{relcut}{extra_tag}"
    dir_path = BASE_DIR / folder_name
    dir_path.mkdir(exist_ok=True)
    
    # Use specified input XYZ or fall back to cleaned XYZ
    xyz_to_use = input_xyz if input_xyz else CLEANED_XYZ
    shutil.copy2(xyz_to_use, dir_path / xyz_to_use.name)
    
    write_input_file(run_type, dir_path, encut, relcut, motion_section)
    subprocess.run([CP2K_CMD, "-i", "input.inp", "-o", "output.out"], cwd=dir_path, check=True)
    
    return dir_path

# ================================================================
# 6. Main Workflow
# ================================================================
def main():
    print("Cleaning XYZ file...")
    clean_xyz_file(ORIGINAL_XYZ, CLEANED_XYZ)

    # Step 1: Cell Optimization (replaces Geometry Optimization)
    print("="*50)
    print("STEP 1: Running Cell Optimization")
    print("="*50)
    cell_opt_dir = run_cp2k("CELL_OPT", CUTOFF, REL_CUTOFF, cell_opt_block(), 
                          extra_tag="_cellopt")
    
    # Get the final optimized structure
    optimized_xyz = cell_opt_dir / "cell_opt-pos-1.xyz"
    if not optimized_xyz.exists():
        print("Warning: Optimized structure not found. Using initial structure for MD.")
        optimized_xyz = CLEANED_XYZ

    # Step 2: NVT Equilibration
    print("="*50)
    print("STEP 2: Running NVT Equilibration")
    print("="*50)
    nvt_temp = 300
    nvt_dir = run_cp2k("MD", CUTOFF, REL_CUTOFF, nvt_md_block(nvt_temp, steps=5000), 
                      extra_tag=f"_NVT_T{nvt_temp}", input_xyz=optimized_xyz)

    # Step 3: NPT Production Run
    print("="*50)
    print("STEP 3: Running NPT Production")
    print("="*50)
    npt_temp = 300
    npt_pressure = 1.0
    run_cp2k("MD", CUTOFF, REL_CUTOFF, npt_md_block(npt_temp, npt_pressure, steps=10000), 
            extra_tag=f"_NPT_T{npt_temp}_P{npt_pressure}", input_xyz=optimized_xyz)

    print("="*50)
    print("All calculations completed successfully!")
    print("="*50)
    print("Output files generated:")
    print(f"  Cell Optimization: {cell_opt_dir}")
    print(f"  NVT Equilibration: {nvt_dir}")
    print("  NPT Production: npt_... folder")
    print("\nEach folder contains:")
    print("  - Position trajectories (pos-*.xyz)")
    print("  - Force trajectories (forces-*.xyz)") 
    print("  - Stress tensor trajectories (stress-*.dat)")
    print("  - Cell parameter trajectories (cell-*.dat)")

if __name__ == "__main__":
    main()
