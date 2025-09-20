#!/usr/bin/env python3
"""
CP2K Simulation Utilities Module

This module contains reusable functions for CP2K simulations including:
- File handling and XYZ processing
- CP2K block generators
- Input file creation
- Job execution utilities
"""

import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional

# ================================================================
# System and File Utilities
# ================================================================

def detect_cp2k() -> str:
    """Detect CP2K executable in PATH"""
    for exe_name in ["cp2k.psmp", "cp2k.popt", "cp2k"]:
        path = shutil.which(exe_name)
        if path:
            return path
    raise FileNotFoundError("No CP2K executable found in PATH — please load CP2K before running.")

def find_xyz_file(mineral_name: str) -> Path:
    """Find XYZ file based on mineral name or path"""
    # If it's already a path to an existing file
    if Path(mineral_name).exists():
        return Path(mineral_name)
    
    # Try common extensions and patterns
    base_name = Path(mineral_name).stem
    possible_files = [
        f"{mineral_name}.xyz",
        f"{base_name}.xyz",
        f"{mineral_name.upper()}.xyz",
        f"{mineral_name.lower()}.xyz",
        mineral_name  # in case it's already the full filename
    ]
    
    for filename in possible_files:
        if Path(filename).exists():
            return Path(filename)
    
    raise FileNotFoundError(f"Could not find XYZ file for mineral '{mineral_name}'. "
                          f"Tried: {', '.join(possible_files)}")

def clean_xyz_file(original_xyz: Path, cleaned_xyz: Path) -> None:
    """Clean XYZ file by removing first two lines (atom count and comment)"""
    lines = original_xyz.read_text().strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ file {original_xyz} is too short (less than 3 lines).")
    
    # Remove atom count and comment lines, keep only coordinates
    coord_lines = lines[2:]
    cleaned_xyz.write_text("\n".join(coord_lines) + "\n")

def parse_cell_parameters(cell_str: str, mineral_name: str) -> Tuple[float, float, float]:
    """Parse cell parameters from string or auto-detect based on mineral"""
    if cell_str.lower() == "auto":
        # Mineral-specific default cell parameters (in Angstroms)
        mineral_defaults = {
            "sio2": "4.916 4.916 5.405",      # alpha-quartz
            "quartz": "4.916 4.916 5.405",
            "caco3": "4.99 4.99 17.06",       # calcite
            "calcite": "4.99 4.99 17.06",
            "mg2sio4": "4.754 10.207 5.980",  # forsterite
            "forsterite": "4.754 10.207 5.980",
            "al2o3": "4.759 4.759 12.991",    # corundum
            "corundum": "4.759 4.759 12.991",
            "fe2o3": "5.034 5.034 13.747",    # hematite
            "hematite": "5.034 5.034 13.747",
        }
        
        mineral_key = mineral_name.lower().replace(".xyz", "")
        if mineral_key in mineral_defaults:
            cell_str = mineral_defaults[mineral_key]
            print(f"Using default cell parameters for {mineral_name}: {cell_str}")
        else:
            # Generic cubic cell
            cell_str = "5.0 5.0 5.0"
            print(f"Warning: Using generic cubic cell (5.0 5.0 5.0) for {mineral_name}")
    
    try:
        a, b, c = map(float, cell_str.split())
        return a, b, c
    except ValueError:
        raise ValueError(f"Invalid cell parameters: {cell_str}. Expected format: 'a b c'")

# ================================================================
# CP2K Block Generators
# ================================================================

def generate_xc_block(functional: str) -> str:
    """Generate XC functional block for CP2K input"""
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
    
    xc_content += "\n    &END XC"
    return xc_content

def generate_cell_block(a: float, b: float, c: float) -> str:
    """Generate cell block for CP2K input"""
    return f"""    &CELL
      A {a:.4f} 0.0    0.0
      B 0.0    {b:.4f} 0.0
      C 0.0    0.0    {c:.4f}
      PERIODIC XYZ
    &END CELL
"""

def generate_kind_blocks(functional: str) -> str:
    """Generate atomic kind blocks for common elements in minerals"""
    pot_func = functional.replace("_D3", "")
    
    # Common elements in minerals with their valence electrons
    elements = {
        "Si": 4, "O": 6, "Ca": 10, "C": 4, "Mg": 10, 
        "Al": 3, "Fe": 8, "Na": 9, "K": 9, "Ti": 12,
        "Mn": 7, "Cr": 6, "Ni": 18, "Cu": 11, "Zn": 12
    }
    
    kind_blocks = ""
    for element, valence in elements.items():
        kind_blocks += f"""    &KIND {element}
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-{pot_func}-q{valence}
    &END KIND
    
"""
    return kind_blocks.rstrip()

def generate_subsys_block(xyz_file: str, functional: str, a: float, b: float, c: float) -> str:
    """Generate complete subsystem block"""
    cell_block = generate_cell_block(a, b, c)
    kind_blocks = generate_kind_blocks(functional)
    
    return f"""{cell_block}
    &COORD
      @INCLUDE {xyz_file}
    &END COORD
    
{kind_blocks}
"""

def generate_motion_block(motion_type: str, **kwargs) -> str:
    """Generate motion blocks for different calculation types"""
    
    if motion_type == "CELL_OPT":
        return """
&MOTION
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
    
    elif motion_type == "NVT":
        temp = kwargs.get('temperature', 300.0)
        steps = kwargs.get('steps', 5000)
        timecon = kwargs.get('timecon', 100.0)
        
        return f"""
&MOTION
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
    
    elif motion_type == "NPT":
        temp = kwargs.get('temperature', 300.0)
        pressure = kwargs.get('pressure', 1.0)
        steps = kwargs.get('steps', 10000)
        timecon = kwargs.get('timecon', 100.0)
        
        return f"""
&MOTION
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
    
    else:
        raise ValueError(f"Unknown motion type: {motion_type}")

# ================================================================
# Input File Generation
# ================================================================

def create_cp2k_input(run_type: str, mineral_name: str, dir_path: Path, 
                     functional: str, cutoff: int, rel_cutoff: int,
                     cell_params: Tuple[float, float, float],
                     motion_section: str = "") -> None:
    """Create complete CP2K input file"""
    
    xc_block = generate_xc_block(functional)
    a, b, c = cell_params
    subsys_block = generate_subsys_block("coords.xyz", functional, a, b, c)
    
    project_name = f"{mineral_name}_{run_type}_{functional}"
    
    input_content = f"""&GLOBAL
  PROJECT {project_name}
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
      CUTOFF {cutoff}
      REL_CUTOFF {rel_cutoff}
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
# Job Execution
# ================================================================

def run_cp2k_calculation(job_name: str, job_dir: Path, cp2k_cmd: str) -> dict:
    """Execute CP2K calculation and return run information"""
    
    print(f"Running {job_name} calculation in {job_dir}")
    start_time = datetime.now()
    
    try:
        result = subprocess.run(
            [cp2k_cmd, "-i", "input.inp", "-o", "output.out"], 
            cwd=job_dir, 
            check=True,
            capture_output=True,
            text=True
        )
        
        end_time = datetime.now()
        runtime = end_time - start_time
        print(f"✓ {job_name} completed successfully in {runtime}")
        
        run_info = {
            'job_name': job_name,
            'start_time': start_time,
            'end_time': end_time,
            'runtime': runtime,
            'success': True,
            'command': ' '.join([cp2k_cmd, '-i', 'input.inp', '-o', 'output.out'])
        }
        
        return run_info
        
    except subprocess.CalledProcessError as e:
        print(f"✗ {job_name} failed!")
        print(f"Error: {e.stderr}")
        
        run_info = {
            'job_name': job_name,
            'start_time': start_time,
            'end_time': datetime.now(),
            'runtime': datetime.now() - start_time,
            'success': False,
            'error': str(e),
            'stderr': e.stderr
        }
        
        return run_info

def write_run_info(job_dir: Path, run_info: dict, mineral_name: str, 
                  temp: float, pressure: float, functional: str, 
                  cutoff: int) -> None:
    """Write run information to file"""
    
    info_file = job_dir / "run_info.txt"
    
    with open(info_file, "w") as f:
        f.write(f"Job Name: {run_info['job_name']}\n")
        f.write(f"Mineral: {mineral_name}\n")
        f.write(f"Temperature: {temp} K\n")
        f.write(f"Pressure: {pressure} bar\n")
        f.write(f"Functional: {functional}\n")
        f.write(f"Cutoff: {cutoff} Ry\n")
        f.write(f"Start Time: {run_info['start_time']}\n")
        f.write(f"End Time: {run_info['end_time']}\n")
        f.write(f"Runtime: {run_info['runtime']}\n")
        f.write(f"Success: {run_info['success']}\n")
        f.write(f"Command: {run_info['command']}\n")
        
        if not run_info['success']:
            f.write(f"Error: {run_info.get('error', 'Unknown error')}\n")

def find_optimized_structure(cell_opt_dir: Path, mineral_name: str, 
                           functional: str, fallback_xyz: Path) -> Path:
    """Find the optimized structure from cell optimization"""
    
    # Try different possible naming patterns
    possible_names = [
        f"{mineral_name}_CELL_OPT_{functional}-pos-1.xyz",
        "cell_opt-pos-1.xyz",
        "*-pos-1.xyz"
    ]
    
    for pattern in possible_names:
        if "*" in pattern:
            files = list(cell_opt_dir.glob(pattern))
            if files:
                return files[0]
        else:
            file_path = cell_opt_dir / pattern
            if file_path.exists():
                return file_path
    
    print("Warning: Optimized structure not found. Using fallback structure for MD.")
    return fallback_xyz

# ================================================================
# Summary and Reporting
# ================================================================

def create_simulation_summary(output_dir: Path, mineral_name: str, 
                            simulation_params: dict, job_dirs: dict, 
                            run_infos: list) -> None:
    """Create comprehensive simulation summary"""
    
    summary_file = output_dir / "simulation_summary.txt"
    
    with open(summary_file, "w") as f:
        f.write("CP2K MINERAL SIMULATION SUMMARY\n")
        f.write("="*50 + "\n\n")
        
        # Simulation parameters
        f.write("SIMULATION PARAMETERS:\n")
        f.write("-"*25 + "\n")
        for key, value in simulation_params.items():
            f.write(f"{key}: {value}\n")
        f.write(f"Simulation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Job directories
        f.write("CALCULATION DIRECTORIES:\n")
        f.write("-"*25 + "\n")
        for job_type, job_dir in job_dirs.items():
            f.write(f"{job_type}: {job_dir.name}\n")
        f.write("\n")
        
        # Runtime information
        f.write("RUNTIME INFORMATION:\n")
        f.write("-"*25 + "\n")
        total_runtime = sum([info['runtime'] for info in run_infos], datetime.timedelta())
        f.write(f"Total Runtime: {total_runtime}\n")
        
        for info in run_infos:
            status = "✓ SUCCESS" if info['success'] else "✗ FAILED"
            f.write(f"{info['job_name']}: {info['runtime']} - {status}\n")
        f.write("\n")
        
        # Output files description
        f.write("OUTPUT FILES:\n")
        f.write("-"*25 + "\n")
        f.write("Each calculation directory contains:\n")
        f.write("  - input.inp (CP2K input file)\n")
        f.write("  - output.out (CP2K output file)\n")
        f.write("  - coords.xyz (coordinate file)\n")
        f.write("  - *-pos-*.xyz (trajectory files)\n")
        f.write("  - *-forces-*.xyz (force files)\n")
        f.write("  - *-stress-*.dat (stress tensor files)\n")
        f.write("  - *-cell-*.dat (cell parameter files)\n")
        f.write("  - run_info.txt (detailed run information)\n")

def print_completion_message(output_dir: Path, run_infos: list) -> None:
    """Print final completion message with status"""
    
    successful_runs = sum(1 for info in run_infos if info['success'])
    total_runs = len(run_infos)
    
    if successful_runs == total_runs:
        print("\n" + "="*60)
        print("🎉 ALL CALCULATIONS COMPLETED SUCCESSFULLY!")
        print("="*60)
    else:
        print("\n" + "="*60)
        print(f"⚠️  SIMULATION COMPLETED WITH {total_runs - successful_runs} FAILURES")
        print(f"Successful: {successful_runs}/{total_runs}")
        print("="*60)
    
    print(f"Results saved to: {output_dir}")
    print(f"Summary file: {output_dir}/simulation_summary.txt")