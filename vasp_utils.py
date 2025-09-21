#!/usr/bin/env python3
"""
VASP Simulation Utilities Module

This module contains reusable functions for VASP simulations including:
- File handling and POSCAR processing
- INCAR generation
- POTCAR creation
- KPOINTS generation
- Job execution utilities
"""

import subprocess
import shutil
import os
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Union
import math

# ================================================================
# System and File Utilities
# ================================================================

def detect_vasp() -> str:
    """Detect VASP executable in PATH"""
    vasp_executables = [
        "vasp_std", "vasp_gam", "vasp_ncl",  # VASP 6.x
        "vasp", "vasp-std", "vasp-gamma",    # Common names
        "mpirun -np 1 vasp_std"              # MPI version fallback
    ]
    
    for exe_name in vasp_executables:
        if ' ' in exe_name:  # Handle MPI commands
            cmd_parts = exe_name.split()
            if shutil.which(cmd_parts[0]):  # Check if mpirun exists
                return exe_name
        else:
            path = shutil.which(exe_name)
            if path:
                return path
    
    raise FileNotFoundError("No VASP executable found in PATH. "
                          "Please ensure VASP is installed and accessible.")

def detect_vasp_potentials() -> str:
    """Detect VASP potential directory"""
    potential_paths = [
        os.environ.get('VASP_PSP_DIR'),
        os.environ.get('VASP_PP_PATH'),
        '/opt/vasp/potentials',
        '/usr/local/vasp/potentials',
        './potentials',
        '../potentials'
    ]
    
    for path in potential_paths:
        if path and Path(path).exists():
            return str(Path(path))
    
    raise FileNotFoundError("VASP potential directory not found. "
                          "Please set VASP_PSP_DIR environment variable.")

def find_poscar_file(mineral_name: str) -> Path:
    """Find POSCAR file based on mineral name or path"""
    # If it's already a path to an existing file
    if Path(mineral_name).exists():
        return Path(mineral_name)
    
    # Try common naming patterns
    base_name = Path(mineral_name).stem
    possible_files = [
        f"{mineral_name}",
        f"{base_name}",
        f"POSCAR_{mineral_name}",
        f"POSCAR_{base_name}",
        f"{mineral_name.upper()}",
        f"{mineral_name.lower()}",
        "POSCAR",
        "CONTCAR"
    ]
    
    for filename in possible_files:
        if Path(filename).exists():
            return Path(filename)
    
    raise FileNotFoundError(f"Could not find POSCAR file for mineral '{mineral_name}'. "
                          f"Tried: {', '.join(possible_files)}")

def parse_poscar(poscar_file: Path) -> Dict:
    """Parse POSCAR file and extract structure information"""
    lines = poscar_file.read_text().strip().split('\n')
    
    if len(lines) < 8:
        raise ValueError(f"POSCAR file {poscar_file} appears to be incomplete")
    
    structure = {
        'comment': lines[0],
        'scaling_factor': float(lines[1]),
        'lattice_vectors': [],
        'element_types': [],
        'element_counts': [],
        'coordinate_type': '',
        'positions': []
    }
    
    # Parse lattice vectors
    for i in range(2, 5):
        vector = [float(x) for x in lines[i].split()]
        structure['lattice_vectors'].append(vector)
    
    # Parse element information
    element_line = lines[5].split()
    count_line = [int(x) for x in lines[6].split()]
    
    structure['element_types'] = element_line
    structure['element_counts'] = count_line
    
    # Parse coordinate type
    structure['coordinate_type'] = lines[7].strip()[0].upper()  # Direct or Cartesian
    
    # Parse atomic positions
    total_atoms = sum(count_line)
    for i in range(8, 8 + total_atoms):
        if i < len(lines):
            pos = [float(x) for x in lines[i].split()[:3]]
            structure['positions'].append(pos)
    
    return structure

def write_poscar(structure: Dict, poscar_file: Path) -> None:
    """Write structure to POSCAR file"""
    lines = []
    
    # Comment line
    lines.append(structure['comment'])
    
    # Scaling factor
    lines.append(f"  {structure['scaling_factor']:.6f}")
    
    # Lattice vectors
    for vector in structure['lattice_vectors']:
        line = "  " + "  ".join(f"{x:15.9f}" for x in vector)
        lines.append(line)
    
    # Element types
    lines.append("  " + "  ".join(structure['element_types']))
    
    # Element counts
    lines.append("  " + "  ".join(str(x) for x in structure['element_counts']))
    
    # Coordinate type
    coord_type = structure['coordinate_type']
    if coord_type.upper() in ['D', 'DIRECT']:
        lines.append("Direct")
    else:
        lines.append("Cartesian")
    
    # Atomic positions
    for pos in structure['positions']:
        line = "  " + "  ".join(f"{x:15.9f}" for x in pos)
        lines.append(line)
    
    poscar_file.write_text('\n'.join(lines) + '\n')

# ================================================================
# INCAR Generation
# ================================================================

def create_incar_content(parameters: Dict, template: str = None) -> str:
    """Generate INCAR file content from parameters dictionary"""
    
    # Start with template if provided
    if template:
        from vasp_config import get_simulation_template
        template_params = get_simulation_template(template)
        # Merge template with user parameters (user parameters override)
        merged_params = template_params.copy()
        merged_params.update(parameters)
        parameters = merged_params
    
    lines = []
    lines.append("# VASP INCAR file generated automatically")
    lines.append(f"# Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    
    # Group parameters by category for better organization
    categories = {
        'System': ['SYSTEM'],
        'Electronic Structure': [
            'ENCUT', 'PREC', 'ALGO', 'EDIFF', 'NELM', 'NELMIN',
            'ISMEAR', 'SIGMA', 'LREAL', 'IALGO'
        ],
        'Exchange-Correlation': [
            'GGA', 'METAGGA', 'LHFCALC', 'HFSCREEN', 'AEXX',
            'IVDW', 'LVDW', 'LUSE_VDW'
        ],
        'Ionic Relaxation': [
            'IBRION', 'NSW', 'ISIF', 'EDIFFG', 'POTIM',
            'ISYM', 'SYMPREC'
        ],
        'Molecular Dynamics': [
            'MDALGO', 'TEBEG', 'TEEND', 'SMASS', 'NBLOCK',
            'KBLOCK', 'NPACO', 'APACO', 'PSTRESS'
        ],
        'Magnetic Properties': [
            'ISPIN', 'MAGMOM', 'NUPDOWN', 'LNONCOLLINEAR',
            'LSORBIT', 'SAXIS'
        ],
        'Output Control': [
            'LWAVE', 'LCHARG', 'LVTOT', 'LORBIT', 'NEDOS',
            'NWRITE', 'LELF'
        ],
        'Parallelization': [
            'NPAR', 'NCORE', 'KPAR', 'NSIM'
        ]
    }
    
    # Write parameters by category
    for category, param_list in categories.items():
        category_params = {k: v for k, v in parameters.items() if k in param_list}
        
        if category_params:
            lines.append(f"# {category}")
            lines.append("-" * (len(category) + 2))
            
            for key, value in category_params.items():
                if isinstance(value, bool):
                    lines.append(f"{key:<12} = .{str(value).upper()}.")
                elif isinstance(value, (int, float)):
                    lines.append(f"{key:<12} = {value}")
                elif isinstance(value, list):
                    value_str = " ".join(str(v) for v in value)
                    lines.append(f"{key:<12} = {value_str}")
                else:
                    lines.append(f"{key:<12} = {value}")
            lines.append("")
    
    # Add any remaining parameters not in categories
    remaining_params = {k: v for k, v in parameters.items() 
                       if not any(k in param_list for param_list in categories.values())}
    
    if remaining_params:
        lines.append("# Other Parameters")
        lines.append("------------------")
        for key, value in remaining_params.items():
            if isinstance(value, bool):
                lines.append(f"{key:<12} = .{str(value).upper()}.")
            elif isinstance(value, (int, float)):
                lines.append(f"{key:<12} = {value}")
            elif isinstance(value, list):
                value_str = " ".join(str(v) for v in value)
                lines.append(f"{key:<12} = {value_str}")
            else:
                lines.append(f"{key:<12} = {value}")
    
    return '\n'.join(lines)

def write_incar(parameters: Dict, incar_file: Path, template: str = None) -> None:
    """Write INCAR file"""
    content = create_incar_content(parameters, template)
    incar_file.write_text(content)

# ================================================================
# KPOINTS Generation
# ================================================================

def generate_kpoints_grid(lattice_vectors: List[List[float]], 
                         kpts_density: float = 1000) -> List[int]:
    """Generate k-point grid based on lattice vectors and density"""
    
    # Calculate reciprocal lattice vector lengths
    reciprocal_lengths = []
    
    for i in range(3):
        a = lattice_vectors[i]
        length = math.sqrt(sum(x*x for x in a))
        reciprocal_lengths.append(2 * math.pi / length)
    
    # Calculate k-point grid
    kpts = []
    for length in reciprocal_lengths:
        k = max(1, int(kpts_density / length))
        kpts.append(k)
    
    return kpts

def create_kpoints_content(kpts_grid: List[int], 
                          kpts_shift: List[float] = None,
                          kpts_type: str = "Gamma") -> str:
    """Generate KPOINTS file content"""
    
    if kpts_shift is None:
        kpts_shift = [0.0, 0.0, 0.0]
    
    lines = []
    lines.append("Automatic mesh")
    lines.append("0")
    lines.append(kpts_type)
    lines.append(f"  {kpts_grid[0]}  {kpts_grid[1]}  {kpts_grid[2]}")
    lines.append(f"  {kpts_shift[0]:.3f}  {kpts_shift[1]:.3f}  {kpts_shift[2]:.3f}")
    
    return '\n'.join(lines)

def write_kpoints(kpts_grid: List[int], kpoints_file: Path,
                 kpts_shift: List[float] = None,
                 kpts_type: str = "Gamma") -> None:
    """Write KPOINTS file"""
    content = create_kpoints_content(kpts_grid, kpts_shift, kpts_type)
    kpoints_file.write_text(content)

def create_kpoints_from_poscar(poscar_file: Path, kpts_density: float = 1000) -> List[int]:
    """Generate k-points grid from POSCAR file"""
    structure = parse_poscar(poscar_file)
    lattice_vectors = structure['lattice_vectors']
    scaling = structure['scaling_factor']
    
    # Scale lattice vectors
    scaled_vectors = [[v * scaling for v in vec] for vec in lattice_vectors]
    
    return generate_kpoints_grid(scaled_vectors, kpts_density)

# ================================================================
# POTCAR Generation
# ================================================================

def create_potcar(elements: List[str], potcar_file: Path, 
                 potcar_dir: str = None) -> None:
    """Create POTCAR file by concatenating individual element POTCARs"""
    
    if potcar_dir is None:
        potcar_dir = detect_vasp_potentials()
    
    potcar_dir = Path(potcar_dir)
    
    # Import POTCAR info
    from vasp_config import POTCAR_INFO
    
    potcar_content = []
    
    for element in elements:
        # Determine POTCAR version to use
        if element in POTCAR_INFO:
            potcar_version = POTCAR_INFO[element]['version']
        else:
            potcar_version = element
            print(f"Warning: Using default POTCAR for {element}")
        
        # Look for POTCAR file
        element_potcar_paths = [
            potcar_dir / potcar_version / "POTCAR",
            potcar_dir / element / "POTCAR",
            potcar_dir / f"{element}_pv" / "POTCAR",
            potcar_dir / f"{element}_sv" / "POTCAR",
        ]
        
        potcar_found = False
        for potcar_path in element_potcar_paths:
            if potcar_path.exists():
                potcar_content.append(potcar_path.read_text())
                print(f"Added POTCAR for {element} from {potcar_path}")
                potcar_found = True
                break
        
        if not potcar_found:
            raise FileNotFoundError(f"POTCAR for element {element} not found. "
                                  f"Searched in: {[str(p) for p in element_potcar_paths]}")
    
    # Write combined POTCAR
    potcar_file.write_text(''.join(potcar_content))

def validate_potcar(potcar_file: Path, expected_elements: List[str]) -> bool:
    """Validate POTCAR file contains expected elements"""
    
    if not potcar_file.exists():
        return False
    
    content = potcar_file.read_text()
    
    # Count POTCAR sections
    potcar_count = content.count('VRHFIN')
    expected_count = len(expected_elements)
    
    if potcar_count != expected_count:
        print(f"Warning: POTCAR contains {potcar_count} elements, "
              f"expected {expected_count}")
        return False
    
    return True

# ================================================================
# Job Execution
# ================================================================

def run_vasp_calculation(job_name: str, job_dir: Path, vasp_cmd: str,
                        mpi_ranks: int = None) -> Dict:
    """Execute VASP calculation and return run information"""
    
    print(f"Running {job_name} calculation in {job_dir}")
    start_time = datetime.now()
    
    # Prepare command
    if mpi_ranks and mpi_ranks > 1:
        if 'mpirun' not in vasp_cmd:
            cmd = f"mpirun -np {mpi_ranks} {vasp_cmd}"
        else:
            cmd = vasp_cmd.replace('-np 1', f'-np {mpi_ranks}')
    else:
        cmd = vasp_cmd
    
    try:
        # Run VASP
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=job_dir,
            check=True,
            capture_output=True,
            text=True
        )
        
        end_time = datetime.now()
        runtime = end_time - start_time
        
        # Check if calculation converged by examining OSZICAR
        converged = check_convergence(job_dir)
        
        if converged:
            print(f"✓ {job_name} completed successfully in {runtime}")
        else:
            print(f"⚠️ {job_name} completed but may not be converged in {runtime}")
        
        run_info = {
            'job_name': job_name,
            'start_time': start_time,
            'end_time': end_time,
            'runtime': runtime,
            'success': True,
            'converged': converged,
            'command': cmd
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
            'converged': False,
            'error': str(e),
            'stderr': e.stderr
        }
        
        return run_info

def check_convergence(job_dir: Path) -> bool:
    """Check if VASP calculation converged by examining OSZICAR"""
    oszicar_file = job_dir / "OSZICAR"
    
    if not oszicar_file.exists():
        return False
    
    try:
        lines = oszicar_file.read_text().strip().split('\n')
        
        # Check for electronic convergence in last few lines
        for line in reversed(lines[-10:]):
            if 'F=' in line and 'E0=' in line:
                return True
        
        return False
    
    except Exception:
        return False

def extract_final_energy(job_dir: Path) -> Optional[float]:
    """Extract final energy from OSZICAR"""
    oszicar_file = job_dir / "OSZICAR"
    
    if not oszicar_file.exists():
        return None
    
    try:
        lines = oszicar_file.read_text().strip().split('\n')
        
        # Look for final energy line
        for line in reversed(lines):
            if 'F=' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == 'F=' and i + 1 < len(parts):
                        return float(parts[i + 1])
        
        return None
    
    except Exception:
        return None

# ================================================================
# File Management
# ================================================================

def backup_previous_calculation(job_dir: Path) -> None:
    """Backup previous calculation files if they exist"""
    
    backup_dir = job_dir / "previous_run"
    
    files_to_backup = [
        "OUTCAR", "OSZICAR", "vasprun.xml", "CONTCAR",
        "CHG", "CHGCAR", "WAVECAR", "EIGENVAL", "DOSCAR"
    ]
    
    files_found = []
    for filename in files_to_backup:
        filepath = job_dir / filename
        if filepath.exists():
            files_found.append(filepath)
    
    if files_found:
        backup_dir.mkdir(exist_ok=True)
        print(f"Backing up {len(files_found)} files to {backup_dir}")
        
        for filepath in files_found:
            backup_path = backup_dir / filepath.name
            shutil.move(str(filepath), str(backup_path))

def copy_structure_files(source_dir: Path, target_dir: Path) -> None:
    """Copy structure files (POSCAR, CONTCAR) between directories"""
    
    # Try CONTCAR first (final structure), then POSCAR
    for filename in ["CONTCAR", "POSCAR"]:
        source_file = source_dir / filename
        if source_file.exists():
            target_file = target_dir / "POSCAR"
            shutil.copy2(source_file, target_file)
            print(f"Copied {filename} from {source_dir} to {target_dir}/POSCAR")
            return
    
    print(f"Warning: No structure file found in {source_dir}")

def write_run_info(job_dir: Path, run_info: Dict, mineral_name: str,
                  simulation_params: Dict) -> None:
    """Write detailed run information to file"""
    
    info_file = job_dir / "run_info.txt"
    
    with open(info_file, "w") as f:
        f.write(f"VASP Calculation Information\n")
        f.write(f"{'='*40}\n\n")
        
        f.write(f"Job Name: {run_info['job_name']}\n")
        f.write(f"Mineral: {mineral_name}\n")
        f.write(f"Start Time: {run_info['start_time']}\n")
        f.write(f"End Time: {run_info['end_time']}\n")
        f.write(f"Runtime: {run_info['runtime']}\n")
        f.write(f"Success: {run_info['success']}\n")
        f.write(f"Converged: {run_info.get('converged', 'Unknown')}\n")
        f.write(f"Command: {run_info['command']}\n\n")
        
        # Extract final energy if available
        final_energy = extract_final_energy(job_dir)
        if final_energy is not None:
            f.write(f"Final Energy: {final_energy:.6f} eV\n\n")
        
        f.write("Simulation Parameters:\n")
        f.write("-" * 20 + "\n")
        for key, value in simulation_params.items():
            f.write(f"{key}: {value}\n")
        
        if not run_info['success']:
            f.write(f"\nError: {run_info.get('error', 'Unknown error')}\n")

# ================================================================
# Summary and Reporting
# ================================================================

def create_simulation_summary(output_dir: Path, mineral_name: str,
                            simulation_params: Dict, job_dirs: Dict,
                            run_infos: List[Dict]) -> None:
    """Create comprehensive simulation summary"""
    
    summary_file = output_dir / "simulation_summary.txt"
    
    with open(summary_file, "w") as f:
        f.write("VASP MINERAL SIMULATION SUMMARY\n")
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
            converged = "CONVERGED" if info.get('converged') else "NOT CONVERGED"
            f.write(f"{info['job_name']}: {info['runtime']} - {status} - {converged}\n")
        f.write("\n")
        
        # Final energies
        f.write("FINAL ENERGIES:\n")
        f.write("-"*25 + "\n")
        for job_type, job_dir in job_dirs.items():
            final_energy = extract_final_energy(job_dir)
            if final_energy is not None:
                f.write(f"{job_type}: {final_energy:.6f} eV\n")
        f.write("\n")
        
        # Output files description
        f.write("OUTPUT FILES:\n")
        f.write("-"*25 + "\n")
        f.write("Each calculation directory contains:\n")
        f.write("  - INCAR (input parameters)\n")
        f.write("  - POSCAR (structure file)\n")
        f.write("  - POTCAR (pseudopotential file)\n")
        f.write("  - KPOINTS (k-point sampling)\n")
        f.write("  - OUTCAR (detailed output)\n")
        f.write("  - OSZICAR (convergence information)\n")
        f.write("  - vasprun.xml (XML output)\n")
        f.write("  - CONTCAR (final structure)\n")
        f.write("  - run_info.txt (job information)\n")

def print_completion_message(output_dir: Path, run_infos: List[Dict]) -> None:
    """Print final completion message with status"""
    
    successful_runs = sum(1 for info in run_infos if info['success'])
    converged_runs = sum(1 for info in run_infos if info.get('converged', False))
    total_runs = len(run_infos)
    
    if successful_runs == total_runs and converged_runs == total_runs:
        print("\n" + "="*60)
        print("🎉 ALL CALCULATIONS COMPLETED AND CONVERGED!")
        print("="*60)
    elif successful_runs == total_runs:
        print("\n" + "="*60)
        print("⚠️  ALL CALCULATIONS COMPLETED BUT SOME MAY NOT BE CONVERGED")
        print(f"Converged: {converged_runs}/{total_runs}")
        print("="*60)
    else:
        print("\n" + "="*60)
        print(f"⚠️  SIMULATION COMPLETED WITH {total_runs - successful_runs} FAILURES")
        print(f"Successful: {successful_runs}/{total_runs}")
        print(f"Converged: {converged_runs}/{total_runs}")
        print("="*60)
    
    print(f"Results saved to: {output_dir}")
    print(f"Summary file: {output_dir}/simulation_summary.txt")