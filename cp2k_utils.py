#!/usr/bin/env python3
"""
CP2K Simulation Utilities Module

This module contains reusable functions for CP2K simulations including:
- File handling and structure processing
- Input file generation  
- Job execution utilities
- Output analysis
"""

import subprocess
import shutil
import os
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Dict, Tuple, Optional, Union
import math
import re
from cp2k_config import BASIS_SET_INFO, PSEUDOPOTENTIAL_INFO


# ================================================================
# System and File Utilities
# ================================================================

def detect_cp2k() -> str:
    """Detect CP2K executable in PATH"""
    cp2k_executables = [
        "cp2k.psmp", "cp2k.popt", "cp2k.ssmp", "cp2k.sopt",  # Standard CP2K
        "cp2k", "cp2k.x",                                      # Generic names
        "mpirun -np 1 cp2k.psmp"                              # MPI version fallback
    ]
    
    for exe_name in cp2k_executables:
        if ' ' in exe_name:  # Handle MPI commands
            cmd_parts = exe_name.split()
            if shutil.which(cmd_parts[0]):  # Check if mpirun exists
                return exe_name
        else:
            path = shutil.which(exe_name)
            if path:
                return path
    
    raise FileNotFoundError("No CP2K executable found in PATH. "
                          "Please ensure CP2K is installed and accessible.")

def detect_cp2k_data() -> str:
    """Detect CP2K data directory for basis sets and potentials"""
    data_paths = [
        os.environ.get('CP2K_DATA_DIR'),
        '/opt/cp2k/data',
        '/usr/local/share/cp2k/data',
        '/usr/share/cp2k/data',
        './data',
        '../data'
    ]
    
    for path in data_paths:
        if path and Path(path).exists():
            # Check for basis sets and potentials subdirectories
            if (Path(path) / 'BASIS_MOLOPT').exists() or (Path(path) / 'POTENTIAL').exists():
                return str(Path(path))
    
    raise FileNotFoundError("CP2K data directory not found. "
                          "Please set CP2K_DATA_DIR environment variable.")

def find_structure_file(mineral_name: str) -> Path:
    """Find structure file based on mineral name or path"""
    # If it's already a path to an existing file
    if Path(mineral_name).exists():
        return Path(mineral_name)
    
    # Try common naming patterns and file extensions
    base_name = Path(mineral_name).stem
    extensions = ['.xyz', '.pdb', '.cif', '.inp', '']
    
    possible_files = []
    for ext in extensions:
        possible_files.extend([
            f"{mineral_name}{ext}",
            f"{base_name}{ext}",
            f"{mineral_name.upper()}{ext}",
            f"{mineral_name.lower()}{ext}",
        ])
    
    # Add common structure file names
    possible_files.extend(['POSCAR', 'CONTCAR', 'structure.xyz', 'geometry.xyz'])
    
    for filename in possible_files:
        if Path(filename).exists():
            return Path(filename)
    
    raise FileNotFoundError(f"Could not find structure file for mineral '{mineral_name}'. "
                          f"Tried: {', '.join(possible_files[:8])}")

def parse_xyz_file(xyz_file: Path) -> Dict:
    """Parse XYZ file and extract structure information"""
    lines = xyz_file.read_text().strip().split('\n')
    
    if len(lines) < 3:
        raise ValueError(f"XYZ file {xyz_file} appears to be incomplete")
    
    # Parse number of atoms
    try:
        n_atoms = int(lines[0])
    except ValueError:
        raise ValueError(f"Invalid number of atoms in first line: {lines[0]}")
    
    # Parse comment line (may contain cell information)
    comment = lines[1] if len(lines) > 1 else ""
    
    # Parse atomic coordinates
    atoms = []
    elements = []
    positions = []
    
    for i in range(2, min(2 + n_atoms, len(lines))):
        parts = lines[i].split()
        if len(parts) < 4:
            continue
        
        element = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
            
        elements.append(element)
        positions.append([x, y, z])
        atoms.append({'element': element, 'position': [x, y, z]})
    
    # Extract unique elements and their counts
    unique_elements = []
    element_counts = []
    
    for element in elements:
        if element not in unique_elements:
            unique_elements.append(element)
            element_counts.append(elements.count(element))
    
    structure = {
        'comment': comment,
        'n_atoms': len(atoms),
        'atoms': atoms,
        'elements': elements,
        'unique_elements': unique_elements,
        'element_counts': element_counts,
        'positions': positions,
        'cell_vectors': None  # Will be set if found in comment or separate
    }
    
    # Try to extract cell vectors from comment if present
    cell_info = extract_cell_from_comment(comment)
    if cell_info:
        structure['cell_vectors'] = cell_info
    
    return structure

def parse_poscar_to_xyz(poscar_file: Path) -> Dict:
    """Convert POSCAR format to XYZ-like structure information"""
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
    structure['coordinate_type'] = lines[7].strip()[0].upper()
    
    # Parse atomic positions
    total_atoms = sum(count_line)
    for i in range(8, 8 + total_atoms):
        if i < len(lines):
            pos = [float(x) for x in lines[i].split()[:3]]
            structure['positions'].append(pos)
    
    # Convert to XYZ-like format
    atoms = []
    elements = []
    positions = []
    
    atom_idx = 0
    for elem_idx, (element, count) in enumerate(zip(element_line, count_line)):
        for _ in range(count):
            pos = structure['positions'][atom_idx]
            atoms.append({'element': element, 'position': pos})
            elements.append(element)
            positions.append(pos)
            atom_idx += 1
    
    xyz_structure = {
        'comment': structure['comment'],
        'n_atoms': total_atoms,
        'atoms': atoms,
        'elements': elements,
        'unique_elements': element_line,
        'element_counts': count_line,
        'positions': positions,
        'cell_vectors': structure['lattice_vectors'],
        'scaling_factor': structure['scaling_factor'],
        'coordinate_type': structure['coordinate_type']
    }
    
    return xyz_structure

def extract_cell_from_comment(comment: str) -> Optional[List[List[float]]]:
    """Extract cell vectors from XYZ comment line"""
    # Look for patterns like: cell="a b c alpha beta gamma" or Lattice="xx yy zz ..."
    patterns = [
        r'cell="([^"]+)"',
        r'Lattice="([^"]+)"',
        r'CELL\s*=\s*"([^"]+)"'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, comment, re.IGNORECASE)
        if match:
            cell_str = match.group(1)
            try:
                values = [float(x) for x in cell_str.split()]
                if len(values) == 6:
                    # Convert from a,b,c,α,β,γ to vectors
                    return cell_params_to_vectors(values)
                elif len(values) == 9:
                    # Direct lattice vectors
                    return [values[0:3], values[3:6], values[6:9]]
            except ValueError:
                continue
    
    return None

def cell_params_to_vectors(params: List[float]) -> List[List[float]]:
    """Convert cell parameters (a,b,c,α,β,γ) to lattice vectors"""
    a, b, c, alpha, beta, gamma = params
    
    # Convert angles to radians
    alpha_rad = math.radians(alpha)
    beta_rad = math.radians(beta)  
    gamma_rad = math.radians(gamma)
    
    # Calculate lattice vectors
    cos_alpha = math.cos(alpha_rad)
    cos_beta = math.cos(beta_rad)
    cos_gamma = math.cos(gamma_rad)
    sin_gamma = math.sin(gamma_rad)
    
    # Volume calculation
    volume = a * b * c * math.sqrt(1 + 2*cos_alpha*cos_beta*cos_gamma 
                                   - cos_alpha**2 - cos_beta**2 - cos_gamma**2)
    
    # Lattice vectors
    vector_a = [a, 0.0, 0.0]
    vector_b = [b * cos_gamma, b * sin_gamma, 0.0]
    vector_c = [c * cos_beta,
                c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
                volume / (a * b * sin_gamma)]
    
    return [vector_a, vector_b, vector_c]

def write_xyz_file(structure: Dict, xyz_file: Path) -> None:
    """Write structure to XYZ file"""
    lines = []
    
    # Number of atoms
    lines.append(str(structure['n_atoms']))
    
    # Comment line
    comment = structure.get('comment', 'Generated by CP2K utilities')
    
    # Add cell information to comment if available
    if structure.get('cell_vectors'):
        cell_vectors = structure['cell_vectors']
        cell_str = ' '.join([f"{v:.6f}" for vec in cell_vectors for v in vec])
        comment += f' Lattice="{cell_str}"'
    
    lines.append(comment)
    
    # Atomic coordinates
    for atom in structure['atoms']:
        element = atom['element']
        x, y, z = atom['position']
        lines.append(f"{element:>2s}  {x:15.9f}  {y:15.9f}  {z:15.9f}")
    
    xyz_file.write_text('\n'.join(lines) + '\n')

# ================================================================
# CP2K Input File Generation
# ================================================================

def create_cp2k_input(parameters: Dict, structure: Dict, input_file: Path) -> None:
    """Generate CP2K input file"""
    
    content = []
    content.append("&GLOBAL")
    content.append(f"  PROJECT_NAME {parameters.get('project_name', 'cp2k_sim')}")
    content.append(f"  RUN_TYPE {parameters.get('run_type', 'GEO_OPT')}")
    content.append(f"  PRINT_LEVEL {parameters.get('print_level', 'MEDIUM')}")
    content.append("&END GLOBAL")
    content.append("")
    
    # Force evaluation section
    content.extend(create_force_eval_section(parameters, structure))
    
    # Motion section (for geometry optimization, MD, etc.)
    if parameters.get('run_type') in ['GEO_OPT', 'CELL_OPT', 'MD']:
        content.extend(create_motion_section(parameters))
    
    input_file.write_text('\n'.join(content))

def create_force_eval_section(parameters: Dict, structure: Dict) -> List[str]:
    """Create FORCE_EVAL section"""
    lines = []
    lines.append("&FORCE_EVAL")
    lines.append(f"  METHOD {parameters.get('method', 'QS')}")
    lines.append("")
    
    # DFT section
    if parameters.get('method', 'QS') == 'QS':
        lines.extend(create_dft_section(parameters))
        lines.append("")
    
    # Subsystem section
    lines.extend(create_subsys_section(parameters, structure))
    
    lines.append("&END FORCE_EVAL")
    return lines

def create_dft_section(parameters: Dict) -> List[str]:
    """Create DFT section"""
    lines = []
    lines.append("  &DFT")
    
    # Basis set file
    lines.append("    BASIS_SET_FILE_NAME BASIS_MOLOPT")
    lines.append("    POTENTIAL_FILE_NAME GTH_POTENTIALS")
    lines.append("")
    
    # QS section (plane waves)
    lines.append("    &QS")
    lines.append(f"      EPS_DEFAULT {parameters.get('scf_eps', 1E-6):.0E}")
    lines.append("    &END QS")
    lines.append("")
    
    # MGRID section
    lines.append("    &MGRID")
    lines.append(f"      CUTOFF {parameters.get('cutoff', 400)}")
    lines.append(f"      REL_CUTOFF {parameters.get('rel_cutoff', 50)}")
    lines.append("    &END MGRID")
    lines.append("")
    
    # Exchange-correlation
    lines.extend(create_xc_section(parameters))
    lines.append("")
    
    # SCF section
    lines.extend(create_scf_section(parameters))
    
    lines.append("  &END DFT")
    return lines
def create_xc_section(parameters: dict) -> list[str]:
    """Create XC section for exchange-correlation functional (CP2K >= 2023.3)"""
    lines = []
    lines.append("    &XC")

    functional = parameters.get('xc_functional', 'PBE').upper()

    # For standard built-in functionals (PBE, LDA, etc.), write simplified syntax
    standard_functionals = {'PBE', 'LDA', 'PADE', 'BLYP', 'SCAN'}  # extend if needed
    if functional in standard_functionals:
        lines.append(f"      &XC_FUNCTIONAL {functional}")
        lines.append("      &END XC_FUNCTIONAL")
    else:
        # Use LibXC for custom functionals
        from cp2k_config import get_functional_info
        func_info = get_functional_info(functional)
        lines.append("      &XC_FUNCTIONAL")
        if 'libxc' in func_info:
            libxc_functionals = func_info['libxc'].split()
            for libxc_func in libxc_functionals:
                lines.append(f"        &{libxc_func}")
                lines.append(f"        &END {libxc_func}")
        lines.append("      &END XC_FUNCTIONAL")

    # Van der Waals section (optional)
    vdw_type = parameters.get('vdw_potential')
    if vdw_type:
        lines.extend(create_vdw_section(vdw_type))

    lines.append("    &END XC")
    return lines


# def create_xc_section(parameters: Dict) -> List[str]:
#     """Create XC section for exchange-correlation functional"""
#     lines = []
#     lines.append("    &XC")
    
#     functional = parameters.get('xc_functional', 'PBE').upper()
    
#     # Import functional information
#     from cp2k_config import get_functional_info
#     func_info = get_functional_info(functional)
    
#     if 'libxc' in func_info:
#         # Use LibXC
#         lines.append("      &XC_FUNCTIONAL")
#         libxc_functionals = func_info['libxc'].split()
#         for libxc_func in libxc_functionals:
#             lines.append(f"        &{libxc_func}")
#             lines.append(f"        &END {libxc_func}")
#         lines.append("      &END XC_FUNCTIONAL")
#     else:
#         # Use built-in functional
#         lines.append("      &XC_FUNCTIONAL")
#         if functional == 'PBE':
#             lines.append("        &PBE")
#             lines.append("        &END PBE")
#         elif functional == 'LDA':
#             lines.append("        &PADE")
#             lines.append("        &END PADE")
#         lines.append("      &END XC_FUNCTIONAL")
    
#     # Van der Waals correction
#     vdw_type = parameters.get('vdw_potential')
#     if vdw_type:
#         lines.extend(create_vdw_section(vdw_type))
    
#     lines.append("    &END XC")
#     return lines

def create_vdw_section(vdw_type: str) -> List[str]:
    """Create van der Waals correction section"""
    lines = []
    
    from cp2k_config import get_vdw_info
    vdw_info = get_vdw_info(vdw_type)
    
    if vdw_info:
        lines.append("      &VDW_POTENTIAL")
        if vdw_type.upper() in ['D2', 'D3', 'D3BJ']:
            lines.append(f"        DISPERSION_FUNCTIONAL PAIR_POTENTIAL")
            lines.append(f"        &PAIR_POTENTIAL")
            lines.append(f"          TYPE DFTD3")
            if vdw_type.upper() == 'D3BJ':
                lines.append(f"          PARAMETER_FILE_NAME dftd3.dat")
                lines.append(f"          REFERENCE_FUNCTIONAL PBE")
            lines.append(f"        &END PAIR_POTENTIAL")
        lines.append("      &END VDW_POTENTIAL")
    
    return lines

def create_scf_section(parameters: Dict) -> List[str]:
    """Create SCF section"""
    lines = []
    lines.append("    &SCF")
    lines.append(f"      MAX_SCF {parameters.get('scf_max_iter', 50)}")
    lines.append(f"      EPS_SCF {parameters.get('scf_eps', 1E-6):.0E}")
    lines.append(f"      SCF_GUESS {parameters.get('scf_guess', 'ATOMIC')}")
    
    # UKS for spin-polarized calculations
    if parameters.get('uks', False):
        lines.append("      UKS .TRUE.")
        if parameters.get('multiplicity'):
            lines.append(f"      MULTIPLICITY {parameters.get('multiplicity')}")
    
    # OT diagonalization (usually faster)
    lines.append("      &OT")
    lines.append("        MINIMIZER DIIS")
    lines.append("        PRECONDITIONER FULL_SINGLE_INVERSE")
    lines.append("      &END OT")
    
    # Mixing for difficult convergence cases
    lines.append("      &MIXING")
    lines.append("        METHOD BROYDEN_MIXING")
    lines.append("        ALPHA 0.4")
    lines.append("      &END MIXING")
    
    lines.append("    &END SCF")
    return lines

def create_subsys_section(parameters: Dict, structure: Dict) -> List[str]:
    """Create SUBSYS section with geometry and cell"""
    lines = []
    lines.append("  &SUBSYS")
    lines.append("")
    
    # Cell information
    if structure.get('cell_vectors'):
        lines.append("    &CELL")
        for i, vector in enumerate(structure['cell_vectors']):
            label = ['A', 'B', 'C'][i]
            lines.append(f"      {label}  {vector[0]:12.6f}  {vector[1]:12.6f}  {vector[2]:12.6f}")
        
        # Periodic boundary conditions
        lines.append("      PERIODIC XYZ")
        lines.append("    &END CELL")
        lines.append("")
    
    # Coordinate section
    lines.append("    &COORD")
    for atom in structure['atoms']:
        element = atom['element']
        x, y, z = atom['position']
        lines.append(f"      {element:>2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
    lines.append("    &END COORD")
    lines.append("")
    
    # Kind sections (basis sets and pseudopotentials)
    for element in structure['unique_elements']:
        lines.extend(create_kind_section(element, parameters))
        lines.append("")
    
    lines.append("  &END SUBSYS")
    return lines

from typing import Dict, List

def create_kind_section(element: str, parameters: Dict) -> List[str]:
    """Create KIND section for an element with correct basis set and GTH-PBE-qN pseudopotential."""
    lines = []
    lines.append(f"    &KIND {element}")

    # Import predefined info
    from cp2k_config import BASIS_SET_INFO, PSEUDOPOTENTIAL_INFO

    # --- BASIS SET ---
    # Check user-defined basis sets first
    basis_set = parameters.get('basis_sets', {}).get(element)
    if not basis_set:
        # Check predefined default
        basis_set = BASIS_SET_INFO.get(element, {}).get('default')
    if not basis_set:
        # Fallback
        basis_set = parameters.get('basis_set_type', 'DZVP-MOLOPT-SR-GTH')
    lines.append(f"      BASIS_SET {basis_set}")

    # --- PSEUDOPOTENTIAL ---
    # Check user-defined pseudopotentials first
    potential = parameters.get('pseudopotentials', {}).get(element)
    if not potential:
        # Use predefined default with -qN if available
        potential = PSEUDOPOTENTIAL_INFO.get(element, {}).get('default')
    if not potential:
        # Fallback: generate GTH-PBE-q{valence} if valence known
        valence = PSEUDOPOTENTIAL_INFO.get(element, {}).get('valence')
        if valence:
            potential = f"GTH-PBE-q{valence}"
        else:
            # Safe ultimate fallback
            potential = "GTH-PBE-q999"
    
    lines.append(f"      POTENTIAL {potential}")
    lines.append(f"    &END KIND")

    return lines



# def create_kind_section(element: str, parameters: Dict) -> List[str]:
#     """Create KIND section for an element"""
#     lines = []
#     lines.append(f"    &KIND {element}")
    
#     # Get basis set and pseudopotential info
#     from cp2k_config import BASIS_SET_INFO, PSEUDOPOTENTIAL_INFO
    
#     # Basis set
#     basis_set = parameters.get('basis_sets', {}).get(element)
#     if not basis_set and element in BASIS_SET_INFO:
#         basis_set = BASIS_SET_INFO[element]['default']
#     if not basis_set:
#         basis_set = parameters.get('basis_set_type', 'DZVP-MOLOPT-SR-GTH')
    
#     lines.append(f"      BASIS_SET {basis_set}")
    
#     # Pseudopotential
#     potential = parameters.get('pseudopotentials', {}).get(element)
#     if not potential and element in PSEUDOPOTENTIAL_INFO:
#         potential = PSEUDOPOTENTIAL_INFO[element]['default']
#     if not potential:
#         potential = f"GTH-{parameters.get('xc_functional', 'PBE')}"
    
#     lines.append(f"      POTENTIAL {potential}")
    
#     lines.append(f"    &END KIND")
#     return lines

# def create_motion_section(parameters: Dict) -> List[str]:
#     """Create MOTION section"""
#     lines = []
#     lines.append("")
#     lines.append("&MOTION")
    
#     run_type = parameters.get('run_type', 'GEO_OPT')
    
#     if run_type in ['GEO_OPT', 'CELL_OPT']:
#         lines.extend(create_geo_opt_section(parameters))
#     elif run_type == 'MD':
#         lines.extend(create_md_section(parameters))
    
#     # Print section for trajectory
#     lines.append("  &PRINT")
#     lines.append("    &TRAJECTORY")
#     lines.append(f"      FORMAT {parameters.get('trajectory_format', 'XYZ')}")
#     lines.append("      &EACH")
#     lines.append("        MD 1")
#     lines.append("      &END EACH")
#     lines.append("    &END TRAJECTORY")
#     lines.append("    &RESTART")
#     lines.append(f"      BACKUP_COPIES {parameters.get('restart_backup', 1)}")
#     lines.append("    &END RESTART")
#     lines.append("  &END PRINT")
    
#     lines.append("&END MOTION")
#     return lines

def create_motion_section(parameters: Dict) -> List[str]:
    """Create MOTION section"""
    lines = []
    lines.append("")
    lines.append("&MOTION")
    
    run_type = parameters.get('run_type', 'GEO_OPT')
    
    if run_type in ['GEO_OPT', 'CELL_OPT']:
        lines.extend(create_geo_opt_section(parameters))
    elif run_type == 'MD':
        lines.extend(create_md_section(parameters))
    
    # Print section
    lines.append("  &PRINT")

    # Trajectory
    lines.append("    &TRAJECTORY")
    lines.append(f"      FORMAT {parameters.get('trajectory_format', 'XYZ')}")
    lines.append(f"      FILENAME {parameters.get('traj_file', 'md-pos.xyz')}")
    lines.append("      &EACH")
    lines.append("        MD 1")
    lines.append("      &END EACH")
    lines.append("    &END TRAJECTORY")

    # Velocities
    lines.append("    &VELOCITIES")
    lines.append(f"      FORMAT {parameters.get('velocity_format', 'XYZ')}")
    lines.append(f"      FILENAME {parameters.get('vel_file', 'md-vel.xyz')}")
    lines.append("      &EACH")
    lines.append("        MD 1")
    lines.append("      &END EACH")
    lines.append("    &END VELOCITIES")

    # Forces
    lines.append("    &FORCES")
    lines.append(f"      FORMAT {parameters.get('forces_format', 'XYZ')}")
    lines.append(f"      FILENAME {parameters.get('forces_file', 'md-forces.xyz')}")
    lines.append("      &EACH")
    lines.append("        MD 1")
    lines.append("      &END EACH")
    lines.append("    &END FORCES")

    # Stress
    lines.append("    &STRESS")
    lines.append(f"      FILENAME {parameters.get('stress_file', 'md-stress.dat')}")
    lines.append("      &EACH")
    lines.append("        MD 1")
    lines.append("      &END EACH")
    lines.append("    &END STRESS")

    # Cell
    lines.append("    &CELL")
    lines.append(f"      FILENAME {parameters.get('cell_file', 'md-cell.dat')}")
    lines.append("      &EACH")
    lines.append("        MD 1")
    lines.append("      &END EACH")
    lines.append("    &END CELL")

    # Restart
    lines.append("    &RESTART")
    lines.append(f"      BACKUP_COPIES {parameters.get('restart_backup', 1)}")
    lines.append("    &END RESTART")

    lines.append("  &END PRINT")
    lines.append("&END MOTION")
    return lines



def create_geo_opt_section(parameters: Dict) -> List[str]:
    """Create geometry optimization section"""
    lines = []
    lines.append("  &GEO_OPT")
    lines.append(f"    MAX_ITER {parameters.get('geo_max_iter', 200)}")
    lines.append(f"    OPTIMIZER {parameters.get('optimizer', 'BFGS')}")
    lines.append("    &BFGS")
    lines.append("    &END BFGS")
    lines.append("  &END GEO_OPT")
    
    # Cell optimization if requested
    if parameters.get('cell_opt', False):
        lines.append("")
        lines.append("  &CELL_OPT")
        lines.append(f"    MAX_ITER {parameters.get('geo_max_iter', 200)}")
        lines.append(f"    OPTIMIZER {parameters.get('optimizer', 'BFGS')}")
        lines.append(f"    PRESSURE_TOLERANCE {parameters.get('pressure_tolerance', 100)}")
        lines.append("  &END CELL_OPT")
    
    return lines

def create_md_section(parameters: Dict) -> List[str]:
    """Create molecular dynamics section"""
    lines = []
    lines.append("  &MD")
    lines.append(f"    ENSEMBLE {parameters.get('md_ensemble', 'NVT')}")
    lines.append(f"    STEPS {parameters.get('md_steps', 1000)}")
    lines.append(f"    TIMESTEP {parameters.get('md_timestep', 0.5)}")
    lines.append(f"    TEMPERATURE {parameters.get('md_temperature', 300.0)}")
    lines.append("")
    
    # Thermostat
    thermostat = parameters.get('thermostat', 'NOSE')
    if thermostat == 'NOSE':
        lines.append("    &THERMOSTAT")
        lines.append("      TYPE NOSE")
        lines.append(f"      &NOSE")
        lines.append(f"        TIMECON {parameters.get('thermostat_tau', 100.0)}")
        lines.append("      &END NOSE")
        lines.append("    &END THERMOSTAT")
    elif thermostat == 'CSVR':
        lines.append("    &THERMOSTAT")
        lines.append("      TYPE CSVR")
        lines.append(f"      &CSVR")
        lines.append(f"        TIMECON {parameters.get('thermostat_tau', 100.0)}")
        lines.append("      &END CSVR")
        lines.append("    &END THERMOSTAT")
    
    # Barostat for NPT
    if parameters.get('md_ensemble') == 'NPT_I':
        lines.append("    &BAROSTAT")
        lines.append("      TYPE NOSE")
        lines.append(f"      PRESSURE {parameters.get('md_pressure', 1.0)}")
        lines.append("      &NOSE")
        lines.append(f"        TIMECON {parameters.get('barostat_tau', 1000.0)}")
        lines.append("      &END NOSE")
        lines.append("    &END BAROSTAT")
    
    lines.append("  &END MD")
    return lines

# ================================================================
# Job Execution
# ================================================================

def run_cp2k_calculation(job_name: str, job_dir: Path, cp2k_cmd: str,
                        input_file: str = "input.inp", mpi_ranks: int = None) -> Dict:
    """Execute CP2K calculation and return run information"""
    
    print(f"Running {job_name} calculation in {job_dir}")
    start_time = datetime.now()
    
    # Prepare command - always suppress CUDA warnings
    mca_flags = "--mca opal_warn_on_missing_libcuda 0 --mca btl ^openib"
    
    if mpi_ranks and mpi_ranks > 1:
        if 'mpirun' not in cp2k_cmd:
            cmd = f"mpirun {mca_flags} -np {mpi_ranks} {cp2k_cmd} -i {input_file}"
        else:
            # Extract the actual CP2K executable from the command
            cp2k_exe = cp2k_cmd.split()[-1] if 'cp2k' in cp2k_cmd else cp2k_cmd
            cmd = f"mpirun {mca_flags} -np {mpi_ranks} {cp2k_exe} -i {input_file}"
    else:
        # For single process, still use mpirun to control the environment
        if 'mpirun' in cp2k_cmd:
            # Extract the actual CP2K executable from the command
            cp2k_exe = cp2k_cmd.split()[-1] if 'cp2k' in cp2k_cmd else cp2k_cmd.replace('mpirun -np 1 ', '')
            cmd = f"mpirun {mca_flags} -np 1 {cp2k_exe} -i {input_file}"
        else:
            # Even for single process, use mpirun to suppress warnings
            cmd = f"mpirun {mca_flags} -np 1 {cp2k_cmd} -i {input_file}"
    
    try:
        # Run CP2K
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
        
        # Check if calculation converged
        converged = check_cp2k_convergence(job_dir)
        
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
            'stderr': e.stderr,
            'command': cmd
        }
        
        return run_info

def check_cp2k_convergence(job_dir: Path) -> bool:
    """Check if CP2K calculation converged by examining output"""
    output_files = [
        job_dir / "cp2k_sim.out",  # Default project name
        job_dir / "output.out",
        job_dir / "cp2k.out"
    ]
    
    # Find output file
    output_file = None
    for f in output_files:
        if f.exists():
            output_file = f
            break
    
    # Also check for project-specific output
    for f in job_dir.glob("*.out"):
        if f.exists():
            output_file = f
            break
    
    if not output_file:
        return False
    
    try:
        content = output_file.read_text()
        
        # Check for successful completion
        success_patterns = [
            "PROGRAM ENDED AT",
            "CP2K                                 1.0",
            "GLOBAL| CP2K",
            "*** SCF run converged in"
        ]
        
        # Check for convergence issues
        error_patterns = [
            "PROGRAM STOPPED",
            "ERROR",
            "SCF has not converged"
        ]
        
        has_success = any(pattern in content for pattern in success_patterns)
        has_error = any(pattern in content for pattern in error_patterns)
        
        return has_success and not has_error
    
    except Exception:
        return False

def extract_final_energy(job_dir: Path) -> Optional[float]:
    """Extract final energy from CP2K output"""
    output_files = [job_dir / f for f in ["cp2k_sim.out", "output.out", "cp2k.out"]]
    
    # Find output file
    output_file = None
    for f in output_files:
        if f.exists():
            output_file = f
            break
    
    # Also check for project-specific output
    for f in job_dir.glob("*.out"):
        if f.exists():
            output_file = f
            break
    
    if not output_file:
        return None
    
    try:
        lines = output_file.read_text().strip().split('\n')
        
        # Look for total energy lines
        energy_patterns = [
            r'Total energy:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',
            r'ENERGY\|\s+Total FORCE_EVAL.*:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',
            r'Total Electronic Energy:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'
        ]
        
        # Search from the end of the file
        for line in reversed(lines[-100:]):
            for pattern in energy_patterns:
                match = re.search(pattern, line)
                if match:
                    return float(match.group(1))
        
        return None
    
    except Exception:
        return None

def extract_forces(job_dir: Path) -> Optional[List[List[float]]]:
    """Extract forces from CP2K output"""
    output_files = [job_dir / f for f in ["cp2k_sim.out", "output.out", "cp2k.out"]]
    
    # Find output file
    output_file = None
    for f in output_files:
        if f.exists():
            output_file = f
            break
    
    if not output_file:
        return None
    
    try:
        content = output_file.read_text()
        
        # Find the last forces section
        force_sections = re.findall(
            r'ATOMIC FORCES.*?\n((?:.*?\n)*?).*?SUM OF ATOMIC FORCES',
            content, re.DOTALL
        )
        
        if not force_sections:
            return None
        
        last_forces = force_sections[-1]
        forces = []
        
        for line in last_forces.split('\n'):
            if re.match(r'\s*\d+', line):
                parts = line.split()
                if len(parts) >= 6:  # atom, kind, element, x, y, z
                    try:
                        fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                        forces.append([fx, fy, fz])
                    except ValueError:
                        continue
        
        return forces if forces else None
    
    except Exception:
        return None

# ================================================================
# File Management
# ================================================================

def backup_previous_calculation(job_dir: Path) -> None:
    """Backup previous calculation files if they exist"""
    
    backup_dir = job_dir / "previous_run"
    
    files_to_backup = [
        "*.out", "*.log", "*.ener", "*.xyz", "*.restart*",
        "*.wfn", "*.cube", "*.pdos", "*.stress"
    ]
    
    files_found = []
    for pattern in files_to_backup:
        files_found.extend(job_dir.glob(pattern))
    
    if files_found:
        backup_dir.mkdir(exist_ok=True)
        print(f"Backing up {len(files_found)} files to {backup_dir}")
        
        for filepath in files_found:
            backup_path = backup_dir / filepath.name
            shutil.move(str(filepath), str(backup_path))

def copy_structure_files(source_dir: Path, target_dir: Path) -> None:
    """Copy structure files between directories"""
    
    # Try restart file first, then XYZ files
    for pattern in ["*-pos-1.xyz", "*.xyz", "input.inp"]:
        source_files = list(source_dir.glob(pattern))
        if source_files:
            source_file = source_files[-1]  # Get the most recent
            target_file = target_dir / "geometry.xyz"
            shutil.copy2(source_file, target_file)
            print(f"Copied {source_file.name} from {source_dir} to {target_file}")
            return
    
    print(f"Warning: No structure file found in {source_dir}")

def write_run_info(job_dir: Path, run_info: Dict, mineral_name: str,
                  simulation_params: Dict) -> None:
    """Write detailed run information to file"""
    
    info_file = job_dir / "run_info.txt"
    
    with open(info_file, "w") as f:
        f.write(f"CP2K Calculation Information\n")
        f.write(f"{'='*40}\n\n")
        
        f.write(f"Job Name: {run_info['job_name']}\n")
        f.write(f"Mineral: {mineral_name}\n")
        f.write(f"Start Time: {run_info.get('start_time', 'N/A')}\n")
        f.write(f"End Time: {run_info.get('end_time', 'N/A')}\n")
        f.write(f"Runtime: {run_info.get('runtime', 'N/A')}\n")
        f.write(f"Success: {run_info['success']}\n")
        f.write(f"Converged: {run_info.get('converged', 'Unknown')}\n")
        f.write(f"Command: {run_info.get('command', 'Not specified')}\n\n")
        
        # Extract final energy if available
        final_energy = extract_final_energy(job_dir)
        if final_energy is not None:
            f.write(f"Final Energy: {final_energy:.6f} Hartree\n\n")
        
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
        total_runtime = sum([info['runtime'] for info in run_infos], timedelta())
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
                f.write(f"{job_type}: {final_energy:.6f} Hartree\n")
        f.write("\n")
        
        # Output files description
        f.write("OUTPUT FILES:\n")
        f.write("-"*25 + "\n")
        f.write("Each calculation directory contains:\n")
        f.write("  - input.inp (CP2K input file)\n")
        f.write("  - *.out (CP2K output file)\n")
        f.write("  - *.xyz (trajectory files)\n")
        f.write("  - *.ener (energy file)\n")
        f.write("  - *.restart (restart files)\n")
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

# ================================================================
# Utility Functions
# ================================================================

def convert_units(value: float, from_unit: str, to_unit: str) -> float:
    """Convert between different units commonly used in CP2K"""
    
    # Energy conversions (to Hartree)
    energy_to_hartree = {
        'hartree': 1.0,
        'ev': 0.03674930814,
        'ry': 0.5,
        'kcal/mol': 0.00159360144,
        'kj/mol': 0.00038087988
    }
    
    # Length conversions (to Bohr)
    length_to_bohr = {
        'bohr': 1.0,
        'angstrom': 1.88972687777,
        'ang': 1.88972687777,
        'pm': 0.0188972687777
    }
    
    # Time conversions (to fs)
    time_to_fs = {
        'fs': 1.0,
        'ps': 1000.0,
        'ns': 1000000.0,
        'au': 0.024188843265857
    }
    
    from_unit = from_unit.lower()
    to_unit = to_unit.lower()
    
    if from_unit in energy_to_hartree and to_unit in energy_to_hartree:
        # Energy conversion
        hartree_value = value * energy_to_hartree[from_unit]
        return hartree_value / energy_to_hartree[to_unit]
    
    elif from_unit in length_to_bohr and to_unit in length_to_bohr:
        # Length conversion
        bohr_value = value * length_to_bohr[from_unit]
        return bohr_value / length_to_bohr[to_unit]
    
    elif from_unit in time_to_fs and to_unit in time_to_fs:
        # Time conversion
        fs_value = value * time_to_fs[from_unit]
        return fs_value / time_to_fs[to_unit]
    
    else:
        raise ValueError(f"Unsupported unit conversion: {from_unit} to {to_unit}")

def validate_cp2k_input(input_file: Path) -> List[str]:
    """Validate CP2K input file and return list of potential issues"""
    
    issues = []
    
    if not input_file.exists():
        issues.append(f"Input file {input_file} does not exist")
        return issues
    
    try:
        content = input_file.read_text().upper()
        
        # Check for required sections
        required_sections = ['&GLOBAL', '&FORCE_EVAL']
        for section in required_sections:
            if section not in content:
                issues.append(f"Missing required section: {section}")
        
        # Check for common issues
        if '&DFT' in content:
            if 'BASIS_SET_FILE_NAME' not in content:
                issues.append("Missing BASIS_SET_FILE_NAME in DFT section")
            if 'POTENTIAL_FILE_NAME' not in content:
                issues.append("Missing POTENTIAL_FILE_NAME in DFT section")
        
        # Check for unmatched sections
        sections = re.findall(r'&(\w+)', content)
        end_sections = re.findall(r'&END (\w+)', content)
        
        for section in sections:
            if section != 'END' and section not in end_sections:
                issues.append(f"Unmatched section: &{section}")
    
    except Exception as e:
        issues.append(f"Error reading input file: {str(e)}")
    
    return issues