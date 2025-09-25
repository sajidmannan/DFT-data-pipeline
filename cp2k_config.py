#!/usr/bin/env python3
"""
CP2K Configuration Module

Contains default parameters, mineral-specific configurations,
and CP2K-specific settings for simulations.
"""

# ================================================================
# Default CP2K Parameters
# ================================================================
import os
def get_functional_tags(functional: str) -> dict:
    """
    Return CP2K tags for the chosen exchange-correlation functional.
    """
    functional_upper = functional.upper()
    
    tags_mapping = {
        "PBE": {"XC_FUNCTIONAL": "PBE", "VDW_POTENTIAL": "D3"},
        "PBE0": {"XC_FUNCTIONAL": "PBE0", "VDW_POTENTIAL": "D3"},
        "B3LYP": {"XC_FUNCTIONAL": "B3LYP", "VDW_POTENTIAL": "D3"},
        "LDA": {"XC_FUNCTIONAL": "LDA", "VDW_POTENTIAL": None},
        "SCAN": {"XC_FUNCTIONAL": "SCAN", "VDW_POTENTIAL": "D3"},
    }
    
    if functional_upper not in tags_mapping:
        raise ValueError(f"Functional '{functional}' not supported. Choose from {list(tags_mapping.keys())}")
    
    return tags_mapping[functional_upper]


DEFAULT_PARAMS = {
    # Global settings
    'project_name': 'mineral_simulation',
    'run_type': 'GEO_OPT',
    'print_level': 'MEDIUM',
    
    # Force evaluation settings
    'method': 'DFT',
    'cutoff': 400,              # Plane wave cutoff (Ry)
    'rel_cutoff': 50,           # Relative cutoff for Gaussian grid (Ry)
    'scf_max_iter': 50,         # Maximum SCF iterations
    'scf_eps': 1E-6,            # SCF convergence criterion
    'scf_guess': 'ATOMIC',      # Initial guess for SCF
    
    # Exchange-correlation
    'xc_functional': 'PBE',     # Exchange-correlation functional
    'vdw_potential': None,      # Van der Waals potential
    
    # Basis sets and pseudopotentials
    'basis_set_type': 'DZVP-MOLOPT-SR-GTH',
    'potential_type': 'GTH-PBE',
    
    # Geometry optimization
    'geo_max_iter': 200,        # Maximum geometry optimization steps
    'geo_max_force': 4.5E-4,    # Maximum force convergence (Hartree/Bohr)
    'geo_max_dr': 3.0E-3,       # Maximum displacement convergence (Bohr)
    'geo_rms_force': 3.0E-4,    # RMS force convergence (Hartree/Bohr)
    'geo_rms_dr': 1.5E-3,       # RMS displacement convergence (Bohr)
    'optimizer': 'BFGS',        # Geometry optimizer
    
    # Cell optimization
    'cell_opt': False,          # Enable cell optimization
    'pressure_tolerance': 100,  # Pressure tolerance (bar)
    
    # Molecular dynamics
    'md_ensemble': 'NVT',       # MD ensemble
    'md_temperature': 300.0,    # Temperature (K)
    'md_timestep': 0.5,         # Time step (fs)
    'md_steps': 1000,           # Number of MD steps
    'thermostat': 'NOSE',       # Thermostat type
    'thermostat_tau': 100.0,    # Thermostat time constant (fs)
    
    # Output control
    'print_forces': True,       # Print forces
    'print_stress': True,       # Print stress tensor
    'trajectory_format': 'XYZ', # Trajectory file format
    'restart_backup': 'ON',     # Backup restart files
}

# ================================================================
# Mineral-Specific Configurations
# ================================================================

MINERAL_DATABASE = {
    'sio2': {
        'name': 'Silicon Dioxide (Quartz)',
        'formula': 'SiO2',
        'space_group': 'P3121',
        'density': 2.65,
        'elements': ['Si', 'O'],
        'cell_params': [4.916, 4.916, 5.405, 90.0, 90.0, 120.0],  # a,b,c,α,β,γ
        'recommended_settings': {
            'cutoff': 500,
            'rel_cutoff': 60,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Si': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Si': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Alpha-quartz structure, hexagonal'
    },
    
    'quartz': {
        'name': 'Quartz',
        'formula': 'SiO2',
        'space_group': 'P3121',
        'density': 2.65,
        'elements': ['Si', 'O'],
        'cell_params': [4.916, 4.916, 5.405, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 500,
            'rel_cutoff': 60,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Si': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Si': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Same as SiO2'
    },
    
    'caco3': {
        'name': 'Calcium Carbonate (Calcite)',
        'formula': 'CaCO3',
        'space_group': 'R-3c',
        'density': 2.71,
        'elements': ['Ca', 'C', 'O'],
        'cell_params': [4.99, 4.99, 17.06, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Ca': 'DZVP-MOLOPT-SR-GTH',
            'C': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Ca': 'GTH-PBE-q10',
            'C': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Calcite structure, rhombohedral'
    },
    
    'calcite': {
        'name': 'Calcite',
        'formula': 'CaCO3',
        'space_group': 'R-3c',
        'density': 2.71,
        'elements': ['Ca', 'C', 'O'],
        'cell_params': [4.99, 4.99, 17.06, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Ca': 'DZVP-MOLOPT-SR-GTH',
            'C': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Ca': 'GTH-PBE-q10',
            'C': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'CaCO3 polymorph'
    },
    
    'mgo': {
        'name': 'Magnesium Oxide (Periclase)',
        'formula': 'MgO',
        'space_group': 'Fm-3m',
        'density': 3.58,
        'elements': ['Mg', 'O'],
        'cell_params': [4.212, 4.212, 4.212, 90.0, 90.0, 90.0],
        'recommended_settings': {
            'cutoff': 450,
            'rel_cutoff': 50,
            'cell_opt': True,
            'scf_eps': 1E-6,
        },
        'basis_sets': {
            'Mg': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Mg': 'GTH-PBE-q10',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Cubic rock salt structure'
    },
    
    'periclase': {
        'name': 'Periclase',
        'formula': 'MgO',
        'space_group': 'Fm-3m',
        'density': 3.58,
        'elements': ['Mg', 'O'],
        'cell_params': [4.212, 4.212, 4.212, 90.0, 90.0, 90.0],
        'recommended_settings': {
            'cutoff': 450,
            'rel_cutoff': 50,
            'cell_opt': True,
            'scf_eps': 1E-6,
        },
        'basis_sets': {
            'Mg': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Mg': 'GTH-PBE-q10',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'MgO cubic structure'
    },
    
    'mg2sio4': {
        'name': 'Magnesium Silicate (Forsterite)',
        'formula': 'Mg2SiO4',
        'space_group': 'Pbnm',
        'density': 3.22,
        'elements': ['Mg', 'Si', 'O'],
        'cell_params': [4.754, 10.207, 5.980, 90.0, 90.0, 90.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Mg': 'DZVP-MOLOPT-SR-GTH',
            'Si': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Mg': 'GTH-PBE-q10',
            'Si': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Olivine group mineral, orthorhombic'
    },
    
    'forsterite': {
        'name': 'Forsterite',
        'formula': 'Mg2SiO4',
        'space_group': 'Pbnm',
        'density': 3.22,
        'elements': ['Mg', 'Si', 'O'],
        'cell_params': [4.754, 10.207, 5.980, 90.0, 90.0, 90.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Mg': 'DZVP-MOLOPT-SR-GTH',
            'Si': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Mg': 'GTH-PBE-q10',
            'Si': 'GTH-PBE-q4',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Mg2SiO4 olivine endmember'
    },
    
    'al2o3': {
        'name': 'Aluminum Oxide (Corundum)',
        'formula': 'Al2O3',
        'space_group': 'R-3c',
        'density': 3.95,
        'elements': ['Al', 'O'],
        'cell_params': [4.759, 4.759, 12.991, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 550,
            'rel_cutoff': 60,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Al': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Al': 'GTH-PBE-q3',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Corundum structure, hexagonal'
    },
    
    'corundum': {
        'name': 'Corundum',
        'formula': 'Al2O3',
        'space_group': 'R-3c',
        'density': 3.95,
        'elements': ['Al', 'O'],
        'cell_params': [4.759, 4.759, 12.991, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 550,
            'rel_cutoff': 60,
            'cell_opt': True,
            'scf_eps': 1E-7,
        },
        'basis_sets': {
            'Al': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Al': 'GTH-PBE-q3',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Al2O3 structure'
    },
    
    'fe2o3': {
        'name': 'Iron Oxide (Hematite)',
        'formula': 'Fe2O3',
        'space_group': 'R-3c',
        'density': 5.26,
        'elements': ['Fe', 'O'],
        'cell_params': [5.034, 5.034, 13.747, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-6,
            'uks': True,  # Unrestricted Kohn-Sham (spin polarized)
            'multiplicity': 1,  # Antiferromagnetic ground state
        },
        'basis_sets': {
            'Fe': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Fe': 'GTH-PBE-q16',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Hematite structure, antiferromagnetic'
    },
    
    'hematite': {
        'name': 'Hematite',
        'formula': 'Fe2O3',
        'space_group': 'R-3c',
        'density': 5.26,
        'elements': ['Fe', 'O'],
        'cell_params': [5.034, 5.034, 13.747, 90.0, 90.0, 120.0],
        'recommended_settings': {
            'cutoff': 600,
            'rel_cutoff': 70,
            'cell_opt': True,
            'scf_eps': 1E-6,
            'uks': True,
            'multiplicity': 1,
        },
        'basis_sets': {
            'Fe': 'DZVP-MOLOPT-SR-GTH',
            'O': 'DZVP-MOLOPT-SR-GTH'
        },
        'pseudopotentials': {
            'Fe': 'GTH-PBE-q16',
            'O': 'GTH-PBE-q6'
        },
        'notes': 'Fe2O3 structure, magnetic'
    },
}

# ================================================================
# Basis Set Information
# ================================================================

BASIS_SET_INFO = {
    # Light elements
    'H': {'default': 'DZVP-MOLOPT-GTH', 'electrons': 1},
    'He': {'default': 'DZVP-MOLOPT-GTH', 'electrons': 2},
    'Li': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 3},
    'Be': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 4},
    'B': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 3},
    'C': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 4},
    'N': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 5},
    'O': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 6},
    'F': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 7},
    'Ne': {'default': 'DZVP-MOLOPT-GTH', 'electrons': 8},
    
    # Main group elements
    'Na': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 9},
    'Mg': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 10},
    'Al': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 3},
    'Si': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 4},
    'P': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 5},
    'S': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 6},
    'Cl': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 7},
    'Ar': {'default': 'DZVP-MOLOPT-GTH', 'electrons': 8},
    'K': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 9},
    'Ca': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 10},
    
    # Transition metals
    'Sc': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 11},
    'Ti': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 12},
    'V': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 13},
    'Cr': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 14},
    'Mn': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 15},
    'Fe': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 16},
    'Co': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 17},
    'Ni': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 18},
    'Cu': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 19},
    'Zn': {'default': 'DZVP-MOLOPT-SR-GTH', 'electrons': 12},
}

# ================================================================
# Pseudopotential Information  
# ================================================================

PSEUDOPOTENTIAL_INFO = {
    # GTH pseudopotentials with valence electron count
    'H': {'default': 'GTH-PBE-q1', 'valence': 1},
    'He': {'default': 'GTH-PBE-q2', 'valence': 2},
    'Li': {'default': 'GTH-PBE-q3', 'valence': 3},
    'Be': {'default': 'GTH-PBE-q4', 'valence': 4},
    'B': {'default': 'GTH-PBE-q3', 'valence': 3},
    'C': {'default': 'GTH-PBE-q4', 'valence': 4},
    'N': {'default': 'GTH-PBE-q5', 'valence': 5},
    'O': {'default': 'GTH-PBE-q6', 'valence': 6},
    'F': {'default': 'GTH-PBE-q7', 'valence': 7},
    'Ne': {'default': 'GTH-PBE-q8', 'valence': 8},
    'Na': {'default': 'GTH-PBE-q9', 'valence': 9},
    'Mg': {'default': 'GTH-PBE-q10', 'valence': 10},
    'Al': {'default': 'GTH-PBE-q3', 'valence': 3},
    'Si': {'default': 'GTH-PBE-q4', 'valence': 4},
    'P': {'default': 'GTH-PBE-q5', 'valence': 5},
    'S': {'default': 'GTH-PBE-q6', 'valence': 6},
    'Cl': {'default': 'GTH-PBE-q7', 'valence': 7},
    'Ar': {'default': 'GTH-PBE-q8', 'valence': 8},
    'K': {'default': 'GTH-PBE-q9', 'valence': 9},
    'Ca': {'default': 'GTH-PBE-q10', 'valence': 10},
    'Sc': {'default': 'GTH-PBE-q11', 'valence': 11},
    'Ti': {'default': 'GTH-PBE-q12', 'valence': 12},
    'V': {'default': 'GTH-PBE-q13', 'valence': 13},
    'Cr': {'default': 'GTH-PBE-q14', 'valence': 14},
    'Mn': {'default': 'GTH-PBE-q15', 'valence': 15},
    'Fe': {'default': 'GTH-PBE-q16', 'valence': 16},
    'Co': {'default': 'GTH-PBE-q17', 'valence': 17},
    'Ni': {'default': 'GTH-PBE-q18', 'valence': 18},
    'Cu': {'default': 'GTH-PBE-q19', 'valence': 19},
    'Zn': {'default': 'GTH-PBE-q12', 'valence': 12},
}

# ================================================================
# Exchange-Correlation Functionals
# ================================================================

XC_FUNCTIONALS = {
    'PBE': {
        'type': 'GGA',
        'description': 'Perdew-Burke-Ernzerhof GGA',
        'libxc': 'XC_GGA_X_PBE XC_GGA_C_PBE',
        'recommended_for': ['general_purpose', 'solids'],
    },
    
    'PBE0': {
        'type': 'Hybrid GGA',
        'description': 'PBE hybrid functional (25% exact exchange)',
        'libxc': 'XC_HYB_GGA_XC_PBEH',
        'recommended_for': ['band_gaps', 'excited_states'],
        'expensive': True,
    },
    
    'B3LYP': {
        'type': 'Hybrid GGA',
        'description': 'Becke 3-parameter Lee-Yang-Parr hybrid functional',
        'libxc': 'XC_HYB_GGA_XC_B3LYP',
        'recommended_for': ['molecules', 'band_gaps'],
        'expensive': True,
    },
    
    'BLYP': {
        'type': 'GGA',
        'description': 'Becke-Lee-Yang-Parr GGA',
        'libxc': 'XC_GGA_X_B88 XC_GGA_C_LYP',
        'recommended_for': ['molecules', 'soft_matter'],
    },
    
    'BP86': {
        'type': 'GGA',
        'description': 'Becke-Perdew 86 GGA',
        'libxc': 'XC_GGA_X_B88 XC_GGA_C_P86',
        'recommended_for': ['molecules', 'transition_metals'],
    },
    
    'LDA': {
        'type': 'LDA',
        'description': 'Local Density Approximation',
        'libxc': 'XC_LDA_X XC_LDA_C_PZ',
        'recommended_for': ['high_pressure', 'simple_systems'],
    },
    
    'SCAN': {
        'type': 'Meta-GGA',
        'description': 'Strongly Constrained and Appropriately Normed',
        'libxc': 'XC_MGGA_X_SCAN XC_MGGA_C_SCAN',
        'recommended_for': ['accurate_energetics', 'phase_stability'],
    },
}

# ================================================================
# Van der Waals Corrections
# ================================================================

VDW_CORRECTIONS = {
    'D2': {
        'description': 'Grimme D2 dispersion correction',
        'potential_type': 'PAIR_POTENTIAL',
        'parameter_file': 'dftd2_parameters',
        'recommended_for': ['layered_materials', 'molecular_crystals'],
    },
    
    'D3': {
        'description': 'Grimme D3 dispersion correction',
        'potential_type': 'PAIR_POTENTIAL', 
        'parameter_file': 'dftd3_parameters',
        'recommended_for': ['layered_materials', 'molecular_crystals', 'general_purpose'],
    },
    
    'D3BJ': {
        'description': 'Grimme D3 with Becke-Johnson damping',
        'potential_type': 'PAIR_POTENTIAL',
        'parameter_file': 'dftd3_parameters',
        'damping': 'BJ',
        'recommended_for': ['layered_materials', 'interfaces'],
    },
    
    'TS': {
        'description': 'Tkatchenko-Scheffler method',
        'potential_type': 'PAIR_POTENTIAL',
        'parameter_file': 'ts_parameters',
        'recommended_for': ['molecules', 'surfaces'],
    },
}

# ================================================================
# Simulation Templates
# ================================================================

SIMULATION_TEMPLATES = {
    'opt': {
        'description': 'Geometry optimization',
        'run_type': 'GEO_OPT',
        'geo_max_iter': 200,
        'optimizer': 'BFGS',
        'cell_opt': False,
    },
    
    'cell_opt': {
        'description': 'Cell and geometry optimization',
        'run_type': 'CELL_OPT',
        'geo_max_iter': 200,
        'optimizer': 'BFGS',
        'cell_opt': True,
        'pressure_tolerance': 100,
    },
    
    'sp': {
        'description': 'Single point energy calculation',
        'run_type': 'ENERGY',
        'scf_max_iter': 50,
    },
    
    'md': {
        'description': 'Molecular dynamics simulation',
        'run_type': 'MD',
        'md_ensemble': 'NVT',
        'md_steps': 1000,
        'md_timestep': 0.5,
        'md_temperature': 300.0,
        'thermostat': 'NOSE',
        'thermostat_tau': 100.0,
    },
    
    'npt': {
        'description': 'NPT molecular dynamics simulation',
        'run_type': 'MD',
        'md_ensemble': 'NPT_I',
        'md_steps': 5000,
        'md_timestep': 0.5,
        'md_temperature': 300.0,
        'thermostat': 'NOSE',
        'thermostat_tau': 100.0,
        'barostat': 'NOSE',
        'barostat_tau': 1000.0,
    },
    
    'vib': {
        'description': 'Vibrational analysis',
        'run_type': 'VIBRATIONAL_ANALYSIS',
        'vibrations_proc': 'PRINT',
        'dx': 0.01,  # Displacement for finite differences
    },
    
    'band': {
        'description': 'Band structure calculation',
        'run_type': 'ENERGY_FORCE',
        'added_mos': 50,
        'scf_max_iter': 50,
        'kpoints': 'band_structure',
    },
}

# Default fallback basis sets and pseudopotentials
DEFAULT_BASIS_SETS = {
    'H': 'DZVP-MOLOPT-SR-GTH',
    'C': 'DZVP-MOLOPT-SR-GTH',
    'N': 'DZVP-MOLOPT-SR-GTH',
    'O': 'DZVP-MOLOPT-SR-GTH',
    'Si': 'DZVP-MOLOPT-SR-GTH',
    'Al': 'DZVP-MOLOPT-SR-GTH',
    'Fe': 'DZVP-MOLOPT-SR-GTH',
}

DEFAULT_PSEUDOPOTENTIALS = {
    'H': 'GTH-PBE',
    'C': 'GTH-PBE-q4',
    'N': 'GTH-PBE-q5',
    'O': 'GTH-PBE-q6',
    'Si': 'GTH-PBE-q4',
    'Al': 'GTH-PBE-q3',
    'Fe': 'GTH-PBE-q8',
}

# ================================================================
# Utility Functions
# ================================================================


def get_mineral_info(mineral_name: str) -> dict:
    """Get information about a specific mineral from the database."""
    base = os.path.basename(mineral_name)  # Strip path
    key = base.lower().replace('.xyz', '')  # Normalize key
    if key in MINERAL_DATABASE:
        return MINERAL_DATABASE[key]
    else:
        # Generic fallback info
        return {
            'name': f'Unknown mineral: {mineral_name}',
            'formula': 'Unknown',
            'space_group': 'Unknown',
            'density': None,
            'elements': [],
            'cell_params': [5.0, 5.0, 5.0, 90.0, 90.0, 90.0],
            'recommended_settings': {
                'cutoff': 400,
                'rel_cutoff': 50,
                'cell_opt': True,
                'scf_eps': 1E-6,
            },
            'basis_sets': {},
            'pseudopotentials': {},
            'notes': 'Generic settings used'
        }

def get_recommended_cutoff(elements: list, margin: float = 1.2) -> int:
    """Get recommended cutoff based on element requirements"""
    if not elements:
        return DEFAULT_PARAMS['cutoff']
    
    # Base cutoffs for different element types
    base_cutoffs = {
        'H': 300, 'He': 400, 'Li': 350, 'Be': 350, 'B': 400, 'C': 400, 'N': 400, 'O': 400, 'F': 400, 'Ne': 400,
        'Na': 350, 'Mg': 400, 'Al': 450, 'Si': 450, 'P': 450, 'S': 450, 'Cl': 450, 'Ar': 450,
        'K': 400, 'Ca': 450, 'Sc': 500, 'Ti': 500, 'V': 500, 'Cr': 500, 'Mn': 500, 'Fe': 500,
        'Co': 500, 'Ni': 500, 'Cu': 500, 'Zn': 450
    }
    
    max_cutoff = 0
    for element in elements:
        cutoff = base_cutoffs.get(element, 450)  # Default for unknown elements
        max_cutoff = max(max_cutoff, cutoff)
    
    return int(max_cutoff * margin)

def get_basis_sets(mineral_name: str) -> dict:
    """Return basis sets for the given mineral, falling back to defaults if missing."""
    mineral_info = get_mineral_info(mineral_name)
    basis = mineral_info.get('basis_sets', {})
    
    # Fill missing elements with default basis sets
    for element in mineral_info.get('elements', []):
        if element not in basis:
            basis[element] = DEFAULT_BASIS_SETS.get(element, 'DZVP-MOLOPT-SR-GTH')
    
    return basis

def get_pseudopotentials(mineral_name: str) -> dict:
    """Return pseudopotentials for the given mineral, falling back to defaults if missing."""
    mineral_info = get_mineral_info(mineral_name)
    pseudo = mineral_info.get('pseudopotentials', {})
    
    # Fill missing elements with default pseudopotentials
    for element in mineral_info.get('elements', []):
        if element not in pseudo:
            pseudo[element] = DEFAULT_PSEUDOPOTENTIALS.get(element, 'GTH-PBE')
    
    return pseudo





# def get_basis_sets(mineral_name: str) -> dict:
#     """Get basis sets for a mineral"""
#     mineral_info = get_mineral_info(mineral_name)
#     print(mineral_info.get('basis_sets', {}))
#     return mineral_info.get('basis_sets', {})

# def get_pseudopotentials(mineral_name: str) -> dict:
#     """Get pseudopotentials for a mineral"""
#     mineral_info = get_mineral_info(mineral_name)
#     return mineral_info.get('pseudopotentials', {})



def validate_functional(functional: str) -> bool:
    """Validate if functional is supported"""
    return functional.upper() in XC_FUNCTIONALS

def get_functional_info(functional: str) -> dict:
    """Get functional information"""
    func_upper = functional.upper()
    if func_upper in XC_FUNCTIONALS:
        return XC_FUNCTIONALS[func_upper]
    else:
        return XC_FUNCTIONALS['PBE']  # Default to PBE

def get_vdw_info(vdw_type: str) -> dict:
    """Get van der Waals correction information"""
    if vdw_type and vdw_type.upper() in VDW_CORRECTIONS:
        return VDW_CORRECTIONS[vdw_type.upper()]
    else:
        return {}

def get_simulation_template(template_name: str) -> dict:
    """Get simulation template parameters"""
    if template_name in SIMULATION_TEMPLATES:
        return SIMULATION_TEMPLATES[template_name].copy()
    else:
        raise ValueError(f"Unknown template: {template_name}. "
                        f"Available: {list(SIMULATION_TEMPLATES.keys())}")

def list_available_minerals():
    """List all minerals in the database"""
    print("Available minerals with predefined parameters:")
    print("=" * 60)
    
    for key, info in MINERAL_DATABASE.items():
        cell = info['cell_params']
        elements = ', '.join(info['elements'])
        print(f"{key:12s} - {info['name']} ({info['formula']})")
        print(f"             Cell: {cell[0]:.3f} {cell[1]:.3f} {cell[2]:.3f} Å")
        print(f"             Elements: {elements}")
        print(f"             Cutoff: {info['recommended_settings']['cutoff']} Ry")
        print()

def list_available_functionals():
    """List available exchange-correlation functionals"""
    print("Available Exchange-Correlation Functionals:")
    print("=" * 50)
    
    for name, info in XC_FUNCTIONALS.items():
        print(f"{name:8s} - {info['description']} ({info['type']})")
        if info.get('expensive'):
            print("           ⚠️  Computationally expensive")
        print(f"           Recommended for: {', '.join(info['recommended_for'])}")
        print()

def list_simulation_templates():
    """List available simulation templates"""
    print("Available Simulation Templates:")
    print("=" * 40)
    
    for name, template in SIMULATION_TEMPLATES.items():
        print(f"{name:12s} - {template['description']}")
        print(f"             Run type: {template['run_type']}")
        if 'md_steps' in template:
            print(f"             Steps: {template['md_steps']}")
        if 'geo_max_iter' in template:
            print(f"             Max iterations: {template['geo_max_iter']}")
        print()

def list_available_basis_sets():
    """List available basis sets"""
    print("Common Basis Sets:")
    print("=" * 30)
    
    basis_sets = [
        'DZVP-MOLOPT-SR-GTH - Double-zeta valence polarized, optimized for GTH pseudopotentials',
        'TZVP-MOLOPT-GTH - Triple-zeta valence polarized',
        'TZV2P-MOLOPT-GTH - Triple-zeta with 2 polarization functions',
        'SZV-MOLOPT-GTH - Single-zeta valence',
        'DZVP-GTH - Standard double-zeta',
        'aug-DZVP-MOLOPT-SR-GTH - Augmented double-zeta'
    ]
    
    for basis_set in basis_sets:
        print(f"  {basis_set}")
    print()

if __name__ == "__main__":
    # Demo the configuration module
    print("CP2K Configuration Module Demo")
    print("=" * 40)
    
    list_available_minerals()
    list_available_functionals()
    list_simulation_templates()
    list_available_basis_sets()
    
    print("\nTesting mineral lookup:")
    test_minerals = ['sio2', 'calcite', 'unknown_mineral']
    for mineral in test_minerals:
        info = get_mineral_info(mineral)
        basis_sets = get_basis_sets(mineral)
        recommended_cutoff = get_recommended_cutoff(info['elements'])
        print(f"{mineral}: {info['name']}")
        print(f"  Elements: {info['elements']}")
        print(f"  Recommended cutoff: {recommended_cutoff} Ry")
        print(f"  Basis sets: {basis_sets}")
        print()