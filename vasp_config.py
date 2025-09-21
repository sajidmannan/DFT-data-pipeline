#!/usr/bin/env python3
"""
VASP Configuration Module

Contains default parameters, mineral-specific configurations,
and VASP-specific settings for simulations.
"""

# ================================================================
# Default VASP Parameters
# ================================================================

DEFAULT_PARAMS = {
    # Electronic structure
    'encut': 500,               # Plane wave cutoff (eV)
    'prec': 'Accurate',         # Precision level
    'algo': 'Fast',             # Electronic minimization algorithm
    'ediff': 1E-6,              # Electronic convergence criterion (eV)
    'nelm': 100,                # Maximum electronic steps
    
    # Ionic relaxation
    'ibrion': 2,                # Ionic relaxation method (2 = conjugate gradient)
    'nsw': 200,                 # Maximum ionic steps
    'ediffg': -0.01,            # Ionic convergence criterion (eV/Å)
    
    # Molecular dynamics
    'potim': 0.5,               # Time step for MD (fs)
    'tebeg': 300,               # Initial temperature (K)
    'teend': 300,               # Final temperature (K)
    'smass': 0,                 # Nose thermostat mass
    
    # K-point sampling
    'kpts_grid': [6, 6, 6],     # Gamma-centered k-point grid
    'kpts_shift': [0.0, 0.0, 0.0],  # K-point shift
    
    # Exchange-correlation
    'xc_functional': 'PBE',     # Exchange-correlation functional
    'luse_vdw': False,          # Van der Waals corrections
    'vdw_type': None,           # Type of vdW correction
    
    # Output control
    'lwave': False,             # Write WAVECAR
    'lcharg': True,             # Write CHGCAR
    'lorbit': 11,               # Orbital projection
    'lreal': 'Auto',            # Real space projection
    
    # Spin settings
    'ispin': 1,                 # Non-spin polarized (1) or spin polarized (2)
    'magmom': None,             # Initial magnetic moments
    
    # Pressure settings (for NPT-like simulations)
    'pstress': 0.0,             # External pressure (kBar)
    'isif': 2,                  # Stress/relaxation flag
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
            'encut': 520,
            'kpts_grid': [8, 8, 6],
            'ediffg': -0.01,
            'isif': 3,  # Relax cell shape and volume
        },
        'potcar_elements': ['Si', 'O'],
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
            'encut': 520,
            'kpts_grid': [8, 8, 6],
            'ediffg': -0.01,
            'isif': 3,
        },
        'potcar_elements': ['Si', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 6, 4],
            'ediffg': -0.015,
            'isif': 3,
        },
        'potcar_elements': ['Ca', 'C', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 6, 4],
            'ediffg': -0.015,
            'isif': 3,
        },
        'potcar_elements': ['Ca', 'C', 'O'],
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
            'encut': 500,
            'kpts_grid': [8, 8, 8],
            'ediffg': -0.01,
            'isif': 3,
        },
        'potcar_elements': ['Mg', 'O'],
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
            'encut': 500,
            'kpts_grid': [8, 8, 8],
            'ediffg': -0.01,
            'isif': 3,
        },
        'potcar_elements': ['Mg', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 4, 6],
            'ediffg': -0.015,
            'isif': 3,
        },
        'potcar_elements': ['Mg', 'Si', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 4, 6],
            'ediffg': -0.015,
            'isif': 3,
        },
        'potcar_elements': ['Mg', 'Si', 'O'],
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
            'encut': 550,
            'kpts_grid': [8, 8, 4],
            'ediffg': -0.01,
            'isif': 3,
        },
        'potcar_elements': ['Al', 'O'],
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
            'encut': 550,
            'kpts_grid': [8, 8, 4],
            'ediffg': -0.01,
            'isif': 3,
        },
        'potcar_elements': ['Al', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 6, 3],
            'ediffg': -0.02,
            'isif': 3,
            'ispin': 2,  # Spin polarized
            'magmom': [5.0, 5.0, 0.0, 0.0, 0.0],  # AFM Fe moments
        },
        'potcar_elements': ['Fe', 'O'],
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
            'encut': 600,
            'kpts_grid': [6, 6, 3],
            'ediffg': -0.02,
            'isif': 3,
            'ispin': 2,
            'magmom': [5.0, 5.0, 0.0, 0.0, 0.0],
        },
        'potcar_elements': ['Fe', 'O'],
        'notes': 'Fe2O3 structure, magnetic'
    },
}

# ================================================================
# POTCAR Information
# ================================================================

POTCAR_INFO = {
    'H': {'enmax': 250, 'version': 'H', 'electrons': 1},
    'He': {'enmax': 479, 'version': 'He', 'electrons': 2},
    'Li': {'enmax': 499, 'version': 'Li_sv', 'electrons': 3},
    'Be': {'enmax': 248, 'version': 'Be_sv', 'electrons': 4},
    'B': {'enmax': 319, 'version': 'B', 'electrons': 3},
    'C': {'enmax': 400, 'version': 'C', 'electrons': 4},
    'N': {'enmax': 400, 'version': 'N', 'electrons': 5},
    'O': {'enmax': 400, 'version': 'O', 'electrons': 6},
    'F': {'enmax': 400, 'version': 'F', 'electrons': 7},
    'Ne': {'enmax': 344, 'version': 'Ne', 'electrons': 8},
    'Na': {'enmax': 102, 'version': 'Na_pv', 'electrons': 9},
    'Mg': {'enmax': 200, 'version': 'Mg_pv', 'electrons': 10},
    'Al': {'enmax': 240, 'version': 'Al', 'electrons': 3},
    'Si': {'enmax': 245, 'version': 'Si', 'electrons': 4},
    'P': {'enmax': 255, 'version': 'P', 'electrons': 5},
    'S': {'enmax': 259, 'version': 'S', 'electrons': 6},
    'Cl': {'enmax': 262, 'version': 'Cl', 'electrons': 7},
    'Ar': {'enmax': 266, 'version': 'Ar', 'electrons': 8},
    'K': {'enmax': 117, 'version': 'K_pv', 'electrons': 9},
    'Ca': {'enmax': 267, 'version': 'Ca_pv', 'electrons': 10},
    'Sc': {'enmax': 155, 'version': 'Sc_sv', 'electrons': 11},
    'Ti': {'enmax': 178, 'version': 'Ti_pv', 'electrons': 10},
    'V': {'enmax': 192, 'version': 'V_pv', 'electrons': 11},
    'Cr': {'enmax': 227, 'version': 'Cr_pv', 'electrons': 12},
    'Mn': {'enmax': 269, 'version': 'Mn_pv', 'electrons': 13},
    'Fe': {'enmax': 267, 'version': 'Fe_pv', 'electrons': 14},
    'Co': {'enmax': 267, 'version': 'Co', 'electrons': 9},
    'Ni': {'enmax': 270, 'version': 'Ni_pv', 'electrons': 16},
    'Cu': {'enmax': 295, 'version': 'Cu_pv', 'electrons': 17},
    'Zn': {'enmax': 277, 'version': 'Zn', 'electrons': 12},
}

# ================================================================
# Exchange-Correlation Functionals
# ================================================================

XC_FUNCTIONALS = {
    'PBE': {
        'type': 'GGA',
        'description': 'Perdew-Burke-Ernzerhof GGA',
        'incar_tags': {'GGA': 'PE'},
        'recommended_for': ['general_purpose', 'solids'],
    },
    
    'PBEsol': {
        'type': 'GGA',
        'description': 'PBE revised for solids',
        'incar_tags': {'GGA': 'PS'},
        'recommended_for': ['solids', 'lattice_parameters'],
    },
    
    'LDA': {
        'type': 'LDA',
        'description': 'Local Density Approximation',
        'incar_tags': {'GGA': '--'},
        'recommended_for': ['molecules', 'high_pressure'],
    },
    
    'SCAN': {
        'type': 'Meta-GGA',
        'description': 'Strongly Constrained and Appropriately Normed',
        'incar_tags': {'METAGGA': 'SCAN'},
        'recommended_for': ['accurate_energetics', 'phase_stability'],
    },
    
    'HSE06': {
        'type': 'Hybrid',
        'description': 'Heyd-Scuseria-Ernzerhof hybrid functional',
        'incar_tags': {'LHFCALC': True, 'HFSCREEN': 0.2},
        'recommended_for': ['band_gaps', 'excited_states'],
        'expensive': True,
    },
    
    'PBE0': {
        'type': 'Hybrid',
        'description': 'PBE hybrid functional',
        'incar_tags': {'LHFCALC': True, 'AEXX': 0.25},
        'recommended_for': ['band_gaps', 'accurate_energetics'],
        'expensive': True,
    },
}

# ================================================================
# Van der Waals Corrections
# ================================================================

VDW_CORRECTIONS = {
    'D2': {
        'description': 'Grimme D2 dispersion correction',
        'incar_tags': {'LVDW': True},
        'recommended_for': ['layered_materials', 'molecular_crystals'],
    },
    
    'D3': {
        'description': 'Grimme D3 dispersion correction',
        'incar_tags': {'IVDW': 11},
        'recommended_for': ['layered_materials', 'molecular_crystals', 'general_purpose'],
    },
    
    'D3BJ': {
        'description': 'Grimme D3 with Becke-Johnson damping',
        'incar_tags': {'IVDW': 12},
        'recommended_for': ['layered_materials', 'interfaces'],
    },
    
    'TS': {
        'description': 'Tkatchenko-Scheffler method',
        'incar_tags': {'IVDW': 20},
        'recommended_for': ['molecules', 'surfaces'],
    },
    
    'BEEF': {
        'description': 'BEEF-vdW functional',
        'incar_tags': {'GGA': 'BF', 'LUSE_VDW': True, 'ZABE_VDW': 'BEEF'},
        'recommended_for': ['surfaces', 'catalysis'],
    },
}

# ================================================================
# Simulation Templates
# ================================================================

SIMULATION_TEMPLATES = {
    'relax': {
        'description': 'Structural optimization',
        'ibrion': 2,
        'nsw': 200,
        'isif': 3,
        'ediffg': -0.01,
        'prec': 'Accurate',
        'algo': 'Fast',
    },
    
    'scf': {
        'description': 'Self-consistent field calculation',
        'ibrion': -1,
        'nsw': 0,
        'prec': 'Accurate',
        'algo': 'Fast',
    },
    
    'dos': {
        'description': 'Density of states calculation',
        'ibrion': -1,
        'nsw': 0,
        'ismear': -5,
        'lorbit': 11,
        'nedos': 3000,
        'prec': 'Accurate',
    },
    
    'band': {
        'description': 'Band structure calculation',
        'ibrion': -1,
        'nsw': 0,
        'icharg': 11,
        'lorbit': 11,
        'prec': 'Accurate',
    },
    
    'md_nvt': {
        'description': 'NVT molecular dynamics',
        'ibrion': 0,
        'nsw': 5000,
        'potim': 1.0,
        'mdalgo': 2,
        'isif': 2,
        'tebeg': 300,
        'teend': 300,
        'smass': 0,
        'prec': 'Normal',
    },
    
    'md_npt': {
        'description': 'NPT molecular dynamics',
        'ibrion': 0,
        'nsw': 10000,
        'potim': 1.0,
        'mdalgo': 3,
        'isif': 3,
        'tebeg': 300,
        'teend': 300,
        'smass': 0,
        'pstress': 0.0,
        'prec': 'Normal',
    },
    
    'phonon': {
        'description': 'Phonon calculation preparation',
        'ibrion': 6,
        'nfree': 2,
        'potim': 0.015,
        'nwrite': 3,
        'prec': 'Accurate',
    },
    
    'elastic': {
        'description': 'Elastic constants calculation',
        'ibrion': 6,
        'nfree': 2,
        'potim': 0.015,
        'isif': 3,
        'nwrite': 3,
        'prec': 'Accurate',
    },
}

# ================================================================
# Utility Functions
# ================================================================

def get_mineral_info(mineral_name: str) -> dict:
    """Get information about a specific mineral"""
    key = mineral_name.lower().replace('.poscar', '').replace('_', '')
    
    if key in MINERAL_DATABASE:
        return MINERAL_DATABASE[key]
    else:
        # Return generic info for unknown minerals
        return {
            'name': f'Unknown mineral: {mineral_name}',
            'formula': 'Unknown',
            'space_group': 'Unknown',
            'density': None,
            'elements': [],
            'cell_params': [5.0, 5.0, 5.0, 90.0, 90.0, 90.0],
            'recommended_settings': {
                'encut': 500,
                'kpts_grid': [6, 6, 6],
                'ediffg': -0.01,
                'isif': 3,
            },
            'potcar_elements': [],
            'notes': 'Generic settings used'
        }

def get_recommended_encut(elements: list, margin: float = 1.3) -> int:
    """Get recommended ENCUT based on POTCAR ENMAX values"""
    if not elements:
        return DEFAULT_PARAMS['encut']
    
    max_enmax = 0
    for element in elements:
        if element in POTCAR_INFO:
            max_enmax = max(max_enmax, POTCAR_INFO[element]['enmax'])
        else:
            print(f"Warning: Unknown element {element}, using default ENMAX=400")
            max_enmax = max(max_enmax, 400)
    
    return int(max_enmax * margin)

def get_potcar_elements(mineral_name: str) -> list:
    """Get POTCAR elements for a mineral"""
    mineral_info = get_mineral_info(mineral_name)
    return mineral_info.get('potcar_elements', mineral_info.get('elements', []))

def validate_functional(functional: str) -> bool:
    """Validate if functional is supported"""
    return functional.upper() in XC_FUNCTIONALS

def get_functional_tags(functional: str) -> dict:
    """Get INCAR tags for a functional"""
    func_upper = functional.upper()
    if func_upper in XC_FUNCTIONALS:
        return XC_FUNCTIONALS[func_upper]['incar_tags']
    else:
        return {'GGA': 'PE'}  # Default to PBE

def get_vdw_tags(vdw_type: str) -> dict:
    """Get INCAR tags for van der Waals corrections"""
    if vdw_type and vdw_type.upper() in VDW_CORRECTIONS:
        return VDW_CORRECTIONS[vdw_type.upper()]['incar_tags']
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
        print(f"             ENCUT: {info['recommended_settings']['encut']} eV")
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
        if 'nsw' in template:
            print(f"             Steps: {template['nsw']}")
        if 'ibrion' in template:
            print(f"             IBRION: {template['ibrion']}")
        print()

if __name__ == "__main__":
    # Demo the configuration module
    print("VASP Configuration Module Demo")
    print("=" * 40)
    
    list_available_minerals()
    list_available_functionals()
    list_simulation_templates()
    
    print("\nTesting mineral lookup:")
    test_minerals = ['sio2', 'calcite', 'unknown_mineral']
    for mineral in test_minerals:
        info = get_mineral_info(mineral)
        potcar_elements = get_potcar_elements(mineral)
        recommended_encut = get_recommended_encut(potcar_elements)
        print(f"{mineral}: {info['name']}")
        print(f"  Elements: {potcar_elements}")
        print(f"  Recommended ENCUT: {recommended_encut} eV")
        print()