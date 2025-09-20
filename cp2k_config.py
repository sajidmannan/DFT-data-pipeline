#!/usr/bin/env python3
"""
CP2K Configuration Module

Contains default parameters and mineral-specific configurations
for CP2K simulations.
"""

# ================================================================
# Default Simulation Parameters
# ================================================================

DEFAULT_PARAMS = {
    'temperature': 300.0,        # Kelvin
    'pressure': 1.0,             # bar
    'functional': 'PBE_D3',      # Exchange-correlation functional
    'cutoff': 600,               # Plane wave cutoff (Ry)
    'rel_cutoff': 40,            # Relative cutoff (Ry)
    'nvt_steps': 5000,           # NVT equilibration steps
    'npt_steps': 10000,          # NPT production steps
    'timecon': 100.0,            # Thermostat time constant (fs)
    'timestep': 0.5,             # MD timestep (fs)
    'scf_convergence': 1.0E-6,   # SCF convergence criterion
    'max_scf': 50,               # Maximum SCF iterations
    'barostat_timecon': 500,     # Barostat time constant (fs)
}

# ================================================================
# Mineral-Specific Configurations
# ================================================================

MINERAL_DATABASE = {
    'sio2': {
        'name': 'Silicon Dioxide (Quartz)',
        'cell_params': [4.916, 4.916, 5.405],  # a, b, c in Angstroms
        'space_group': 'P3121',
        'density': 2.65,  # g/cm³
        'common_elements': ['Si', 'O'],
        'recommended_cutoff': 600,
        'notes': 'Alpha-quartz structure'
    },
    
    'quartz': {
        'name': 'Quartz',
        'cell_params': [4.916, 4.916, 5.405],
        'space_group': 'P3121', 
        'density': 2.65,
        'common_elements': ['Si', 'O'],
        'recommended_cutoff': 600,
        'notes': 'Same as SiO2'
    },
    
    'caco3': {
        'name': 'Calcium Carbonate (Calcite)',
        'cell_params': [4.99, 4.99, 17.06],
        'space_group': 'R-3c',
        'density': 2.71,
        'common_elements': ['Ca', 'C', 'O'],
        'recommended_cutoff': 650,
        'notes': 'Calcite structure'
    },
    
    'calcite': {
        'name': 'Calcite',
        'cell_params': [4.99, 4.99, 17.06],
        'space_group': 'R-3c',
        'density': 2.71,
        'common_elements': ['Ca', 'C', 'O'],
        'recommended_cutoff': 650,
        'notes': 'CaCO3 polymorph'
    },
    
    'mg2sio4': {
        'name': 'Magnesium Silicate (Forsterite)',
        'cell_params': [4.754, 10.207, 5.980],
        'space_group': 'Pbnm',
        'density': 3.22,
        'common_elements': ['Mg', 'Si', 'O'],
        'recommended_cutoff': 700,
        'notes': 'Olivine group mineral'
    },
    
    'forsterite': {
        'name': 'Forsterite',
        'cell_params': [4.754, 10.207, 5.980],
        'space_group': 'Pbnm',
        'density': 3.22,
        'common_elements': ['Mg', 'Si', 'O'],
        'recommended_cutoff': 700,
        'notes': 'Mg2SiO4 olivine endmember'
    },
    
    'al2o3': {
        'name': 'Aluminum Oxide (Corundum)',
        'cell_params': [4.759, 4.759, 12.991],
        'space_group': 'R-3c',
        'density': 3.95,
        'common_elements': ['Al', 'O'],
        'recommended_cutoff': 650,
        'notes': 'Corundum structure'
    },
    
    'corundum': {
        'name': 'Corundum',
        'cell_params': [4.759, 4.759, 12.991],
        'space_group': 'R-3c',
        'density': 3.95,
        'common_elements': ['Al', 'O'],
        'recommended_cutoff': 650,
        'notes': 'Al2O3 structure'
    },
    
    'fe2o3': {
        'name': 'Iron Oxide (Hematite)',
        'cell_params': [5.034, 5.034, 13.747],
        'space_group': 'R-3c',
        'density': 5.26,
        'common_elements': ['Fe', 'O'],
        'recommended_cutoff': 700,
        'notes': 'Hematite structure'
    },
    
    'hematite': {
        'name': 'Hematite',
        'cell_params': [5.034, 5.034, 13.747],
        'space_group': 'R-3c',
        'density': 5.26,
        'common_elements': ['Fe', 'O'],
        'recommended_cutoff': 700,
        'notes': 'Fe2O3 structure'
    },
    
    'mgo': {
        'name': 'Magnesium Oxide (Periclase)',
        'cell_params': [4.212, 4.212, 4.212],
        'space_group': 'Fm-3m',
        'density': 3.58,
        'common_elements': ['Mg', 'O'],
        'recommended_cutoff': 600,
        'notes': 'Cubic halite structure'
    },
    
    'periclase': {
        'name': 'Periclase',
        'cell_params': [4.212, 4.212, 4.212],
        'space_group': 'Fm-3m',
        'density': 3.58,
        'common_elements': ['Mg', 'O'],
        'recommended_cutoff': 600,
        'notes': 'MgO cubic structure'
    }
}

# ================================================================
# Element Properties for Basis Sets
# ================================================================

ELEMENT_PROPERTIES = {
    'H':  {'valence': 1,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 1.008},
    'He': {'valence': 2,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 4.003},
    'Li': {'valence': 3,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 6.941},
    'Be': {'valence': 4,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 9.012},
    'B':  {'valence': 3,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 10.811},
    'C':  {'valence': 4,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 12.011},
    'N':  {'valence': 5,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 14.007},
    'O':  {'valence': 6,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 15.999},
    'F':  {'valence': 7,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 18.998},
    'Ne': {'valence': 8,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 20.180},
    'Na': {'valence': 9,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 22.990},
    'Mg': {'valence': 10, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 24.305},
    'Al': {'valence': 3,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 26.982},
    'Si': {'valence': 4,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 28.086},
    'P':  {'valence': 5,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 30.974},
    'S':  {'valence': 6,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 32.065},
    'Cl': {'valence': 7,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 35.453},
    'Ar': {'valence': 8,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 39.948},
    'K':  {'valence': 9,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 39.098},
    'Ca': {'valence': 10, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 40.078},
    'Sc': {'valence': 11, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 44.956},
    'Ti': {'valence': 12, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 47.867},
    'V':  {'valence': 13, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 50.942},
    'Cr': {'valence': 6,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 51.996},
    'Mn': {'valence': 7,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 54.938},
    'Fe': {'valence': 8,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 55.845},
    'Co': {'valence': 9,  'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 58.933},
    'Ni': {'valence': 18, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 58.693},
    'Cu': {'valence': 11, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 63.546},
    'Zn': {'valence': 12, 'basis': 'DZVP-MOLOPT-SR-GTH', 'mass': 65.409},
}

# ================================================================
# Supported Functionals
# ================================================================

FUNCTIONALS = {
    'PBE': {
        'type': 'GGA',
        'dispersion': False,
        'recommended_for': ['general_purpose', 'solids'],
        'notes': 'Standard PBE functional'
    },
    
    'PBE_D3': {
        'type': 'GGA+D3',
        'dispersion': True,
        'recommended_for': ['general_purpose', 'solids', 'layered_materials'],
        'notes': 'PBE with Grimme D3 dispersion correction'
    },
    
    'BLYP': {
        'type': 'GGA',
        'dispersion': False,
        'recommended_for': ['molecules', 'soft_materials'],
        'notes': 'Becke-Lee-Yang-Parr functional'
    },
    
    'BLYP_D3': {
        'type': 'GGA+D3',
        'dispersion': True,
        'recommended_for': ['molecules', 'soft_materials', 'layered_materials'],
        'notes': 'BLYP with D3 dispersion correction'
    },
    
    'PBE0': {
        'type': 'Hybrid',
        'dispersion': False,
        'recommended_for': ['accurate_energetics', 'band_gaps'],
        'notes': 'Hybrid functional, computationally expensive'
    }
}

# ================================================================
# Simulation Templates
# ================================================================

SIMULATION_TEMPLATES = {
    'quick_test': {
        'nvt_steps': 1000,
        'npt_steps': 2000,
        'cutoff': 400,
        'rel_cutoff': 30,
        'description': 'Quick test run with minimal steps'
    },
    
    'standard': {
        'nvt_steps': 5000,
        'npt_steps': 10000,
        'cutoff': 600,
        'rel_cutoff': 40,
        'description': 'Standard production run'
    },
    
    'high_quality': {
        'nvt_steps': 10000,
        'npt_steps': 50000,
        'cutoff': 800,
        'rel_cutoff': 50,
        'description': 'High-quality run with extensive sampling'
    },
    
    'convergence_test': {
        'nvt_steps': 2000,
        'npt_steps': 5000,
        'cutoff': [400, 600, 800],  # Multiple cutoffs for testing
        'rel_cutoff': 40,
        'description': 'Test multiple cutoffs for convergence'
    }
}

# ================================================================
# Utility Functions
# ================================================================

def get_mineral_info(mineral_name: str) -> dict:
    """Get information about a specific mineral"""
    key = mineral_name.lower().replace('.xyz', '')
    
    if key in MINERAL_DATABASE:
        return MINERAL_DATABASE[key]
    else:
        # Return generic info for unknown minerals
        return {
            'name': f'Unknown mineral: {mineral_name}',
            'cell_params': [5.0, 5.0, 5.0],
            'space_group': 'Unknown',
            'density': None,
            'common_elements': [],
            'recommended_cutoff': 600,
            'notes': 'Generic cubic cell used'
        }

def get_recommended_cutoff(mineral_name: str, functional: str = 'PBE_D3') -> int:
    """Get recommended cutoff for a mineral"""
    mineral_info = get_mineral_info(mineral_name)
    base_cutoff = mineral_info['recommended_cutoff']
    
    # Adjust for functional type
    if 'PBE0' in functional:
        return int(base_cutoff * 1.2)  # Hybrids need higher cutoff
    elif 'BLYP' in functional:
        return int(base_cutoff * 0.9)   # BLYP can use slightly lower
    else:
        return base_cutoff

def list_available_minerals():
    """List all minerals in the database"""
    print("Available minerals with predefined parameters:")
    print("=" * 50)
    
    for key, info in MINERAL_DATABASE.items():
        cell = info['cell_params']
        print(f"{key:12s} - {info['name']}")
        print(f"             Cell: {cell[0]:.3f} {cell[1]:.3f} {cell[2]:.3f} Å")
        print(f"             Elements: {', '.join(info['common_elements'])}")
        print()

def validate_functional(functional: str) -> bool:
    """Validate if functional is supported"""
    return functional in FUNCTIONALS

def get_simulation_template(template_name: str) -> dict:
    """Get simulation template parameters"""
    if template_name in SIMULATION_TEMPLATES:
        return SIMULATION_TEMPLATES[template_name].copy()
    else:
        raise ValueError(f"Unknown template: {template_name}. "
                        f"Available: {list(SIMULATION_TEMPLATES.keys())}")

if __name__ == "__main__":
    # Demo the configuration module
    print("CP2K Configuration Module Demo")
    print("=" * 40)
    list_available_minerals()
    
    print("\nTesting mineral lookup:")
    test_minerals = ['sio2', 'calcite', 'unknown_mineral']
    for mineral in test_minerals:
        info = get_mineral_info(mineral)
        print(f"{mineral}: {info['name']}, cutoff: {info['recommended_cutoff']}")