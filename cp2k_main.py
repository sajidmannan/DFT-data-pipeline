#!/usr/bin/env python3
"""
CP2K Mineral Simulation Script

A clean, command-line driven script for running CP2K simulations on minerals.
Performs geometry optimization, electronic structure calculations, and molecular dynamics.

Usage:
    python cp2k_sim.py mineral_name --template opt --functional PBE --output results/
"""

import argparse
import sys
import shutil
from pathlib import Path

# Import our utility functions
from cp2k_utils import (
    # System utilities
    detect_cp2k, detect_cp2k_data, find_structure_file, parse_xyz_file, 
    parse_poscar_to_xyz, write_xyz_file,
    
    # File generators
    create_cp2k_input,
    
    # Job execution
    run_cp2k_calculation, backup_previous_calculation, copy_structure_files,
    write_run_info, validate_cp2k_input,
    
    # Reporting
    create_simulation_summary, print_completion_message
)

from cp2k_config import (
    get_mineral_info, get_recommended_cutoff, get_basis_sets, get_pseudopotentials,
    get_functional_info, get_vdw_info, get_simulation_template,
    DEFAULT_PARAMS, list_available_minerals, list_available_functionals,
    list_simulation_templates, list_available_basis_sets
)

# ================================================================
# Command Line Interface
# ================================================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Run CP2K simulations for minerals with geometry optimization, electronic structure, and MD capabilities.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cp2k_sim.py SiO2.xyz --template opt --output results/
  python cp2k_sim.py quartz --functional PBE0 --cutoff 600 --output high_accuracy/
  python cp2k_sim.py CaCO3 --template md --temp 500 --steps 5000
  python cp2k_sim.py mineral.xyz --template sp --basis TZVP-MOLOPT-GTH --functional B3LYP

Available templates: opt, cell_opt, sp, md, npt, vib, band

Supported minerals with auto parameters:
  SiO2, quartz, CaCO3, calcite, MgO, periclase, Mg2SiO4, forsterite, 
  Al2O3, corundum, Fe2O3, hematite
        """
    )
    
    # Required arguments
    parser.add_argument("mineral", 
                       help="Mineral name or path to structure file (XYZ, POSCAR)")
    
    # Calculation type
    parser.add_argument("--template", "-t",
                       type=str, default="opt",
                       choices=["opt", "cell_opt", "sp", "md", "npt", "vib", "band"],
                       help="Simulation template (default: opt)")
    
    # Output settings
    parser.add_argument("--output", "-o",
                       type=str, default="cp2k_results",
                       help="Output directory path (default: cp2k_results)")
    
    # Electronic structure settings
    parser.add_argument("--functional",
                       type=str, default="PBE",
                       choices=["PBE", "PBE0", "B3LYP", "BLYP", "BP86", "LDA", "SCAN"],
                       help="Exchange-correlation functional (default: PBE)")
    
    parser.add_argument("--cutoff",
                       type=int, default=None,
                       help="Plane wave cutoff in Ry (default: auto from elements)")
    
    parser.add_argument("--rel-cutoff",
                       type=int, default=50,
                       help="Relative cutoff for Gaussian grid in Ry (default: 50)")
    
    parser.add_argument("--scf-eps",
                       type=float, default=1E-6,
                       help="SCF convergence criterion (default: 1E-6)")
    
    parser.add_argument("--scf-guess",
                       type=str, default="ATOMIC",
                       choices=["ATOMIC", "RANDOM", "RESTART"],
                       help="SCF initial guess (default: ATOMIC)")
    
    # Basis sets and pseudopotentials
    parser.add_argument("--basis",
                       type=str, default=None,
                       help="Basis set type (default: DZVP-MOLOPT-SR-GTH)")
    
    parser.add_argument("--potential",
                       type=str, default=None,
                       help="Pseudopotential type (default: GTH-PBE)")
    
    # Geometry optimization
    parser.add_argument("--max-iter",
                       type=int, default=200,
                       help="Maximum optimization iterations (default: 200)")
    
    parser.add_argument("--optimizer",
                       type=str, default="BFGS",
                       choices=["BFGS", "CG", "LBFGS"],
                       help="Geometry optimizer (default: BFGS)")
    
    parser.add_argument("--cell-opt",
                       action="store_true",
                       help="Enable cell optimization")
    
    parser.add_argument("--pressure-tol",
                       type=float, default=100.0,
                       help="Pressure tolerance for cell opt in bar (default: 100)")
    
    # Molecular dynamics
    parser.add_argument("--temp",
                       type=float, default=300.0,
                       help="Temperature for MD in Kelvin (default: 300.0)")
    
    parser.add_argument("--pressure",
                       type=float, default=1.0,
                       help="External pressure in bar for NPT (default: 1.0)")
    
    parser.add_argument("--steps",
                       type=int, default=1000,
                       help="Number of MD steps (default: 1000)")
    
    parser.add_argument("--timestep",
                       type=float, default=0.5,
                       help="MD timestep in fs (default: 0.5)")
    
    parser.add_argument("--ensemble",
                       type=str, default="NVT",
                       choices=["NVT", "NPT_I", "NVE"],
                       help="MD ensemble (default: NVT)")
    
    parser.add_argument("--thermostat",
                       type=str, default="NOSE",
                       choices=["NOSE", "CSVR", "GLE"],
                       help="Thermostat type (default: NOSE)")
    
    # Van der Waals corrections
    parser.add_argument("--vdw",
                       type=str, default=None,
                       choices=["D2", "D3", "D3BJ", "TS"],
                       help="Van der Waals correction method")
    
    # Magnetic properties
    parser.add_argument("--uks",
                       action="store_true",
                       help="Enable unrestricted Kohn-Sham (spin polarization)")
    
    parser.add_argument("--multiplicity",
                       type=int, default=None,
                       help="Spin multiplicity (2S+1)")
    
    # Parallelization
    parser.add_argument("--mpi",
                       type=int, default=1,
                       help="Number of MPI ranks (default: 1)")
    
    # Job control
    parser.add_argument("--dry-run",
                       action="store_true",
                       help="Generate input files only, don't run CP2K")
    
    parser.add_argument("--continue-run",
                       action="store_true",
                       help="Continue previous calculation from restart file")
    
    # Utility options
    parser.add_argument("--list-minerals",
                       action="store_true",
                       help="List available minerals and exit")
    
    parser.add_argument("--list-functionals",
                       action="store_true",
                       help="List available functionals and exit")
    
    parser.add_argument("--list-templates",
                       action="store_true",
                       help="List available templates and exit")
    
    parser.add_argument("--list-basis",
                       action="store_true",
                       help="List available basis sets and exit")

    parser.add_argument("--cp2k-exe", default="cp2k.psmp",
                    help="Path or name of CP2K executable (default: cp2k.psmp)")
    
    return parser.parse_args()

# ================================================================
# Job Setup Functions
# ================================================================

def handle_utility_options(args):
    """Handle utility options like listing minerals, functionals, etc."""
    if args.list_minerals:
        list_available_minerals()
        sys.exit(0)
    
    if args.list_functionals:
        list_available_functionals()
        sys.exit(0)
    
    if args.list_templates:
        list_simulation_templates()
        sys.exit(0)
    
    if args.list_basis:
        list_available_basis_sets()
        sys.exit(0)

def setup_simulation(args):
    """Setup simulation environment and return configuration"""
    
    # Print header
    print("="*60)
    print("CP2K MINERAL SIMULATION PIPELINE")
    print("="*60)
    print(f"Mineral: {args.mineral}")
    print(f"Template: {args.template}")
    print(f"Functional: {args.functional}")
    print(f"Output directory: {args.output}")
    print("="*60)
    
    # Detect CP2K and data directory
    cp2k_cmd = detect_cp2k()
    try:
        data_dir = detect_cp2k_data()
        print(f"Using CP2K data directory: {data_dir}")
    except FileNotFoundError as e:
        print(f"Warning: {e}")
        data_dir = None
    
    print(f"Using CP2K executable: {cp2k_cmd}")
    
    # Find and parse structure file
    structure_file = find_structure_file(args.mineral)
    mineral_name = Path(args.mineral).stem
    
    print(f"Found structure file: {structure_file}")
    
    # Parse structure based on file type
    if structure_file.suffix.lower() in ['.xyz']:
        structure = parse_xyz_file(structure_file)
    elif structure_file.name in ['POSCAR', 'CONTCAR'] or 'POSCAR' in structure_file.name:
        structure = parse_poscar_to_xyz(structure_file)
    else:
        # Try to parse as XYZ first
        try:
            structure = parse_xyz_file(structure_file)
        except:
            try:
                structure = parse_poscar_to_xyz(structure_file)
            except:
                raise ValueError(f"Unable to parse structure file: {structure_file}")
    
    print(f"Structure: {structure['n_atoms']} atoms, "
          f"Elements: {' '.join(structure['unique_elements'])}")
    
    # Create output directory
    output_dir = Path(args.output) / f"{mineral_name}_{args.template}_{args.functional}"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Get mineral information and recommendations
    mineral_info = get_mineral_info(mineral_name)
    elements = structure['unique_elements']
    
    print(f"Elements: {elements}")
    
    # Setup simulation parameters
    sim_params = create_simulation_parameters(args, mineral_info, elements, structure)
    
    return {
        'cp2k_cmd': cp2k_cmd,
        'data_dir': data_dir,
        'structure_file': structure_file,
        'structure': structure,
        'mineral_name': mineral_name,
        'output_dir': output_dir,
        'mineral_info': mineral_info,
        'elements': elements,
        'sim_params': sim_params
    }

def create_simulation_parameters(args, mineral_info, elements, structure):
    """Create simulation parameters dictionary"""
    
    # Start with template parameters
    template_params = get_simulation_template(args.template)
    
    # Get recommended cutoff
    recommended_cutoff = get_recommended_cutoff(elements)
    cutoff = args.cutoff if args.cutoff else recommended_cutoff
    
    # Get functional info
    functional_info = get_functional_info(args.functional)
    
    # Get VdW info if specified
    vdw_info = get_vdw_info(args.vdw) if args.vdw else {}
    
    # Get basis sets and pseudopotentials
    mineral_basis_sets = get_basis_sets(args.mineral)
    mineral_pseudopotentials = get_pseudopotentials(args.mineral)
    
    # Override with user-specified or mineral-specific values
    basis_sets = {}
    pseudopotentials = {}
    
    for element in elements:
        # Basis set
        if args.basis:
            basis_sets[element] = args.basis
        elif element in mineral_basis_sets:
            basis_sets[element] = mineral_basis_sets[element]
        else:
            basis_sets[element] = DEFAULT_PARAMS['basis_set_type']
        
        # Pseudopotential
        if args.potential:
            pseudopotentials[element] = args.potential
        elif element in mineral_pseudopotentials:
            pseudopotentials[element] = mineral_pseudopotentials[element]
        else:
            functional_name = args.functional.upper()
            pseudopotentials[element] = f"GTH-{functional_name}"
    
    # Merge all parameters
    params = {}
    params.update(template_params)  # Template base
    
    # Override with user-specified parameters
    params.update({
        'project_name': f"{mineral_info['name'].replace(' ', '_').replace('(', '').replace(')', '')}_{args.template}",
        'cutoff': cutoff,
        'rel_cutoff': args.rel_cutoff,
        'scf_eps': args.scf_eps,
        'scf_guess': args.scf_guess,
        'xc_functional': args.functional,
        'basis_sets': basis_sets,
        'pseudopotentials': pseudopotentials,
    })
    
    # Add optional parameters
    if args.max_iter:
        params['geo_max_iter'] = args.max_iter
    
    if args.optimizer:
        params['optimizer'] = args.optimizer
    
    if args.cell_opt or args.template == 'cell_opt':
        params['cell_opt'] = True
        params['pressure_tolerance'] = args.pressure_tol
    
    if args.uks:
        params['uks'] = True
        if args.multiplicity:
            params['multiplicity'] = args.multiplicity
    
    if args.template in ['md', 'npt']:
        params['md_temperature'] = args.temp
        params['md_steps'] = args.steps
        params['md_timestep'] = args.timestep
        params['md_ensemble'] = args.ensemble
        params['thermostat'] = args.thermostat
        
        if args.template == 'npt' or args.ensemble == 'NPT_I':
            params['md_pressure'] = args.pressure
    
    if args.vdw:
        params['vdw_potential'] = args.vdw
    
    # Create simulation info dictionary
    sim_info = {
        'Mineral': mineral_info['name'],
        'Formula': mineral_info.get('formula', 'Unknown'),
        'Template': args.template,
        'Functional': args.functional,
        'VdW Correction': args.vdw or 'None',
        'Cutoff': f"{cutoff} Ry",
        'Rel Cutoff': f"{args.rel_cutoff} Ry",
        'Basis Set': args.basis or 'Auto',
        'Pseudopotential': args.potential or 'Auto',
        'SCF Convergence': f"{args.scf_eps:.0E}",
    }
    
    if args.template in ['md', 'npt']:
        sim_info['Temperature'] = f"{args.temp} K"
        sim_info['MD Steps'] = str(args.steps)
        sim_info['Timestep'] = f"{args.timestep} fs"
        if args.template == 'npt' or args.ensemble == 'NPT_I':
            sim_info['Pressure'] = f"{args.pressure} bar"
    
    if args.cell_opt or args.template == 'cell_opt':
        sim_info['Cell Optimization'] = 'Enabled'
        sim_info['Pressure Tolerance'] = f"{args.pressure_tol} bar"
    
    return {
        'cp2k_params': params,
        'sim_info': sim_info
    }

def create_job_directory(output_dir: Path, job_name: str, setup_config: dict, args):
    """Create and setup a job directory with all CP2K input files"""
    
    job_dir = output_dir / job_name
    job_dir.mkdir(exist_ok=True)
    
    # Backup previous calculation if continuing
    if args.continue_run:
        backup_previous_calculation(job_dir)
    
    # Copy structure file or use restart geometry if continuing
    if args.continue_run:
        restart_files = list(job_dir.glob("*-pos-1.xyz"))
        if restart_files:
            latest_restart = sorted(restart_files)[-1]
            restart_structure = parse_xyz_file(latest_restart)
            write_xyz_file(restart_structure, job_dir / "geometry.xyz")
            print("Using restart geometry as starting structure")
        else:
            write_xyz_file(setup_config['structure'], job_dir / "geometry.xyz")
    else:
        write_xyz_file(setup_config['structure'], job_dir / "geometry.xyz")
    
    # Generate CP2K input file
    cp2k_params = setup_config['sim_params']['cp2k_params']
    create_cp2k_input(cp2k_params, setup_config['structure'], job_dir / "input.inp")
    print(f"Generated input.inp with cutoff {cp2k_params['cutoff']} Ry")
    
    # Validate input file
    issues = validate_cp2k_input(job_dir / "input.inp")
    if issues:
        print("Warning: Input file validation issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("Input file validation: OK")
    
    return job_dir

# ================================================================
# Main Simulation Workflow
# ================================================================

def run_cp2k_job(setup_config: dict, args) -> tuple:
    """Run main CP2K calculation"""
    
    job_name = f"{args.template}_{args.functional}"
    if args.vdw:
        job_name += f"_{args.vdw}"
    
    print(f"\n{'='*50}")
    print(f"RUNNING {args.template.upper()} CALCULATION")
    print(f"{'='*50}")
    
    # Create job directory and input files
    job_dir = create_job_directory(setup_config['output_dir'], job_name, 
                                  setup_config, args)
    
    if args.dry_run:
        print(f"Dry run: Input files generated in {job_dir}")
        print("Skipping CP2K execution")
        return job_dir, {
            'job_name': args.template,
            'success': True,
            'converged': None,
            'runtime': None
        }
    
    # Run CP2K calculation
    run_info = run_cp2k_calculation(
        job_name=args.template,
        job_dir=job_dir,
        cp2k_cmd=setup_config['cp2k_cmd'],
        mpi_ranks=args.mpi
    )
    
    # Write run information
    write_run_info(job_dir, run_info, setup_config['mineral_name'], 
                  setup_config['sim_params']['sim_info'])
    
    if not run_info['success']:
        raise RuntimeError(f"{args.template} calculation failed")
    
    return job_dir, run_info

def main():
    """Main simulation workflow"""
    
    # Parse arguments
    args = parse_arguments()
    
    # Handle utility options
    handle_utility_options(args)
    
    try:
        # Setup simulation environment
        setup_config = setup_simulation(args)
        
        # Run main calculation
        job_dir, run_info = run_cp2k_job(setup_config, args)
        
        # Create summary
        job_dirs = {args.template: job_dir}
        run_infos = [run_info]
        
        if not args.dry_run:
            create_simulation_summary(
                output_dir=setup_config['output_dir'],
                mineral_name=setup_config['mineral_name'],
                simulation_params=setup_config['sim_params']['sim_info'],
                job_dirs=job_dirs,
                run_infos=run_infos
            )
        
        # Print completion message
        if args.dry_run:
            print(f"\n{'='*60}")
            print("📁 INPUT FILES GENERATED SUCCESSFULLY!")
            print(f"{'='*60}")
            print(f"Files created in: {setup_config['output_dir']}")
        else:
            print_completion_message(setup_config['output_dir'], run_infos)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Simulation interrupted by user")
        sys.exit(1)
        
    except Exception as e:
        print(f"\n❌ Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()