#!/usr/bin/env python3
"""
VASP Mineral Simulation Script

A clean, command-line driven script for running VASP simulations on minerals.
Performs structural optimization, electronic structure calculations, and molecular dynamics.

Usage:
    python vasp_sim.py mineral_name --template relax --functional PBE --output results/
"""

import argparse
import sys
import shutil
from pathlib import Path

# Import our utility functions
from vasp_utils import (
    # System utilities
    detect_vasp, detect_vasp_potentials, find_poscar_file, parse_poscar,
    
    # File generators
    write_incar, write_kpoints, create_potcar, create_kpoints_from_poscar,
    
    # Job execution
    run_vasp_calculation, backup_previous_calculation, copy_structure_files,
    write_run_info,
    
    # Reporting
    create_simulation_summary, print_completion_message
)

from vasp_config import (
    get_mineral_info, get_recommended_encut, get_potcar_elements,
    get_functional_tags, get_vdw_tags, get_simulation_template,
    DEFAULT_PARAMS, list_available_minerals, list_available_functionals,
    list_simulation_templates
)

# ================================================================
# Command Line Interface
# ================================================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Run VASP simulations for minerals with structural optimization, electronic structure, and MD capabilities.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python vasp_sim.py SiO2.POSCAR --template relax --output results/
  python vasp_sim.py quartz --functional PBEsol --encut 600 --output high_accuracy/
  python vasp_sim.py CaCO3 --template md_nvt --temp 500 --steps 10000
  python vasp_sim.py mineral.POSCAR --template scf --kpts 12,12,8 --functional HSE06

Available templates: relax, scf, dos, band, md_nvt, md_npt, phonon, elastic

Supported minerals with auto parameters:
  SiO2, quartz, CaCO3, calcite, MgO, periclase, Mg2SiO4, forsterite, 
  Al2O3, corundum, Fe2O3, hematite
        """
    )
    
    # Required arguments
    parser.add_argument("mineral", 
                       help="Mineral name or path to POSCAR file")
    
    # Calculation type
    parser.add_argument("--template", "-t",
                       type=str, default="relax",
                       choices=["relax", "scf", "dos", "band", "md_nvt", "md_npt", 
                               "phonon", "elastic"],
                       help="Simulation template (default: relax)")
    
    # Output settings
    parser.add_argument("--output", "-o",
                       type=str, default="vasp_results",
                       help="Output directory path (default: vasp_results)")
    
    # Electronic structure settings
    parser.add_argument("--functional",
                       type=str, default="PBE",
                       choices=["PBE", "PBEsol", "LDA", "SCAN", "HSE06", "PBE0"],
                       help="Exchange-correlation functional (default: PBE)")
    
    parser.add_argument("--encut",
                       type=int, default=None,
                       help="Plane wave cutoff in eV (default: auto from POTCAR)")
    
    parser.add_argument("--prec",
                       type=str, default="Accurate",
                       choices=["Low", "Medium", "High", "Normal", "Accurate"],
                       help="Precision level (default: Accurate)")
    
    parser.add_argument("--algo",
                       type=str, default="Fast",
                       choices=["Normal", "VeryFast", "Fast", "Conjugate", "All"],
                       help="Electronic minimization algorithm (default: Fast)")
    
    # K-points
    parser.add_argument("--kpts",
                       type=str, default="auto",
                       help="K-points grid 'nx,ny,nz' or 'auto' or density (default: auto)")
    
    parser.add_argument("--kpts-shift",
                       type=str, default="0,0,0",
                       help="K-points shift 'x,y,z' (default: 0,0,0)")
    
    # Structural optimization
    parser.add_argument("--ediffg",
                       type=float, default=-0.01,
                       help="Ionic convergence criterion eV/Å (default: -0.01)")
    
    parser.add_argument("--isif",
                       type=int, default=None,
                       help="Stress/relaxation flag (default: from template)")
    
    parser.add_argument("--steps",
                       type=int, default=None,
                       help="Maximum number of ionic steps (default: from template)")
    
    # Molecular dynamics
    parser.add_argument("--temp",
                       type=float, default=300.0,
                       help="Temperature for MD in Kelvin (default: 300.0)")
    
    parser.add_argument("--pressure",
                       type=float, default=0.0,
                       help="External pressure in kBar (default: 0.0)")
    
    # Van der Waals corrections
    parser.add_argument("--vdw",
                       type=str, default=None,
                       choices=["D2", "D3", "D3BJ", "TS", "BEEF"],
                       help="Van der Waals correction method")
    
    # Magnetic properties
    parser.add_argument("--spin",
                       action="store_true",
                       help="Enable spin polarization")
    
    parser.add_argument("--magmom",
                       type=str, default=None,
                       help="Initial magnetic moments 'atom1,atom2,...'")
    
    # Parallelization
    parser.add_argument("--mpi",
                       type=int, default=1,
                       help="Number of MPI ranks (default: 1)")
    
    parser.add_argument("--ncore",
                       type=int, default=None,
                       help="NCORE parameter for parallelization")
    
    # Job control
    parser.add_argument("--dry-run",
                       action="store_true",
                       help="Generate input files only, don't run VASP")
    
    parser.add_argument("--continue-run",
                       action="store_true",
                       help="Continue previous calculation from CONTCAR")
    
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

def setup_simulation(args):
    """Setup simulation environment and return configuration"""
    
    # Print header
    print("="*60)
    print("VASP MINERAL SIMULATION PIPELINE")
    print("="*60)
    print(f"Mineral: {args.mineral}")
    print(f"Template: {args.template}")
    print(f"Functional: {args.functional}")
    print(f"Output directory: {args.output}")
    print("="*60)
    
    # Detect VASP and potentials
    vasp_cmd = detect_vasp()
    try:
        potcar_dir = detect_vasp_potentials()
        print(f"Using POTCAR directory: {potcar_dir}")
    except FileNotFoundError as e:
        print(f"Warning: {e}")
        potcar_dir = None
    
    print(f"Using VASP executable: {vasp_cmd}")
    
    # Find and parse POSCAR file
    poscar_file = find_poscar_file(args.mineral)
    structure = parse_poscar(poscar_file)
    mineral_name = Path(args.mineral).stem
    
    print(f"Found POSCAR file: {poscar_file}")
    print(f"Structure: {len(structure['positions'])} atoms, "
          f"Elements: {' '.join(structure['element_types'])}")
    
    # Create output directory
    output_dir = Path(args.output) / f"{mineral_name}_{args.template}_{args.functional}"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Get mineral information and recommendations
    mineral_info = get_mineral_info(mineral_name)
    potcar_elements = structure['element_types']  # Use elements from POSCAR
    
    if not potcar_elements:
        potcar_elements = get_potcar_elements(mineral_name)
    
    print(f"POTCAR elements: {potcar_elements}")
    
    # Setup simulation parameters
    sim_params = create_simulation_parameters(args, mineral_info, potcar_elements, structure)
    
    return {
        'vasp_cmd': vasp_cmd,
        'potcar_dir': potcar_dir,
        'poscar_file': poscar_file,
        'structure': structure,
        'mineral_name': mineral_name,
        'output_dir': output_dir,
        'mineral_info': mineral_info,
        'potcar_elements': potcar_elements,
        'sim_params': sim_params
    }

def create_simulation_parameters(args, mineral_info, potcar_elements, structure):
    """Create simulation parameters dictionary"""
    
    # Start with template parameters
    template_params = get_simulation_template(args.template)
    
    # Get recommended ENCUT
    recommended_encut = get_recommended_encut(potcar_elements)
    encut = args.encut if args.encut else recommended_encut
    
    # Get functional tags
    functional_tags = get_functional_tags(args.functional)
    
    # Get VdW tags if specified
    vdw_tags = get_vdw_tags(args.vdw) if args.vdw else {}
    
    # Parse k-points
    kpts_grid, kpts_shift = parse_kpoints_input(args, structure)
    
    # Parse magnetic moments
    magmom = parse_magnetic_moments(args.magmom, structure) if args.magmom else None
    
    # Merge all parameters
    params = {}
    params.update(template_params)  # Template base
    params.update(functional_tags)  # XC functional
    params.update(vdw_tags)         # VdW corrections
    
    # Override with user-specified parameters
    params.update({
        'SYSTEM': f"{mineral_info['name']} - {args.template}",
        'ENCUT': encut,
        'PREC': args.prec,
        'ALGO': args.algo,
        'EDIFFG': args.ediffg,
    })
    
    # Add optional parameters
    if args.isif is not None:
        params['ISIF'] = args.isif
    
    if args.steps is not None:
        params['NSW'] = args.steps
    
    if args.spin or magmom:
        params['ISPIN'] = 2
        if magmom:
            params['MAGMOM'] = magmom
    
    if args.template in ['md_nvt', 'md_npt']:
        params['TEBEG'] = args.temp
        params['TEEND'] = args.temp
        
        if args.template == 'md_npt':
            params['PSTRESS'] = args.pressure
    
    if args.ncore:
        params['NCORE'] = args.ncore
    
    # Create simulation info dictionary
    sim_info = {
        'Mineral': mineral_info['name'],
        'Formula': mineral_info.get('formula', 'Unknown'),
        'Template': args.template,
        'Functional': args.functional,
        'VdW Correction': args.vdw or 'None',
        'ENCUT': f"{encut} eV",
        'K-points': f"{kpts_grid[0]}×{kpts_grid[1]}×{kpts_grid[2]}",
        'Precision': args.prec,
        'Algorithm': args.algo,
    }
    
    if args.template in ['md_nvt', 'md_npt']:
        sim_info['Temperature'] = f"{args.temp} K"
        if args.template == 'md_npt':
            sim_info['Pressure'] = f"{args.pressure} kBar"
    
    return {
        'incar_params': params,
        'kpts_grid': kpts_grid,
        'kpts_shift': kpts_shift,
        'sim_info': sim_info
    }

def parse_kpoints_input(args, structure):
    """Parse k-points input from command line"""
    
    if args.kpts == "auto":
        # Use mineral-specific recommendation or auto-generate
        mineral_info = get_mineral_info(Path(args.mineral).stem)
        if 'kpts_grid' in mineral_info.get('recommended_settings', {}):
            kpts_grid = mineral_info['recommended_settings']['kpts_grid']
        else:
            # Auto-generate from POSCAR
            kpts_grid = create_kpoints_from_poscar(Path(args.mineral))
    
    elif args.kpts.replace('.', '').isdigit():
        # K-points density specified
        density = float(args.kpts)
        kpts_grid = create_kpoints_from_poscar(Path(args.mineral), density)
    
    else:
        # Explicit grid specified
        try:
            kpts_grid = [int(x) for x in args.kpts.split(',')]
            if len(kpts_grid) != 3:
                raise ValueError
        except ValueError:
            raise ValueError(f"Invalid k-points format: {args.kpts}. Use 'nx,ny,nz' or 'auto'")
    
    # Parse k-points shift
    try:
        kpts_shift = [float(x) for x in args.kpts_shift.split(',')]
        if len(kpts_shift) != 3:
            raise ValueError
    except ValueError:
        raise ValueError(f"Invalid k-points shift format: {args.kpts_shift}. Use 'x,y,z'")
    
    return kpts_grid, kpts_shift

def parse_magnetic_moments(magmom_str, structure):
    """Parse magnetic moments from command line input"""
    
    try:
        magmom_list = [float(x) for x in magmom_str.split(',')]
        total_atoms = sum(structure['element_counts'])
        
        if len(magmom_list) == 1:
            # Same moment for all atoms
            return [magmom_list[0]] * total_atoms
        elif len(magmom_list) == len(structure['element_types']):
            # One moment per element type
            result = []
            for i, count in enumerate(structure['element_counts']):
                result.extend([magmom_list[i]] * count)
            return result
        elif len(magmom_list) == total_atoms:
            # Individual moments for each atom
            return magmom_list
        else:
            raise ValueError(f"Magnetic moments count mismatch. Provided {len(magmom_list)}, "
                           f"expected 1, {len(structure['element_types'])}, or {total_atoms}")
    
    except ValueError as e:
        raise ValueError(f"Invalid magnetic moments format: {magmom_str}. {str(e)}")

def create_job_directory(output_dir: Path, job_name: str, setup_config: dict, args):
    """Create and setup a job directory with all VASP input files"""
    
    job_dir = output_dir / job_name
    job_dir.mkdir(exist_ok=True)
    
    # Backup previous calculation if continuing
    if args.continue_run:
        backup_previous_calculation(job_dir)
    
    # Copy POSCAR (or use CONTCAR if continuing)
    if args.continue_run and (job_dir / "CONTCAR").exists():
        shutil.copy2(job_dir / "CONTCAR", job_dir / "POSCAR")
        print("Using CONTCAR as starting structure")
    else:
        shutil.copy2(setup_config['poscar_file'], job_dir / "POSCAR")
    
    # Write INCAR
    incar_params = setup_config['sim_params']['incar_params']
    write_incar(incar_params, job_dir / "INCAR", args.template)
    print(f"Generated INCAR with {len(incar_params)} parameters")
    
    # Write KPOINTS
    kpts_grid = setup_config['sim_params']['kpts_grid']
    kpts_shift = setup_config['sim_params']['kpts_shift']
    write_kpoints(kpts_grid, job_dir / "KPOINTS", kpts_shift)
    print(f"Generated KPOINTS: {kpts_grid[0]}×{kpts_grid[1]}×{kpts_grid[2]}")
    
    # Create POTCAR
    if setup_config['potcar_dir']:
        try:
            create_potcar(setup_config['potcar_elements'], 
                         job_dir / "POTCAR", 
                         setup_config['potcar_dir'])
            print(f"Generated POTCAR for elements: {setup_config['potcar_elements']}")
        except FileNotFoundError as e:
            print(f"Warning: Could not create POTCAR: {e}")
            print("You may need to create POTCAR manually")
    else:
        print("Warning: POTCAR directory not found. Create POTCAR manually.")
    
    return job_dir

# ================================================================
# Main Simulation Workflow
# ================================================================

def run_vasp_job(setup_config: dict, args) -> tuple:
    """Run main VASP calculation"""
    
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
        print("Skipping VASP execution")
        return job_dir, {
            'job_name': args.template,
            'success': True,
            'converged': None,
            'runtime': None
        }
    
    # Run VASP calculation
    run_info = run_vasp_calculation(
        job_name=args.template,
        job_dir=job_dir,
        vasp_cmd=setup_config['vasp_cmd'],
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
        job_dir, run_info = run_vasp_job(setup_config, args)
        
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