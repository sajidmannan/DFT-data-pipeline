#!/usr/bin/env python3
"""
CP2K Mineral Simulation Script

A clean, command-line driven script for running CP2K simulations on minerals.
Performs cell optimization, NVT equilibration, and NPT production runs.

Usage:
    python cp2k_sim.py mineral_name --temp 300 --pressure 1.0 --output results/
"""

import argparse
import sys
import shutil
from pathlib import Path

# Import our utility functions
from cp2k_utils import (
    # System utilities
    detect_cp2k, find_xyz_file, clean_xyz_file, parse_cell_parameters,
    
    # CP2K generators
    generate_motion_block, create_cp2k_input,
    
    # Job execution
    run_cp2k_calculation, write_run_info, find_optimized_structure,
    
    # Reporting
    create_simulation_summary, print_completion_message
)

# ================================================================
# Command Line Interface
# ================================================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Run CP2K simulations for minerals with cell optimization, NVT equilibration, and NPT production runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cp2k_sim.py SiO2 --temp 300 --pressure 1.0 --output results/
  python cp2k_sim.py quartz.xyz --temp 500 --pressure 2.5 --output /path/to/results
  python cp2k_sim.py CaCO3 --temp 300 --pressure 1.0 --functional PBE_D3 --cutoff 800

Supported minerals with auto cell parameters:
  SiO2, quartz, CaCO3, calcite, Mg2SiO4, forsterite, Al2O3, corundum, Fe2O3, hematite
        """
    )
    
    # Required arguments
    parser.add_argument("mineral", 
                       help="Mineral name or path to XYZ file")
    
    # Thermodynamic conditions
    parser.add_argument("--temp", "-T", 
                       type=float, default=300.0,
                       help="Temperature in Kelvin (default: 300.0)")
    
    parser.add_argument("--pressure", "-P", 
                       type=float, default=1.0,
                       help="Pressure in bar (default: 1.0)")
    
    # Output settings
    parser.add_argument("--output", "-o", 
                       type=str, default="cp2k_results",
                       help="Output directory path (default: cp2k_results)")
    
    # DFT settings
    parser.add_argument("--functional", 
                       type=str, default="PBE_D3",
                       choices=["PBE", "PBE_D3", "BLYP", "BLYP_D3", "PBE0"],
                       help="Exchange-correlation functional (default: PBE_D3)")
    
    parser.add_argument("--cutoff", 
                       type=int, default=600,
                       help="Plane wave cutoff in Ry (default: 600)")
    
    parser.add_argument("--rel-cutoff", 
                       type=int, default=40,
                       help="Relative cutoff in Ry (default: 40)")
    
    # MD settings
    parser.add_argument("--nvt-steps", 
                       type=int, default=5000,
                       help="Number of NVT equilibration steps (default: 5000)")
    
    parser.add_argument("--npt-steps", 
                       type=int, default=10000,
                       help="Number of NPT production steps (default: 10000)")
    
    parser.add_argument("--timecon", 
                       type=float, default=100.0,
                       help="Thermostat time constant in fs (default: 100.0)")
    
    # Cell parameters
    parser.add_argument("--cell", 
                       type=str, default="auto",
                       help="Cell parameters 'a b c' or 'auto' for defaults (default: auto)")
    
    # Job control
    parser.add_argument("--skip-cellopt", 
                       action="store_true",
                       help="Skip cell optimization (use initial structure)")
    
    parser.add_argument("--skip-nvt", 
                       action="store_true",
                       help="Skip NVT equilibration")
    
    parser.add_argument("--only-cellopt", 
                       action="store_true",
                       help="Only run cell optimization")
    
    return parser.parse_args()

# ================================================================
# Job Setup Functions
# ================================================================

def setup_simulation(args):
    """Setup simulation environment and return configuration"""
    
    # Print header
    print("="*60)
    print("CP2K MINERAL SIMULATION PIPELINE")
    print("="*60)
    print(f"Mineral: {args.mineral}")
    print(f"Temperature: {args.temp} K")
    print(f"Pressure: {args.pressure} bar")
    print(f"Output directory: {args.output}")
    print(f"Functional: {args.functional}")
    print("="*60)
    
    # Detect CP2K and find XYZ file
    cp2k_cmd = detect_cp2k()
    original_xyz = find_xyz_file(args.mineral)
    mineral_name = Path(original_xyz).stem
    
    print(f"Using CP2K executable: {cp2k_cmd}")
    print(f"Found XYZ file: {original_xyz}")
    
    # Create output directory structure
    output_dir = Path(args.output) / f"{mineral_name}_T{args.temp}_P{args.pressure}"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Prepare coordinate file
    cleaned_xyz = output_dir / "coords.xyz"
    clean_xyz_file(original_xyz, cleaned_xyz)
    shutil.copy2(original_xyz, output_dir / f"original_{original_xyz.name}")
    
    print(f"Prepared coordinate file: {cleaned_xyz}")
    
    # Parse cell parameters
    cell_params = parse_cell_parameters(args.cell, mineral_name)
    
    # Create simulation parameters dictionary
    sim_params = {
        'Mineral': mineral_name,
        'Temperature': f"{args.temp} K",
        'Pressure': f"{args.pressure} bar",
        'Functional': args.functional,
        'Cutoff': f"{args.cutoff} Ry",
        'Relative Cutoff': f"{args.rel_cutoff} Ry",
        'NVT Steps': args.nvt_steps,
        'NPT Steps': args.npt_steps,
        'Thermostat Time Constant': f"{args.timecon} fs",
        'Cell Parameters': f"{cell_params[0]:.3f} {cell_params[1]:.3f} {cell_params[2]:.3f}"
    }
    
    return {
        'cp2k_cmd': cp2k_cmd,
        'original_xyz': original_xyz,
        'mineral_name': mineral_name,
        'output_dir': output_dir,
        'cleaned_xyz': cleaned_xyz,
        'cell_params': cell_params,
        'sim_params': sim_params
    }

def create_job_directory(output_dir: Path, job_name: str, coords_xyz: Path, 
                        mineral_name: str, run_type: str, args, cell_params, 
                        motion_section: str):
    """Create and setup a job directory"""
    
    job_dir = output_dir / job_name
    job_dir.mkdir(exist_ok=True, parents=True)
    
    # Copy coordinate file
    shutil.copy2(coords_xyz, job_dir / "coords.xyz")
    
    # Create CP2K input file
    create_cp2k_input(
        run_type=run_type,
        mineral_name=mineral_name,
        dir_path=job_dir,
        functional=args.functional,
        cutoff=args.cutoff,
        rel_cutoff=args.rel_cutoff,
        cell_params=cell_params,
        motion_section=motion_section
    )
    
    return job_dir

# ================================================================
# Main Simulation Workflow
# ================================================================

def run_cell_optimization(setup_config, args):
    """Run cell optimization step"""
    
    print("\n" + "="*50)
    print("STEP 1: CELL OPTIMIZATION")
    print("="*50)
    
    motion_section = generate_motion_block("CELL_OPT")
    job_name = f"01_cellopt_{args.functional}_T{args.temp}_P{args.pressure}"
    
    job_dir = create_job_directory(
        output_dir=setup_config['output_dir'],
        job_name=job_name,
        coords_xyz=setup_config['cleaned_xyz'],
        mineral_name=setup_config['mineral_name'],
        run_type="CELL_OPT",
        args=args,
        cell_params=setup_config['cell_params'],
        motion_section=motion_section
    )
    
    # Run calculation
    run_info = run_cp2k_calculation("Cell Optimization", job_dir, setup_config['cp2k_cmd'])
    
    # Write run info
    write_run_info(job_dir, run_info, setup_config['mineral_name'], 
                  args.temp, args.pressure, args.functional, args.cutoff)
    
    if not run_info['success']:
        raise RuntimeError("Cell optimization failed")
    
    return job_dir, run_info

def run_nvt_equilibration(setup_config, args, input_structure):
    """Run NVT equilibration step"""
    
    print("\n" + "="*50)
    print("STEP 2: NVT EQUILIBRATION")
    print("="*50)
    
    motion_section = generate_motion_block(
        "NVT",
        temperature=args.temp,
        steps=args.nvt_steps,
        timecon=args.timecon
    )
    
    job_name = f"02_nvt_{args.functional}_T{args.temp}_P{args.pressure}"
    
    job_dir = create_job_directory(
        output_dir=setup_config['output_dir'],
        job_name=job_name,
        coords_xyz=input_structure,
        mineral_name=setup_config['mineral_name'],
        run_type="MD",
        args=args,
        cell_params=setup_config['cell_params'],
        motion_section=motion_section
    )
    
    # Run calculation
    run_info = run_cp2k_calculation("NVT Equilibration", job_dir, setup_config['cp2k_cmd'])
    
    # Write run info
    write_run_info(job_dir, run_info, setup_config['mineral_name'], 
                  args.temp, args.pressure, args.functional, args.cutoff)
    
    if not run_info['success']:
        raise RuntimeError("NVT equilibration failed")
    
    return job_dir, run_info

def run_npt_production(setup_config, args, input_structure):
    """Run NPT production step"""
    
    print("\n" + "="*50)
    print("STEP 3: NPT PRODUCTION")
    print("="*50)
    
    motion_section = generate_motion_block(
        "NPT",
        temperature=args.temp,
        pressure=args.pressure,
        steps=args.npt_steps,
        timecon=args.timecon
    )
    
    job_name = f"03_npt_{args.functional}_T{args.temp}_P{args.pressure}"
    
    job_dir = create_job_directory(
        output_dir=setup_config['output_dir'],
        job_name=job_name,
        coords_xyz=input_structure,
        mineral_name=setup_config['mineral_name'],
        run_type="MD",
        args=args,
        cell_params=setup_config['cell_params'],
        motion_section=motion_section
    )
    
    # Run calculation
    run_info = run_cp2k_calculation("NPT Production", job_dir, setup_config['cp2k_cmd'])
    
    # Write run info
    write_run_info(job_dir, run_info, setup_config['mineral_name'], 
                  args.temp, args.pressure, args.functional, args.cutoff)
    
    if not run_info['success']:
        raise RuntimeError("NPT production failed")
    
    return job_dir, run_info

def main():
    """Main simulation workflow"""
    
    # Parse arguments and setup simulation
    args = parse_arguments()
    
    try:
        # Setup simulation environment
        setup_config = setup_simulation(args)
        
        # Initialize tracking variables
        run_infos = []
        job_dirs = {}
        current_structure = setup_config['cleaned_xyz']
        
        # Step 1: Cell Optimization (unless skipped)
        if not args.skip_cellopt:
            cellopt_dir, cellopt_info = run_cell_optimization(setup_config, args)
            run_infos.append(cellopt_info)
            job_dirs['Cell Optimization'] = cellopt_dir
            
            # Find optimized structure for next step
            current_structure = find_optimized_structure(
                cellopt_dir, 
                setup_config['mineral_name'], 
                args.functional, 
                setup_config['cleaned_xyz']
            )
            
            if args.only_cellopt:
                print("\n" + "="*60)
                print("🎉 CELL OPTIMIZATION COMPLETED!")
                print("="*60)
                print(f"Results saved to: {setup_config['output_dir']}")
                return
        
        # Step 2: NVT Equilibration (unless skipped)
        if not args.skip_nvt:
            nvt_dir, nvt_info = run_nvt_equilibration(setup_config, args, current_structure)
            run_infos.append(nvt_info)
            job_dirs['NVT Equilibration'] = nvt_dir
            
            # Update structure for NPT (use same optimized structure)
            # In practice, you might want to use the final NVT structure
        
        # Step 3: NPT Production
        npt_dir, npt_info = run_npt_production(setup_config, args, current_structure)
        run_infos.append(npt_info)
        job_dirs['NPT Production'] = npt_dir
        
        # Create comprehensive summary
        create_simulation_summary(
            output_dir=setup_config['output_dir'],
            mineral_name=setup_config['mineral_name'],
            simulation_params=setup_config['sim_params'],
            job_dirs=job_dirs,
            run_infos=run_infos
        )
        
        # Print completion message
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