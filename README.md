# CP2K Simulation System - Setup and Usage Guide

## 🛠️ Prerequisites

1. **CP2K Installation**: Make sure CP2K is installed and available in your PATH
   ```bash
   # Test if CP2K is available
   which cp2k.psmp  # or cp2k.popt or cp2k
   ```

2. **Python Environment**: Python 3.6+ with standard libraries (no additional packages needed)

3. **XYZ Files**: Have your mineral structure files ready (e.g., `SiO2.xyz`, `quartz.xyz`)

## 🚀 Basic Usage

### 1. Simple Run
```bash
python cp2k_sim.py SiO2.xyz --temp 300 --pressure 1.0 --output results/
```

### 2. With Custom Parameters
```bash
python cp2k_sim.py quartz.xyz \
    --temp 500 \
    --pressure 2.5 \
    --functional PBE_D3 \
    --cutoff 800 \
    --npt-steps 20000 \
    --output high_temp_quartz/
```

### 3. Quick Test Run
```bash
python cp2k_sim.py CaCO3.xyz \
    --temp 300 \
    --nvt-steps 1000 \
    --npt-steps 2000 \
    --cutoff 400 \
    --output test_run/
```

## 📋 Command Line Options

### Required Arguments
- `mineral` - Mineral name or path to XYZ file

### Optional Arguments
- `--temp`, `-T` - Temperature in Kelvin (default: 300.0)
- `--pressure`, `-P` - Pressure in bar (default: 1.0)
- `--output`, `-o` - Output directory (default: cp2k_results)
- `--functional` - XC functional: PBE, PBE_D3, BLYP, BLYP_D3, PBE0 (default: PBE_D3)
- `--cutoff` - Plane wave cutoff in Ry (default: 600)
- `--rel-cutoff` - Relative cutoff in Ry (default: 40)
- `--nvt-steps` - NVT equilibration steps (default: 5000)
- `--npt-steps` - NPT production steps (default: 10000)
- `--timecon` - Thermostat time constant in fs (default: 100.0)
- `--cell` - Cell parameters 'a b c' or 'auto' (default: auto)

### Job Control Options
- `--skip-cellopt` - Skip cell optimization
- `--skip-nvt` - Skip NVT equilibration  
- `--only-cellopt` - Only run cell optimization


## 🎯 Step-by-Step Example

### 1. Prepare Your System
```bash
# Create working directory
mkdir my_simulation
cd my_simulation

# Copy the Python files
cp /path/to/cp2k_*.py .

# Create or copy your XYZ file
# (example: SiO2.xyz)
```

### 2. Run a Test Simulation
```bash
# Quick test with minimal steps
python cp2k_sim.py SiO2.xyz \
    --temp 300 \
    --nvt-steps 1000 \
    --npt-steps 2000 \
    --output test_SiO2/
```

### 3. Check the Results
```bash
# Look at the output structure
ls -la test_SiO2/
cat test_SiO2/SiO2_T300.0_P1.0/simulation_summary.txt
```

## 📊 Expected Output Structure

After running, you'll get:
```
output_directory/
└── mineral_T300_P1.0/
    ├── simulation_summary.txt           # Overall summary
    ├── original_mineral.xyz             # Your original file
    ├── coords.xyz                       # Cleaned coordinates
    ├── 01_cellopt_PBE_D3_T300_P1.0/    # Cell optimization
    │   ├── input.inp
    │   ├── output.out
    │   ├── coords.xyz
    │   └── *trajectory files*
    ├── 02_nvt_PBE_D3_T300_P1.0/        # NVT equilibration
    └── 03_npt_PBE_D3_T300_P1.0/        # NPT production
```

## ⚡ Quick Start Commands

### For SiO2 (Quartz):
```bash
python cp2k_sim.py SiO2.xyz --temp 300 --output SiO2_ambient/
```

### For CaCO3 (Calcite):
```bash
python cp2k_sim.py CaCO3.xyz --temp 500 --pressure 5.0 --output CaCO3_heated/
```


# VASP Simulation System - Setup and Usage Guide

## 🛠️ Prerequisites

### 1. VASP Installation
Make sure VASP is installed and available in your PATH:
```bash
# Test if VASP is available
which vasp_std  # or vasp_gam, vasp_ncl, vasp
```

### 2. VASP Pseudopotentials
Set up the VASP pseudopotential directory:
```bash
# Set environment variable (add to ~/.bashrc)
export VASP_PSP_DIR=/path/to/vasp/potentials

# Or use alternative names
export VASP_PP_PATH=/path/to/vasp/potentials
```

The potentials directory should have structure like:
```
potentials/
├── H/POTCAR
├── C/POTCAR
├── O/POTCAR
├── Si/POTCAR
├── Ca_pv/POTCAR
├── Mg_pv/POTCAR
└── ...
```

### 3. POSCAR Files
Have your mineral structure files ready in VASP POSCAR format

## 🚀 Basic Usage Examples

### 1. Structural Optimization
```bash
# Basic relaxation
python vasp_sim.py SiO2.POSCAR --template relax --output results/

# High-precision relaxation
python vasp_sim.py quartz --template relax --functional PBEsol --encut 600 --prec High
```

### 2. Electronic Structure Calculations
```bash
# Self-consistent field calculation
python vasp_sim.py CaCO3.POSCAR --template scf --kpts 8,8,6

# Density of states
python vasp_sim.py mineral.POSCAR --template dos --functional HSE06

# Band structure calculation
python vasp_sim.py structure.POSCAR --template band --kpts 12,12,12
```

### 3. Molecular Dynamics
```bash
# NVT molecular dynamics at 500K
python vasp_sim.py forsterite.POSCAR --template md_nvt --temp 500 --steps 10000

# NPT molecular dynamics with pressure
python vasp_sim.py MgO.POSCAR --template md_npt --temp 1000 --pressure 10 --steps 20000
```

### 4. Advanced Calculations
```bash
# Phonon calculation preparation
python vasp_sim.py structure.POSCAR --template phonon --encut 600

# Elastic constants
python vasp_sim.py crystal.POSCAR --template elastic --prec Accurate

# Spin-polarized calculation
python vasp_sim.py Fe2O3.POSCAR --template relax --spin --magmom "5,5,-5,-5,0,0,0"
```

## 📋 Command Line Options

### Required Arguments
- `mineral` - Mineral name or path to POSCAR file

### Calculation Templates
- `--template`, `-t` - Simulation type:
  - `relax` - Structural optimization (default)
  - `scf` - Self-consistent field calculation
  - `dos` - Density of states
  - `band` - Band structure
  - `md_nvt` - NVT molecular dynamics
  - `md_npt` - NPT molecular dynamics
  - `phonon` - Phonon calculation prep
  - `elastic` - Elastic constants

### Electronic Structure
- `--functional` - XC functional: PBE, PBEsol, LDA, SCAN, HSE06, PBE0
- `--encut` - Plane wave cutoff in eV (auto from POTCAR if not specified)
- `--prec` - Precision: Low, Medium, High, Normal, Accurate
- `--algo` - Algorithm: Normal, VeryFast, Fast, Conjugate, All

### K-points
- `--kpts` - K-point grid: 'nx,ny,nz', 'auto', or density value
- `--kpts-shift` - K-point shift: 'x,y,z'

### Optimization & MD
- `--ediffg` - Ionic convergence (eV/Å, default: -0.01)
- `--isif` - Stress/relaxation flag
- `--steps` - Maximum ionic steps
- `--temp` - Temperature for MD (K)
- `--pressure` - External pressure (kBar)

### Advanced Options
- `--vdw` - Van der Waals: D2, D3, D3BJ, TS, BEEF
- `--spin` - Enable spin polarization
- `--magmom` - Initial magnetic moments
- `--mpi` - Number of MPI ranks
- `--ncore` - NCORE parallelization parameter

### Utility Options
- `--dry-run` - Generate input files only
- `--continue-run` - Continue from CONTCAR
- `--list-minerals` - Show available minerals
- `--list-functionals` - Show available functionals  
- `--list-templates` - Show available templates
<!-- 
## 📝 POSCAR File Format Example

Create a file called `SiO2.POSCAR`:
```
Si2O4 unit cell
   1.00000000
     4.9160000000    0.0000000000    0.0000000000
     0.0000000000    4.9160000000    0.0000000000
     0.0000000000    0.0000000000    5.4050000000
   Si   O
     3     6
Direct
  0.4701  0.0000  0.0000
  0.0000  0.4701  0.6667
  0.4701  0.4701  0.3333
  0.4139  0.2669  0.1188
  0.7331  0.4139  0.2146
  0.2669  0.7331  0.5479
  0.5861  0.9701  0.4521
  0.0299  0.5861  0.7854
  0.9701  0.0299  0.8812
```

## 🎯 Step-by-Step Examples

### Example 1: Basic Structural Optimization
```bash
# 1. Prepare your system
mkdir quartz_optimization
cd quartz_optimization

# 2. Copy the Python files
cp /path/to/vasp_*.py .

# 3. Create or copy your POSCAR file
# (example: quartz.POSCAR)

# 4. Run optimization
python vasp_sim.py quartz.POSCAR --template relax --output results/

# 5. Check results
ls -la results/
cat results/quartz_relax_PBE/simulation_summary.txt
```

### Example 2: High-Accuracy Electronic Structure
```bash
# High-precision SCF with hybrid functional
python vasp_sim.py structure.POSCAR \
    --template scf \
    --functional HSE06 \
    --encut 600 \
    --kpts 12,12,8 \
    --prec Accurate \
    --output high_accuracy/
```

### Example 3: Temperature Series MD
```bash
# Create temperature series
for temp in 300 500 700 1000; do
    python vasp_sim.py mineral.POSCAR \
        --template md_nvt \
        --temp $temp \
        --steps 5000 \
        --output MD_T${temp}/
done
``` -->

## 📊 Expected Output Structure

After running a calculation:
```
output_directory/
└── mineral_template_functional/
    ├── simulation_summary.txt    # Overall summary
    ├── template_functional/      # Main calculation directory
    │   ├── INCAR                # Input parameters
    │   ├── POSCAR               # Structure file
    │   ├── POTCAR               # Pseudopotentials
    │   ├── KPOINTS              # K-point sampling
    │   ├── OUTCAR               # Detailed output
    │   ├── OSZICAR              # Convergence info
    │   ├── vasprun.xml          # XML output
    │   ├── CONTCAR              # Final structure
    │   ├── CHGCAR               # Charge density
    │   ├── run_info.txt         # Job information
    │   └── previous_run/        # Backup (if continuing)
    └── ...
```
