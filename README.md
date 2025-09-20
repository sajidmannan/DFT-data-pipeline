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
