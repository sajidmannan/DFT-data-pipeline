# DFT-data-pipeline
Framework for DFT data generation


cp2k_pipeline/
│── main.py              # entry point (parsing + workflow)
│── utils.py             # helper functions (xyz, input writer, runners, etc.)
│── config.py            # default settings, mineral-specific cells, functional maps
│── init.py          # optional

# CP2K Simulation System - Setup and Usage Guide

## 📁 File Structure

First, create a directory and save all the files:

```bash
mkdir cp2k_simulations
cd cp2k_simulations
```

Save these files in your directory:
- `cp2k_utils.py` - Utility functions (from artifact 1)
- `cp2k_sim.py` - Main simulation script (from artifact 2) 
- `cp2k_config.py` - Configuration module (from artifact 3)

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

### Get Help
```bash
python cp2k_sim.py --help
```

## 📂 Example XYZ File Format

Create a file called `SiO2.xyz`:
```
9
Silicon dioxide unit cell
Si    0.4701    0.0000    0.0000
Si    0.0000    0.4701    0.6667
Si    0.4701    0.4701    0.3333
O     0.4139    0.2669    0.1188
O     0.7331    0.4139    0.2146
O     0.2669    0.7331    0.5479
O     0.5861    0.9701    0.4521
O     0.0299    0.5861    0.7854
O     0.9701    0.0299    0.8812
```

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

## 🔧 Troubleshooting

### Common Issues:

1. **CP2K not found**:
   ```bash
   # Load CP2K module (on clusters)
   module load cp2k
   
   # Or set PATH manually
   export PATH=/path/to/cp2k/exe:$PATH
   ```

2. **XYZ file not found**:
   ```bash
   # Make sure file exists and has correct format
   ls -la *.xyz
   head -5 your_mineral.xyz
   ```

3. **Permission errors**:
   ```bash
   # Make sure Python files are executable
   chmod +x cp2k_sim.py
   
   # Check output directory permissions
   mkdir -p results && ls -ld results/
   ```

4. **Import errors**:
   ```bash
   # Make sure all Python files are in same directory
   ls -la cp2k_*.py
   
   # Run from the directory containing the scripts
   python cp2k_sim.py --help
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

### High-Quality Run:
```bash
python cp2k_sim.py mineral.xyz \
    --cutoff 800 \
    --nvt-steps 10000 \
    --npt-steps 50000 \
    --output high_quality/
```

### Only Cell Optimization:
```bash
python cp2k_sim.py mineral.xyz --only-cellopt --output cellopt_only/
```

## 📈 Batch Processing Example

Create a bash script (`run_batch.sh`):
```bash
#!/bin/bash

# Temperature series
for temp in 300 500 700 1000; do
    python cp2k_sim.py SiO2.xyz --temp $temp --output SiO2_T${temp}/
done

# Pressure series  
for pressure in 1 10 50 100; do
    python cp2k_sim.py CaCO3.xyz --pressure $pressure --temp 300 --output CaCO3_P${pressure}/
done
```

Run it:
```bash
chmod +x run_batch.sh
./run_batch.sh
```

## 🔍 Monitoring Progress

The script provides real-time feedback:
```
CP2K MINERAL SIMULATION PIPELINE
============================================================
Mineral: SiO2
Temperature: 300.0 K
Pressure: 1.0 bar
Output directory: results
Functional: PBE_D3
============================================================
Using CP2K executable: /usr/local/bin/cp2k.psmp
Found XYZ file: SiO2.xyz

==================================================
STEP 1: CELL OPTIMIZATION
==================================================
Running Cell Optimization calculation in results/SiO2_T300.0_P1.0/01_cellopt_PBE_D3_T300_P1.0
✓ Cell Optimization completed successfully in 0:05:23
```

That's it! You now have a complete, modular CP2K simulation system ready to use.