import numpy as np

# Conversion factor: 1 Hartree = 27.211386245988 eV
HARTREE_TO_EV = 27.211386245988

# Input files
xyz_file = "SiO2.xyz"
pos_file = "SiO2_position"
force_file = "SiO2_force"
stress_file = "SiO2.stress"

# Output file
out_file = "SiO2.extxyz"

# --- Read lattice and atom types from xyz (first file) ---
with open(xyz_file, "r") as f:
    lines = f.readlines()

natoms = int(lines[0].strip())
lattice = lines[1].split('Lattice="')[1].split('"')[0]

species = []
for i in range(2, 2 + natoms):
    species.append(lines[i].split()[0])

# --- Read positions and energies ---
positions_by_step = []
energies = []

with open(pos_file, "r") as f:
    step_positions = []
    for line in f:
        if "E =" in line:
            if step_positions:
                positions_by_step.append(step_positions)
                step_positions = []
            # Convert energy from Hartree to eV
            energy_hartree = float(line.split("=")[-1])
            energy_ev = energy_hartree * HARTREE_TO_EV
            energies.append(energy_ev)
        elif line.strip().startswith(("Si", "O")):
            parts = line.split()
            step_positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if step_positions:
        positions_by_step.append(step_positions)

# --- Read forces ---
forces_by_step = []
with open(force_file, "r") as f:
    step_forces = []
    for line in f:
        if line.strip().startswith("i ="):
            if step_forces:
                forces_by_step.append(step_forces)
                step_forces = []
        elif line.strip().startswith(("Si", "O")):
            parts = line.split()
            step_forces.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if step_forces:
        forces_by_step.append(step_forces)

# --- Read stresses ---
stresses = []
with open(stress_file, "r") as f:
    for line in f:
        parts = line.split()
        if len(parts) == 11 and parts[0].isdigit():
            stress_vals = parts[2:11]  # xx xy xz yx yy yz zx zy zz
            stresses.append(" ".join(stress_vals))

# --- Write EXTXYZ trajectory ---
with open(out_file, "w") as f:
    nsteps = min(len(positions_by_step), len(forces_by_step), len(stresses), len(energies))
    for step in range(nsteps):
        f.write(f"{natoms}\n")
        f.write(f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:forces:R:3 '
                f'energy={energies[step]} stress="{stresses[step]}" free_energy={energies[step]} pbc="T T T"\n')
        for i in range(natoms):
            sx = species[i]
            px, py, pz = positions_by_step[step][i]
            fx, fy, fz = forces_by_step[step][i]
            f.write(f"{sx:<2} {px:15.8f} {py:15.8f} {pz:15.8f} {fx:15.8e} {fy:15.8e} {fz:15.8e}\n")

print(f"✅ Written EXTXYZ trajectory with {nsteps} steps: {out_file}")
