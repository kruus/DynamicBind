from openbabel import openbabel as ob
from pathlib import Path
datadir = Path("/datasets/data/PDBBind_processed")
fname = datadir / "4op1" / "4op1_protein_processed.pdb"

obconv = ob.OBConversion()
obconv.SetInAndOutFormats("pdb", "fasta")
mol = ob.OBMol()
obconv.ReadFile(mol, str(fname))
fasta = obconv.WriteString(mol)
print(f"{fname} -> fast string using python openbabel:\n{fasta}", flush=True)

# convert input pdb vi openbabel into openmm format.
obconv.SetOutFormat("pdb")
str_pdb = obconv.WriteString(mol)
from openmm.app import PDBFile, Topology, Modeller, ForceField
from io import StringIO
om_pdb = PDBFile(StringIO(str_pdb))
topo = om_pdb.getTopology()
print(f"{topo=}", flush=True)

# openmm add hydrogens and minimize
# ff = ForceField('amber14-all.xml')
# modeller = Modeller(om_pdb.topology, om_pdb.positions)
# modeller.addHydrogens(ff)       # error terminal group probably missing?
from pdbfixer import PDBFixer
# fix_pdb = PDBFixer(StringIO(str_pdb))
sio = StringIO(str_pdb)
fixer = PDBFixer(pdbfile=sio)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens()
PDBFile.writeFile(fixer.topology,
                  fixer.positions,
                  open("fixer.pdb","w"))

# we can also optimize w/ openmm and our own force field
output = StringIO()
PDBFile.writeFile(fixer.topology, fixer.positions, output)
str_fix = output.getvalue() #.replace('\n', "\n")
# output.seek(0)
# str_fix = output.read()
output.close()
# Note: f-string will escape things.  I want unescaped newlines here
print("str_fix[0:100]", str_fix[0:100], flush=True)
#
ff = ForceField('amber14-all.xml')
modeller = Modeller(fixer.topology, fixer.positions)
modeller.addHydrogens(ff)       # error terminal group probably missing?
PDBFile.writeFile(modeller.topology,
                  modeller.positions,
                  open("foo.pdb","w"))

from openmm.app import (PME, HBonds, Simulation)
from openmm import (LangevinMiddleIntegrator, )
from openmm.unit import (nanometer, kelvin, picosecond, picoseconds, angstrom)
system = ff.createSystem(fixer.topology,
                         # nonbondedMethod=PME,
                         # nonbondedCutoff = 1 * nanometer,
                         constraints=HBonds,
                         )
for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)

state = simulation.context.getState(getEnergy=True)
print(f"Initial energy {state.getPotentialEnergy()}", flush=True)
for i, f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    print(f"{f.getName():30s} energy {state.getPotentialEnergy()}", flush=True)

pos0 = fixer.positions.copy()
print(f"{type(pos0)=}", flush=True)

simulation.minimizeEnergy()
state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"{ff=}\nminimized energy {state.getPotentialEnergy()}", flush=True)
pos1p = state.getPositions(asNumpy=False)               # positions (w/ .x .y .z)
pos1q = state.getPositions(asNumpy=True)                # bunch of Quantity
pos1 = state.getPositions(asNumpy=True) / angstrom      # plain numpy array
print(f"{type(pos1)=}", flush=True)
import numpy as np
print(f"mean( |pos1-pos0| ) = {np.sqrt(((((pos1 - pos0)** 2))*3).mean())} angstrom", flush=True)

from rdkit import Chem
from rdkit.Chem import AllChem
def mk_rdkit(topology, positions):
    rdkit_mol = Chem.RWMol()
    # Add atoms to the RDKit molecule
    for atom in topology.atoms():
        rdkit_atom = Chem.Atom(atom.element.symbol)
        rdkit_mol.AddAtom(rdkit_atom)
    # Add bonds to the RDKit molecule
    for bond in topology.bonds():
        rdkit_mol.AddBond(bond.atom1.index, bond.atom2.index, Chem.BondType.SINGLE)
    # Set atomic coordinates
    conf = Chem.Conformer(rdkit_mol.GetNumAtoms())
    for i, pos in enumerate(positions):
        conf.SetAtomPosition(i, (pos.x, pos.y, pos.z))
    rdkit_mol.AddConformer(conf)
    return rdkit_mol


rdkit0 = mk_rdkit(fixer.topology, fixer.positions)
print("OK, have rdkit0", flush=True)
rdkit1 = mk_rdkit(fixer.topology, pos1p)
# Optionally, you can sanitize the molecule
# Chem.SanitizeMol(rdkit_mol)
# Now you have an RDKit molecule
# print(Chem.MolToSmiles(rdkit_mol)
al_rmsd = Chem.rdMolAlign.AlignMol(rdkit0, rdkit1)
print(f"Aligned rmsd, from rdkit: {al_rmsd=}", flush=True)

#from openmmml import MLPotential
#potential = MLPotential('ani2x')
# done
