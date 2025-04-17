from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import writers


def crest_input(args):
    # Generate minimum conformer
    molecule = find_min_conformer(
        args.SmilesString,
        num_conf=args.NumberConformers,
        max_opt_iters=args.MaximumIterations
    )

    # Write CREST input files
    writers.CRESTWriter(
        molecule,
        args.FileName
    )


def find_min_conformer(smiles, num_conf: int = 100, max_opt_iters: int = 1000):
    print("\n----------------------------------------------------------------")
    print("Generating initial CREST input structure:\n")

    # Generate structure of a random unoptimized conformer from input SMILES 
    # string
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    # Generate conformers and optimize them
    print(
        f"Generating {num_conf} conformers for initial sorting and "
        "optimisation using RDKit."
    )
    rdDistGeom.EmbedMultipleConfs(
        molecule_h,
        num_conf,
        params=rdDistGeom.ETKDGv3()
    )
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
        molecule_h,
        maxIters=max_opt_iters,
        ignoreInterfragInteractions=False
    )

    # Check convergence of optimization
    if all(conf_set[0] == 0 for conf_set in conf_energy):
        print("All conformers converged.")
    else:
        print("WARNING! Not all conformers converged.")

    # Find minimum energy conformer
    min_energy = float("inf")
    min_index = 0
    for index, (_, energy) in enumerate(conf_energy):
        if energy < min_energy:
            min_energy = energy
            min_index = index
    mol_min = Chem.Mol(molecule_h, False, min_index)

    return mol_min
