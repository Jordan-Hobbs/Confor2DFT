from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import writers
from exceptions import InvalidSmilesError

def conformer_gen(args):
    try:
        molecule = find_min_conformer(
            args.SmilesString,
            num_conf=args.NumberConformers,
            max_opt_iters=args.MaximumIterations
        )
    except InvalidSmilesError as e:
        print(f"Error: {e}. Please provide a valid SMILES code")
        return
    
    if args.ConformerMode == "CREST":
        writers.CRESTWriter(molecule, args.FileName)
    elif args.ConformerMode == "ORCA":
        print("stuff")

def find_min_conformer(smiles, num_conf: int = 100, max_opt_iters: int = 1000):
    print("\n----------------------------------------------------------------")
    print("Generating initial CREST input structure:\n")

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise InvalidSmilesError(f"Invalid SMILES string: '{smiles}'")

    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    print(
        f"Generating {num_conf} conformers for initial sorting and "
        "optimization using RDKit."
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

    if all(conf_set[0] == 0 for conf_set in conf_energy):
        print("All conformers converged.")
    else:
        print("WARNING! Not all conformers converged.")

    min_energy = float("inf")
    min_index = 0
    for index, (_, energy) in enumerate(conf_energy):
        if energy < min_energy:
            min_energy = energy
            min_index = index
    mol_min = Chem.Mol(molecule_h, False, min_index)

    return mol_min