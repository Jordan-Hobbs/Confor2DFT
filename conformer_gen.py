from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import writers
from exceptions import InvalidSmilesError

def gen(args):
    try:
        molecule = find_min_conformer(args)
    except InvalidSmilesError as e:
        print(f"Error: {e}. Please provide a valid SMILES code")
        return
    
    if args.ConformerMode == "crest":
        crest_files = writers.CRESTWriter(molecule, args)
        crest_files.write_xyz()
        crest_files.write_toml()
        crest_files.write_sh()
    elif args.ConformerMode == "orca":
        orca_files = writers.ORCAWriter(molecule, args)
        orca_files.write_inp()
        orca_files.write_sh()

def find_min_conformer(args):
    print("\n----------------------------------------------------------------")
    print("Generating initial CREST input structure:\n")

    molecule = Chem.MolFromSmiles(args.SmilesString)
    if molecule is None:
        raise InvalidSmilesError(f"Invalid SMILES string: '{args.SmilesString}'")

    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    print(
        f"Generating {args.RDNumConf} conformers for initial sorting and "
        "optimization using RDKit."
    )
    rdDistGeom.EmbedMultipleConfs(
        molecule_h,
        args.RDNumConf,
        params=rdDistGeom.ETKDGv3()
    )
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
        molecule_h,
        maxIters=args.RDMaxIter,
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