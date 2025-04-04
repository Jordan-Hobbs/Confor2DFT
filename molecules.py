from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, Draw


def find_min_confromer(smiles, num_conf: int = 100, max_opt_iters: int = 1000): 


    ## Generate structure of a random unoptimised conformer from input smiles string
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)


    ## Generate conformers and optimise them
    rdDistGeom.EmbedMultipleConfs(molecule_h, num_conf, params=rdDistGeom.ETKDGv3())
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(molecule_h, maxIters=max_opt_iters, ignoreInterfragInteractions=False)
    ## Check convergence of optimisation
    if all(conf_set[0] == 0 for conf_set in conf_energy):
        print("All conformers converged")
    else:
        print("Not all conformers converged")


    ## Find minimum energy conformer
    min_energy_MMFF = 10000
    for index, energy in enumerate(conf_energy):
        if min_energy_MMFF > energy[1]:
            min_energy_MMFF = energy[1]
            min_energy_index_MMFF = index
    mol_min = Chem.Mol(molecule_h, False, min_energy_index_MMFF)

    return mol_min