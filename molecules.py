from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, Draw


def find_min_confromer(smiles, no_conf: int = 1, max_opt_iter: int = 1): 

    ## Generate the minimimum structure of the a random conformer of input simles string
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)

    ## Generate conformers and optimise them
    rdDistGeom.EmbedMultipleConfs(molecule_h, no_conf, params=rdDistGeom.ETKDGv3())
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(molecule_h, max_opt_iter)
    print(conf_energy)
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

    Draw.MolToImage(mol_min)
    return mol_min

    

    