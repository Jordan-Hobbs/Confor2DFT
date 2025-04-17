from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import writers


def write_crest_input(args):

    ## Gen minimum conformer
    molecule = find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F", 
                                num_conf=2, max_opt_iters=1000)
    
    ## Write CREST input files
    writers.CRESTWriter("test", molecule)



def find_min_confromer(smiles, num_conf: int = 100, max_opt_iters: int = 1000): 

    #molecule = molecules.find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F", 
    #                            num_conf=2, max_opt_iters=1000)

    print("----------------------------------------------------------------")
    print("Generating initial CREST input structure:\n")

    ## Generate structure of a random unoptimised conformer from input smiles string
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    ## Generate conformers and optimise them
    print(f"Generating {num_conf} conformers for initial sorting "
          "and optimisation using RDKit")
    rdDistGeom.EmbedMultipleConfs(molecule_h, num_conf, 
        params = rdDistGeom.ETKDGv3())
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(molecule_h, 
        maxIters=max_opt_iters, ignoreInterfragInteractions=False)
    ## Check convergence of optimisation
    if all(conf_set[0] == 0 for conf_set in conf_energy):
        print("All conformers converged")
    else:
        print("WARNING! Not all conformers converged")

    ## Find minimum energy conformer
    min_energy_MMFF = 10000
    for index, energy in enumerate(conf_energy):
        if min_energy_MMFF > energy[1]:
            min_energy_MMFF = energy[1]
            min_energy_index_MMFF = index
    mol_min = Chem.Mol(molecule_h, False, min_energy_index_MMFF)

    print("")
    return mol_min