from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, Draw
import py3Dmol

class RDMolGen:
    def __init__(self, smiles):
        self.smiles = smiles
      
    def molecule_init(self):
        self.molecule = Chem.MolFromSmiles(self.smiles)
        self.molecule_h = Chem.AddHs(self.molecule)
        Chem.rdCoordGen.AddCoords(self.molecule_h)
        rdDistGeom.EmbedMolecule(self.molecule_h)
        rdForceFieldHelpers.MMFFOptimizeMolecule(self.molecule_h, maxIters=1000, ignoreInterfragInteractions=False)

    def conformer_gen(self):
        rdDistGeom.EmbedMultipleConfs(self.molecule_h, 10, params=rdDistGeom.ETKDGv3())
        self.conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(self.molecule_h, maxIters=1000, ignoreInterfragInteractions=False)
        if all(conf_set[0] == 0 for conf_set in self.conf_energy):
            print("All conformers converged")
        else:
            print("Not all conformers converged")

    def find_min_confromer(self):
        min_energy_MMFF = 10000
        for index, energy in enumerate(self.conf_energy):
            if min_energy_MMFF > energy[1]:
                min_energy_MMFF = energy[1]
                min_energy_index_MMFF = index
        self.mol_min = Chem.Mol(self.molecule_h, False, min_energy_index_MMFF)
        Draw.MolToImage(self.mol_min) # This doesnt work?

    def return_RD_molecule(self):
        return self.molecule_h

    def threeD_view(self):
        view = py3Dmol.view(data=Chem.MolToMolBlock(self.molecule_h), style={"stick": {}, "sphere": {"scale": 0.3}})
        view.zoomTo()
        return view
    

    