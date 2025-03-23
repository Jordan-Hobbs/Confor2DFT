from rdkit import Chem
from rdkit.Chem import rdDistGeom
import py3Dmol

class molec:
    def __init__(self, smiles):
        self.smiles = smiles
      
    def molecule_init(self):
        self.molecule = Chem.MolFromSmiles(self.smiles)
        self.molecule_h = Chem.AddHs(self.molecule)
        Chem.rdCoordGen.AddCoords(self.molecule_h)
        rdDistGeom.EmbedMolecule(self.molecule_h)

    def threeD_view(self):
        view = py3Dmol.view(data=Chem.MolToMolBlock(self.molecule_h), style={"stick": {}, "sphere": {"scale": 0.3}})
        view.zoomTo()
        return view
    

    