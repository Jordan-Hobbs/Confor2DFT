from rdkit import Chem

def molecule_init(smiles):
        molecule = Chem.MolFromSmiles(smiles)
        molecule = Chem.AddHs(molecule)
        molecule = Chem.rdCoordGen.AddCoords(molecule)
        return molecule