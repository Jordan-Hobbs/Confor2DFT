from molecules import RDMolGen

def main():
    mol = RDMolGen("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F")
    mol.molecule_init()
    mol.conformer_gen()
    mol.find_min_confromer()

if __name__ == "__main__":
    main()
