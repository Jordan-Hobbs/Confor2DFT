import molecules
import crest_writer

def main():

    ## gen minimum conformer 
    molecule = molecules.find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F", no_conf=10)

    CREST = crest_writer.CRESTGen("test", molecule)
    CREST.write_xyz()
    CREST.write_toml("ancopt")
    CREST.write_sh()

if __name__ == "__main__":
    main()
