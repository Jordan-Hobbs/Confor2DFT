import molecules
import writers

def main():

    print("   ___                 __               ___   ___    ___   _____ ")
    print("  / __|  ___   _ _    / _|  ___   _ _  |_  ) |   \\  | __| |_   _|")
    print(" | (__  / _ \\ | ' \\  |  _| / _ \\ | '_|  / /  | |) | | _|    | |  ")
    print("  \\___| \\___/ |_||_| |_|   \\___/ |_|   /___| |___/  |_|     |_| ")
    print("Author: Dr Jordan Hobbs, University of Leeds\n")


    ## gen minimum conformer 
    molecule = molecules.find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F", 
                                num_conf=2, max_opt_iters=1000)

    CREST = writers.CRESTWriter("test", molecule)
    print("Writing input files:")
    CREST.write_xyz()
    CREST.write_toml()
    CREST.write_sh()
    print("")

if __name__ == "__main__":
    main()
