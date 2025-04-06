import cmd_parser
import molecules
import writers

def main():
    print("")
    print("   ___                 __               ___   ___    ___   _____ ")
    print("  / __|  ___   _ _    / _|  ___   _ _  |_  ) |   \\  | __| |_   _|")
    print(" | (__  / _ \\ | ' \\  |  _| / _ \\ | '_|  / /  | |) | | _|    | |  ")
    print("  \\___| \\___/ |_||_| |_|   \\___/ |_|   /___| |___/  |_|     |_| ")
    print("Author: Dr Jordan Hobbs, University of Leeds")
    print("https://github.com/Jordan-Hobbs/")
    print("")

    ## Load program settings from cmd line and check for input config file
    parser = cmd_parser.ProgramController()
    args = parser.get_args()
    if args.Config:
        print("Config file found. Loading parameters from config file")
        parser.load_from_config()

    ## Gen minimum conformer
    #molecule = molecules.find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C(OC3=CC(F)=C(F)C(F)=C3)(F)F)C(F)=C2)C(F)=C1)C4=C(F)C=C(C5OCC(CCC)CO5)C=C4F", 
    #                            num_conf=2, max_opt_iters=1000)

    ## Write CREST input files
    #writers.CRESTWriter("test", molecule)

if __name__ == "__main__":
    main()
