import cmd_parser


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
        print("Config input detected. Attempting to load parameters from config file.")
        parser.load_from_config()

    




if __name__ == "__main__":
    main()
