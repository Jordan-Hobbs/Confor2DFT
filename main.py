import cmd_parser
import crest_input


def main():
    """Main function to run the program."""
    print("")
    print("   ___                 __               ___   ___    ___   _____ ")
    print("  / __|  ___   _ _    / _|  ___   _ _  |_  ) |   \\  | __| |_   _|")
    print(
        " | (__  / _ \\ | ' \\  |  _| / _ \\ | '_|  / /  | |) | | _|    | |  "
    )
    print("  \\___| \\___/ |_||_| |_|   \\___/ |_|   /___| |___/  |_|     |_| ")
    print("Author: Dr. Jordan Hobbs, University of Leeds")
    print("https://github.com/Jordan-Hobbs/")
    print("")

    # Load program settings from command line and check for input config file
    commands = {
        "crest_input": crest_input.crest_input
    }

    # Initialize parser and load config if provided
    print("\n----------------------------------------------------------------")
    print("Setting up config:\n")
    parser = cmd_parser.ProgramController(commands)
    args = parser.get_args()

    if args.Config:
        print(
            "Config input detected. Attempting to load parameters from "
            "config file."
        )
        parser.load_from_config()
    else:
        print(
            "No config file provided. Using command line and defaults "
            "instead."
        )

    # Execute the function corresponding to the selected command
    if hasattr(args, 'func'):
        args.func(args)
        print("")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
