import cmd_parser
import conformer_gen


def main():
    """Main function to run the program."""

    # Sets command functions
    commands = {
        "conformer_gen": conformer_gen.conformer_gen
    }

    # Initialize parser
    parser = cmd_parser.ProgramController(commands)
    args = parser.get_args()

    print("")
    print("   ____             __            ____  _____ _____ ")
    print("  / ___|___  _ __  / _| ___  _ __|  _ \\|  ___|_   _|")
    print(" | |   / _ \\| '_ \\| |_ / _ \\| '__| | | | |_    | |  ")
    print(" | |__| (_) | | | |  _| (_) | |  | |_| |  _|   | |  ")
    print("  \\____\\___/|_| |_|_|  \\___/|_|  |____/|_|     |_|  ")
    print("")                                                                                        
    print("Author: Dr. Jordan Hobbs, University of Leeds")
    print("https://github.com/Jordan-Hobbs/")

    print("\n----------------------------------------------------------------")
    print("Setting up config:\n")

    # Loads config file if config file provided
    if args.Config:
        print(
            "Config input detected. Attempting to load parameters from "
            f"\"{args.Config}\""
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
