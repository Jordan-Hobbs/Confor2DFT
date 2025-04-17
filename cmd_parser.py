import argparse
import tomllib
import sys

import utils


class ProgramController:
    def __init__(self, commands):
        self.commands = commands

        # Create parent parser for arguments shared across subparsers
        self.parent_parser = argparse.ArgumentParser(add_help=False)
        self.parent_parser.add_argument(
            "-c", "--Config",
            type=str,
            help="Name of the config file used to input required settings."
        )

        # Create parser and subparsers for input arguments
        self.parser = argparse.ArgumentParser(
            prog="Confor2DFT",
            parents=[self.parent_parser],
            description=(
                "Combines CREST and some DFT package to obtain optimised "
                "molecules."
            )
        )
        self.subparsers = self.parser.add_subparsers(
            required=True,
            dest="Command",
            help=(
                "Select operating mode: 'crest_input' or 'orca_write'. "
                "In a typical workflow, 'crest_input' is the first stage, "
                "followed by 'orca_write'."
            )
        )

        self.parse_crest_input()
        self.parse_orca_input()

    def parse_crest_input(self):
        """Create parser for CREST input mode."""
        crest_write = self.subparsers.add_parser(
            "crest_input",
            parents=[self.parent_parser],
            help="Generate input files for CREST conformer search."
        )

        crest_write.add_argument(
            "SmilesString",
            type=str,
            help="SMILES string of the compound to be optimised."
        )
        crest_write.add_argument(
            "-f", "--FileName",
            type=str,
            default="conf2dft",
            help="File name to use for CREST input files."
        )
        crest_write.add_argument(
            "-n", "--NumberConformers",
            type=int,
            default=100,
            help="Number of conformers for the initial RDKit search."
        )
        crest_write.add_argument(
            "-m", "--MaximumIterations",
            type=int,
            default=1000,
            help=(
                "Number of iterations for the initial RDKit optimisation "
                "step."
            )
        )
        crest_write.add_argument(
            "-ol", "--OptimisationLevel",
            type=str,
            default="extreme",
            choices=[
                "crude", "vloose", "loose", "normal",
                "tight", "vtight", "extreme"
            ],
            help="Convergence conditions for the CREST optimisation steps."
        )
        crest_write.add_argument(
            "-xtb", "--GFNMethod",
            type=str,
            default="gfn2",
            choices=["gfnff", "gfn0", "gfn1", "gfn2"],
            help=(
                "Tight binding method used by CREST for optimisation and "
                "conformer searching."
            )
        )
        crest_write.add_argument(
            "-nc", "--NumCPUs",
            type=int,
            default=4,
            help="Number of CPUs allocated to the HPC calculation."
        )
        crest_write.add_argument(
            "-rt", "--RunTime",
            type=str,
            default="24:00:00",
            help="Maximum runtime allowed on the HPC system."
        )
        crest_write.add_argument(
            "-e", "--Email",
            type=str,
            help=(
                "Email address to receive updates on job start, finish, "
                "and failure."
            )
        )

        crest_write.set_defaults(func=self.commands["crest_input"])

    def parse_orca_input(self):
        """Create parser for ORCA input mode."""
        orca_write = self.subparsers.add_parser(
            "orca_write",
            parents=[self.parent_parser],
            help="Generate ORCA input files from CREST outputs."
        )
        orca_write.add_argument(
            "-f",
            type=str,
            default="conf2dft",
            help="File name to use for ORCA input files."
        )

    def load_from_config(self):
        try:
            with open(self.args.Config, "rb") as file:
                config = tomllib.load(file)
            print("Config file found.")
        except FileNotFoundError:
            print("Config file not found. Using default parameters.")
            return

        # Flatten the config file to remove nested structure
        flat_config = utils.flatten_dict(config.get(self.args.Command))

        # Load input args from command line except for -c / --Config
        cmdline_args = [
            arg for arg in sys.argv[2:]
            if arg.startswith(("-", "--")) and arg not in ("-c", "--Config")
        ]

        # Convert any short arg names in cmdline_args to the long form
        subparser = self.subparsers.choices[self.args.Command]
        for action in subparser._actions[2:]:
            arg_names = action.option_strings
            if len(arg_names) > 1:
                cmdline_args = [
                    arg_names[1].strip("--") if item == arg_names[0] else item
                    for item in cmdline_args
                ]

        # Remove keys from config that were already provided via CLI
        for key in cmdline_args:
            flat_config.pop(key, None)

        # Update args with remaining config values
        for key, value in flat_config.items():
            setattr(self.args, key, value)

        print("Config parameters loaded successfully.")

    def get_args(self):
        """Parse and return command-line arguments."""
        self.args = self.parser.parse_args()
        return self.args
