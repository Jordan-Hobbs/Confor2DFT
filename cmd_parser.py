import argparse
import tomllib
import sys
from typing import Any, Dict

import utils


class ProgramController:
    def __init__(self, commands: Dict[str, Any]) -> None:
        """Initialize the ProgramController with available commands."""
        self.commands = commands

        # Create parent parser for arguments shared across all subparsers
        self.top_parent_parser = argparse.ArgumentParser(add_help=False)

        self.top_parent_parser.add_argument(
            "-c", "--Config",
            type=str,
            help="Name of the config file used to input required settings."
        )

        # Create parser and subparsers for input arguments
        self.parser = argparse.ArgumentParser(
            prog="ConforDFT",
            parents=[self.top_parent_parser],
            description=(
                "Combines CREST and a DFT package to obtain optimised "
                "molecules."
            )
        )
        self.top_subparsers = self.parser.add_subparsers(
            required=True,
            dest="Command",
            help=(
                "Select operating mode: \"crest_input\" or \"orca_write\". "
                "Typical workflow: \"conformation_gen\" first, then "
                "\"orca_write\"."
            )
        )

        self.parse_conformer_input()
        self.parse_orca_input()

    def parse_conformer_input(self) -> None:
        """Create parser for conformer input mode."""
        self.conformer_write = self.top_subparsers.add_parser(
            "conformer_gen",
            parents=[self.top_parent_parser],
            help="Generate input files for a conformer search."
        )

        self.conformer_subparsers = self.conformer_write.add_subparsers(
            required=True,
            dest="ConformerMode",
            help="Specify which package to write input files for."
        )

        # Create shared args parser for conformer mode
        self.conformer_common_parser = argparse.ArgumentParser(add_help=False)

        self.conformer_common_parser.add_argument(
            "SmilesString",
            type=str,
            help="SMILES string of the compound to be optimised."
        )
        self.conformer_common_parser.add_argument(
            "-f", "--FileName",
            type=str,
            default="confdft",
            help="File name to use for conformer_gen input files."
        )
        self.conformer_common_parser.add_argument(
            "-nc", "--NumCPUs",
            type=int,
            default=4,
            help="Number of CPUs allocated to the HPC calculation."
        )
        self.conformer_common_parser.add_argument(
            "-rt", "--RunTime",
            type=str,
            default="24:00:00",
            help="Maximum runtime allowed on the HPC system."
        )
        self.conformer_common_parser.add_argument(
            "-e", "--Email",
            type=str,
            help=(
                "Email address to receive updates on job start, finish, "
                "and failure."
            )
        )
        self.conformer_common_parser.add_argument(
            "-n", "--NumberConformers",
            type=int,
            default=100,
            help="Number of conformers for the initial RDKit search."
        )
        self.conformer_common_parser.add_argument(
            "-m", "--MaximumIterations",
            type=int,
            default=1000,
            help=(
                "Number of iterations for the initial RDKit optimisation "
                "step."
            )
        )

        self.parse_conformer_crest_input()
        self.parse_conformer_orca_input()

    def parse_conformer_crest_input(self) -> None:
        """Create parser for CREST conformer input mode."""
        crest_write = self.conformer_subparsers.add_parser(
            "crest",
            parents=[self.top_parent_parser, self.conformer_common_parser],
            help="CREST input mode."
        )
        crest_write.add_argument(
            "-ol", "--OptimisationLevel",
            type=str,
            default="extreme",
            choices=[
                "crude", "vloose", "loose", "normal",
                "tight", "vtight", "extreme"
            ],
            help="Convergence conditions for CREST optimisation steps."
        )
        crest_write.add_argument(
            "-me", "--CRESTMethod",
            type=str,
            default="gfn2",
            choices=["gfnff", "gfn0", "gfn1", "gfn2"],
            help=(
                "Tight binding method used by CREST for optimisation and "
                "conformer searching."
            )
        )
        crest_write.set_defaults(func=self.commands["conformer_gen"])

    def parse_conformer_orca_input(self) -> None:
        """Create parser for ORCA conformer input mode."""
        orca_write = self.conformer_subparsers.add_parser(
            "orca",
            parents=[self.top_parent_parser, self.conformer_common_parser],
            help="ORCA input mode using ORCA GOAT."
        )
        orca_write.add_argument(
            "-mo", "--GOATMode",
            type=str,
            default="GOAT",
            choices=["GOAT", "GOAT-ENTROPY", "GOAT-EXPLORE"],
            help="Specific mode for ORCA GOAT to run."
        )
        orca_write.add_argument(
            "-me", "--GOATMethod",
            type=str,
            default="gfn2",
            choices=["gfnff", "gfn0", "gfn1", "gfn2"],
            help=(
                "Tight binding method used by ORCA GOAT for optimisation "
                "and conformer searching."
            )
        )
        orca_write.set_defaults(func=self.commands["conformer_gen"])

    def parse_orca_input(self) -> None:
        """Create parser for ORCA input mode."""
        self.top_subparsers.add_parser(
            "orca_write",
            parents=[self.top_parent_parser],
            help="Generate ORCA input files from CREST outputs."
        )

    def load_from_config(self) -> None:
        """Load configuration parameters from a TOML file if provided."""
        try:
            with open(self.args.Config, "rb") as file:
                config = tomllib.load(file)
            print("Config file found.")
        except FileNotFoundError:
            print(
                f"Config file not found. Is \"{self.args.Config}\" a valid "
                "path? Using default parameters instead."
            )
            return

        flat_config = utils.flatten_dict(config.get(self.args.Command))

        # Load input args from command line except -c / --Config
        cmdline_args = [
            arg for arg in sys.argv[2:]
            if arg.startswith(("-", "--")) and arg not in ("-c", "--Config")
        ]

        is_nested = getattr(
            self.args, "Command", "conformer_gen"
        ) == "conformer_gen"

        parser_source = (
            self.conformer_subparsers if is_nested else self.top_subparsers
        )
        mode_key = self.args.ConformerMode if is_nested else self.args.Command
        subparser = parser_source.choices[mode_key]

        # Convert short arg names to long form
        for action in subparser._actions[2:]:
            arg_names = action.option_strings
            if len(arg_names) > 1:
                primary, secondary = arg_names[0], arg_names[1].strip("--")
                cmdline_args = [
                    secondary if item == primary else item
                    for item in cmdline_args
                ]

        # Remove keys from config that were already provided via CLI
        for key in cmdline_args:
            flat_config.pop(key, None)

        # Update args with remaining config values
        for key, value in flat_config.items():
            setattr(self.args, key, value)

        print("Config parameters loaded successfully.")

    def get_args(self) -> argparse.Namespace:
        """Parse and return command-line arguments."""
        self.args = self.parser.parse_args()
        return self.args
