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
                "Select operating mode: \"conformation_gen\" or \"conformation_sort\". "
                "Typical workflow: \"conformation_gen\" first, then "
                "\"conformation_sort\"."
            )
        )

        self.parse_conf_gen_input()
        self.parse_sort_input()

    def parse_conf_gen_input(self) -> None:
        """Create parser for conformer input mode."""
        self.conf_gen = self.top_subparsers.add_parser(
            "conformer_gen",
            parents=[self.top_parent_parser],
            help="Generate input files for a conformer search."
        )
        self.conf_gen_subparsers = self.conf_gen.add_subparsers(
            required=True,
            dest="ConformerMode",
            help="Specify which package to write input files for."
        )
        # Create shared args parser for conformer mode
        self.conf_gen_common_parser = argparse.ArgumentParser(add_help=False)

        self.conf_gen_common_parser.add_argument(
            "SmilesString",
            type=str,
            help="SMILES string of the compound to be optimised."
        )
        self.conf_gen_common_parser.add_argument(
            "-f", "--FileName",
            type=str,
            default="confdft",
            help="File name to use for conformer_gen input files."
        )
        self.conf_gen_common_parser.add_argument(
            "-nc", "--NumCPUs",
            type=int,
            default=4,
            help="Number of CPUs allocated to the HPC calculation."
        )
        self.conf_gen_common_parser.add_argument(
            "-t", "--RunTime",
            type=str,
            default="24:00:00",
            help="Maximum runtime allowed on the HPC system."
        )
        self.conf_gen_common_parser.add_argument(
            "-e", "--Email",
            type=str,
            help=(
                "Email address to receive updates on job start, finish, "
                "and failure."
            )
        )
        self.conf_gen_common_parser.add_argument(
            "-rn", "--RDNumConf",
            type=int,
            default=100,
            help="Number of conformers for the initial RDKit search."
        )
        self.conf_gen_common_parser.add_argument(
            "-rm", "--RDMaxIter",
            type=int,
            default=1000,
            help=(
                "Number of iterations for the initial RDKit optimisation "
                "step."
            )
        )
        self.parse_conf_gen_crest_input()
        self.parse_conformer_orca_input()

    def parse_conf_gen_crest_input(self) -> None:
        """Create parser for CREST conformer input mode."""
        crest_write = self.conf_gen_subparsers.add_parser(
            "crest",
            parents=[self.top_parent_parser, self.conf_gen_common_parser],
            help="CREST output mode."
        )
        crest_write.add_argument(
            "-cl", "--CRESTOptLevel",
            type=str,
            default="extreme",
            choices=[
                "crude", "vloose", "loose", "normal",
                "tight", "vtight", "extreme"
            ],
            help="Convergence conditions for CREST optimisation steps."
        )
        crest_write.add_argument(
            "-cm", "--CRESTMethod",
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
        orca_write = self.conf_gen_subparsers.add_parser(
            "orca",
            parents=[self.top_parent_parser, self.conf_gen_common_parser],
            help="ORCA output mode for using ORCA GOAT."
        )
        orca_write.add_argument(
            "-omo", "--GOATMode",
            type=str,
            default="GOAT",
            choices=["GOAT", "GOAT-ENTROPY", "GOAT-EXPLORE"],
            help="Specific mode for ORCA GOAT to run."
        )
        orca_write.add_argument(
            "-om", "--GOATMethod",
            type=str,
            default="gfn2",
            choices=["gfnff", "gfn0", "gfn1", "gfn2"],
            help=(
                "Tight binding method used by ORCA GOAT for optimisation "
                "and conformer searching."
            )
        )
        orca_write.add_argument(
            "-cf", "--ConvForced",
            type=str2bool,
            choices=[True, False],
            help=(
                "Specifies whether ORCA can proceed with the GOAT calulation "
                "if the initial scf optimisation fails. \"false\" specifies it can "
                "continue while \"true\" specifies it cant."
            )
        )
        orca_write.add_argument(
            "-ol", "--GOATOptLevel",
            type=str,
            default="Default",
            choices=[
                "Sloppy", "Loose",  "Medium", "Default", "Strong",
                "Tight", "VeryTight", "Extreme"
            ],
            help=(
                "Convergence conditions for ORCA optimisation steps. More "
                "details provided in the ORCA manual under \"SCF Convergence\""
            )
        )
        orca_write.set_defaults(func=self.commands["conformer_gen"])

    def parse_sort_input(self) -> None:
        """Create parser for ORCA input mode."""
        self.conf_sort = self.top_subparsers.add_parser(
            "conformer_sort",
            parents=[self.top_parent_parser],
            help=(
                "Generate optimisation input files for ORCA based on the "
                "results of an CREST or ORCA GOAT conformation analysis "
                "calculation."
            )
        )
        self.conf_sort_subparsers = self.conf_sort.add_subparsers(
            required=True,
            dest="InputMode",
            help="Specify which package to read input files from."
        )
        self.conf_sort_common_parser = argparse.ArgumentParser(add_help=False)
        self.conf_sort_common_parser.add_argument(
            "InputFile",
            type=str,
            help="Name of the input file to be optimised and refined."
        )
        self.conf_sort_common_parser.add_argument(
            "-f", "--FileName",
            type=str,
            default="confdft",
            help="File name to use for conformer_gen input files."
        )
        self.conf_sort_common_parser.add_argument(
            "-nc", "--NumCPUs",
            type=int,
            default=4,
            help="Number of CPUs allocated to the HPC calculation."
        )
        self.conf_sort_common_parser.add_argument(
            "-rt", "--RunTime",
            type=str,
            default="24:00:00",
            help="Maximum runtime allowed on the HPC system."
        )
        self.conf_sort_common_parser.add_argument(
            "-e", "--Email",
            type=str,
            help=(
                "Email address to receive updates on job start, finish, "
                "and failure."
            )
        )
        self.parse_conf_sort_crest_input()

    def parse_conf_sort_crest_input(self) -> None:
        crest_write = self.conf_sort_subparsers.add_parser(
            "crest",
            parents=[self.top_parent_parser, self.conf_sort_common_parser],
            help="CREST input mode."
        )
        crest_write.add_argument







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
            self.conf_gen_subparsers if is_nested else self.top_subparsers
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

def str2bool(string):
    if isinstance(string, bool):
        return string
    if string.lower() in ('true', 't', 'yes', '1'):
        return True
    elif string.lower() in ('false', 'f', 'no', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected (true/false).")