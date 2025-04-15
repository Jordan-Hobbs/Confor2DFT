import argparse
import tomllib
import utils
import sys

class ProgramController:
    def __init__(self):

        ## Create parent parser for arguments shared across subparsers
        self.parent_parser = argparse.ArgumentParser(add_help = False)
        self.parent_parser.add_argument("-c", "--Config", type = str,
            help = "Name of the config file used to input required settings"
        )
        ## Create parser and subparsers for input arguments
        self.parser = argparse.ArgumentParser(
            prog = "Confor2DFT",
            parents = [self.parent_parser],
            description = "Combines CREST and some DFT package to obtain " 
            "optimised molecules."
        )
        self.subparsers = self.parser.add_subparsers(
            required = True, dest = "Command",
            help = "Select operating mode, options are crest_write and "
            "orca_write. In a normal workflow operation crest_write acts as "
            "the first stage while orca_write acts as the second stage."
        )
        ## Create parameters for subparsers. Doesnt need to be in seperate
        ## functions but it helps keep it tidy. Idk if this a "good" thing to 
        ## do but it works for me.
        self.parse_crest_input()
        self.parse_orca_input()

        self.subparser_map = self.subparsers.choices


    def parse_crest_input(self):
        """
        Create parser for CREST input mode
        """
        crest_write = self.subparsers.add_parser("crest_write",
            parents = [self.parent_parser],
            help = "Needs adding"
        )
        crest_write.add_argument("SmilesString", type = str, 
            help = "Smiles string of the compound to be optimised."                   
        )
        crest_write.add_argument("-f", "--FileName", type = str, 
            default = "conf2dft", 
            help = "File name to use for CREST input files."
        )
        crest_write.add_argument("-n", "--NumberConformers", type = int, 
            default = 100, 
            help = "Number of conformers for the initial RDKit search to "
            "generate."
        )
        crest_write.add_argument("-m", "--MaximumIterations", type = int, 
            default = 1000,
            help = "Numbers of iterations for the initial RDKit otimisation "
            "step for the RDKit generator conformers."
        )
        crest_write.add_argument("-ol", "--OptimisationLevel", type = str, 
            default = "extreme", choices = ["crude", "vloose", "loose", 
            "normal", "tight", "vtight", "extreme"],
            help = "Sets convergence conditions for the CREST optimisation "
            "steps."
        )
        crest_write.add_argument("-xtb", "--GFNMethod", type = str, 
            default = "gfn2", choices = ["gfnff", "gfn0", "gfn1", "gfn2"],
            help = "Select tight binding method used by CREST optimisation and "
            "conformer searching."
        )
        crest_write.add_argument("-nc", "--NumCPUs", type = int, default = 4,
            help = "Number of CPUs to be allocated to the the HPC calulation."
        )        
        crest_write.add_argument("-rt", "--RunTime", type = str, 
            default = "24:00:00",
            help = "Maximum time allowed for running the reuslting set of "
            "files on the HPC system."
        )
        crest_write.add_argument("-e", "--Email", type = str,
            help = "Provide an email to be kept updated of start, finish and "
            "fail times."                         
        )
    
    def parse_orca_input(self):
        """
        Create parser for ORCA input mode
        """ 
        orca_write = self.subparsers.add_parser("orca_write",
            parents = [self.parent_parser],
            help = "Needs adding."
        )
        orca_write.add_argument("-f", type = str, default = "conf2dft", 
            help = "File name to use for CREST input files."
        )

    def load_from_config(self):
        """
        """

        with open(self.args.Config, "rb") as file:
            config = tomllib.load(file)

        flat_config = utils.flatten_dict(
            config.get(self.args.Command))

        ## Replace values in args with config file values
        print(self.args)
        for key, value in flat_config.items():
            setattr(self.args, key, value)

        print(self.args)
        ## Load args specified on cmd line
        cmd_line = sys.argv[2:]
        cmdline_args = [arg for arg in cmd_line if arg.startswith(("-", "--")) and arg not in ("-c", "--Config")]
        

        ## Change cmdline_ars from - values to -- values
        subparser = self.subparsers.choices[self.args.Command]
        for action in subparser._actions[2:]:
            codes = (action.option_strings)
            if len(codes) > 1:
                cmdline_args = [codes[1].strip("--") if item == codes[0] else item for item in cmdline_args]

        ## Load values specified on cmd line
        cmdline_vals = [cmd_line[i + 1] for i in range(len(cmd_line) - 1) if cmd_line[i].startswith(("-", "--")) and cmd_line[i] not in ("-c", "--Config")]

        ## Replace values in args with cmd line values
        cmd_args = dict(zip(cmdline_args, cmdline_vals))

        for key, value in cmd_args.items():
            setattr(self.args, key, value)
        print(self.args)

    def get_args(self):
        """
        Small function for returning arges from command line
        """
        self.args = self.parser.parse_args()
        return self.args
