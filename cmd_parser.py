import argparse
import configparser

class ProgramController:
    def __init__(self):

        config = configparser.ConfigParser()
        config.read("input_config.ini")

        ## Creat parent parser for arguments shared across subparsers
        self.parent_parser = argparse.ArgumentParser(add_help = False)
        self.parent_parser.add_argument("-c", "--config", type = str,
            default = "config_input.ini",
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
            required = True,
            help = "Select operating mode, options are crest_write and "
            "orca_write. In normal operation crest_write acts as the first "
            "stage while orca_write acts as second stage."
        )

        ## Create parameters for subparsers. Doesnt need to be in seperate
        ## functions but it helps keep it tidy
        self.parse_crest_input()
        self.parse_orca_input()

    def parse_crest_input(self):
        """
        Create parser for CREST input mode
        """
        crest_write = self.subparsers.add_parser("crest_write",
            parents = [self.parent_parser],
            help = "Needs adding"
        )
        crest_write.add_argument("s", type = str, 
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

    def load_config(self, config_file):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)



    def get_args(self):
        """
        Small function for returning arges from command line
        """
        return self.parser.parse_args()