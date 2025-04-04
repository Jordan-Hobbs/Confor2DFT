import argparse


class ProgramController:
    def __init__(self):
        ## Create parser for input arguments
        self.parser = argparse.ArgumentParser(
            prog = "Confor2DFT",
            description = """
            Combines CREST and some DFT package to obtain optimised molecules
            """
        )
        self.subparsers = self.parser.add_subparsers(
            help = """
            Operating mode, options are:\n
            CREST_Input - takes a smiles string and generators and initial 
            input structure for optimisation, and then conformational searching 
            in the CREST program.
            """
        )
        ## Create parser for CREST_Input mode
        
        

