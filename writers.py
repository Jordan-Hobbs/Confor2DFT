import numpy as np
import textwrap


class CRESTWriter:
    def __init__(self, file_name, molecule):
        self.molecule = molecule
        self.file_name = file_name
    

    def write_xyz(self,
        comment: str = "Confor2DFT"
        ):

        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 3)
        no_atoms = self.molecule.GetConformer().GetNumAtoms()

        with open(f"{self.file_name}.xyz", "w", newline="\n") as file:
            file.write(f"{no_atoms}\n{comment}\n")
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x,y,z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")
        
        print(f"{self.file_name}.xyz file written successfully")


    def write_toml(self, 
        threads = 4, 
        optlevel = "extreme", 
        gfn_method = "gfn2"
        ):

        run_types = ("opt", "conf")

        for run in run_types:
            if run == "opt":
                run_type = "ancopt"
            elif run == "conf":
                run_type = "imtd-gc"
            toml_text = f"""                # CREST 3 input file
                input = \"{self.file_name}.xyz\"
                runtype = \"{run_type}\"
                threads = {threads}

                [calculation]
                optlev = \"{optlevel}\"

                [[calculation.level]]
                method = \"{gfn_method}\""""

            toml_name = f"{self.file_name}_{run}.toml"

            with open(toml_name, "w", newline="\n") as file:
                file.write(textwrap.dedent(toml_text))

            print(f"{toml_name} file written successfully")
    

    def write_sh(self,
        run_time = "48:00:00",
        num_cpus = 4,
        email = "j.l.hobbs@leeds.ac.uk"
        ):
        
        sh_text = f"""        #!/bin/bash
        #SBATCH --job-name={self.file_name}
        #SBATCH --output={self.file_name}.out
        #SBATCH --time={run_time}
        #SBATCH --ntasks=1
        #SBATCH --mem-per-cpu=1G
        #SBATCH --cpus-per-task={num_cpus}
        #SBATCH --mail-type=BEGIN,END,FAIL
        #SBATCH --mail-user={email}

        module load crest

        # =============================================
        # Optmise initial structure
        # =============================================

        crest {self.file_name}_opt.toml
        wait

        # =============================================
        # Manage files
        # =============================================

        find . -type f ! -name "{self.file_name}.sh" ! -name "crestopt.xyz" ! -name "{self.file_name}_opt.toml" ! -name "{self.file_name}_conf.toml" ! -name "{self.file_name}.out" -delete
        mv "crestopt.xyz" "{self.file_name}.xyz"

        # =============================================
        # Generate and optimise conformer ensemble
        # =============================================

        crest {self.file_name}_conf.toml
        wait"""

        with open(f"{self.file_name}.sh", "w", newline="\n") as file:
            file.write(textwrap.dedent(sh_text))

        print(f"{self.file_name}.sh file written successfully")