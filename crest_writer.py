import numpy as np


class CRESTGen:
    def __init__(self, file_name, molecule):
        self.molecule = molecule
        self.xyz_file = f"{file_name}.xyz"
        self.toml_file = f"{file_name}.toml"
        self.sh_file = f"{file_name}.sh"
    

    def write_xyz(self,
        comment: str = ""
        ):
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 3)
        no_atoms = self.molecule.GetConformer().GetNumAtoms()

        with open(self.xyz_file, "w", newline="\n") as file:
            file.write(f"{no_atoms}\n")
            file.write(f"{comment}\n")
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x,y,z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")
        
        print(f"{self.xyz_file} file written successfully")


    def write_toml(self, 
        runtype: str, 
        threads = 4, 
        optlevel = "extreme", 
        gfn_method = "gfn2"
        ):

        with open(self.toml_file, "w", newline="\n") as file:
            file.write("# CREST 3 input file\n")
            file.write(f"input = \"{self.xyz_file}\"\n")
            file.write(f"runtype = \"{runtype}\"\n")
            file.write(f"threads = {threads}\n")
            file.write("\n")
            file.write("[calculation]\n")
            file.write(f"optlev = \"{optlevel}\"\n")
            file.write("\n")
            file.write("[[calculation.level]]\n")
            file.write(f"method = \"{gfn_method}\"")

        print(f"{self.toml_file} file written successfully")
    

    def write_sh(self,
        run_time = "48:00:00",
        num_cpus = 4,
        email = "j.l.hobbs@leeds.ac.uk"
        ):
        with open(self.sh_file, "w", newline="\n") as file:
            file.write("#!\\bin\\bash\n")
            file.write(f"#SBATCH --time={run_time}\n")
            file.write("#SBATCH --ntasks=1\n")
            file.write("#SBATCH --mem-per-cpu=1G\n")
            file.write(f"#SBATCH --cpus-per-task={num_cpus}\n")
            file.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
            file.write(f"#SBATCH --mail-user={email}\n")
            file.write("\n")
            file.write("module load crest\n")
            file.write(f"crest {self.toml_file}")

        print(f"{self.sh_file} file written successfully")