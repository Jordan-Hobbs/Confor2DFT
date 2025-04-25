import numpy as np
import textwrap


class CRESTWriter:
    def __init__(self, molecule, file_name):
        print("\n----------------------------------------------------------------")
        print("Writing input files:\n")

        self.molecule = molecule
        self.file_name = file_name

        self.write_xyz()
        self.write_toml()
        self.write_sh()

    def write_xyz(self):
        """Write the XYZ file containing molecular coordinates."""
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 4)
        no_atoms = self.molecule.GetConformer().GetNumAtoms()

        with open(f"{self.file_name}.xyz", "w", newline="\n") as file:
            file.write(f"{no_atoms}\n{self.file_name}\n")
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x, y, z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")

        print(f"{self.file_name}.xyz file written successfully")

    def write_toml(self, num_cpus=4, optlevel="extreme", gfn_method="gfn2"):
        """Write TOML files for CREST runs (opt and conf)."""
        run_types = ("opt", "conf")

        for run in run_types:
            run_type = "ancopt" if run == "opt" else "imtd-gc"

            toml_text = f"""\
                # CREST 3 input file
                input = "{self.file_name}.xyz"
                runtype = "{run_type}"
                threads = {num_cpus}

                [calculation]
                optlev = "{optlevel}"

                [[calculation.level]]
                method = "{gfn_method}"
            """

            toml_name = f"{self.file_name}_{run}.toml"
            with open(toml_name, "w", newline="\n") as file:
                file.write(textwrap.dedent(toml_text))

            print(f"{toml_name} file written successfully")

    def write_sh(self, run_time="24:00:00", num_cpus=4, email="j.l.hobbs@leeds.ac.uk"):
        """Write the SLURM shell script to run CREST jobs."""
        sh_text =(
            "#!/bin/bash\n"
            f"#SBATCH --job-name={self.file_name}\n"
            f"#SBATCH --output={self.file_name}.out\n"
            f"#SBATCH --time={run_time}\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --mem-per-cpu=1G\n"
            f"#SBATCH --cpus-per-task={num_cpus}\n"
            "#SBATCH --mail-type=BEGIN,END,FAIL\n"
            f"#SBATCH --mail-user={email}\n"
            "\n"
            "module load crest\n"
            "\n"
            "# =============================================\n"
            "# Optimise initial structure\n"
            "# =============================================\n"
            "\n"
            f"crest {self.file_name}_opt.toml\n"
            "wait\n"
            "\n"
            "# =============================================\n"
            "# Delete files unneeded for Conf gen\n"
            "# =============================================\n"
            "\n"
            f"find . -type f ! -name \"{self.file_name}.sh\" ! -name \"crestopt.xyz\" ! -name \"{self.file_name}_opt.toml\" ! -name \"{self.file_name}_conf.toml\" ! -name \"{self.file_name}.out\" -delete\n"
            "\n"
            f"mv \"crestopt.xyz\" \"{self.file_name}.xyz\"\n"
            "\n"
            "# =============================================\n"
            "# Generate and optimise conformer ensemble\n"
            "# =============================================\n"
            "\n"
            f"crest {self.file_name}_conf.toml\n"
            "wait\n"
        )

        with open(f"{self.file_name}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.file_name}.sh file written successfully")