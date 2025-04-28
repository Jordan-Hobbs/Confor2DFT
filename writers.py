import numpy as np
import textwrap



class ORCAWriter:
    def __init__(self, molecule, file_name):
        print("\n----------------------------------------------------------------")
        print("Writing ORCA input files:\n")
        self.molecule = molecule
        self.file_name = file_name

    def write_inp(self, mode, method):
        """Write the inp file containing molecular coordinates."""
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 4)

        methods = {
            "gfnff": "XTBFF", 
            "gfn0": "XTB0", 
            "gfn1": "XTB1", 
            "gfn2": "XTB2"
        }

        header_text = (
            f"!{mode} {methods.get(method)}\n"
            "* xyz 0 1\n"
        )

        with open(f"{self.file_name}.inp", "w", newline="\n") as file:
            file.write(header_text)
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x, y, z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")

        print(f"{self.file_name}.inp file written successfully")

    def write_sh(self, run_time="24:00:00", num_cpus=4, email="j.l.hobbs@leeds.ac.uk"):

        sh_text = (
            "#!/bin/bash\n"
            f"#SBATCH --job-name={self.file_name}\n"
            f"#SBATCH --time={run_time}\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --mem-per-cpu=1G\n"
            f"#SBATCH --cpus-per-task={num_cpus}\n"
            "#SBATCH --mail-type=BEGIN,END,FAIL\n"
            f"#SBATCH --mail-user={email}\n"
            "\n"
            "module load orca\n"
            "\n"
            "# =============================================\n"
            "# Write folders and set variables for script\n"
            "# =============================================\n"
            "\n"
            f"INPUT_FILE =\"{self.file_name}\"\n"
            "LAUNCH_DIR=$PWD\n"
            "WORKING_DIR=\"/mnt/flash/tmp/job.$SLURM_JOB_ID\"\n"
            "mkdir -p \"${WORKING_DIR}\"\n"
            "cp \"${LAUNCH_DIR}/${INPUT_FILE}.inp\" \"${WORKING_DIR}/${INPUT_FILE}.inp\"\n"
            "\n"
            "# =============================================\n"
            "# Run ORCA calulations\n"
            "# =============================================\n"
            "\n"
            "orca \"${WORKING_DIR}/${INPUT_FILE}.inp\" > \"${LAUNCH_DIR}/${INPUT_FILE}.out\"\n"
            "wait\n"
            "\n"
            "# =============================================\n"
            "# Copy required files from the working directory\n"
            "# =============================================\n"
            "\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.gbw\" \"${LAUNCH_DIR}/${INPUT_FILE}.gbw\"\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.globalminimum.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.globalminimum.xyz\"\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.finalensemble.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.finalensemble.xyz\"\n"
            "\n"
            "rm -rf \"${WORKING_DIR}\""
        )

        with open(f"{self.file_name}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.file_name}.sh file written successfully")


class CRESTWriter:
    def __init__(self, molecule, file_name):
        print("\n----------------------------------------------------------------")
        print("Writing CREST input files:\n")
        self.molecule = molecule
        self.file_name = file_name

    def write_xyz(self):
        """Write the XYZ file containing molecular coordinates."""
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 4)
        no_atoms = self.molecule.GetConformer().GetNumAtoms()

        header_text = (
            f"{no_atoms}\n"
            f"{self.file_name}\n"
        )

        with open(f"{self.file_name}.xyz", "w", newline="\n") as file:
            file.write(header_text)
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