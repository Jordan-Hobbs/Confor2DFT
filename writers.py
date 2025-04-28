import numpy as np
import textwrap
from rdkit.Chem import Mol


class ORCAWriter:
    """Class to generate ORCA input and SLURM shell files."""

    def __init__(self, molecule: Mol, file_name: str) -> None:
        """Initialize ORCAWriter with a molecule and a file name."""
        print("\n----------------------------------------------------------------")
        print("Writing ORCA input files:\n")
        self.molecule = molecule
        self.file_name = file_name

    def write_inp(self, mode: str, method: str) -> None:
        """Write the .inp file containing molecular coordinates."""
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

    def write_sh(self, run_time: str = "24:00:00", num_cpus: int = 4, email: str = "j.l.hobbs@leeds.ac.uk") -> None:
        """Write the SLURM shell script to run ORCA jobs."""
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
            "INPUT_FILE=\"{self.file_name}\"\n"
            "LAUNCH_DIR=$PWD\n"
            "WORKING_DIR=\"/mnt/flash/tmp/job.$SLURM_JOB_ID\"\n"
            "mkdir -p \"${WORKING_DIR}\"\n"
            "cp \"${LAUNCH_DIR}/${INPUT_FILE}.inp\" \"${WORKING_DIR}/${INPUT_FILE}.inp\"\n"
            "\n"
            "orca \"${WORKING_DIR}/${INPUT_FILE}.inp\" > \"${LAUNCH_DIR}/${INPUT_FILE}.out\"\n"
            "wait\n"
            "\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.gbw\" \"${LAUNCH_DIR}/${INPUT_FILE}.gbw\"\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.globalminimum.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.globalminimum.xyz\"\n"
            "mv \"${WORKING_DIR}/${INPUT_FILE}.finalensemble.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.finalensemble.xyz\"\n"
            "\n"
            "rm -rf \"${WORKING_DIR}\"\n"
        )

        with open(f"{self.file_name}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.file_name}.sh file written successfully")


class CRESTWriter:
    """Class to generate CREST input, TOML, and SLURM shell files."""

    def __init__(self, molecule: Mol, file_name: str) -> None:
        """Initialize CRESTWriter with a molecule and a file name."""
        print("\n----------------------------------------------------------------")
        print("Writing CREST input files:\n")
        self.molecule = molecule
        self.file_name = file_name

    def write_xyz(self) -> None:
        """Write the .xyz file containing molecular coordinates."""
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

    def write_toml(self, num_cpus: int = 4, optlevel: str = "extreme", gfn_method: str = "gfn2") -> None:
        """Write TOML files for CREST runs (optimization and conformer generation)."""
        run_types = ("opt", "conf")

        for run in run_types:
            run_type = "ancopt" if run == "opt" else "imtd-gc"

            toml_text = (
                "# CREST 3 input file"
                f"input = \"{self.file_name}.xyz\""
                f"runtype = \"{run_type}\""
                f"threads = {num_cpus}"
                ""
                "[calculation]"
                f"optlev = \"{optlevel}\""
                ""
                "[[calculation.level]]"
                f"method = \"{gfn_method}\""
            )

            toml_name = f"{self.file_name}_{run}.toml"
            with open(toml_name, "w", newline="\n") as file:
                file.write(textwrap.dedent(toml_text))

            print(f"{toml_name} file written successfully")

    def write_sh(self, run_time: str = "24:00:00", num_cpus: int = 4, email: str = "j.l.hobbs@leeds.ac.uk") -> None:
        """Write the SLURM shell script to run CREST jobs."""
        sh_text = (
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
            f"crest {self.file_name}_opt.toml\n"
            "wait\n"
            "\n"
            f"find . -type f ! -name \"{self.file_name}.sh\" "
            f"! -name \"crestopt.xyz\" "
            f"! -name \"{self.file_name}_opt.toml\" "
            f"! -name \"{self.file_name}_conf.toml\" "
            f"! -name \"{self.file_name}.out\" -delete\n"
            "\n"
            f"mv \"crestopt.xyz\" \"{self.file_name}.xyz\"\n"
            "\n"
            f"crest {self.file_name}_conf.toml\n"
            "wait\n"
        )

        with open(f"{self.file_name}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.file_name}.sh file written successfully")