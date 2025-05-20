import numpy as np
from rdkit.Chem import Mol


class ORCAWriter:
    """Class to generate ORCA input and SLURM shell files."""

    def __init__(self, molecule: Mol, args) -> None:
        """Initialize ORCAWriter with a molecule and a file name."""
        print("\n----------------------------------------------------------------")
        print("Writing ORCA input files:\n")
        self.molecule = molecule
        self.args = args

    def write_inp(self) -> None:
        """Write the .inp file containing molecular coordinates."""
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 4)

        pal_block = f"%PAL\n nprocs {self.args.NumCPUs}\nEND\n"

        methods = {
            "gfnff": "XTBFF",
            "gfn0": "XTB0",
            "gfn1": "XTB1",
            "gfn2": "XTB2"
        }
        methods_text = (
            f"!{self.args.GOATMode} {methods.get(self.args.GOATMethod)}\n"
        )

        goat_block = "%GOAT\n ALIGN TRUE\n"
        if self.args.GOATUphill:
            uphill = {
                "gfnff": "GFNFF",
                "gfn0": "GFN0XTB",
                "gfn1": "GFN1XTB",
                "gfn2": "GFN2XTB"
            }
            goat_block += f" GFNUPHILL {uphill.get(self.args.GOATUphill)}\n"
        goat_block += "END\n"
        
        scf_lines = []
        include_GOATOptLevel = (
            self.args.GOATOptLevel is not None 
            and self.args.GOATOptLevel.lower() != "default"
        )

        include_convforced = self.args.ConvForced is not None
        if include_GOATOptLevel or include_convforced:
            scf_lines.append("%scf")
            if include_GOATOptLevel:
                scf_lines.append(f" Convergence {self.args.GOATOptLevel}")
            if include_convforced:
                scf_lines.append(f" ConvForced {1 if self.args.ConvForced is True else 0}")
            scf_lines.append("end")
        scf_block = "\n".join(scf_lines) + "\n"

        with open(f"{self.args.FileName}.inp", "w", newline="\n") as file:
            file.write(methods_text)
            file.write(pal_block)
            file.write(goat_block)
            if scf_block:
                file.write(scf_block)
            file.write("* xyz 0 1\n")
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x, y, z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")
            file.write("*")

        print(f"{self.args.FileName}.inp file written successfully")

    def write_sh(self) -> None:
        """Write the SLURM shell script to run ORCA jobs."""
        sh_text = (
            "#!/bin/bash\n"
            f"#SBATCH --job-name={self.args.FileName}\n"
            f"#SBATCH --time={self.args.RunTime}\n"
            f"#SBATCH --ntasks={self.args.NumCPUs}\n"
            "#SBATCH --mem-per-cpu=1G\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --mail-type=BEGIN,END,FAIL\n"
            f"#SBATCH --mail-user={self.args.Email}\n"
            "\n"
            "module load orca\n"
            "module load crest\n"
            "module load openmpi\n"
            "export RSH_COMMAND=\"/usr/bin/ssh\"\n"
            "\n"
            f"INPUT_FILE=\"{self.args.FileName}\"\n"
            "LAUNCH_DIR=$PWD\n"
            "WORKING_DIR=\"/mnt/flash/tmp/job.$SLURM_JOB_ID\"\n"
            "mkdir -p \"${WORKING_DIR}\"\n"
            "cp \"${INPUT_FILE}.inp\" \"${WORKING_DIR}/${INPUT_FILE}.inp\"\n"
            "\n"
            "cd ${WORKING_DIR}\n"
            "\n"
            "/opt/apps/pkg/applications/orca/orca_6_0_1_linux_x86-64_shared_openmpi416/orca \"${WORKING_DIR}/${INPUT_FILE}.inp\" > \"${LAUNCH_DIR}/${INPUT_FILE}.out\"\n"
            "wait\n"
            "\n"
            "sed -i \"s/ converged=true//\" \"${INPUT_FILE}.finalensemble.xyz\"\n"
            "\n"
            "crest ${INPUT_FILE}.globalminimum.xyz --cregen ${INPUT_FILE}.finalensemble.xyz\n"
            "wait\n"
            "mv \"${INPUT_FILE}.globalminimum.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.globalminimum.xyz\"\n"
            "mv \"${INPUT_FILE}.finalensemble.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.finalensemble.xyz\"\n"
            "\n"
            "mv \"${INPUT_FILE}.finalensemble.xyz.sorted\" \"${LAUNCH_DIR}/${INPUT_FILE}.sortedensemble.xyz\"\n"
            "mv \"crest_ensemble.xyz\" \"${LAUNCH_DIR}/${INPUT_FILE}.totalensemble.xyz\"\n"
            "mv \"crest.energies\" \"${LAUNCH_DIR}/${INPUT_FILE}.energies\""
        )

        with open(f"{self.args.FileName}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.args.FileName}.sh file written successfully")


class CRESTWriter:
    """Class to generate CREST input, TOML, and SLURM shell files."""

    def __init__(self, molecule: Mol, args) -> None:
        """Initialize CRESTWriter with a molecule and a file name."""
        print("\n----------------------------------------------------------------")
        print("Writing CREST input files:\n")
        self.molecule = molecule
        self.args = args

    def write_xyz(self) -> None:
        """Write the .xyz file containing molecular coordinates."""
        xyz_pos = self.molecule.GetConformer().GetPositions()
        xyz_pos = np.round(xyz_pos, 4)
        no_atoms = self.molecule.GetConformer().GetNumAtoms()

        header_text = (
            f"{no_atoms}\n"
            f"{self.file_name}\n"
        )

        with open(f"{self.args.FileName}.xyz", "w", newline="\n") as file:
            file.write(header_text)
            for index, atom in enumerate(self.molecule.GetAtoms()):
                atom_symbol = atom.GetSymbol()
                x, y, z = xyz_pos[index]
                file.write(f"{atom_symbol} {x} {y} {z}\n")

        print(f"{self.args.FileName}.xyz file written successfully")

    def write_toml(self) -> None:
        """Write TOML files for CREST runs (optimization and conformer generation)."""
        run_types = ("opt", "conf")

        for run in run_types:
            run_type = "ancopt" if run == "opt" else "imtd-gc"

            toml_text = (
                "# CREST 3 input file"
                f"input = \"{self.args.FileName}.xyz\""
                f"runtype = \"{run_type}\""
                f"threads = {self.args.NumCPUs}"
                ""
                "[calculation]"
                f"optlev = \"{self.args.CRESTOptLevel}\""
                ""
                "[[calculation.level]]"
                f"method = \"{self.args.CRESTMethod}\""
            )

            toml_name = f"{self.args.FileName}_{run}.toml"
            with open(toml_name, "w", newline="\n") as file:
                file.write(toml_text)

            print(f"{toml_name} file written successfully")

    def write_sh(self) -> None:
        """Write the SLURM shell script to run CREST jobs."""
        sh_text = (
            "#!/bin/bash\n"
            f"#SBATCH --job-name={self.args.FileName}\n"
            f"#SBATCH --output={self.args.FileName}.out\n"
            f"#SBATCH --time={self.args.RunTime}\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --mem-per-cpu=1G\n"
            f"#SBATCH --cpus-per-task={self.args.NumCPUs}\n"
            "#SBATCH --mail-type=BEGIN,END,FAIL\n"
            f"#SBATCH --mail-user={self.args.Email}\n"
            "\n"
            "module load crest\n"
            "\n"
            f"crest {self.args.FileName}_opt.toml\n"
            "wait\n"
            "\n"
            f"find . -type f ! -name \"{self.args.FileName}.sh\" "
            f"! -name \"crestopt.xyz\" "
            f"! -name \"{self.args.FileName}_opt.toml\" "
            f"! -name \"{self.args.FileName}_conf.toml\" "
            f"! -name \"{self.args.FileName}.out\" -delete\n"
            "\n"
            f"mv \"crestopt.xyz\" \"{self.args.FileName}.xyz\"\n"
            "\n"
            f"crest {self.args.FileName}_conf.toml\n"
            "wait\n"
        )

        with open(f"{self.args.FileName}.sh", "w", newline="\n") as file:
            file.write(sh_text)

        print(f"{self.args.FileName}.sh file written successfully")