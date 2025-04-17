# Confor2DFT

Confor2DFT is a Python tool that integrates CREST and DFT methods to optimize molecular structures. It generates conformers using RDKit and performs optimizations with CREST, writing all necessary input files for the process.

### Author: Dr. Jordan Hobbs  
University of Leeds  
[GitHub Repository](https://github.com/Jordan-Hobbs/)

---

## Features

- Generate multiple conformers from a molecule using RDKit.
- Optimize conformers with CREST.
- Create input files for CREST (`.xyz`, `.toml`, `.sh`).
- Optionally, prepare ORCA input files for quantum mechanical calculations.

---

## Requirements

To use Confor2DFT, you need:

- **Python 3.x**
- **RDKit** (for molecular manipulation and conformer generation)
- **CREST** (for conformer optimization)
- **ORCA** (for quantum mechanical calculations, optional)
- **tomllib** (for reading TOML configuration files)

You can install the required Python dependencies via `pip`:

```bash
pip install rdkit toml
