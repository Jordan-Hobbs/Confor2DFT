# Confor2DFT

Confor2DFT is a Python tool that integrates CREST and DFT methods to optimize molecular structures. It generates an initiale structure guess using RDKit to generate the structure. Various options for then optimising and checking the conformational landscape are available. This program is generally set up to work on a SLURM based HPC system, specifically AIRE at the University of Leeds, UK.

### Author: Dr. Jordan Hobbs  
University of Leeds  
[GitHub Repository](https://github.com/Jordan-Hobbs/)

---

## Requirements

To use Confor2DFT, you need:

- **Python 3.x**
- **RDKit** (for molecular manipulation and conformer generation)
- Some ability to run the calulations generated for CREST or ORCA
