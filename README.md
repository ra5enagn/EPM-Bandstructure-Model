# EPM Bandstructure Model

## Overview
This repository contains a MATLAB implementation of the **Empirical Pseudopotential Method (EPM)** for computing the electronic bandstructure of crystalline semiconductors. The code evaluates energy–momentum (E–k) relationships by solving the Schrödinger equation in reciprocal space using empirically fitted pseudopotentials.

---

## Methodology
The Empirical Pseudopotential Method replaces the true crystal potential with a simplified, empirically fitted pseudopotential that reproduces experimental bandstructure data. In this implementation:

- The crystal potential is represented in reciprocal space
- Plane-wave basis functions are used to expand the electron wavefunction
- The Hamiltonian matrix is constructed using empirical form factors
- Eigenvalues are computed to obtain the electronic energy bands

This approach enables efficient bandstructure calculation while attempting to retain physical accuracy.

---

## Implementation Details
- Language: MATLAB
- Numerical approach: Reciprocal-space Schrödinger equation
- Output: Energy bands as a function of crystal momentum

The code is structured as a standalone script and can be easily extended to support multiple materials, additional k-points, or alternative pseudopotential parameter sets.

---

## Usage
1. Open MATLAB
2. Navigate to the repository directory
3. Run the main script:
   ```matlab
   epm_bandstructure
