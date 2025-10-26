# SCFEnergyComputations

Jupyter notebooks using **PySCF** to perform **Self-Consistent Field (SCF)** energy calculations on molecules. Includes quick demos of **Koopmans’ theorem** for rough ionization potential (IP) and electron affinity (EA) estimates from orbital energies.

## Overview

* Run HF/DFT SCF and report total electronic energies.
* Compare computed energies vs. reference/experimental values.
* Use Koopmans (IP ≈ −ε_HOMO, EA ≈ −ε_LUMO) and note its limitations.
* Explore how **geometry changes** affect energies.

## Install

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install numpy matplotlib pyscf jupyter
```

## Quickstart

```bash
git clone https://github.com/yourusername/SCFEnergyComputations.git
cd SCFEnergyComputations
jupyter notebook
```

Open `UE06.ipynb` and run all cells.

## Minimal PySCF Example

```python
import numpy as np
from pyscf import gto, scf

mol = gto.Mole(atom="O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587", basis="sto-3g")
mol.build()

mf = scf.RHF(mol).run()                 # SCF energy
E_tot = mf.e_tot

eps = np.array(mf.mo_energy)            # orbital energies (Hartree)
nocc = mol.nelec[0] + mol.nelec[1]      # closed-shell: equals number of electrons / 2
eps_homo, eps_lumo = eps[nocc-1], eps[nocc]

hartree_to_ev = 27.211386245988
IP_koop = -eps_homo * hartree_to_ev     # eV
EA_koop = -eps_lumo * hartree_to_ev     # eV

print("E_total (Hartree):", E_tot)
print("Koopmans IP (eV):", IP_koop)
print("Koopmans EA (eV):", EA_koop)
```

## Contents

* **UE06.ipynb** — SCF energy workflow + Koopmans estimates
* (Optional) `data/` — geometries or reference values
* (Optional) `scripts/` — small helpers for scans/plots

## Notes & Caveats

* Koopmans ignores orbital relaxation/correlation → **IP/EA are rough**.
* Use consistent basis sets when comparing energies.
* Geometry matters: small bond changes can shift energies noticeably.

## License

MIT. 
