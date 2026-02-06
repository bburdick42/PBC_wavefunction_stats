# PBC_stats.py

Calculates the center of mass, and normalized spread of a wavefunction under periodic boundary conditions.

## Overview
To quantify localization, wavefunction spread needs to be computed. However, periodic boundary conditions can cause the naive calculation of wavefunction center and spread to be incorrect. This package implements a set of tools used to calculate the center and spread of wavefunctions even when the wavefunction intersects with the periodic boundary.

## Usage
Extract your wavefunction as a numpy array, and specify the physical dimensions (lengths) of the array. Then, call the relevant functions in "PBC_stats.py"

```python
import numpy as np
from PBC_stats_v3 import get_periodic_com, get_wf_spread

com = get_periodic_com(my_wavefunction, my_wavefunction_lengths)
spread = get_wf_spread(my_wavefunction, my_wavefunction_lengths)
```
## Requirements
- Python >= 3.12.4
- numpy >= 2.4.2
- scipy >= 1.17.0

