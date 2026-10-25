# Wave_induced_breakup_of_sea_ice

# Modified DualSPHysics Code and MATLAB Scripts

## Overview

This repository contains modified DualSPHysics source code incorporating a discrete element method (DEM)-based bonding model, together with MATLAB scripts for analyzing wave-sea ice simulation results. 

## Repository Contents

### Modified DualSPHysics Code

The original DualSPHysics source code has been modified to incorporate a DEM-based bonding model that represents connections between adjacent discrete elements (DEs). The model enables bond failure and the subsequent formation of disconnected floes during the simulation.

### MATLAB Scripts

- **`elevation.m`** — Processes SPH simulation data to calculate surface elevation.
- **`stress.m`** — Calculates local stress using surface curvature derived from the elevation data.
- **`reference.m`** — Records the initial configuration of intact bonds at $t=0$ to establish a reference for tracking bond breakage.
- **`analysis.m`** — Compares the final bond configuration with the initial reference to identify broken bonds and determine the resulting floes.
- **`L_eff.m`** — Performs floe-size analysis and determines a characteristic floe size using the effective length scale, $L_{\mathrm{eff}}$.
- **`Figure4.m`** — Generates Figure 4 presented in the associated manuscript.

## License

Refer to the original DualSPHysics licensing terms for conditions governing the use and redistribution of the modified code.
