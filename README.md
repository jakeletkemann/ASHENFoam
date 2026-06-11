# ASHENFoam
Ablative Surface High-Enthalpy Non-equilibrium reaction package for OpenFOAM's Hy2Foam. This package allows surface fluxes from ablation models to be included in Hy2Foam.

## What does this package do?

ASHENFoam produces wall-boundary mass fractions from a mass balance on the surface, where source terms from ablation models are included. The mixture-averaged sum of the source terms is used to calculate the flux on the surface, from which the boundary-value mass fractions are transported into the flow. 

The complete version of this package also solves the electron energy equation and includes collisional radiative chemical kinetics modelling for N$_2$ and NO, as well as relevant plasma chemical kinetics modelling. 

To install, Hy2Foam must be first installed into your OpenFOAM installation (see https://github.com/hystrath/hyStrath/tree/master to install and https://hystrath.github.io/solvers/fleming/hy2foam/ for guidance). When this has been done, unzip ASHENFoam.zip into "/hyStrath/applications/solvers/compressible/hy2Foam/" which is located wherever you chose to install Hy2Foam. Run Allwmake to complete the installation.

## How to use:

A new dictionary must be included in /constant in the project directory, called ablationProperties. It allows for the following options:

ablationType takes a string "surfaceReactionPark" or "surfaceReactionACA." If none of these strings are present, no ablation model will be used.

- surfaceReactionPark utilises Park's phenomonological model of surface chemistry. 
- surfaceReactionACA utilises the ACA model of finite-rate surface chemistry. 

For ablation to be implemented properly, mass fraction entries and U in the 0 folder **must** have wall boundary conditions of fixedValue. If another boundary condition is selected, ablation calculations will be overwritten.  


