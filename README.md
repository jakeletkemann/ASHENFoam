# ASHENFOAM
Ablative Surface High-Enthalpy Non-equilibrium reaction package for OpenFOAM's Hy2Foam. This package allows surface flux to be included in Hy2Foam either through constant flux distribution or surface reactions.


To install, unzip ASHENFOAM.zip into /hyStrath-OF-v1612-/applications/solvers/compressible/hy2Foam/.

## How to use:

A new dictionary must be included in /constant in the project directory, called ablationProperties. It allows for the following options:

ablationType takes a string "constantFlux" or "surfaceReactions"

{
constantFlux detects the wall boundary and creates an exponentially decaying distribution of ablation flux with a maximum ablationMagnitude at the leftmost x-value. Currently, the flux is only distributed over x.

surfaceReactions is unimplemented.
}

ablationComposition is unimplemented but will take a subdictionary of species and their mass (mole?) fractions. This information will be written to the wall boundary composition.


ablationMagnitude takes a value in kg/m2s


distributionEquation is currently unimplemented but will take an expression for the spatial decay of ablation, for example, "exp(-93.3*x)"


