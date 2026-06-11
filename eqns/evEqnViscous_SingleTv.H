


// --- Solve vibrational energy equation, VISCOUS
if (!inviscid)
{
    fvScalarMatrix evEqnViscous
    (
        fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
      - fvm::laplacian(transport().kappave()/thermo.CvvelMix(), ev)
      ==
        fvOptions(rho, ev)
    );
    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        evEqnViscous +=
            fvc::div
            (
                speciesDiffusion().multiSpeciesVEHeatSource(), //Electron energy must go in correct equation
                //speciesDiffusion().multiSpeciesVEHeatSource() - thermo_.composition().hevel(posElectrons)*Jcorrected(posElectrons), //Electron energy must go in correct equation
                "div(diffusiveFluxes)"
            );
    }

    evEqnViscous.solve();

    ev.correctBoundaryConditions();
    rhoEv = rho*ev;
}

