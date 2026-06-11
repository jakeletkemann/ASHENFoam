// #include <fstream>


if (solveSpeciesEqns)
{
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        speciesDiffusion().correct();
    }
    
    reaction->correct();
    dQ = reaction->dQ();

    const volScalarField& Wmix = composition.Wmix();
    
    forAll(Y, speciei) 
    {
        volScalarField& Yi = Y[speciei];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + fvm::Sp(Sp[speciei], Y[speciei])
          ==
          fvOptions(rho, Yi)
          + Su[speciei]
        );
        
        if (reaction->active())
        {
            YiEqn -= reaction->R(Yi);
            // YiEqn -= ablationSource[speciei];
            // Info << "Species: " << speciei << ", max flux term: " << max(ablationSource[speciei]) << ", min flux term: " << min(ablationSource[speciei]) <<endl;
            // Info << "Species: " << speciei << ", reaction term: " << max(reaction->R(Yi)) << endl;
        }

        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            if (containsElectrons and posElectrons == speciei)
            {
                YiEqn += 
                    fvc::div
                    (
                        speciesDiffusion().J(posElectrons),
                        "div(diffusiveFluxes)"
                    );
            }
            else
            {
                YiEqn += 
                    fvm::laplacian
                    (
                       -speciesDiffusion().rhoD(speciei), // TEST FOR DIFFUSION
                        Yi,
                        "laplacian(rhoD,Yi)"
                    )
                  + fvc::div
                    (
                        speciesDiffusion().JGradp(speciei)
                      + speciesDiffusion().JGradT(speciei),
                        "div(diffusiveFluxes)"
                    );
                    
                if (not speciesDiffusion().useNonCorrected())
                {
                    YiEqn -=
                        fvc::div
                        (
                            speciesDiffusion().sumDiffusiveFluxes()*Yi,
                            "div(diffusiveFluxes)"
                        );
                }      
                    
                if (speciesDiffusion().solvingForX())
                {
                    YiEqn -=
                        fvc::div
                        (
                            speciesDiffusion().rhoD(speciei)*Yi/Wmix
                          * fvc::grad(Wmix),
                            "div(diffusiveFluxes)"
                        );
                }
            }
        }
        
        YiEqn.relax();
        
        fvOptions.constrain(YiEqn);
        
        YiEqn.solve(mesh.solver("Yi"));
        
        fvOptions.correct(Yi);
    }
}

// Info<<U.boundaryFieldRef<<nl;

thermo.correctVolFractions();

if (solveSpeciesEqns)
{
    forAll(Y, speciei) 
    {
        Y[speciei].correctBoundaryConditions();
    }
}

thermo.correctBdrFractions();

if (hyLight)
{
    thermo.correct2T_Light();
}
else
{
    thermo.correct2T();
}

transport().correct();


// std::ofstream outFile("U_boundaryFieldRef.txt", std::ios_base::app);

// // forAll(U.boundaryFieldRef(), patchi)
// //     {
// //         const auto& patchField = U.boundaryFieldRef()[patchi];
// //         outFile << "Patch " << patchi << ":\n";
// //         forAll(patchField, facei)
// //         {
// //             outFile << "  " << patchField[facei] << "\n";
// //         }
// //     }

//     outFile << "test" << nl;

//     // Close the file stream
//     outFile.close();
//     Info << "Boundary field values have been written to U_boundaryFieldRef.txt" << nl;


