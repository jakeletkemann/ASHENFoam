// ------------------------------------------------------------
// Ablation patch
// ------------------------------------------------------------
const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

List<label> wallPatchIDs;

forAll(bMesh, patchI)
{
    const polyPatch& pp = bMesh[patchI];

    if (isA<wallPolyPatch>(pp))
    {
        wallPatchIDs.append(patchI);

        Info << "Found wall patch for ablation: "
             << pp.name()
             << " (ID = " << patchI << ")"
             << nl;
    }
}

if (wallPatchIDs.empty())
{
    FatalErrorInFunction
        << "No wall patches found. Ablation must occur on wall patch." << exit(FatalError);
}

label patchID = wallPatchIDs[0]; // Must eventually update so this works on all wall patches in sequence, currently defaults to first one

IOdictionary ablationProperties
(
    IOobject
    (
        "ablationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Access parameters from customDict
string ablationType = ablationProperties.lookup("ablationType");

// ------------------------------------------------------------
// Constants (Park sublimation, hard-coded as requested)
// ------------------------------------------------------------
scalar Av = 6.02e23;
scalar B  = 3.8240917e19 / Av;   // [kg/m2/s] Park prefactor
scalar Es = 8.71e5;              // [J/mol] Park sublimation energy
scalar Rg = 8.314462618;         // [J/mol/K]

// ------------------------------------------------------------
// Diagnostics fields
// ------------------------------------------------------------

label nSpecies = composition.Y().size();

autoPtr<volScalarField> ablationFluxTotalPtr;
autoPtr<volScalarField> ablationSourcePtr;

if (ablationType == "surfaceReactionACA" || ablationType == "surfaceReactionPark")
{
    Info << "Initializing surface flux fields" << nl;

    ablationFluxTotalPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "BlowingRate",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/(dimArea*dimTime), 0.0)
        )
    );

    // ablationSourcePtr.reset
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "ablationSource",
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         mesh,
    //         dimensionedScalar("zero", dimMass/(dimVolume*dimTime), 0.0)
    //     )
    // );
}

// volScalarField mass_flux_ablation
// (
//     IOobject
//     (
//         "mass_flux_ablation",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0)
// );

volScalarField D_Output
(
    IOobject
    (
        "D",                      // field name
        runTime.timeName(),       // current timestep folder
        mesh,
        IOobject::NO_READ,        // don't read from 0/
        IOobject::AUTO_WRITE      // ensure it's written each timestep
    ),
    mesh,
    dimensionedScalar("D", dimensionSet(0,2,-1,0,0,0,0), 0.0) // init
);

// ------------------------------------------------------------
// Species-resolved surface fluxes
// ------------------------------------------------------------

autoPtr<PtrList<volScalarField>> ablationFluxPtr;


if (ablationType == "surfaceReactionACA" || ablationType == "surfaceReactionPark")
{
    Info << "Initializing species ablation flux fields" << nl;

    ablationFluxPtr.reset(new PtrList<volScalarField>(nSpecies));

    PtrList<volScalarField>& ablationFlux = ablationFluxPtr();

    for (label Yi = 0; Yi < nSpecies; Yi++)
    {
        ablationFlux.set
        (
            Yi,
            new volScalarField
            (
                IOobject
                (
                    "ablationFlux_" + composition.species()[Yi],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0)
            )
        );
    }
}

// label nSpecies = composition.Y().size();

// Foam::PtrList<Foam::volScalarField> ablationFlux(nSpecies);

// for (label Yi = 0; Yi < nSpecies; Yi++)
// {
//     ablationFlux.set
//     (
//         Yi,
//         new Foam::volScalarField
//         (
//             IOobject
//             (
//                 "ablationFlux_" + composition.species()[Yi],
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0)
//         )
//     );
// }

// ------------------------------------------------------------
// Volumetric source terms (this feeds YEqn.H)
// ------------------------------------------------------------
Foam::PtrList<Foam::volScalarField> ablationSource(nSpecies);

for (label Yi = 0; Yi < nSpecies; Yi++)
{
    ablationSource.set
    (
        Yi,
        new Foam::volScalarField
        (
            IOobject
            (
                "ablationSource_" + composition.species()[Yi],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}

volScalarField X_Os
    (
        IOobject
        (
            "Os",               // Field name
            runTime.timeName(),   // Time folder (e.g., "0")
            mesh,
            IOobject::MUST_READ,  // Read field from file
            IOobject::AUTO_WRITE  // Automatically write to new time folders
        ),
        mesh
    );

volScalarField X_Ns
    (
        IOobject
        (
            "Ns",               // Field name
            runTime.timeName(),   // Time folder (e.g., "0")
            mesh,
            IOobject::MUST_READ,  // Read field from file
            IOobject::AUTO_WRITE  // Automatically write to new time folders
        ),
        mesh
    );
    
volScalarField X_Os_star
    (
        IOobject
        (
            "Os_star",               // Field name
            runTime.timeName(),   // Time folder (e.g., "0")
            mesh,
            IOobject::MUST_READ,  // Read field from file
            IOobject::AUTO_WRITE  // Automatically write to new time folders
        ),
        mesh
    );
    
volScalarField X_s
    (
        IOobject
        (
            "s",               // Field name
            runTime.timeName(),   // Time folder (e.g., "0")
            mesh,
            IOobject::MUST_READ,  // Read field from file
            IOobject::AUTO_WRITE  // Automatically write to new time folders
        ),
        mesh
    );

scalarField X_O(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_N(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_O2(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_N2(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_C3(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_CO(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_CN(mesh.C().boundaryField()[patchID].size(), 0.0);
scalarField X_CO2(mesh.C().boundaryField()[patchID].size(), 0.0);

scalarField blowing_rate(mesh.C().boundaryField()[patchID].size(), 0.0);

const scalar pi = 3.14159265;
const scalar kb = 1.38e-23;
const scalar h = 6.62607e-34;
const scalar R = 8.314;


Foam::PtrList<Foam::volScalarField> Sp(nSpecies);

for (label Yi = 0; Yi < nSpecies; Yi++)
{
    Sp.set
    (
        Yi,
        new Foam::volScalarField
        (
            IOobject
            (
                "Sp_" + composition.species()[Yi],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}

Foam::PtrList<Foam::volScalarField> Su(nSpecies);

for (label Yi = 0; Yi < nSpecies; Yi++)
{
    Su.set
    (
        Yi,
        new Foam::volScalarField
        (
            IOobject
            (
                "Su_" + composition.species()[Yi],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );
}












IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
//bool includeTeEqn = thermophysicalProperties.lookup("includeTeEqn");
bool sourceActive = true;