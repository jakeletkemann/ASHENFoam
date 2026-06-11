// --- Re-read run-time editable dictionaries
#include "runTimeEditing/hTCPropertiesDictModification.H"
#include "runTimeEditing/transportDictModification.H"
#include "runTimeEditing/twoTemperatureDictModification.H"

// --- Upwind interpolation of primitive fields on faces
#include "numerics/upwindInterpolation.H"
// --- Compute numerical fluxes
#include "numerics/fluxesCalculation.H"

// --- Time control and local time stepping
#include "numerics/compressibleCourantNo.H"
#include "readTimeControls.H"
#include "LTS/readLTSTimeControls.H"
#include "setDeltaT.H"

// label patchID =mesh.boundaryMesh().findPatchID("wedge");
const fvPatch& patch = mesh.boundary()[patchID];

runTime++;
Info<< "Time = " << runTime.timeName() << nl << endl;

volScalarField muEff(turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

/*if (hyDev)
{
    Zrot = ZrotInf/(1.0 + sqr(pi)/2.0*sqrt(TZ/thermo.T())
        + (sqr(pi)/4.0 + pi)*TZ/thermo.T());
    Z = 3.0*pi/(4.0*(thermo.zetar()+3.0))*Zrot;

    mub = mu*(2.0*thermo.zetar()*Z)/(3.0*(thermo.zetar()+3.0));
    bulkToShearRatio = mub/mu;

    tauMC = tauMC + mub*tr(fvc::grad(U))*I;
}*/

// --- Local time stepping (LTS)
if (activateLTS)
{
    #include "LTS/setrDeltaT.H"
}

// --- Re-set the switch/counter that serves as a warning if the
//     temperature goes unbounded
thermo.resetTemperatureBoundingInfo();

if (ablationType=="surfaceReactionACA")
{
//X_Os.correctBoundaryConditions();
//X_Ns.correctBoundaryConditions();
//X_Os_star.correctBoundaryConditions();
// Solve ablation stuff
#include "ablation/surfaceReactionsACA_boundaryCondition_NewtonFirst.H"
}
else if (ablationType=="surfaceReactionPark") {
#include "ablation/surfaceReactionsPark_boundaryCondition_NewtonFirst.H"
}

//Info << "Check2" << endl;

// --- Solve continuity equation
#include "eqns/rhoEqn.H"


// --- Solve species transport and reaction equations
#include "eqns/YEqn.H"

//#include "eqns/YEqn_old.H"

// --- Solve momentum equations
#include HY2FOAM_EXTERNAL_FILE_MOMENTUM_EQNS

//#include "eqns/UEqn_old.H"



//// Set thermal conductivity and enthalpy flux to zero on the surface for ablation
//if (ablationType == "constantFlux" || ablationType == "surfaceReactionPark" || ablationType == "surfaceReactionZA" || ablationType == "surfaceReactionACA")
//{
//tmp<volVectorField> tmpHeatSource = speciesDiffusion().multiSpeciesHeatSource();
//tmp<volVectorField> tmpVEHeatSource = speciesDiffusion().multiSpeciesVEHeatSource();
//
//volVectorField& heatSource = const_cast<volVectorField&>(tmpHeatSource.ref());
//volVectorField& veHeatSource = const_cast<volVectorField&>(tmpVEHeatSource.ref());
//
//heatSource.boundaryFieldRef()[patchID] = vector::zero;
//veHeatSource.boundaryFieldRef()[patchID] = vector::zero;
//
//}


if (downgradeSingleT)
{
    // --- Solve the total energy equation
    //- inviscid
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else if (downgradeSingleTv)
{
    // --- Solve the vibrational energy equation
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    
    //#include "eqns/evEqnInviscid_SingleTv_old.H"
    //#include "eqns/eEqnInviscid_old.H"
    
    //- viscous
    #include "eqns/evEqnViscous_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
    
    //#include "eqns/evEqnViscous_SingleTv_old.H"
    //#include "eqns/eEqnViscous_old.H"
}
else if (downgradeSingleVibMode)
{
    // --- Solve the vibrational energy equations
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include "eqns/evEqnViscous.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else
{
    // ABORTIVE WORK
    // --- Solve the vibrational energy equations (one per vib. mode)
    // --- Solve the total energy equation
    //- inviscid
    //#include "eqns/evEqn_MultiModes.H"
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}


// --- Electro-magnetic fields calculation
#include "eqns/BEqn.H"

// --- Pressure field calculation
#include "eqns/pEqn.H"

// --- Enable input from / coupling with another application
#ifdef HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#include HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#endif

turbulence->correct();

mhd->update(U);

rarefactionParameters().correct(U);

scalar ablationTime = runTime.elapsedCpuTime();

scalar postAblationTime = runTime.elapsedCpuTime();

//forAll(mesh.C(), celli)
//{
//    gamma[celli] = thermo.gamma()()[celli]; // unwrap const tmp on the fly
//}



//tmp<volScalarField> gammaTmp = thermo.gamma();
//gamma = gammaTmp(); // move the data into gamma

//gamma.boundaryFieldRef() = gamma(thermo.p(), thermo.T(), thermo.Tv(), mesh.boundaryMesh().findPatchID("wedge"))


// --- Print a report in the log file if temperature had to be bounded
thermo.temperatureBoundingReport();

if (runTime.outputTime())
{
    runTime.write();
    #include "write/write.H"

    #ifdef HY2FOAM_EXTERNAL_FILE_OUTPUT
    #include HY2FOAM_EXTERNAL_FILE_OUTPUT
    #endif
}

previousIterationTime =
    max(runTime.elapsedCpuTime()-currentIterationTime, 1e-3);
     
//Pout << endl << "Percentage of ablation module of total runtime: " << (postAblationTime - ablationTime) / previousIterationTime << endl << endl;
//Pout << "Ablation time: " << (postAblationTime - ablationTime) << endl;
//Info << "Pre ablation time: " << ablationTime << endl;
//Info << "Post ablation time: " << postAblationTime << endl;

Info<< "Phase no " << noRestart << "." << noSubRestart
    << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << "  Iteration no " << noIteration<<" (" << previousIterationTime
    << " s)" << nl << endl;

currentIterationTime = runTime.elapsedCpuTime();
noIteration += 1;
