
//==============================================================
// ACA with sublimation surface finite-rate chemistry
// Implemented as conservative volumetric source terms
//==============================================================

PtrList<volScalarField>& ablationFlux = ablationFluxPtr();
volScalarField& mass_flux_ablation = ablationFluxTotalPtr();

// -----------------------------
// Species indices
// -----------------------------
label idO   = composition.species()["O"];
label idO2  = composition.species()["O2"];
label idN   = composition.species()["N"];
label idN2  = composition.species()["N2"];
label idC  = composition.species()["C"];
label idCO  = composition.species()["CO"];
label idCO2 = composition.species()["CO2"];
label idCN  = composition.species()["CN"];
label idC3  = composition.species()["C3"];
label idNO  = composition.species()["NO"];
label idNCO  = composition.species()["NCO"];
label idC2  = composition.species()["C2"];

// -----------------------------
// Zero fields (EVERY timestep)
// -----------------------------
for (label Yi = 0; Yi < nSpecies; Yi++)
{
    ablationFlux[Yi] =
        dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0);

    ablationSource[Yi] =
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0);

}

mass_flux_ablation =
    dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0);

scalar w1 = 0; scalar w2 = 0; scalar w3 = 0; scalar w4 = 0; scalar w5 = 0; scalar w6 = 0; scalar w7 = 0; scalar w8 = 0; scalar w9 = 0; scalar w10 = 0; 
scalar w11 = 0; scalar w12 = 0; scalar w13 = 0; scalar w14 = 0; scalar w15 = 0; scalar w16 = 0; scalar w17 = 0; scalar w18 = 0; scalar w19 = 0; scalar w20 = 0; 
scalar w21 = 0; scalar w22 = 0;

scalar sublimationEfficiency = 1.0;
scalar Le = 1.0;

scalar localMin = 0.0;
scalar maxIter = 0;
scalar iterPrev = 0;
speciesDiffusion().correct();

// -----------------------------
// Patch access
// -----------------------------
//const fvPatch& patch = mesh.boundary()[patchID];
const scalarField& magSf = patch.magSf();
const labelUList& faceCells = patch.faceCells();
const scalarField& Tw = thermo.T().boundaryField()[patchID];
scalar dtSurf = runTime.deltaTValue();
scalar epsilon_timestep = 1e-8;
scalar num_divisions = 0;
//const volScalarField& D = speciesDiffusion().Di(0);
const scalarField& Dfield =
    speciesDiffusion().rhoD(0);
scalar Dprev = 0.00013;

label indepI = 0;
List<label> indepSpecies;
List<label> mapSpeciesID(nSpecies);

for (label Yi = 0; Yi < nSpecies; Yi++)
{
  if (composition.particleCharge(Yi) == 0 && Yi != idN2)
  {
    indepSpecies.append(Yi);
    mapSpeciesID[Yi] = indepI;
    indepI++;
  }
}


forAll(Tw, facei)
{

    scalar alpha = 1.0; // Damping
    label cellI = faceCells[facei];
    scalar Vc = mesh.V()[cellI];
    scalar AfVc = magSf[facei]/Vc;
    scalar wallTemp = Tw[facei];
    const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
    
    // Initial guess: gas-phase mole fractions at the wall
    scalarField Yold(indepSpecies.size());
    
    forAll(indepSpecies, Yi)
    {
        Yold[Yi] = composition.Y(indepSpecies[Yi]).boundaryFieldRef()[patchID][facei];
    }
    scalar Yindep = composition.Y(idN2).boundaryField()[patchID][facei];
    scalar Yindep_perturb = Yindep;

    // Newton iteration
    label iter = 0;
    scalar residual = 1.0;
    const scalar tol = 1e-6;


    scalar v_O    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idO)/1000.0)/Av))/4.0; // Thermal velocity
    scalar v_O2   = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idO2)/1000.0)/Av))/4.0;
    scalar v_N    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idN)/1000.0)/Av))/4.0;
    scalar v_C    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idC)/1000.0)/Av))/4.0;
    scalar v_C3   = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idC3)/1000.0)/Av))/4.0;
    scalar v_O_2D = Foam::sqrt(pi*kb*wallTemp / 2.0 / ((composition.W(idO)/1000.0)/Av));
    scalar v_N_2D = Foam::sqrt(pi*kb*wallTemp / 2.0 / ((composition.W(idN)/1000.0)/Av));
    

    // ----------------------------------------------------------
    // RATE CONSTANTS (ACA / Prata + Park)
    // ----------------------------------------------------------
    scalar gamma1  = v_O / B * 0.3;
    scalar k1 = gamma1;

    scalar gamma2 = (2.0*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k2 = gamma2 * Foam::exp(-44277.0/wallTemp);

    scalar gamma3 = (v_O/B)*100.0;
    scalar k3 = gamma3 * Foam::exp(-4000.0/wallTemp);

    scalar gamma4 = (v_O/B);
    scalar k4 = gamma4 * Foam::exp(-500.0/wallTemp);

    scalar gamma5 = (v_O/B)*0.7;
    scalar k5 = gamma5;

    scalar gamma6 = (2.0*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k6 = gamma6 * Foam::exp(-96500.0/wallTemp);

    scalar gamma7 = (v_O/B)*1000.0;
    scalar k7 = gamma7 * Foam::exp(-4000.0/wallTemp);

    scalar gamma8 = Foam::sqrt(Av/B)*v_O_2D*1e-3;
    scalar k8 = gamma8 * Foam::exp(-15000.0/wallTemp);

    scalar gamma9 = Foam::sqrt(Av/B)*v_O_2D*5e-5;
    scalar k9 = gamma9 * Foam::exp(-15000.0/wallTemp);

    scalar gamma10 = (v_N/B);
    scalar k10 = gamma10 * Foam::exp(-2500.0/wallTemp);

    scalar gamma11 = (2.0*pi*((composition.W(idN)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                     /(Av*B*h*h*h);
    scalar k11 = gamma11 * Foam::exp(-73971.0/wallTemp);

    scalar gamma12 = (v_N/B)*1.5;
    scalar k12 = gamma12 * Foam::exp(-7000.0/wallTemp);

    scalar gamma13 = (v_N/B)*0.5;
    scalar k13 = gamma13 * Foam::exp(-2000.0/wallTemp);

    scalar gamma14 = Foam::sqrt(Av/B)*v_N_2D*0.1;
    scalar k14 = gamma14 * Foam::exp(-21000.0/wallTemp);

    scalar gamma15 = 1e8;
    scalar k15 = gamma15 * Foam::exp(-20676.0/wallTemp);

    scalar gamma16 = (v_O2/(B*B));
    scalar k16 = gamma16 * Foam::exp(-8000.0/wallTemp);

    scalar gamma17 = (v_O2/B)*100.0;
    scalar k17 = gamma17 * Foam::exp(-4000.0/wallTemp);

    scalar gamma18 = (v_O2/B);
    scalar k18 = gamma18 * Foam::exp(-500.0/wallTemp);

    scalar gamma19 = (v_O2/(B*B));
    scalar k19 = gamma19 * Foam::exp(-8000.0/wallTemp);

    scalar gamma20 = (v_O2/B)*1000.0;
    scalar k20 = gamma20 * Foam::exp(-4000.0/wallTemp);

    scalar gamma21 = 5.19e14;
    scalar k21 = v_C3/(R*wallTemp) * gamma21 * Foam::exp(-90845.0/wallTemp);

    scalar k22 = v_C3;

    X_O[facei]  = thermo.rho().boundaryField()[patchID][facei] * Yold[mapSpeciesID[idO]]  / (composition.W(idO)/1000.0);
    X_N[facei]  = thermo.rho().boundaryField()[patchID][facei] * Yold[mapSpeciesID[idN]]  / (composition.W(idN)/1000.0);
    X_O2[facei] = thermo.rho().boundaryField()[patchID][facei] * Yold[mapSpeciesID[idO2]]  / (composition.W(idO2)/1000.0);
    X_N2[facei] = thermo.rho().boundaryField()[patchID][facei] * Yindep   / (composition.W(idN2)/1000.0);
    X_C3[facei] = thermo.rho().boundaryField()[patchID][facei] * Yold[mapSpeciesID[idC3]] / (composition.W(idC3)/1000.0);

    updateSurfaceCoverageExplicit(
        patchID,
        facei,
        dtSurf,
        B,
        X_O, X_O2, X_N,
        X_Os, X_Os_star, X_Ns,
        k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,
        k11,k12,k13,k14,k15,k16,k17,k18,k19,k20
    );
    
    
    X_s.boundaryFieldRef()[patchID][facei] = max(0.0, B - X_Os.boundaryFieldRef()[patchID][facei] - X_Os_star.boundaryFieldRef()[patchID][facei] - X_Ns.boundaryFieldRef()[patchID][facei]);
    
    //X_Os.boundaryFieldRef()[patchID].updateCoeffs();
    //X_Os_star.boundaryFieldRef()[patchID].updateCoeffs();
    //X_Ns.boundaryFieldRef()[patchID].updateCoeffs();
    

    //Pout << "Ns: " << X_Ns.boundaryField()[patchID][facei] << ", x-location: " << mesh.C()[cellI] << endl;   
}


// Number of independent species at wall
const label nIndep = indepSpecies.size();

// Per-face unknowns
List<scalarField> Ywall(patch.size(), scalarField(nIndep, 0.0));
scalarField YN2(patch.size(), 0.0);

// Residuals and updates
List<scalarField> Rwall(patch.size(), scalarField(nIndep+3, 0.0));
List<scalarField> dYwall(patch.size(), scalarField(nIndep+3, 0.0));
List<scalarField> J_diff(patch.size(), scalarField(nIndep, 0.0));
List<scalarField> Ygrad(patch.size(), scalarField(nIndep, 0.0));
scalarField YgradN2(patch.size(), 0.0);

// Newton control
const scalar newtonTol = 1e-8;
const label maxNewtonIter = 300;
const scalar relax = 1.0;


forAll(Tw, facei)
{
    forAll(indepSpecies, k)
    {
        Ywall[facei][k] =
            composition.Y(indepSpecies[k])
            .boundaryField()[patchID][facei];
    }

    scalar sumY = 0.0;
    forAll(indepSpecies, k) sumY += Ywall[facei][k];
    YN2[facei] = max(1.0 - sumY, SMALL);
}


scalarField W(nSpecies, 0.0);
    forAll(composition.Y(), k)
    {
        W[k] = composition.W(k);
    }

label newtonIter = 0;
scalar globalRes = GREAT;
scalar trueIterationNumber = 0.0;

while (globalRes > newtonTol && newtonIter < maxNewtonIter)
{

    // --------------------------------------------------
    // Write ALL face values (MPI-safe)
    // --------------------------------------------------
    forAll(patch, facei)
    {
        forAll(indepSpecies, k)
        {
            composition.Y(indepSpecies[k])
                .boundaryFieldRef()[patchID][facei]
                = Ywall[facei][k];
        }

        composition.Y(idN2)
            .boundaryFieldRef()[patchID][facei]
            = YN2[facei];
    }

    
    forAll(indepSpecies, k)
    {
//    composition.Y(indepSpecies[k]).correctBoundaryConditions();
composition.Y(indepSpecies[k]).boundaryFieldRef().updateCoeffs(); // local patch only
    }




    // --------------------------------------------------
    // ONE collective diffusion update
    // --------------------------------------------------
//    speciesDiffusion().correct();
//    
//    forAll(Tw, facei)
//    {
//    
//    forAll(indepSpecies,k)
//        {
//          J_diff[facei][k] = speciesDiffusion().Jcorr(indepSpecies[k])
//              .boundaryField()[patchID][facei]
//            & (-mesh.Sf().boundaryField()[patchID][facei]
//               / mag(mesh.Sf().boundaryField()[patchID][facei]));
//        }
//
//    }

    // --------------------------------------------------
    // ONE collective snGrad term update
    // --------------------------------------------------
    

    forAll(indepSpecies, Yi)
    {
    volVectorField YgradField = fvc::grad(composition.Y(indepSpecies[Yi]));
    //Info << "Ygradboundary: " << YgradField.boundaryField()[patchID] << endl;
      forAll(Tw, facei)
      {
        const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
        Ygrad[facei][Yi] = YgradField.boundaryField()[patchID][facei] & norm;
      }
    }
    volVectorField YgradField = fvc::grad(composition.Y("N2"));

    forAll(Tw, facei)
    {
      const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
      YgradN2[facei] = YgradField.boundaryField()[patchID][facei] & norm;
    }
    
    
    // --------------------------------------------------
    // Residual evaluation
    // --------------------------------------------------
    globalRes = 0.0;

    forAll(Tw, facei)
    {
//    sublimationEfficiency = 1.0;
//    if (Tw[facei] < 3000)
//    {
//      sublimationEfficiency = 0.0;
//    }
    
    label cellI = faceCells[facei];
    scalar Vc = mesh.V()[cellI];
    scalar AfVc = magSf[facei]/Vc;
    scalar wallTemp = Tw[facei];
    const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
    
    scalar v_O    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idO)/1000.0)/Av))/4; // Thermal velocity
    scalar v_O2   = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idO2)/1000.0)/Av))/4;
    scalar v_N    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idN)/1000.0)/Av))/4;
    scalar v_C    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idC)/1000.0)/Av))/4;
    scalar v_C3   = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idC3)/1000.0)/Av))/4;
    scalar v_O_2D = Foam::sqrt(pi*kb*wallTemp / 2 / ((composition.W(idO)/1000.0)/Av));
    scalar v_N_2D = Foam::sqrt(pi*kb*wallTemp / 2 / ((composition.W(idN)/1000.0)/Av));
    

    // ----------------------------------------------------------
    // RATE CONSTANTS (ACA / Prata + Park)
    // ----------------------------------------------------------
    scalar gamma1  = v_O / B * 0.3;
    scalar k1 = gamma1;

    scalar gamma2 = (2*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k2 = gamma2 * Foam::exp(-44277/wallTemp);

    scalar gamma3 = (v_O/B)*100;
    scalar k3 = gamma3 * Foam::exp(-4000/wallTemp);

    scalar gamma4 = (v_O/B);
    scalar k4 = gamma4 * Foam::exp(-500/wallTemp);

    scalar gamma5 = (v_O/B)*0.7;
    scalar k5 = gamma5;

    scalar gamma6 = (2*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k6 = gamma6 * Foam::exp(-96500/wallTemp);

    scalar gamma7 = (v_O/B)*1000;
    scalar k7 = gamma7 * Foam::exp(-4000/wallTemp);

    scalar gamma8 = Foam::sqrt(Av/B)*v_O_2D*1e-3;
    scalar k8 = gamma8 * Foam::exp(-15000/wallTemp);

    scalar gamma9 = Foam::sqrt(Av/B)*v_O_2D*5e-5;
    scalar k9 = gamma9 * Foam::exp(-15000/wallTemp);

    scalar gamma10 = (v_N/B);
    scalar k10 = gamma10 * Foam::exp(-2500/wallTemp);

    scalar gamma11 = (2*pi*((composition.W(idN)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                     /(Av*B*h*h*h);
    scalar k11 = gamma11 * Foam::exp(-73971/wallTemp);

    scalar gamma12 = (v_N/B)*1.5;
    scalar k12 = gamma12 * Foam::exp(-7000/wallTemp);

    scalar gamma13 = (v_N/B)*0.5;
    scalar k13 = gamma13 * Foam::exp(-2000/wallTemp);

    scalar gamma14 = Foam::sqrt(Av/B)*v_N_2D*0.1;
    scalar k14 = gamma14 * Foam::exp(-21000/wallTemp);

    scalar gamma15 = 1e8;
    scalar k15 = gamma15 * Foam::exp(-20676/wallTemp);

    scalar gamma16 = (v_O2/(B*B));
    scalar k16 = gamma16 * Foam::exp(-8000/wallTemp);

    scalar gamma17 = (v_O2/B)*100;
    scalar k17 = gamma17 * Foam::exp(-4000/wallTemp);

    scalar gamma18 = (v_O2/B);
    scalar k18 = gamma18 * Foam::exp(-500/wallTemp);

    scalar gamma19 = (v_O2/(B*B));
    scalar k19 = gamma19 * Foam::exp(-8000/wallTemp);

    scalar gamma20 = (v_O2/B)*1000;
    scalar k20 = gamma20 * Foam::exp(-4000/wallTemp);

    scalar gamma21 = 5.19e14;
    scalar k21 = v_C3/(R*wallTemp) * gamma21 * Foam::exp(-90845.0/wallTemp);

    scalar k22 = v_C3;
    
    scalar PO = Ywall[facei][mapSpeciesID[idO]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO)/1000.0)/Av*kb*wallTemp));
    scalar PO2 = Ywall[facei][mapSpeciesID[idO2]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO2)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO2)/1000.0)/Av*kb*wallTemp));
    scalar PN = Ywall[facei][mapSpeciesID[idN]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idN)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idN)/1000.0)/Av*kb*wallTemp));
    
    
    scalarField Ycell(nSpecies,0.0);
    
    forAll(composition.Y(), Yi)
    {
      Ycell[Yi] = composition.Y(Yi)[cellI];
    }
    
    
//    vector Cf = mesh.Cf().boundaryField()[patchID][facei];  // face center
//    vector Cc = mesh.C()[cellI]; 
//    scalar delta = Foam::mag(Cf-Cc);

      scalar delta = 1.0 / (mesh.deltaCoeffs().boundaryField()[patchID][facei]);
//      scalarField D(nIndep,0.0);
      
//      forAll(indepSpecies, Yi)
//      {
//      D[Yi] = speciesDiffusion().Di(indepSpecies[Yi])[cellI];
//      }
        //D = thermo.alphatr().boundaryField()[patchID][facei] / thermo.rho().boundaryField()[patchID][facei];
//        if (Foam::mag(D-Dprev)>3e-4)
//        {
//        D = Dprev;
//        }
//        scalar Dprev = D;
//const scalar& D = speciesDiffusion().Di(0)[cellI];
scalar D = Dfield[cellI]/thermo.rho()[cellI] / Le;

    
        computeResiduals
        (
        patchID, facei,
        k1,    k2,    k3,    k4,
        k5,    k6,    k7,    k8,
        k9,    k10,   k11,   k12,
        k13,   k14,   k15,   k16,
        k17,   k18,   k19,   k20,
        k21,   k22,   sublimationEfficiency,
        idO,   idO2,   idN,   idN2,
        idCO,   idCN,   idCO2,   idC3,
        PO, PO2, PN, 
        nSpecies, thermo.rho().boundaryField()[patchID][facei], W, B, delta, D,
        X_Os.boundaryField()[patchID][facei], X_Os_star.boundaryField()[patchID][facei], X_Ns.boundaryField()[patchID][facei],
        indepSpecies, Ygrad[facei], Ycell, Ywall[facei], YN2, YgradN2, J_diff, mapSpeciesID, Rwall[facei]
        );
        
        if (max(mag(Rwall[facei])) > 1.0e-9)
        {
        //Info << "|R|: " << max(mag(Rwall[facei])) << endl;
        //Info << "face index: " << facei << endl;        
        //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;
        }
        //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;
        scalar maxRes = 0.0;
        forAll(Rwall[facei], k)
        {
            if (k != 2)  // only consider species with meaningful concentration
            {
                maxRes = max(maxRes, mag(Rwall[facei][k]));
            }
            //Info << "species " << k << ", Y: " << Ywall[facei][k] << ", maxRes: " << maxRes << endl;
        }
        globalRes = max(globalRes, maxRes);
    }
//Info << "Next iteration! " << endl << endl;
//Info << "|R|: " << globalRes << endl;
//Info << "R: " << Rwall[200] << endl;


if (globalRes > newtonTol)
{
trueIterationNumber++;
}

globalRes = Foam::returnReduce(globalRes, maxOp<scalar>());
trueIterationNumber = Foam::returnReduce(trueIterationNumber, maxOp<scalar>());

    // --------------------------------------------------
    // Jacobian + Newton update (face-local)
    // --------------------------------------------------
    forAll(Tw, facei)
    {
//    sublimationEfficiency = 1.0;
//    if (Tw[facei] < 3000)
//    {
//      sublimationEfficiency = 0.0;
//    }
    label cellI = faceCells[facei];
    scalar Vc = mesh.V()[cellI];
    scalar AfVc = magSf[facei]/Vc;
    scalar wallTemp = Tw[facei];
    const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
    
    scalar v_O    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idO)/1000.0)/Av))/4; // Thermal velocity
    scalar v_O2   = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idO2)/1000.0)/Av))/4;
    scalar v_N    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idN)/1000.0)/Av))/4;
    scalar v_C    = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idC)/1000.0)/Av))/4;
    scalar v_C3   = Foam::sqrt(8*kb*wallTemp / pi / ((composition.W(idC3)/1000.0)/Av))/4;
    scalar v_O_2D = Foam::sqrt(pi*kb*wallTemp / 2 / ((composition.W(idO)/1000.0)/Av));
    scalar v_N_2D = Foam::sqrt(pi*kb*wallTemp / 2 / ((composition.W(idN)/1000.0)/Av));
    
    // ----------------------------------------------------------
    // RATE CONSTANTS (ACA / Prata + Park)
    // ----------------------------------------------------------
    scalar gamma1  = v_O / B * 0.3;
    scalar k1 = gamma1;

    scalar gamma2 = (2*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k2 = gamma2 * Foam::exp(-44277/wallTemp);

    scalar gamma3 = (v_O/B)*100;
    scalar k3 = gamma3 * Foam::exp(-4000/wallTemp);

    scalar gamma4 = (v_O/B);
    scalar k4 = gamma4 * Foam::exp(-500/wallTemp);

    scalar gamma5 = (v_O/B)*0.7;
    scalar k5 = gamma5;

    scalar gamma6 = (2*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k6 = gamma6 * Foam::exp(-96500/wallTemp);

    scalar gamma7 = (v_O/B)*1000;
    scalar k7 = gamma7 * Foam::exp(-4000/wallTemp);

    scalar gamma8 = Foam::sqrt(Av/B)*v_O_2D*1e-3;
    scalar k8 = gamma8 * Foam::exp(-15000/wallTemp);

    scalar gamma9 = Foam::sqrt(Av/B)*v_O_2D*5e-5;
    scalar k9 = gamma9 * Foam::exp(-15000/wallTemp);

    scalar gamma10 = (v_N/B);
    scalar k10 = gamma10 * Foam::exp(-2500/wallTemp);

    scalar gamma11 = (2*pi*((composition.W(idN)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                     /(Av*B*h*h*h);
    scalar k11 = gamma11 * Foam::exp(-73971/wallTemp);

    scalar gamma12 = (v_N/B)*1.5;
    scalar k12 = gamma12 * Foam::exp(-7000/wallTemp);

    scalar gamma13 = (v_N/B)*0.5;
    scalar k13 = gamma13 * Foam::exp(-2000/wallTemp);

    scalar gamma14 = Foam::sqrt(Av/B)*v_N_2D*0.1;
    scalar k14 = gamma14 * Foam::exp(-21000/wallTemp);

    scalar gamma15 = 1e8;
    scalar k15 = gamma15 * Foam::exp(-20676/wallTemp);

    scalar gamma16 = (v_O2/(B*B));
    scalar k16 = gamma16 * Foam::exp(-8000/wallTemp);

    scalar gamma17 = (v_O2/B)*100;
    scalar k17 = gamma17 * Foam::exp(-4000/wallTemp);

    scalar gamma18 = (v_O2/B);
    scalar k18 = gamma18 * Foam::exp(-500/wallTemp);

    scalar gamma19 = (v_O2/(B*B));
    scalar k19 = gamma19 * Foam::exp(-8000/wallTemp);

    scalar gamma20 = (v_O2/B)*1000;
    scalar k20 = gamma20 * Foam::exp(-4000/wallTemp);

    scalar gamma21 = 5.19e14;
    scalar k21 = v_C3/(R*wallTemp) * gamma21 * Foam::exp(-90845.0/wallTemp);

    scalar k22 = v_C3;

    
//    vector Cf = mesh.Cf().boundaryField()[patchID][facei];  // face center
//    vector Cc = mesh.C()[cellI]; 
//    scalar delta = Foam::mag(Cf-Cc);

    scalar delta = 1.0 / (mesh.deltaCoeffs().boundaryField()[patchID][facei]);
    
    //Info << "delta: " << delta << endl;
    
    scalarField Ycell(nSpecies,0.0);
    forAll(composition.Y(), Yi)
    {
      Ycell[Yi] = composition.Y(Yi)[cellI];
    }

        List<scalarField> J(nIndep+3, scalarField(nIndep+3, 0.0));
        
//        forAll(indepSpecies, Yi)
//        {
//        D[Yi] = speciesDiffusion().Di(indepSpecies[Yi])[cellI];
//        }

        //D = thermo.alphatr().boundaryField()[patchID][facei] / thermo.rho().boundaryField()[patchID][facei];
//        if (Foam::mag(D-Dprev)>3e-4)
//        {
//        D = Dprev;
//        }
//        scalar Dprev = D;

//const scalar& D = speciesDiffusion().Di(0)[cellI];
scalar D = Dfield[cellI]/thermo.rho()[cellI] / Le;

        forAll(indepSpecies, j)
        {
                
            

            scalarField Ypert = Ywall[facei];
            scalarField Ygradpert = Ygrad[facei];
            scalar eps = 1e-8;
            Ypert[j] += eps;
            Ygradpert[j] -= eps/delta;
            
            scalar PO = Ypert[mapSpeciesID[idO]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO)/1000.0)/Av*kb*wallTemp));
            scalar PO2 = Ypert[mapSpeciesID[idO2]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO2)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO2)/1000.0)/Av*kb*wallTemp));
            scalar PN = Ypert[mapSpeciesID[idN]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idN)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idN)/1000.0)/Av*kb*wallTemp));
    

            scalarField Rpert(nIndep+3, 0.0);
            computeResiduals
            (
            patchID, facei,
            k1,    k2,    k3,    k4,
            k5,    k6,    k7,    k8,
            k9,    k10,   k11,   k12,
            k13,   k14,   k15,   k16,
            k17,   k18,   k19,   k20,
            k21,   k22,   sublimationEfficiency,
            idO,   idO2,   idN,   idN2,
            idCO,   idCN,   idCO2,   idC3,
            PO, PO2, PN, 
            nSpecies, thermo.rho().boundaryField()[patchID][facei], W, B, delta, D,
            X_Os.boundaryField()[patchID][facei], X_Os_star.boundaryField()[patchID][facei], X_Ns.boundaryField()[patchID][facei],
            indepSpecies, Ygradpert, Ycell, Ypert, YN2, YgradN2, J_diff, mapSpeciesID, Rpert
            );

            forAll(indepSpecies, i)
            {
                J[i][j] = (Rpert[i] - Rwall[facei][i]) / eps;
            }
            

            J[nIndep+0][j] = (Rpert[nIndep+0] - Rwall[facei][nIndep+0]) / eps;
            J[nIndep+1][j] = (Rpert[nIndep+1] - Rwall[facei][nIndep+1]) / eps;
            J[nIndep+2][j] = (Rpert[nIndep+2] - Rwall[facei][nIndep+2]) / eps;
            
        }
        
            scalarField Rpert(nIndep+3, 0.0);
            scalar eps = 1e-8;
            
            scalar PO = Ywall[facei][mapSpeciesID[idO]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO)/1000.0)/Av*kb*wallTemp));
            scalar PO2 = Ywall[facei][mapSpeciesID[idO2]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idO2)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idO2)/1000.0)/Av*kb*wallTemp));
            scalar PN = Ywall[facei][mapSpeciesID[idN]] * R * wallTemp * thermo.rho().boundaryField()[patchID][facei] / (composition.W(idN)/1000.0) / (Av * Foam::sqrt(2.0*pi*(composition.W(idN)/1000.0)/Av*kb*wallTemp));
    
            

            computeResiduals
            (
            patchID, facei,
            k1,    k2,    k3,    k4,
            k5,    k6,    k7,    k8,
            k9,    k10,   k11,   k12,
            k13,   k14,   k15,   k16,
            k17,   k18,   k19,   k20,
            k21,   k22,   sublimationEfficiency,
            idO,   idO2,   idN,   idN2,
            idCO,   idCN,   idCO2,   idC3,
            PO, PO2, PN, 
            nSpecies, thermo.rho().boundaryField()[patchID][facei], W, B, delta, D,
            X_Os.boundaryField()[patchID][facei] + eps, X_Os_star.boundaryField()[patchID][facei], X_Ns.boundaryField()[patchID][facei],
            indepSpecies, Ygrad[facei], Ycell, Ywall[facei], YN2, YgradN2, J_diff, mapSpeciesID, Rpert
            );
        
            forAll(indepSpecies, i)
            {
                J[i][nIndep+0] = (Rpert[i] - Rwall[facei][i]) / eps;
            }
            J[nIndep+0][nIndep+0] = (Rpert[nIndep+0] - Rwall[facei][nIndep+0]) / eps;
            J[nIndep+1][nIndep+0] = (Rpert[nIndep+1] - Rwall[facei][nIndep+1]) / eps;
            J[nIndep+2][nIndep+0] = (Rpert[nIndep+2] - Rwall[facei][nIndep+2]) / eps;

            Rpert = Rpert * 0;

            computeResiduals
            (
            patchID, facei,
            k1,    k2,    k3,    k4,
            k5,    k6,    k7,    k8,
            k9,    k10,   k11,   k12,
            k13,   k14,   k15,   k16,
            k17,   k18,   k19,   k20,
            k21,   k22,   sublimationEfficiency,
            idO,   idO2,   idN,   idN2,
            idCO,   idCN,   idCO2,   idC3,
            PO, PO2, PN, 
            nSpecies, thermo.rho().boundaryField()[patchID][facei], W, B, delta, D,
            X_Os.boundaryField()[patchID][facei], X_Os_star.boundaryField()[patchID][facei] + eps, X_Ns.boundaryField()[patchID][facei],
            indepSpecies, Ygrad[facei], Ycell, Ywall[facei], YN2, YgradN2, J_diff, mapSpeciesID, Rpert
            );

            forAll(indepSpecies, i)
            {
                J[i][nIndep+1] = (Rpert[i] - Rwall[facei][i]) / eps;
            }
            J[nIndep+0][nIndep+1] = (Rpert[nIndep+0] - Rwall[facei][nIndep+0]) / eps;
            J[nIndep+1][nIndep+1] = (Rpert[nIndep+1] - Rwall[facei][nIndep+1]) / eps;
            J[nIndep+2][nIndep+1] = (Rpert[nIndep+2] - Rwall[facei][nIndep+2]) / eps;

            Rpert = Rpert * 0;
            
            computeResiduals
            (
            patchID, facei,
            k1,    k2,    k3,    k4,
            k5,    k6,    k7,    k8,
            k9,    k10,   k11,   k12,
            k13,   k14,   k15,   k16,
            k17,   k18,   k19,   k20,
            k21,   k22,   sublimationEfficiency,
            idO,   idO2,   idN,   idN2,
            idCO,   idCN,   idCO2,   idC3,
            PO, PO2, PN, 
            nSpecies, thermo.rho().boundaryField()[patchID][facei], W, B, delta, D,
            X_Os.boundaryField()[patchID][facei], X_Os_star.boundaryField()[patchID][facei], X_Ns.boundaryField()[patchID][facei] + eps,
            indepSpecies, Ygrad[facei], Ycell, Ywall[facei], YN2, YgradN2, J_diff, mapSpeciesID, Rpert
            );

            forAll(indepSpecies, i)
            {
                J[i][nIndep+2] = (Rpert[i] - Rwall[facei][i]) / eps;
            }
            J[nIndep+0][nIndep+2] = (Rpert[nIndep+0] - Rwall[facei][nIndep+0]) / eps;
            J[nIndep+1][nIndep+2] = (Rpert[nIndep+1] - Rwall[facei][nIndep+1]) / eps;
            J[nIndep+2][nIndep+2] = (Rpert[nIndep+2] - Rwall[facei][nIndep+2]) / eps;

        //Info << J << endl;
        solveLinearSystem(J, Rwall[facei], dYwall[facei]);
        //Info << "dY: " << dYwall[facei] << endl;
        
    }

    // --------------------------------------------------
    // Update unknowns
    // --------------------------------------------------
    forAll(Tw, facei)
    {
        forAll(indepSpecies, k)
        {
               
            Ywall[facei][k] -= relax * dYwall[facei][k];
            //Ywall[facei][k] = max(Ywall[facei][k], SMALL);
        }

        YN2[facei] = max(1.0 - sum(Ywall[facei]), SMALL);
        
    X_Os.boundaryFieldRef()[patchID][facei] -= relax * dYwall[facei][nIndep+0];
    X_Os_star.boundaryFieldRef()[patchID][facei] -= relax * dYwall[facei][nIndep+1];
    X_Ns.boundaryFieldRef()[patchID][facei] -= relax * dYwall[facei][nIndep+2];
    X_s.boundaryFieldRef()[patchID][facei] = B - X_Os.boundaryFieldRef()[patchID][facei] - X_Os_star.boundaryFieldRef()[patchID][facei] - X_Ns.boundaryFieldRef()[patchID][facei];

        
    //Info << "x-location: " << mesh.Cf().boundaryField()[patchID][facei] << endl;
    }
    

    newtonIter++;
}

    forAll(Tw, facei)
    {
    //Info << "Y: " << Ywall[facei] << endl;
    }
    
    
//Info << "Y: " << Ywall[0] << endl;
//Info << "Coverages: " << X_Os.boundaryFieldRef()[patchID][0] << ", " << X_Os_star.boundaryFieldRef()[patchID][0] << ", " << X_Ns.boundaryFieldRef()[patchID][0] << endl;
    
Info << "Iterations: " << trueIterationNumber << endl;
//Pout << "N mass fraction: " << max(Ywall[:][8] << endl;
//Pout << "CN mass fraction: " << max(Ywall[:][3] << endl;



scalar omega = 0.5;
label nIter = 5;

//smoothPatchY(patch, Ywall, omega, nIter);
forAll(Tw, facei)
{
  scalar totalComposition = 0.0;
  forAll(Ywall[facei],k)
  {
    totalComposition += Ywall[facei][k];
  }
  totalComposition += YN2[facei];
  forAll(Ywall[facei],Yi)
  {
    composition.Y(indepSpecies[Yi]).boundaryFieldRef()[patchID][facei] = Ywall[facei][Yi]/totalComposition;
  }
  composition.Y("N2").boundaryFieldRef()[patchID][facei] = YN2[facei] / totalComposition;
}




forAll(Tw, facei)
{

    const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];

    label cellI = faceCells[facei];
    scalar wallTemp = Tw[facei];
    
    // ----------------------------------------------------------
    // Thermal velocities (m/s)
    // ----------------------------------------------------------

    scalar v_O    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idO)/1000.0)/Av))/4.0; // Thermal velocity
    scalar v_O2   = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idO2)/1000.0)/Av))/4.0;
    scalar v_N    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idN)/1000.0)/Av))/4.0;
    scalar v_C    = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idC)/1000.0)/Av))/4.0;
    scalar v_C3   = Foam::sqrt(8.0*kb*wallTemp / pi / ((composition.W(idC3)/1000.0)/Av))/4.0;
    scalar v_O_2D = Foam::sqrt(pi*kb*wallTemp / 2.0 / ((composition.W(idO)/1000.0)/Av));
    scalar v_N_2D = Foam::sqrt(pi*kb*wallTemp / 2.0 / ((composition.W(idN)/1000.0)/Av));
    
    // ----------------------------------------------------------
    // RATE CONSTANTS (ACA / Prata + Park)
    // ----------------------------------------------------------
    scalar gamma1  = v_O / B * 0.3;
    scalar k1 = gamma1;

    scalar gamma2 = (2.0*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k2 = gamma2 * Foam::exp(-44277.0/wallTemp);

    scalar gamma3 = (v_O/B)*100.0;
    scalar k3 = gamma3 * Foam::exp(-4000.0/wallTemp);

    scalar gamma4 = (v_O/B);
    scalar k4 = gamma4 * Foam::exp(-500.0/wallTemp);

    scalar gamma5 = (v_O/B)*0.7;
    scalar k5 = gamma5;

    scalar gamma6 = (2.0*pi*((composition.W(idO)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                    /(Av*B*h*h*h);
    scalar k6 = gamma6 * Foam::exp(-96500.0/wallTemp);

    scalar gamma7 = (v_O/B)*1000.0;
    scalar k7 = gamma7 * Foam::exp(-4000.0/wallTemp);

    scalar gamma8 = Foam::sqrt(Av/B)*v_O_2D*1e-3;
    scalar k8 = gamma8 * Foam::exp(-15000.0/wallTemp);

    scalar gamma9 = Foam::sqrt(Av/B)*v_O_2D*5e-5;
    scalar k9 = gamma9 * Foam::exp(-15000.0/wallTemp);

    scalar gamma10 = (v_N/B);
    scalar k10 = gamma10 * Foam::exp(-2500.0/wallTemp);

    scalar gamma11 = (2.0*pi*((composition.W(idN)/1000.0)/Av)*kb*kb*wallTemp*wallTemp)
                     /(Av*B*h*h*h);
    scalar k11 = gamma11 * Foam::exp(-73971.0/wallTemp);

    scalar gamma12 = (v_N/B)*1.5;
    scalar k12 = gamma12 * Foam::exp(-7000.0/wallTemp);

    scalar gamma13 = (v_N/B)*0.5;
    scalar k13 = gamma13 * Foam::exp(-2000.0/wallTemp);

    scalar gamma14 = Foam::sqrt(Av/B)*v_N_2D*0.1;
    scalar k14 = gamma14 * Foam::exp(-21000.0/wallTemp);

    scalar gamma15 = 1e8;
    scalar k15 = gamma15 * Foam::exp(-20676.0/wallTemp);

    scalar gamma16 = (v_O2/(B*B));
    scalar k16 = gamma16 * Foam::exp(-8000.0/wallTemp);

    scalar gamma17 = (v_O2/B)*100.0;
    scalar k17 = gamma17 * Foam::exp(-4000.0/wallTemp);

    scalar gamma18 = (v_O2/B);
    scalar k18 = gamma18 * Foam::exp(-500.0/wallTemp);

    scalar gamma19 = (v_O2/(B*B));
    scalar k19 = gamma19 * Foam::exp(-8000.0/wallTemp);

    scalar gamma20 = (v_O2/B)*1000.0;
    scalar k20 = gamma20 * Foam::exp(-4000.0/wallTemp);

    scalar gamma21 = 5.19e14;
    scalar k21 = v_C3/(R*wallTemp) * gamma21 * Foam::exp(-90845.0/wallTemp);

    scalar k22 = v_C3;


    X_O[facei]  = thermo.rho().boundaryField()[patchID][facei] * composition.Y(idO).boundaryField()[patchID][facei]  / (composition.W(idO)/1000.0);
    X_N[facei]  = thermo.rho().boundaryField()[patchID][facei] * composition.Y(idN).boundaryField()[patchID][facei]  / (composition.W(idN)/1000.0);
    X_O2[facei] = thermo.rho().boundaryField()[patchID][facei] * composition.Y(idO2).boundaryField()[patchID][facei] / (composition.W(idO2)/1000.0);
    X_N2[facei] = thermo.rho().boundaryField()[patchID][facei] * composition.Y(idN2).boundaryField()[patchID][facei] / (composition.W(idN2)/1000.0);
    X_C3[facei] = thermo.rho().boundaryField()[patchID][facei] * composition.Y(idC3).boundaryField()[patchID][facei] / (composition.W(idC3)/1000.0);

    w1  = k1  * X_O[facei]  * X_s.boundaryField()[patchID][facei];
    w2  = k2  * X_Os.boundaryField()[patchID][facei];
    w3  = k3  * X_O[facei]  * X_Os.boundaryField()[patchID][facei];
    w4  = k4  * X_O[facei]  * X_Os.boundaryField()[patchID][facei];
    w5  = k5  * X_O[facei]  * X_s.boundaryField()[patchID][facei];
    w6  = k6  * X_Os_star.boundaryField()[patchID][facei];
    w7  = k7  * X_O[facei]  * X_Os_star.boundaryField()[patchID][facei];
    w8  = k8  * sqr(X_Os_star.boundaryField()[patchID][facei]);
    w9  = k9  * sqr(X_Os.boundaryField()[patchID][facei]);
    w10 = k10 * X_N[facei]  * X_s.boundaryField()[patchID][facei];
    w11 = k11 * X_Ns.boundaryField()[patchID][facei];
    w12 = k12 * X_N[facei]  * X_Ns.boundaryField()[patchID][facei];
    w13 = k13 * X_N[facei]  * X_Ns.boundaryField()[patchID][facei];
    w14 = k14 * sqr(X_Ns.boundaryField()[patchID][facei]);
    w15 = k15 * X_Ns.boundaryField()[patchID][facei];
    w16 = k16 * X_O2[facei] * sqr(X_s.boundaryField()[patchID][facei]);
    w17 = k17 * X_O2[facei] * X_Os.boundaryField()[patchID][facei];
    w18 = k18 * X_O2[facei] * X_Os.boundaryField()[patchID][facei];
    w19 = k19 * X_O2[facei] * sqr(X_s.boundaryField()[patchID][facei]);
    w20 = k20 * X_O2[facei] * X_Os_star.boundaryField()[patchID][facei];
    w21 = k21 * sublimationEfficiency;
    w22 = k22 * X_C3[facei] * sublimationEfficiency;

    // ----------------------------------------------------------
    // Gas-phase production (kg/m2/s)
    // ----------------------------------------------------------
    scalar J_CO  = (w3 + w7 + w17 + w20) * composition.W(idC)/1000.0;
    scalar J_CO2 = (w4 + w18) * composition.W(idC)/1000.0;
    scalar J_CN  = (w12 + w15) * composition.W(idC)/1000.0;
    scalar J_C3  = (w21 - w22) * composition.W(idC3)/1000.0;
    scalar J_O   = (-w1 + w2 - w4 - w5 + w6 + w18) * composition.W(idO)/1000.0;
    scalar J_O2  = (w8 + w9 - w16 - w18 - w19) * composition.W(idO2)/1000.0;
    scalar J_N   = (-w10 + w11 - w13) * composition.W(idN)/1000.0;
    scalar J_N2  = (w13 + w14) * composition.W(idN2)/1000.0;

    // Info << "J_CO: " << J_CO << endl;
    // Info << "J_C3: " << J_C3 << endl;

    // Surface diagnostics
    ablationFlux[idCO ].boundaryFieldRef()[patchID][facei] = J_CO;
    ablationFlux[idCO2].boundaryFieldRef()[patchID][facei] = J_CO2;
    ablationFlux[idCN ].boundaryFieldRef()[patchID][facei] = J_CN;
    ablationFlux[idC3 ].boundaryFieldRef()[patchID][facei] = J_C3;
    ablationFlux[idO  ].boundaryFieldRef()[patchID][facei] = J_O;
    ablationFlux[idO2 ].boundaryFieldRef()[patchID][facei] = J_O2;
    ablationFlux[idN  ].boundaryFieldRef()[patchID][facei] = J_N;
    ablationFlux[idN2 ].boundaryFieldRef()[patchID][facei] = J_N2;

    blowing_rate[facei] = -w3 - w4 - w7 - w12 - w15 - w17 - w18 - w20 - 3*w21 + 3*w22;
    
    mass_flux_ablation.boundaryFieldRef()[patchID][facei] = -blowing_rate[facei] * composition.W(idC)/1000.0;
    //Info << "Blowing rate: " << blowing_rate[facei] << endl;
    
    //Info << "Phi: " << phi.boundaryFieldRef()[patchID][facei] << endl;
    //phi.boundaryFieldRef()[patchID][facei] -= mass_flux_ablation.boundaryField()[patchID][facei]*magSf[facei]; // I believe + produces negative velocity? Hard to say currently.
    U.boundaryFieldRef()[patchID][facei] = blowing_rate[facei] * composition.W(idC)/1000.0 / thermo.rho().boundaryField()[patchID][facei] * -norm;
    //Info << "Phi2: " << phi.boundaryFieldRef()[patchID][facei] << endl;

   // Info << "Source term C3: " << ablationFlux[idC3 ].boundaryFieldRef()[patchID][facei] << endl;
    //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;

    
}

//Pout << thermo.kappa().boundaryFieldRef()[patchID][1] / thermo.CpMix().boundaryFieldRef()[patchID][1] << endl;
//Pout << thermo.alphatr().boundaryFieldRef()[patchID][1] << endl;

// TRY: thermo_.kappa()*Le_ / this->thermo_.CpMix() INSTEAD OF THERMO.ALPHA()

