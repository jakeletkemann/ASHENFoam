
//==============================================================
// Park ablation phenomonological model
// Implemented as conservative mass flux source terms
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
label idC2d  = composition.species()["C2d"];
label idCNA  = composition.species()["CNA"];
label idCNB  = composition.species()["CNB"];
label idNOA  = composition.species()["NOA"];
label idNOB  = composition.species()["NOB"];
label idNOC  = composition.species()["NOC"];
label idNOD  = composition.species()["NOD"];
label idN2A  = composition.species()["N2A"];
label idN2B  = composition.species()["N2B"];
label idN2C  = composition.species()["N2C"];

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


scalar sublimationEfficiency = 1.0;
scalar Le = 1.0;

//Pout << "sublimation efficiency: " << sublimationEfficiency << endl;
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
D_Output = speciesDiffusion().rhoD(0) / thermo.rho(); // CURRENTLY USES N2 DIFFUSION COEFFICIENT
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

Info << "mapfields test: " << mapSpeciesID[idO] << ", " << mapSpeciesID[idO2] << ", " << mapSpeciesID[idN] << ", " << mapSpeciesID[idC3] << endl;


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
    scalar gamma1  = thermo.rho().boundaryField()[patchID][facei] * v_O * composition.W("C") / composition.W("O");
    scalar k1 = gamma1 * 0.63 * Foam::exp(-1160.0/wallTemp);

    scalar gamma2 = 2.0 * thermo.rho().boundaryField()[patchID][facei] * v_O2 * composition.W("C") / composition.W("O2");
    scalar k2 = gamma2 * 0.5;

    scalar gamma3 = thermo.rho().boundaryField()[patchID][facei] * v_N * composition.W("C") / composition.W("N");
    scalar k3 = gamma3 * 0.3;

    scalar gamma4 = composition.W("C3")/1000.0 * 5.19e14;
    scalar k4 = v_C3/(R*wallTemp) * gamma4 * Foam::exp(-90845.0/wallTemp);

    scalar k5 = thermo.rho().boundaryField()[patchID][facei] * v_C3;
 
}


// Number of independent species at wall
const label nIndep = indepSpecies.size();

// Per-face unknowns
List<scalarField> Ywall(patch.size(), scalarField(nIndep, 0.0));
scalarField YN2(patch.size(), 0.0);

// Residuals and updates
List<scalarField> Rwall(patch.size(), scalarField(nIndep, 0.0));
List<scalarField> dYwall(patch.size(), scalarField(nIndep, 0.0));
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

//
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
        label cellI = faceCells[facei];
        const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
        //Ygrad[facei][Yi] = YgradField.boundaryField()[patchID][facei] & norm;
        Ygrad[facei][Yi] = YgradField.boundaryField()[patchID][facei] & norm;
        //Info << "Ygrad: "<< Ygrad[facei][Yi] << endl;
      }

      //Info << Ygrad[0] << endl;
    }
    volVectorField YgradField = fvc::grad(composition.Y("N2"));

    forAll(Tw, facei)
    {
        scalar sumGrad = 0.0;
        forAll(indepSpecies, Yi)
        {
            sumGrad += Ygrad[facei][Yi];
        }
      label cellI = faceCells[facei];
      const vector norm = -mesh.Sf().boundaryField()[patchID][facei] / magSf[facei];
      //YgradN2[facei] = YgradField.boundaryField()[patchID][facei] & norm;
      YgradN2[facei] = YgradField.boundaryField()[patchID][facei] & norm;
        //YgradN2[facei] = -sumGrad;
        //Info << "YgradN2: " << YgradN2[facei] << endl;
    }
    

    

    // --------------------------------------------------
    // Residual evaluation
    // --------------------------------------------------
    globalRes = 0.0;

    forAll(Tw, facei)
    {

    
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
    scalar gamma1  = thermo.rho().boundaryField()[patchID][facei] * v_O * composition.W("C") / composition.W("O");
    scalar k1 = gamma1 * 0.63 * Foam::exp(-1160.0/wallTemp);

    scalar gamma2 = 2.0 * thermo.rho().boundaryField()[patchID][facei] * v_O2 * composition.W("C") / composition.W("O2");
    scalar k2 = gamma2 * 0.5;

    scalar gamma3 = thermo.rho().boundaryField()[patchID][facei] * v_N * composition.W("C") / composition.W("N");
    scalar k3 = gamma3 * 0.3;

    scalar gamma4 = composition.W("C3")/1000.0 * 5.19e14;
    scalar k4 = v_C3/(R*wallTemp) * gamma4 * Foam::exp(-90845.0/wallTemp) * sublimationEfficiency;

    scalar k5 = thermo.rho().boundaryField()[patchID][facei] * v_C3 * sublimationEfficiency;
    
    
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
scalar D = Dfield[cellI]/thermo.rho()[cellI]; // * 0.225 / Le;

    
        computeResidualsPark
        (
        patchID, facei,
        k1,    k2,    k3,    k4,    k5,
        idO,   idO2,   idN,   idN2,
        idCO,   idCN,   idCO2,   idC3,    idC,
        nSpecies, thermo.rho().boundaryField()[patchID][facei], W, delta, D,
        indepSpecies, Ygrad[facei], Ycell, Ywall[facei], YN2, YgradN2, J_diff, mapSpeciesID, Rwall[facei]
        );
        
        if (max(mag(Rwall[facei])) > 1.0e-9)
        {
        //Info << "|R|: " << max(mag(Rwall[facei])) << endl;
        //Info << "face index: " << facei << endl;        
        //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;
        }
        //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;

        globalRes = max(globalRes, max(mag(Rwall[facei])));
    }
//Info << "Next iteration! " << endl;
//Info << "|R|: " << globalRes << endl;

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
    scalar gamma1  = thermo.rho().boundaryField()[patchID][facei] * v_O * composition.W("C") / composition.W("O");
    scalar k1 = gamma1 * 0.63 * Foam::exp(-1160.0/wallTemp);

    scalar gamma2 = 2.0 * thermo.rho().boundaryField()[patchID][facei] * v_O2 * composition.W("C") / composition.W("O2");
    scalar k2 = gamma2 * 0.5;

    scalar gamma3 = thermo.rho().boundaryField()[patchID][facei] * v_N * composition.W("C") / composition.W("N");
    scalar k3 = gamma3 * 0.3;

    scalar gamma4 = composition.W("C3")/1000.0 * 5.19e14;
    scalar k4 = v_C3/(R*wallTemp) * gamma4 * Foam::exp(-90845.0/wallTemp) * sublimationEfficiency;

    scalar k5 = thermo.rho().boundaryField()[patchID][facei] * v_C3 * sublimationEfficiency;
    
    
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

        List<scalarField> J(nIndep, scalarField(nIndep, 0.0));


//const scalar& D = speciesDiffusion().Di(0)[cellI];
scalar D = Dfield[cellI]/thermo.rho()[cellI]; // * 0.225 / Le;

        forAll(indepSpecies, j)
        {
            scalarField Ypert = Ywall[facei];
            scalarField Ygradpert = Ygrad[facei];
            scalar eps = 1.0e-12;
            Ypert[j] += eps;
            Ygradpert[j] -= eps/delta;
            
            
            
//                forAll(indepSpecies, k)
//                {
//                    composition.Y(indepSpecies[k])
//                        .boundaryFieldRef()[patchID][facei]
//                        = Ypert[k];
//                }
//            Info << "Correcting diffusion..." << endl;
//            speciesDiffusion().correct();
//            
//            forAll(Tw, facei)
//            {
//            
//            forAll(indepSpecies,k)
//                {
//                  J_diff[facei][k] = speciesDiffusion().Jcorr(indepSpecies[k])
//                      .boundaryField()[patchID][facei]
//                    & (-mesh.Sf().boundaryField()[patchID][facei]
//                       / mag(mesh.Sf().boundaryField()[patchID][facei]));
//                }
//        
//            }
            scalarField Rpert(nIndep, 0.0);
            computeResidualsPark
            (
            patchID, facei,
            k1,    k2,    k3,    k4,     k5,
            idO,   idO2,   idN,   idN2,
            idCO,   idCN,   idCO2,   idC3,   idC, 
            nSpecies, thermo.rho().boundaryField()[patchID][facei], W, delta, D,
            indepSpecies, Ygradpert, Ycell, Ypert, YN2, YgradN2, J_diff, mapSpeciesID, Rpert
            );

            forAll(indepSpecies, i)
            {
                J[i][j] = (Rpert[i] - Rwall[facei][i]) / eps;
            }
            
        }
        

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
        
    //Info << "x-location: " << mesh.Cf().boundaryField()[patchID][facei] << endl;
    }
    

    newtonIter++;
}
    forAll(Tw, facei)
    {
    //Info << "Y: " << Ywall[facei] << endl;
    }
    
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
//    if (facei > 146)
//    {
//      Ywall[facei][k] = Ywall[facei-1][k];
//      YN2[facei] = YN2[facei-1];
//    }
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
    //Info << "Y: " << Ywall[facei] << endl;
    }

//Info << "Y: " << Ywall[0] << endl;

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
    scalar gamma1  = thermo.rho().boundaryField()[patchID][facei] * v_O * composition.W("C") / composition.W("O");
    scalar k1 = gamma1 * 0.63 * Foam::exp(-1160.0/wallTemp);

    scalar gamma2 = 2.0 * thermo.rho().boundaryField()[patchID][facei] * v_O2 * composition.W("C") / composition.W("O2");
    scalar k2 = gamma2 * 0.5;

    scalar gamma3 = thermo.rho().boundaryField()[patchID][facei] * v_N * composition.W("C") / composition.W("N");
    scalar k3 = gamma3 * 0.3;

    scalar gamma4 = composition.W("C3")/1000.0 * 5.19e14;
    scalar k4 = v_C3/(R*wallTemp) * gamma4 * Foam::exp(-90845.0/wallTemp);

    scalar k5 = thermo.rho().boundaryField()[patchID][facei] * v_C3;

    scalar mdot1  = k1  * Ywall[facei][mapSpeciesID[idO]];
    scalar mdot2  = k2  * Ywall[facei][mapSpeciesID[idO2]];
    scalar mdot3  = k3  * Ywall[facei][mapSpeciesID[idN]];
    scalar mdot4  = k4 * sublimationEfficiency;
    scalar mdot5  = k5 * Ywall[facei][mapSpeciesID[idC3]] * sublimationEfficiency;
    
    //mdot3 = 0;

    // ----------------------------------------------------------
    // Gas-phase production (kg/m2/s)
    // ----------------------------------------------------------
    scalar J_CO  = (mdot1*W[idCO]/W[idC] + mdot2*W[idCO]/W[idC]);
    scalar J_CO2 = 0.0;
    scalar J_CN  = (mdot3*W[idCN]/W[idC]);
    scalar J_C3  = (mdot4 - mdot5);
    scalar J_O   = (-mdot1*W[idO]/W[idC]); 
    scalar J_O2  = (-mdot2*W[idO2]/W[idC]);
    scalar J_N   = (-mdot3*W[idN]/W[idC]);
    scalar J_N2  = 0.0;
    
    J_CO = mdot1 + mdot2;
    J_CN = mdot3;

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

    blowing_rate[facei] = mdot1 + mdot2 + mdot3 + mdot4 - mdot5;

    mass_flux_ablation.boundaryFieldRef()[patchID][facei] = blowing_rate[facei];
    //Info << "Blowing rate: " << blowing_rate[facei] << endl;
    
    //Info << "Phi: " << phi.boundaryFieldRef()[patchID][facei] << endl;
    //phi.boundaryFieldRef()[patchID][facei] -= mass_flux_ablation.boundaryField()[patchID][facei]*magSf[facei]; // I believe + produces negative velocity? Hard to say currently.
    U.boundaryFieldRef()[patchID][facei] = -blowing_rate[facei] / thermo.rho().boundaryField()[patchID][facei] * -norm;
    //Info << "Phi2: " << phi.boundaryFieldRef()[patchID][facei] << endl;

    //Info << "Source term: " << ablationSource[idO ][cellI] << endl;
    //Info << "x-location: " << mesh.C().boundaryField()[patchID][cellI] << endl;

}

//Info << thermo.kappatr().boundaryFieldRef()[patchID][1] << endl;
//Info << thermo.CpMix().boundaryFieldRef()[patchID][1] << endl;
//Info << thermo.mu().boundaryFieldRef()[patchID][1] << endl;
//Pout << thermo.alphatr().boundaryFieldRef()[patchID][1] << endl;

// TRY: thermo_.kappa()*Le_ / this->thermo_.CpMix() INSTEAD OF THERMO.ALPHA()

