
void updateSurfaceCoverageExplicit
(
    const label patchID,
    const label facei,
    const scalar dtSurf,
    const scalar B,
    const scalarField& X_O,
    const scalarField& X_O2,
    const scalarField& X_N,
    volScalarField& X_Os,
    volScalarField& X_Os_star,
    volScalarField& X_Ns,
    const scalar k1,  const scalar k2,  const scalar k3,  const scalar k4,
    const scalar k5,  const scalar k6,  const scalar k7,  const scalar k8,
    const scalar k9,  const scalar k10, const scalar k11, const scalar k12,
    const scalar k13, const scalar k14, const scalar k15, const scalar k16,
    const scalar k17, const scalar k18, const scalar k19, const scalar k20
)
{
    // --- Initial values ---
    scalar Os  = X_Os.boundaryField()[patchID][facei];
    scalar Oss = X_Os_star.boundaryField()[patchID][facei];
    scalar Ns  = X_Ns.boundaryField()[patchID][facei];

    const scalar XO  = X_O[facei];
    const scalar XO2 = X_O2[facei];
    const scalar XN  = X_N[facei];

    scalar Os_old = 0.0;
    scalar Oss_old = 0.0;
    scalar Ns_old = 0.0;
    
    const label maxIter = 10000;
    const scalar dtPseudo = dtSurf/maxIter;
    const scalar tol = 1e-10;

    for (label iter = 0; iter < maxIter; ++iter)
    {
        const scalar Xs = max(B - Os - Oss - Ns, scalar(0));

        // --- Reaction rates ---
        const scalar w1  = k1  * XO  * Xs;
        const scalar w2  = k2  * Os;
        const scalar w3  = k3  * XO  * Os;
        const scalar w4  = k4  * XO  * Os;
        const scalar w5  = k5  * XO  * Xs;
        const scalar w6  = k6  * Oss;
        const scalar w7  = k7  * XO  * Oss;
        const scalar w8  = k8  * sqr(Oss);
        const scalar w9  = k9  * sqr(Os);
        const scalar w10 = k10 * XN  * Xs;
        const scalar w11 = k11 * Ns;
        const scalar w12 = k12 * XN  * Ns;
        const scalar w13 = k13 * XN  * Ns;
        const scalar w14 = k14 * sqr(Ns);
        const scalar w15 = k15 * Ns;
        const scalar w16 = k16 * XO2 * sqr(Xs);
        const scalar w17 = k17 * XO2 * Os;
        const scalar w18 = k18 * XO2 * Os;
        const scalar w19 = k19 * XO2 * sqr(Xs);
        const scalar w20 = k20 * XO2 * Oss;

        // --- Net rates ---
        const scalar RO =
              w1 - w2 - w3 - w4 - 2.0*w9
            + 2.0*w16 - w17 - w18;

        const scalar ROs =
              w5 - w6 - w7 - 2.0*w8
            + 2.0*w19 - w20;

        const scalar RN =
              w10 - w11 - w12 - w13
            - 2.0*w14 - w15;

        // --- Forward Euler update ---
        const scalar dOs  = dtPseudo * RO;
        const scalar dOss = dtPseudo * ROs;
        const scalar dNs  = dtPseudo * RN;

        Os  += dOs;
        Oss += dOss;
        Ns  += dNs;

    }

    // --- Write back ---
    X_Os.boundaryFieldRef()[patchID][facei]      = Os;
    X_Os_star.boundaryFieldRef()[patchID][facei] = Oss;
    X_Ns.boundaryFieldRef()[patchID][facei]      = Ns;
    
}

void computeResiduals
(
    const label patchID,
    const label facei,
    const scalar k1,  const scalar k2,  const scalar k3,  const scalar k4,
    const scalar k5,  const scalar k6,  const scalar k7,  const scalar k8,
    const scalar k9,  const scalar k10, const scalar k11, const scalar k12,
    const scalar k13, const scalar k14, const scalar k15, const scalar k16,
    const scalar k17, const scalar k18, const scalar k19, const scalar k20,
    const scalar k21, const scalar k22, const scalar sublimationEfficiency,
    const scalar idO, const scalar idO2, const scalar idN, const scalar idN2,
    const scalar idCO, const scalar idCN, const scalar idCO2, const scalar idC3,  
    const scalar PO, const scalar PO2, const scalar PN, 
    const scalar nspecies,
    const scalar rhoW,
    const scalarField W, 
    const scalar B,
    const scalar delta,
    scalar D,
    const scalar X_Os,
    const scalar X_Os_star,
    const scalar X_Ns,
    const List<label> indepSpecies,
    const scalarField Ygrad,
    const scalarField Ycell, 
    const scalarField& Ywall,
    const scalarField Yindep,
    const scalarField Ygradindep,
    const List<scalarField> J_diff, 
    const List<label> mapSpeciesID,
    scalarField& R
)
{
    scalar WN2 = 28.0;

    scalar X_O  = rhoW * Ywall[mapSpeciesID[idO]]  / (W[idO]/1000.0);
    scalar X_N  = rhoW * Ywall[mapSpeciesID[idN]]  / (W[idN]/1000.0);
    scalar X_O2 = rhoW * Ywall[mapSpeciesID[idO2]]  / (W[idO2]/1000.0);
    scalar X_N2 = rhoW * Yindep[facei]   / (WN2/1000.0);
    scalar X_C3 = rhoW * Ywall[mapSpeciesID[idC3]]  / (W[idC3]/1000.0);
    
    scalar X_s = B - X_Os - X_Os_star - X_Ns;

    scalar w1  = k1  * X_O  * X_s;
    scalar w2  = k2  * X_Os;
    scalar w3  = k3  * X_O  * X_Os;
    scalar w4  = k4  * X_O  * X_Os;
    scalar w5  = k5  * X_O  * X_s;
    scalar w6  = k6  * X_Os_star;
    scalar w7  = k7  * X_O  * X_Os_star;
    scalar w8  = k8  * sqr(X_Os_star);
    scalar w9  = k9  * sqr(X_Os);
    scalar w10 = k10 * X_N  * X_s;
    scalar w11 = k11 * X_Ns;
    scalar w12 = k12 * X_N  * X_Ns;
    scalar w13 = k13 * X_N  * X_Ns;
    scalar w14 = k14 * sqr(X_Ns);
    scalar w15 = k15 * X_Ns;
    scalar w16 = k16 * X_O2 * sqr(X_s);
    scalar w17 = k17 * X_O2 * X_Os;
    scalar w18 = k18 * X_O2 * X_Os;
    scalar w19 = k19 * X_O2 * sqr(X_s);
    scalar w20 = k20 * X_O2 * X_Os_star;
    scalar w21 = k21 * sublimationEfficiency;
    scalar w22 = k22 * X_C3 * sublimationEfficiency;
    
    
    // ------------------------------
    // omega_dot contribution
    // ------------------------------ 
    
    scalar dOdt       = (-w1 + w2 - w4 - w5 + w6 + w18); 
    scalar dO2dt      = (w8 + w9 - w16 - w18 - w19);
    scalar dNdt       = (-w10 + w11 - w13);
    scalar dN2dt      = (w13 + w14);
    scalar dCOdt      = (w3 + w7 + w17 + w20);
    scalar dCNdt      = (w12 + w15);
    scalar dC3dt      = (w21 - w22);
    scalar dCO2dt     = (w4 + w18);
    
    scalarField ddtList(nspecies, 0.0);
    ddtList[idO]   = dOdt;
    ddtList[idO2]  = dO2dt;
    ddtList[idN]   = dNdt;
    ddtList[idN2]  = dN2dt;
    ddtList[idCO]  = dCOdt;
    ddtList[idCN]  = dCNdt;
    ddtList[idC3]  = dC3dt;
    ddtList[idCO2] = dCO2dt;
    
    scalar mDot = 0.0;

    // ------------------------------
    // m_dot contribution
    // ------------------------------ 
    
    forAll(ddtList, Yi)
      {
        mDot += (W[Yi]/1000.0 * ddtList[Yi]);
      }

    scalarField J_species(indepSpecies.size(), 0.0);
    scalar sumGradY = 0.0;
    
        forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        sumGradY += rhoW*D*(Ywall[m] - Ycell[Ym])/delta;
    }
    sumGradY += rhoW*D*(Yindep[facei] - Ycell[idN2])/delta;
    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        J_species[m] = rhoW*D*(Ywall[m] - Ycell[Ym])/delta - Ywall[m]*sumGradY;
    }
    
    scalarField J_species2(indepSpecies.size(), 0.0);
    scalar sumGradY2 = 0.0;

    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        sumGradY2 += rhoW*D*(-Ygrad[m]);
    }
    sumGradY2 += rhoW*D*(-Ygradindep[facei]);
    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        J_species2[m] = rhoW*D*(-Ygrad[m]); // - Ywall[m]*sumGradY2;
    }
    
    forAll(Ywall, k)
    {
        R[k] = Ywall[k] * mDot - W[indepSpecies[k]]/1000.0 * ddtList[indepSpecies[k]] + J_species2[k];// + J_diff[facei][k];
    }
    R[indepSpecies.size()+0] =  w1 - w2 - w3 - w4 - 2.0*w9 + 2.0*w16 - w17 - w18; // Os
    R[indepSpecies.size()+1] = w5 - w6 - w7 - 2.0*w8 + 2.0*w19 - w20;             // Os*
    R[indepSpecies.size()+2] = w10 - w11 - w12 - w13 - 2.0*w14 - w15;             // Ns
    
}

void computeResidualsPark
(
    const label patchID,
    const label facei,
    const scalar k1,  const scalar k2,  const scalar k3,  const scalar k4, const scalar k5,
    const scalar idO, const scalar idO2, const scalar idN, const scalar idN2,
    const scalar idCO, const scalar idCN, const scalar idCO2, const scalar idC3, const scalar idC,
    const scalar nspecies,
    const scalar rhoW,
    const scalarField W, 
    const scalar delta,
    scalar D,
    const List<label> indepSpecies,
    const scalarField Ygrad,
    const scalarField Ycell, 
    const scalarField& Ywall,
    const scalarField Yindep,
    const scalarField Ygradindep,
    const List<scalarField> J_diff, 
    const List<label> mapSpeciesID,
    scalarField& R
)
{
    scalar WN2 = 28.0;

    scalar mdot1  = k1  * Ywall[mapSpeciesID[idO]];
    scalar mdot2  = k2  * Ywall[mapSpeciesID[idO2]];
    scalar mdot3  = k3  * Ywall[mapSpeciesID[idN]];
    scalar mdot4  = k4;
    scalar mdot5  = k5 * Ywall[mapSpeciesID[idC3]];
    
    //mdot3 = 0; // THIS IS FOR NO NITRIDATION        

    // ------------------------------
    // omega_dot contribution
    // ------------------------------ 
    
    scalar dOdt       = (-mdot1*W[idO]/W[idC]); 
    scalar dO2dt      = (-mdot2*W[idO2]/2.0/W[idC]);
    scalar dNdt       = (-mdot3*W[idN]/W[idC]);
    scalar dCOdt      = (mdot1*W[idCO]/W[idC] + mdot2*W[idCO]/W[idC]);
    scalar dCNdt      = (mdot3*W[idCN]/W[idC]);
    scalar dC3dt      = (mdot4 - mdot5);
    
    scalarField ddtList(nspecies, 0.0);
    ddtList[idO]   = dOdt;
    ddtList[idO2]  = dO2dt;
    ddtList[idN]   = dNdt;
    ddtList[idCO]  = dCOdt;
    ddtList[idCN]  = dCNdt;
    ddtList[idC3]  = dC3dt;
    
    scalar mDot = 0.0;
    // ------------------------------
    // m_dot contribution
    // ------------------------------ 
    
    forAll(ddtList, Yi)
      {
        mDot += ddtList[Yi];
      }
    
    scalarField J_species(indepSpecies.size(), 0.0);
    scalar sumGradY = 0.0;
    
        forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        sumGradY += rhoW*D*(1.0-Ywall[m])*(Ywall[m] - Ycell[Ym])/delta;
    }
    sumGradY += rhoW*D*(1.0-Yindep[facei])*(Yindep[facei] - Ycell[idN2])/delta;
    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        J_species[m] = rhoW*D*(1.0-Ywall[m])*(Ywall[m] - Ycell[Ym])/delta; // - Ywall[m]*sumGradY;
    }

    scalarField J_species2(indepSpecies.size(), 0.0);
    scalar sumGradY2 = 0.0;

    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        sumGradY2 += rhoW*D*(-Ygrad[m]);
    }
    //sumGradY2 += rhoW*D*(-Ygradindep[facei]);
    forAll(indepSpecies, m)
    {
        label Ym = indepSpecies[m];
        J_species2[m] = rhoW*D*(-Ygrad[m]); // - Ywall[m]*sumGradY2;
    }

    forAll(Ywall, k)
    {
        R[k] = Ywall[k] * mDot - ddtList[indepSpecies[k]] + J_species2[k]; // + J_species2[k];// + J_diff[facei][k];  
    }

}

void smoothPatchY
(
    const fvPatch& patch,
    List<scalarField>& Ywall,   // [face][species]
    const scalar omega,         // smoothing strength (0 < omega <= 0.5)
    const label nIter           // smoothing passes
)
{

for (label iter = 0; iter < nIter; ++iter)
{
    for (label f = 1; f < patch.size()-1; ++f)
    {
        forAll(Ywall[f], k)
        {
          if(k!=6) // DO NOT SMOOTH C3
          {
            Ywall[f][k] =
                (1.0 - omega)*Ywall[f][k]
              + omega*0.5*(Ywall[f-1][k] + Ywall[f+1][k]);
          }
        }
    }
}

}


void solveLinearSystem(const List<scalarField>& J,
                       const scalarField& R,
                       scalarField& dY)
{

    const label N = R.size();
    dY = scalarField(N, 0.0);

    scalarField Rcopy = R;
    List<scalarField> Jcopy = J;

    // Forward elimination with partial pivoting
    for (label i = 0; i < N; ++i)
    {
        label p = i;
        scalar maxVal = Foam::mag(Jcopy[i][i]);
        for (label r = i+1; r < N; ++r)
        {
            if (Foam::mag(Jcopy[r][i]) > maxVal)
            {
                maxVal = Foam::mag(Jcopy[r][i]);
                p = r;
            }
        }


        if (maxVal < SMALL)
            FatalErrorInFunction << "Singular matrix";

        if (p != i)
        {

            Foam::Swap(Jcopy[p], Jcopy[i]);
            Foam::Swap(Rcopy[p], Rcopy[i]);

        }

        scalar pivot = Jcopy[i][i];

        for (label j = i+1; j < N; ++j)
        {
            scalar factor = Jcopy[j][i] / pivot;
            for (label k = i; k < N; ++k)
                Jcopy[j][k] -= factor * Jcopy[i][k];
            Rcopy[j] -= factor * Rcopy[i];
        }
    }

    // Back substitution
    for (label i = N-1; i >= 0; --i)
    {
        scalar sum = 0.0;
        for (label j = i+1; j < N; ++j)
            sum += Jcopy[i][j] * dY[j];

        dY[i] = (Rcopy[i] - sum) / Jcopy[i][i];
    }
}
