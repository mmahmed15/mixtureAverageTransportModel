{
    combustion->correct();

    // dQ is only for diagnostic purposes
    dQ = combustion->Sh();
//        + mesh.lookupObject<volScalarField>("sootOxidationRate") * singleMixture.qFuel()
//        + mesh.lookupObject<volScalarField>("sootOxidationRate") * (h0_CO - h0_F);

    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i == fuelIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + fvm::div(phi, Yi, "div(phi,Yi)")
              + transport->Yflux(Yi) //MAT model
              ==
//                parcels.SYi(i, Yi)
//              + surfaceFilm.Srho(i)
               combustion->R(Yi)
//              - mesh.lookupObject<volScalarField>("sootFormationRate")
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }

        else if (i != fuelIndex && i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + fvm::div(phi, Yi, "div(phi,Yi)")
              + transport->Yflux(Yi) //MAT model
              ==
//                parcels.SYi(i, Yi)
//              + surfaceFilm.Srho(i)
               combustion->R(Yi)
//              + singleMixture.specieStoichCoeffs()[i] * mesh.lookupObject<volScalarField>("sootOxidationRate")
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }

    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    forAll(Y, i){
        Info << setw(8) << Y[i].name() << "\tmin/ave/max "<< " =\t"
             << setw(12) << setprecision(5) << gMin(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gAverage(Y[i]) << "\t"
             << setw(12) << setprecision(5) << gMax(Y[i]) << nl;
    }



    radiation->correct();

    const volScalarField::Internal T4 = pow(thermo.T(),4).ref();
    const volScalarField::Internal rad = radiation->Ru() - radiation->Rp().ref()*T4;
    volScalarField::Internal & radField = dQrad.ref();
    radField = rad;

 
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he, "div(phi,he)")
      + (
            he.name() == "e"
          ? (
              fvc::div
              (
                  fvc::absolute(phi/fvc::interpolate(rho), U),
                  p,
                "div(phiv,p)"
              )
            - fvm::laplacian(transport->alphaE(), he)
            + transport->Econduction() 
            )
          : (
               - fvm::laplacian(transport->alpha(), he)
//               + transport->Hconduction()
            )
        )
//      + transport->JHs() 
     ==
        combustion->Sh()
//      + mesh.lookupObject<volScalarField>("sootOxidationRate") * singleMixture.qFuel()
//      + mesh.lookupObject<volScalarField>("sootOxidationRate") * (h0_CO - h0_F)
      + radiation->Sh(thermo)
//      + parcels.Sh(he)
//      + surfaceFilm.Sh()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

}
