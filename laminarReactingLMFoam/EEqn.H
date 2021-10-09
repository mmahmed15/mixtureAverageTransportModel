{

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
            - fvm::laplacian(transport.alphaE(), he)
            + transport.Econduction() 
            )
          : (
               - fvm::laplacian(transport.alpha(), he)
               + transport.Hconduction()
            )
        )
      + transport.JHs()      
     ==
        reaction->Sh()
      + radiation->Sh(thermo)
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
