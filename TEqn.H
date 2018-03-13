// Solving the energy equation

// Create surface field Tsurf which is T interpolated at the face centers of the grid
surfaceScalarField Tsurf(fvc::interpolate(T));

// Interpolate the pflag and alpha fields to face centers
surfaceScalarField pflagSurf(fvc::interpolate(pflag));
surfaceScalarField alphaSurf(fvc::interpolate(alpha));

// Assign conductivity based on pflagSurf - default value is solid phase conductivity
surfaceScalarField lambdaSurf(a0 + a1*Tsurf);
forAll(lambdaSurf,I)
{
  // If pflagSurf = 1 then it is powder-powder interface and should have reduced
  // conductivity of kfactorPP
  if(pflagSurf[I] == 1)
    lambdaSurf[I] = kfactorPP.value()*lambdaSurf[I];
  // This handles pflagSurf > 0 and < 1 since == 1 is already tested
  // Indicates the interface is powder-fluid or powder-solid
  else if(pflagSurf[I] > 0)
  {
    // Reduce the conductivity if interpolation temperature is less than
    if(Tsurf[I] < Tl.value())
      lambdaSurf[I] = kfactorPS.value()*lambdaSurf[I];

  }
}

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      // Melting phase change
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*fvm::ddt(T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*(U & fvc::grad(T))
      // Vaporization Phase change
      + hv*4.0*exp(-pow(4.0*(T-Tvap)/(Tv2-Tv1),2))/Foam::sqrt(pi)/(Tv2-Tv1)*fvm::ddt(T)
      + hv*4.0*exp(-pow(4.0*(T-Tvap)/(Tv2-Tv1),2))/Foam::sqrt(pi)/(Tv2-Tv1)*(U & fvc::grad(T))
      // Conduction
      - fvm::laplacian(lambdaSurf/rho, T)
      // Laser heating source
      - laserSource
    );

    TEqn.relax();
    TEqn.solve();

    // Update liquid fraction
    alpha = 0.5*Foam::erf(4.0*(T-Tmelt)/(Tl-Ts))+scalar(0.5);

    // Update Tvisc
    Tvisc = Foam::max(Tl, T);

    // Update pflag to account for melted powder
    forAll(pflag,I)
    {
      if(T[I] > Tl.value())
        pflag[I] = 0;
    };

    // Update heat capacity
    cp = b0 + b1*T;

    // Update kinematic viscosity - temperature Tvisc > Tl
    nu = mu0*exp(A/Tvisc)/rho;

    // Update dynamic viscosity as well
    mu = nu*rho;

    // Update density for Boussinesq approximation
    rhok = 1.0 - Foam::max(beta*(T - Tl),0.0);

    // Upcate momentum sink coefficient
    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
}
