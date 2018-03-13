// Solving the energy equation

// Create surface field Tsurf which is T interpolated at the face centers of the
// grid
surfaceScalarField Tsurf(fvc::interpolate(T));
//surfaceScalarField Tsurf = fvc::interpolate(T);

// Interpolate the pflag and alpha fields to face centers
surfaceScalarField pflagSurf(fvc::interpolate(pflag));
surfaceScalarField alphaSurf(fvc::interpolate(alpha));
//surfaceScalarField pflagSurf = fvc::interpolate(pflag);
//surfaceScalarField alphaSurf = fvc::interpolate(alpha);

// Assign conductivity based on pflagSurf - default to material conductivity
surfaceScalarField lambdaSurf(a0 + a1*Tsurf);
//surfaceScalarField lambdaSurf = a0 + a1*Tsurf;
//surfaceScalarField lambdaSurf = pflagSurf;
forAll(lambdaSurf,I)
{
  // If pflagSurf = 1 then it is powder-powder interface and should have reduced
  // conductivity of kfactorPP
  if(pflagSurf[I] == 1)
    lambdaSurf[I] = kfactorPP.value()*lambdaSurf[I];
  // This handles pflagSurf > 0 and < 1 since == 1 is already tested
  // This case is powder-fluid or powder-solid
  else if(pflagSurf[I] > 0)
  {
    // Reduce the conductivity if interpolation temperature is less than
    // liquidus
    if(Tsurf[I] < Tl.value())
      lambdaSurf[I] = kfactorPS.value()*lambdaSurf[I];

    // Otherwise, Fully liquid-powder contact should take on the higher melt
    // conductivity because the liquid contacts it fully and fills in gaps

    /*
    // Conduction facilitated if interface is very liquidy
    if(alphaSurf[I] >= alphaCondLim)
      lambdaSurf[I] = (kfactorPS.value() + (1.0-kfactorPS.value())*alphaSurf[I])*lambdaSurf[I];
    // Otherwise, the conduction should retain the powder-solid interaction
    // value
    else
      lambdaSurf[I] = kfactorPS.value()*lambdaSurf[I];
    */
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
      //- fvm::laplacian(lambda/rho, T)
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
    // Loop over cells to determine if heat capacity should be modified for a
    // vaporized phase change region
    /*
    forAll(T, I)
    {
      if(T[I] >= 3300.0 && T[I] <= 3400.0)
        //continue;
        cp[I] = 3.0e7; // deltaH/(Tv2 - Tv1)
        //cp[I] = 3.0e10; // deltaH/(Tv2 - Tv1)
        //cp[I] = 3000000.0e4/(3400.0 - 3300.0); // deltaH/(Tv2 - Tv1)
    }
    */

    
    // Update conductivity
    //lambda = pflag*kfactor*(a0 + a1*T) + (1-pflag)*(a0 + a1*T);
    //lambda = pflag*0.1*(a0 + a1*T) + (1-pflag)*(a0 + a1*T);
    //lambda = a0 + a1*T;
    //lambda = alpha*lambdaL+(1.0-alpha)*lambdaS;

    // Update kinematic viscosity
    nu = mu0*exp(A/Tvisc)/rho;
    //nu = mu0*exp(A/T)/rho;
    //nu = alpha*nuL+(1.0-alpha)*nuS;

    // Update dynamic viscosity as well
    mu = nu*rho;

    // Update density for Boussinesq approximation
    rhok = 1.0 - Foam::max(beta*(T - Tl),0.0);

    // Upcate momentum sink coefficient
    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
}
