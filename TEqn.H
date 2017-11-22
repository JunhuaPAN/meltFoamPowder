// Solving the energy equation

// Create surface field Tsurf which is T interpolated at the face centers of the
// grid
surfaceScalarField Tsurf = fvc::interpolate(T);

// Interpolate the pflag and alpha fields to face centers
surfaceScalarField pflagSurf = fvc::interpolate(pflag);
surfaceScalarField alphaSurf = fvc::interpolate(alpha);

// Assign conductivity based on pflagSurf - default to material conductivity
surfaceScalarField lambdaSurf = a0 + a1*Tsurf;
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
    // Let conductivity be in between powder-solid (worst) and powder-fluid
    // (best, = 1.0) defined as a mixture governed by alpha
    lambdaSurf[I] = (kfactorPS.value() + (1.0-kfactorPS.value())*alphaSurf[I])*lambdaSurf[I];
  }
}

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*fvm::ddt(T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*(U & fvc::grad(T))
      - fvm::laplacian(lambdaSurf/rho, T)
      //- fvm::laplacian(lambda/rho, T)
      - laserSource
      //- I0*(1/th)/rho*exp(-2*R2/(w*w))
      //- I0*(1/th)/rho
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
    //cp = alpha*cpL+(1.0-alpha)*cpS;
    
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
