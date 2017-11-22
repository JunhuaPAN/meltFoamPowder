   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	dimensionedScalar pi = constant::mathematical::pi;

    // solid -> phase 1
    // liquid -> phase 2

    // Reading density rho
    dimensionedScalar rho(transportProperties.lookup("rho"));

    // Reading thermal conductivity lambda
    dimensionedScalar a0(transportProperties.lookup("a0"));
    dimensionedScalar a1(transportProperties.lookup("a1"));
    //dimensionedScalar lambdaS(transportProperties.lookup("lambdaS"));
    //dimensionedScalar lambdaL(transportProperties.lookup("lambdaL"));

    // Thermal conductivity reduction for powder
    dimensionedScalar kfactorPP(transportProperties.lookup("kfactorPP")); // Powder-powder
    dimensionedScalar kfactorPS(transportProperties.lookup("kfactorPS")); // Powder-solid

    // Reading heat capacity cp
    dimensionedScalar b0(transportProperties.lookup("b0"));
    dimensionedScalar b1(transportProperties.lookup("b1"));
    //dimensionedScalar cpS(transportProperties.lookup("cpS"));
    //dimensionedScalar cpL(transportProperties.lookup("cpL"));

    // Reading kinematic viscosity
    dimensionedScalar mu0(transportProperties.lookup("mu0"));
    dimensionedScalar A(transportProperties.lookup("A"));
    //dimensionedScalar nuS(transportProperties.lookup("nuS"));
    //dimensionedScalar nuL(transportProperties.lookup("nuL"));

    // Reading latent heat of fusion hs
    dimensionedScalar hs(transportProperties.lookup("hs"));

    // Reading solid bound of melting temperature Ts
    dimensionedScalar Ts(transportProperties.lookup("Ts"));

    // Reading liquid bound of melting temperature Tl
    dimensionedScalar Tl(transportProperties.lookup("Tl"));

    // Reading volume expansion factor beta
    dimensionedScalar beta(transportProperties.lookup("beta"));

    // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));

    // Reading small D'arcy-type source term constant DCs
    dimensionedScalar DCs(transportProperties.lookup("DCs"));
