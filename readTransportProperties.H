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

    // Thermal conductivity reduction for powder
    dimensionedScalar kfactorPP(transportProperties.lookup("kfactorPP")); // Powder-powder
    dimensionedScalar kfactorPS(transportProperties.lookup("kfactorPS")); // Powder-solid

    // Reading heat capacity cp
    dimensionedScalar b0(transportProperties.lookup("b0"));
    dimensionedScalar b1(transportProperties.lookup("b1"));

    // Reading kinematic viscosity
    dimensionedScalar mu0(transportProperties.lookup("mu0"));
    dimensionedScalar A(transportProperties.lookup("A"));

    // Reading latent heat of fusion hs
    dimensionedScalar hs(transportProperties.lookup("hs"));

    // Reading latent heat of vaporization
    dimensionedScalar hv(transportProperties.lookup("hv"));

    // Reading solid bound of melting temperature Ts
    dimensionedScalar Ts(transportProperties.lookup("Ts"));

    // Reading liquid bound of melting temperature Tl
    dimensionedScalar Tl(transportProperties.lookup("Tl"));

    // Reading vaporization temperature ranges
    dimensionedScalar Tv1(transportProperties.lookup("Tv1"));
    dimensionedScalar Tv2(transportProperties.lookup("Tv2"));

    // Reading volume expansion factor beta
    dimensionedScalar beta(transportProperties.lookup("beta"));

    // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));

    // Reading small D'arcy-type source term constant DCs
    dimensionedScalar DCs(transportProperties.lookup("DCs"));
