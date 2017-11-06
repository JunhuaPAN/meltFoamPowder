   IOdictionary laserSettings
    (
        IOobject
        (
            "laserSettings",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Reading power P
    dimensionedScalar P(laserSettings.lookup("P"));

    // Reading effective laser absorptivity beta
    dimensionedScalar Efrac(laserSettings.lookup("Efrac"));

    // Weibull parameters for power penetration
    dimensionedScalar c1(laserSettings.lookup("c1"));
    dimensionedScalar c2(laserSettings.lookup("c2"));
    dimensionedScalar c3(laserSettings.lookup("c3"));

    // Reading spot size w
    dimensionedScalar w(laserSettings.lookup("w"));

    // Layer thickness
    dimensionedScalar th(laserSettings.lookup("th"));

