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

    // Reading spot size w
    dimensionedScalar w(laserSettings.lookup("w"));

    // Layer thickness
    dimensionedScalar th(laserSettings.lookup("th"));


