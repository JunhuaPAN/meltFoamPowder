// Melt pool dimension placeholders
scalar Lx, Ly, Lz;

// Define the scan direction by angle ang using the velocity vectors for the
// current time step
scalar ang = Foam::atan2(Vlaser[1],Vlaser[0]);

Info << "Scan angle = " << (180.0/3.14)*ang << endl;

// Mesh cell coordinates for local processor
scalarField cx(mesh.C().component(0));
scalarField cy(mesh.C().component(1));
scalarField cz(mesh.C().component(2));
//scalarField cx = mesh.C().component(0);
//scalarField cy = mesh.C().component(1);
//scalarField cz = mesh.C().component(2);

// Rotate the x and y mesh coordinate fields by NEGATIVE ang
ang = -ang;
scalarField cxrot(Foam::cos(ang)*cx - Foam::sin(ang)*cy);
//scalarField cxrot = Foam::cos(ang)*cx - Foam::sin(ang)*cy;
scalarField cyrot(Foam::sin(ang)*cx + Foam::cos(ang)*cy);
//scalarField cyrot = Foam::sin(ang)*cx + Foam::cos(ang)*cy;

// Default min/max coordinate values for local processor
scalar xmin = max(cxrot);
scalar xmax = min(cxrot);
scalar ymin = max(cyrot);
scalar ymax = min(cyrot);
scalar zmin = max(cz);
scalar zmax = min(cz);

// Make these coordinate values global
// The mixing of min/max is NOT a typo or bug - you want default xmin coordinate
// to bet set to the max X coordinate so that any deviation is registered 
reduce(xmin, maxOp<scalar>());
reduce(xmax, minOp<scalar>());
reduce(ymin, maxOp<scalar>());
reduce(ymax, minOp<scalar>());
reduce(zmin, maxOp<scalar>());
reduce(zmax, minOp<scalar>());

// Set minimum value to qualify as part of melt pool
scalar minfrac = 0.99;

// Looping on alpha loops on cells
scalar alphamax = 0;
forAll(alpha, I)
{
  if(alpha[I] > minfrac)
  {
    if(cxrot[I] < xmin)
      xmin = cxrot[I];
    if(cxrot[I] > xmax)
      xmax = cxrot[I];
    if(cyrot[I] < ymin)
      ymin = cyrot[I];
    if(cyrot[I] > ymax)
      ymax = cyrot[I];
    if(cz[I] < zmin)
      zmin = cz[I];
    if(cz[I] > zmax)
      zmax = cz[I];
  }
  if(alpha[I] > alphamax)
    alphamax = alpha[I];
}

// For parallel case, combine min/max values over all procs
reduce(xmin, minOp<scalar>());
reduce(xmax, maxOp<scalar>());
reduce(ymin, minOp<scalar>());
reduce(ymax, maxOp<scalar>());
reduce(zmin, minOp<scalar>());
reduce(zmax, maxOp<scalar>());
reduce(alphamax, maxOp<scalar>());

// If melt pool does not exist
if(alphamax < 1.0)
{
  Lx = 0;
  Ly = 0;
  Lz = 0;
}
else
{
  // Melt pool dimensions in microns
  Lx = 1e4*(xmax - xmin);
  Ly = 1e4*(ymax - ymin);
  Lz = 1e4*(zmax - zmin);
}

// setprecision(x) uses x decimal points with 'fixed format'
//outfile << std::setprecision(8);

// 'fixed' format writes without scientific notation (more readable)
//outfile << std::fixed;

// For writing in parallel
// Format: <time, xlaser, ylaser, Lx, Ly, Lz>
if(Pstream::master())
{
  OFstream& outfile = *outfilePtr;
  outfile << runTime.timeName() << tab
          << Xlaser << tab
          << Ylaser << tab
          << Lx << tab
          << Ly << tab
          << Lz << endl;
}


