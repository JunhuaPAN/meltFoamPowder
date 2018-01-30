// Melt pool dimension placeholders
scalar Lx, Ly, Lz;

// Mesh cell coordinates for local processor
scalarField cx = mesh.C().component(0);
scalarField cy = mesh.C().component(1);
scalarField cz = mesh.C().component(2);

// Default min/max coordinate values for local processor
scalar xmin = max(cx);
scalar xmax = min(cx);
scalar ymin = max(cy);
scalar ymax = min(cy);
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

// Set minimum value to quality as part of melt pool
scalar minfrac = 0.1;

// Looping on alpha loops on cells
scalar alphamax = 0;
forAll(alpha, I)
{
  if(alpha[I] > minfrac)
  {
    if(cx[I] < xmin)
      xmin = cx[I];
    if(cx[I] > xmax)
      xmax = cx[I];
    if(cy[I] < ymin)
      ymin = cy[I];
    if(cy[I] > ymax)
      ymax = cy[I];
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


