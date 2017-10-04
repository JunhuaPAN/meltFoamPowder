// Interpolation function
scalar interp(scalar tval, std::vector<scalar> tvec, std::vector<scalar> yvec) {

  // First check that query time is not outside table bounds
  if(tval < tvec[0])
  {
    Info << "Below table range" << endl;
    return yvec[0];
  }
  if(tval > tvec[tvec.size()-1])
  {
    Info << "Above table range" << endl;
    return yvec[yvec.size()-1];
  }

  // If you make it through the above:
  // Loop on time values until you find one that is greater than or equal to the
  // query time, then perform linear interpolation between the two values
  for(label i = 0; i < label(tvec.size()-1); i++)
  {
    if(tval >= tvec[i] && tval < tvec[i+1])
    {
      Info << "Interpolating laser coordinate" << endl;
      scalar dt = tvec[i+1] - tvec[i];
      scalar dy = yvec[i+1] - yvec[i];
      return (yvec[i] + (tval-tvec[i])*(dy/dt));
    }
  }

  // If none of these work out, then return 0
  Info << "Didn't find a valid interpolation" << endl;
  return 0;

}
