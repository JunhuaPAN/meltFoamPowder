
// Calculate max/min variable values over simulation duration
forAll(T, I)
{
  // Max heating rate
  if(max_dTdt[I] < dTdt[I])
    max_dTdt[I] = dTdt[I];
  // Max cooling rate (min heating rate)
  if(min_dTdt[I] > dTdt[I])
    min_dTdt[I] = dTdt[I];
  // Temperature gradient
  if(maxMagGradT[I] < mag(gradT[I]))
    maxMagGradT[I] = mag(gradT[I]);
  if(minMagGradT[I] > mag(gradT[I]))
    minMagGradT[I] = mag(gradT[I]);

  /* Calculate the values of the following durring solidification:
     - Temperature gradient magnitude Gn = mag(gradT)
     - Solidification rate Rdot = Vlaser & (gradT/Gn) ( = 0 if Gn = 0)
     - Sratio = Gn/Rdot
     Calculation performed only if dTdt is NEGATIVE and if alpha is between 0 and 1
  */
  // Gn, Rdot, Sratio = 0 if fluid
  if(alpha[I] >= 0.999)
  {
    Gn[I]     = 0;
    Rdot[I]   = 0;
    Sratio[I] = 0;
  }
  // Calculation is made only if region is fluid and is cooling down
  else if(dTdt[I] < 0 && alpha[I] > 0.001 && alpha[I] < 0.999)
  {
    // Temperature gradient magnitude
    Gn[I] = mag(gradT[I]);
    // Solidification rate
    if(mag(gradT[I]) == 0)
      Rdot[I] = 0;
    else
    {
      Rdot[I] = Vlaser[0]*gradT[I].x() + Vlaser[1]*gradT[I].y();
      Rdot[I] = std::abs(Rdot[I]/mag(gradT[I])); // Should be positive
    }
    // Sratio = magGradT/Rdot
    if(Rdot[I] == 0)
      Sratio[I] = 0;
    else
      Sratio[I] = Gn[I]/Rdot[I]; // Positive since Rdot > 0
  }
}