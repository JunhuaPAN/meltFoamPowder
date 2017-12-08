
// Read alpha from output created by solidificationMeltingSource
const volScalarField& alpha = mesh.lookupObject<volScalarField>("alpha");

// Calculate dT/dt in case you need it
//volScalarField dTdt = fvc::ddt(thermo.T());

// Loop on cells
forAll(alpha, i)
{
  // If region is solid, then retain whatever the previous value of gg was
  if(alpha[i] < 0.001)
  {
    //Info << "Solid zone for cell " << i << endl;
    continue; // Do nothing
    
  }

  // If cell is mushy, then compute gg for the cell
  else if(alpha[i] >= 0.001 && alpha[i] <= 0.999)
  {
    scalar Tdot = gradT[i].x()*Vlaser[0]
                + gradT[i].y()*Vlaser[1];
    

    // Calculate radius of the grains and volume fraction of grains
    if(abs(Tdot)>1e-3)
    {
	   
      rad[i]=1e-4*((1.689e3-T[i])*(1.689e3-T[i])*(1.689e3-T[i]))/(3*abs(Tdot));


      gg[i] =(100*4*22*(rad[i])*(rad[i])*(rad[i]))/(3*7);
     

    }
    if(((alpha.oldTime())[i]<=0.001)||(((alpha.oldTime())[i])>=0.999))
    {
      if((Foam::sqrt((gradT[i].x())*(gradT[i].x())+(gradT[i].y())*(gradT[i].y())+(gradT[i].z())*(gradT[i].z()))*Foam::sqrt((Vlaser[0]*Vlaser[0])+(Vlaser[1]*Vlaser[1])))>1e-6)
      {
      //calculate the angle between the [100] (direcion of gradT) and laser direction
     
      theta[i]=Foam::acos(Tdot/(Foam::sqrt((gradT[i].x())*(gradT[i].x())+(gradT[i].y())*(gradT[i].y())+(gradT[i].z())*(gradT[i].z()))*Foam::sqrt((Vlaser[0]*Vlaser[0])+(Vlaser[1]*Vlaser[1]))));
      theta[i]=(theta[i]*180)/(2*Foam::acos(0.0));
    
      }
      else
      {
        theta[i]=0;
      }
    }
    
  }
  //If cell is liquid, set its value to given
  else if(alpha[i] > 0.999)
  {
    gg[i] = 0;
    theta[i]=180;
  }
	}
//const volScalarField& xCoord = mesh.C().component(vector::X);
//forAll(xCoord, ii)
//{
//if((xCoord[ii]<0.04)||(xCoord[ii]>0.16))
//{
//rad[ii]=0.0001;
//gg[ii]=1;

//}
//} 
