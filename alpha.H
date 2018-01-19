// Read alpha from output created by solidificationMeltingSource

const volScalarField& alpha = mesh.lookupObject<volScalarField>("alpha");

// Calculate dT/dt in case you need it

//volScalarField dTdt = fvc::ddt(thermo.T());
// Loop on cells

//scalar VV=Foam::sqrt(Vlaser[0]*Vlaser[0]+Vlaser[1]*Vlaser[1]);


forAll(alpha, i)
{

  //scalar GG=Foam::sqrt((gradT[i].x())*(gradT[i].x())+(gradT[i].y())*(gradT[i].y())+(gradT[i].z())*(gradT[i].z()));

  scalar Tdot = gradT[i].x()*Vlaser[0]+ gradT[i].y()*Vlaser[1];

  // If region is solid, then retain whatever the previous value of gg was
  if(alpha[i] < 0.00001)
  {

    //Info << "Solid zone for cell " << i << endl;

    //continue; 

    if((((alpha.oldTime())[i])>=0.00001)&&(((alpha.oldTime())[i])<=0.99999))
    {
      //GradT_2[i]=GG;
      Tdot_2[i]=Tdot;				//cooling rate at end
      End_time[i]=runTime.value(); 		//end time for solidification

      //calculate radius of equiaxed grains depends on solidification time
      rad[i]=8.02*Foam::sqrt(2.58e-4*(End_time[i]-Start_time[i]));

      //Calculate the volume fraction of equiaxed grains
      scalar phi=abs((3.4e6)*4*22*rad[i]*rad[i]*rad[i]/(3*7));

      gg[i]=1-Foam::exp(-phi);

      tf[i]=End_time[i]-Start_time[i];   	//solidification time
      SDAS[i]=165.0*Foam::pow(tf[i],0.33333); 	//SDAS in micrometers depends on solidification time

      theta[i]=(theta[i]*180)/(2*Foam::acos(0.0));	//convert theta to degrees

    }

  }
  else if(alpha[i] >= 0.00001 && alpha[i] <= 0.99999)
  {

    if((((alpha.oldTime())[i])>0.99999))    	//||(((alpha.oldTime())[i])<0.00001))
    {
      Tdot_1[i]=Tdot;
      Start_time[i]=runTime.value(); //start time for solidification
    }

    if((Foam::sqrt((gradT[i].x())*(gradT[i].x())+(gradT[i].y())*(gradT[i].y())+(gradT[i].z())*(gradT[i].z()))*Foam::sqrt((Vlaser[0]*Vlaser[0])+(Vlaser[1]*Vlaser[1])))>1e-6)
    {

     //calculate the angle between the [001] (direcion of gradT) and laser direction

     theta[i]=Foam::acos(Tdot/(Foam::sqrt((gradT[i].x())*(gradT[i].x())
		+(gradT[i].y())*(gradT[i].y())+(gradT[i].z())*(gradT[i].z()))*Foam::sqrt((Vlaser[0]*Vlaser[0])
		+(Vlaser[1]*Vlaser[1]))));

    }
    else
    {
       theta[i]=0;
    }
  }
  //If cell is liquid, set its value to given
  else if(alpha[i] > 0.999)
  {
    gg[i] = 0;
    theta[i]=180;
  }

}

