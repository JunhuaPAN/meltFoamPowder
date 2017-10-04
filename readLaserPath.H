
    // Read the array of x and y laser coordinates from file
    std::ifstream xfile("xlaser.txt");
    std::ifstream yfile("ylaser.txt");

    // Store t, x, y in vectors
    std::vector<scalar> tx, ty, x, y;
    scalar tdata, xdata, ydata;
    // Read x laser data
    while(xfile >> tdata >> xdata)
    {
      tx.push_back(tdata);
      x.push_back(xdata);
    }
    // Read y laser data
    while(yfile >> tdata >> ydata)
    {
      ty.push_back(tdata);
      y.push_back(ydata);
    }

    /*
    // Read the laser velocity data
    std::ifstream vxfile("vxlaser.txt");
    std::ifstream vyfile("vylaser.txt");

    // Store t, vx, vy in vectors
    std::vector<scalar> tvx, tvy, vx, vy;
    scalar vxdata, vydata;
    // Read x laser data
    while(vxfile >> tdata >> vxdata)
    {
      tvx.push_back(tdata);
      vx.push_back(xdata);
    }
    // Read y laser data
    while(vyfile >> tdata >> vydata)
    {
      tvy.push_back(tdata);
      vy.push_back(ydata);
    }
    */
