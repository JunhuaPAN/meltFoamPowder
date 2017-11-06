
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

