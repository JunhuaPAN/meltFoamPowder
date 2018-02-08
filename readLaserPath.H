
    // Read the array of x and y laser coordinates from file
    std::ifstream xyfile("xylaser.txt");

    // Store t, x, y in vectors
    std::vector<scalar> t, x, y;
    scalar tdata, xdata, ydata;
    // Read xy laser data
    while(xyfile >> tdata >> xdata >> ydata)
    {
      t.push_back(tdata);
      x.push_back(xdata);
      y.push_back(ydata);
    }

