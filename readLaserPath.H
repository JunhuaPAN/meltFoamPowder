
    // Read the array of x and y laser coordinates from file
    std::ifstream xyfile("xylaser.txt");

    // Store t, x, y in vectors
    std::vector<scalar> t, x, y;
    
    // If file is empty, 't' vector will contain one entry (start time) and x/y
    // vectors will contain arbitrary coordinate
    bool is_empty = (xyfile.peek()==std::ifstream::traits_type::eof());
    if(is_empty)
    //if(xyfile.peek() == std::ifstream::traits_type::eof())
    {
      t.push_back(runTime.value());
      x.push_back(0.0);
      y.push_back(0.0);
    }
    else
    {
      // Read xy laser data
      scalar tdata, xdata, ydata;
      while(xyfile >> tdata >> xdata >> ydata)
      {
        t.push_back(tdata);
        x.push_back(xdata);
        y.push_back(ydata);
      }
    }

