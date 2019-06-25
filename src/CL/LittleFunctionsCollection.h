




// Collection of all the little functions used by the main functions

// Function to check whether a point [PointX, PointY, PointZ] is within a cylider of a given radius
// and height from the center point of the top circle [TopX, TopY, TopZ]
bool InCylinder(const double PointX, const double PointY, const double PointZ, const double Radius, const double Height, const double CtrX, const double CtrY, const double CtrZ){
  if ((pow((PointX - CtrX), 2.0) + pow((PointY - CtrY), 2.0) <= pow(Radius, 2.0)) && (PointZ >= (CtrZ - (0.5*Height))) && (PointZ <= (CtrZ + (0.5*Height))) == true) {
    return true;
  }	else {
    return false;
  }
}

// Help functions for vertical filter
double VerticalDistance(const double Height, const double CtrZ, const double PointZ){
  const double BottomDistance = (double) fabs((CtrZ-Height/4-PointZ)/(3*Height/8));
  const double TopDistance = (double) fabs((CtrZ+Height/2-PointZ)/(3*Height/8));
  double MinDistance = fmin(BottomDistance, TopDistance);
  return MinDistance;
}

//Equivalent R code
//distx <- function(h, CtrZ, PointZ){
  //  bottomdist <- abs((CtrZ-h/4-PointZ)/(3*h/8))
  //  topdist <- abs((CtrZ+h/2-PointZ)/(3*h/8))
  //  mindist <- pmin(bottomdist, topdist)
  //  return(mindist)
  //}

double VerticalMask(const double Height, const double CtrZ, const double PointZ){
  if((PointZ >= CtrZ-Height/4) && (PointZ <= CtrZ+Height/2)){
    return 1;
  }
  else
  {
    return 0;
  }
}

//Equivalent R code
//maskx <- function(h, CtrZ, PointZ){
  //  maskvec <- ifelse(PointZ >= CtrZ-h/4 & PointZ <= CtrZ+h/2, 1, 0)
  //  return(maskvec)
  //}

// Epanechnikov function for vertical filter
double EpanechnikovFunction(const double Height, const double CtrZ, const double PointZ){
  double Result = VerticalMask(Height, CtrZ, PointZ)*(1-(pow((1-VerticalDistance(Height, CtrZ, PointZ)), 2.0)));
  return Result;
}

//Equivalent R code
//Epanechnikov <- function(h, CtrZ, PointZ){
  //  output <- maskx(h, CtrZ, PointZ)*(1-(1-distx(h, CtrZ, PointZ))^2)
  //  return(output)
  //}

// Gauss function for horizontal filter
double GaussFunction(const double Width, const double CtrX, const double CtrY, const double PointX, const double PointY){
  const double Distance = pow((pow((PointX-CtrX), 2.0)+pow((PointY-CtrY), 2.0)), 0.5);
  const double NormDistance = Distance/Width;
  double Result = exp(-5.0*pow(NormDistance, 2.0));
  return Result;
}

//Equivalent R code
//gauss <- function(w, CtrX, CtrY, PointX, PointY){
  //  distance <- ((PointX-CtrX)^2+(PointY-CtrY)^2)^0.5
  //  norm.distance <- distance/w
  //  output <- exp(-5*norm.distance^2)
  //  return(output)
  //}
