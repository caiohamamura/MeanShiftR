#include <Rcpp.h>
#include <cmath>
#include "LittleFunctionsCollection.h"
#include "Progress.h"
#include <unordered_map>
using namespace Rcpp;
using namespace std;


vector<int> getAround(unordered_map<uint64_t, vector<int>> *hashMap, const double x, const double y, const double dist, const uint8_t bitShift) {
  vector<int> myVec;
  uint64_t xMin = (uint64_t)(x - dist);
  uint64_t xMax = (uint64_t)(x + dist + 0.5);
  uint64_t yMin = (uint64_t)(y - dist);
  uint64_t yMax = (uint64_t)(y + dist + 0.5);
  for (uint64_t curY = yMin; curY <= yMax; curY++) {
    for (uint64_t curX = xMin; curX <= xMax; curX++) {
      uint64_t idx = (curX << bitShift) + curY;
      for (int it : (*hashMap)[idx]) {
        myVec.push_back(it);
      }
    }
  }
  return myVec;
}

//' Mean shift clustering
//'
//' @title Mean shift clustering
//' @description
//' Adaptive mean shift clustering to delineate tree crowns from lidar point clouds
//' @param pc Point cloud has to be in matrix format with 3-columns representing X, Y and Z and each row representing one point
//' @param H2CW_fac Factor for the ratio of height to crown width. Determines kernel diameter based on its height above ground.
//' @param H2CL_fac Factor for the ratio of height to crown length. Determines kernel height based on its height above ground.
//' @param UniformKernel Boolean to enable the application of a simple uniform kernel without distance weighting (Default False)
//' @param MaxIter Maximum number of iterations, i.e. steps that the kernel can move for each point. If centroid is not found after all iteration, the last position is assigned as centroid and the processing jumps to the next point
//' @return data.frame with X, Y and Z coordinates of each point in the point cloud and  X, Y and Z coordinates of the centroid to which the point belongs
// [[Rcpp::export]]
DataFrame C_MeanShift_Classical(NumericMatrix pc, const double H2CW_fac, double H2CL_fac, const bool UniformKernel=false, const int MaxIter=20){

  // Create three vectors that all have the length of the incoming point cloud.
  // In these vectors the coodinates of the centroids will be stored
  const int nrows = pc.nrow();
  NumericVector centroidx(nrows);
  NumericVector centroidy(nrows);
  NumericVector centroidz(nrows);

  Progress pb(nrows, "Meanshifting: ");

  //////////////////////////
  // Create an index tree
  /////////////////////////
  std::unordered_map<uint64_t, std::vector<int>> mapIndex;

  // Remove min
  NumericVector X = pc( _, 0 );
  NumericVector Y = pc( _, 1 );
  NumericVector Z = pc( _, 2 );
  NumericVector C1 = pc( _, 3 );
  NumericVector C2 = pc( _, 4 );
  NumericVector C3 = pc( _, 5 );
  // NumericVector Z = pc( _, 3 );
  NumericVector rngX = Rcpp::range(X);
  NumericVector rngY = Rcpp::range(Y);
  const uint8_t bitShift = log2(rngY[1]) + 1;
  // std::cout << "Max y:" << rngY[1] << std::endl;
  for (int i = 0; i < nrows; i++) {
    X[i] -= rngX[0];
    Y[i] -= rngY[0];
    uint64_t idx = ((uint64_t)X[i] << bitShift) + Y[i];
    mapIndex[idx].push_back(i);
  }


  // Loop through all points to process one after the other
  for(int i=0; i<nrows; i++){
    pb.increment();

    // Initialize variables to store the mean coordinates of all neigbors with the
    // actual coordinates of the focal point from where the kernel starts moving
    double meanx = (double) X[i];
    double meany = (double) Y[i];
    double meanz = (double) pc(i,2);
    double meanc1 = (double) C1[i];
    double meanc2 = (double) C2[i];
    double meanc3 = (double) C3[i];


    // Initialize variables to store the old coodinates with unrealistic values of -100
    double oldx = (double) meanx;
    double oldy = (double) meany;
    double oldz = (double) meanz;
    double oldc1 = (double) meanc1;
    double oldc2 = (double) meanc2;
    double oldc3 = (double) meanc3;



    // Keep iterating as long as the centroid (or the maximum number of iterations) is not reached
    int IterCounter = 0;
    do {

      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      double sumc1 = 0.0;
      double sumc2 = 0.0;
      double sumc3 = 0.0;
      double sump = 0.0;

      // Increase the iteration counter
      IterCounter++;

      // Calculate cylinder dimensions based on point height
      const double r = H2CW_fac * oldz * 0.5;
      const double d = H2CW_fac * oldz;
      const double h = H2CL_fac * oldz;

      //
      vector<int> nn_idx = getAround(&mapIndex, oldx, oldy, r, bitShift);
      const int n_points = nn_idx.size();

      // Loop through all points to identify the neighbors of the focal point
      for(int j=0; j < n_points; j++) {
        // for(int j=0; j<nrows; j++){

        const int idx = nn_idx[j];
        const double jx = (double) X[idx];
        const double jy = (double) Y[idx];
        const double jz = (double) pc(idx, 2);
        const double jc1 = (double) C1[idx];
        const double jc2 = (double) C2[idx];
        const double jc3 = (double) C3[idx];

        // double jx = (double) pc(j, 0);
        // double jy = (double) pc(j, 1);
        // double jz = (double) pc(j, 2);

        if(InCylinder(jx, jy, jz, r, h, oldx, oldy, oldz)){
          // if(abs(jz - meanz) <= (h/2.0)){

          // If the option of a uniform kernel is set to false calculate the centroid
          // by multiplying all coodinates by their weights, depending on their relative
          // position within the cylinder, summing up the products and dividing by the
          // sum of all weights
          if(UniformKernel == false){
            const double verticalweight = EpanechnikovFunction(h, oldz, jz);
            const double horizontalweight = GaussFunction(d, oldx, oldy, jx, jy);
            const double weight = verticalweight * horizontalweight;
            sumx += weight * jx;
            sumy += weight * jy;
            sumz += weight * jz;
            sumc1 += weight * jc1;
            sumc2 += weight * jc2;
            sumc3 += weight * jc3;
            sump += weight;
          }
          // If the option of a uniform kernel is set to true calculate the centroid
          // by summing up all coodinates and dividing by the number of points
          else
          {
            sumx = sumx + jx;
            sumy = sumy + jy;
            sumz = sumz + jz;
            sumc1 = sumc1 + jc1;
            sumc2 = sumc2 + jc2;
            sumc3 = sumc3 + jc3;
            sump = sump + 1.0;
          }
        }
      }
      meanx = sumx / sump;
      meany = sumy / sump;
      meanz = sumz / sump;
      meanc1 = sumc1 / sump;
      meanc2 = sumc2 / sump;
      meanc3 = sumc3 / sump;

      if (meanx == oldx &&
          meany == oldy &&
          meanz == oldz &&
          meanc1 == oldc1 &&
          meanc2 == oldc2 &&
          meanc3 == oldc3
          )
        break;

      // Remember the coordinate means (kernel position) of previous iteration
      oldx = meanx;
      oldy = meany;
      oldz = meanz;
      oldc1 = meanc1;
      oldc2 = meanc2;
      oldc3 = meanc3;

      // If the new position equals the previous position (kernel stopped moving), or if the
      // maximum number of iterations is reached, stop the iterations
    }  while(IterCounter < MaxIter);

    // Store the found position as the centroid position for the focal point
    centroidx[i] = meanx + rngX[0];
    centroidy[i] = meany + rngY[0];
    centroidz[i] = meanz;

    if ((i % 1000) == 999) {
      try
      {
        Rcpp::checkUserInterrupt();
      }
      catch(Rcpp::internal::InterruptedException e)
      {
        return DataFrame::create();
      }
    }
  }



  // Return the result as a data.frame with XYZ-coordinates of all points and their corresponding centroids
  return DataFrame::create(_["X"]= pc(_,0),_["Y"]= pc(_,1),_["Z"]= pc(_,2),_["CtrX"]= centroidx,_["CtrY"]= centroidy,_["CtrZ"]= centroidz);
}

