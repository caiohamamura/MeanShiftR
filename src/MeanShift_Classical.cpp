#include <Rcpp.h>
#include <cmath>
#include "LittleFunctionsCollection.h"
#include "Progress.h"
#include <unordered_map>
#include "hashtable.h"

using namespace Rcpp;
using namespace std;

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
  // NumericVector Z = pc( _, 3 );
  NumericVector rngX = Rcpp::range(X);
  NumericVector rngY = Rcpp::range(Y);
  const int mult = rngY[1] + 0.5;
  // std::cout << "Max y:" << rngY[1] << std::endl;
  for (int i = 0; i < nrows; i++) {
    uint64_t idx = ((uint64_t)X[i] * mult) + Y[i];
    mapIndex[idx].push_back(i);
  }

  int maxSize = mapIndex.size();
  hashtable *ht = ht_create(maxSize*1.1, nrows);


  for (auto& it: mapIndex) {
    ht_insert(ht, it.first, (int*)&(it.second[0]), it.second.size());
  }


  // uint64_t idx = ((uint64_t)X[1] * mult) + Y[1];
  // int count = -1;
  // int* vals = ht_get(ht, idx, &count);
  // Rcout << "idx: " << idx << endl;
  // Rcout << "count: " << count << endl;
  // for (int i = 0; i < count; i++) {
  //       Rcout << "id: " << vals[i] << endl;
  // }
  // Rcout << "Build OK" << endl;
  // return DataFrame::create();



  // for (auto& it: mapIndex) {
  //   int size = -1;
  //   int *ptr = ht_get(ht, it.first, &size);
  //   for (int i = 0; i < size; i++) {
  //     Rcout << ptr[i] << endl;
  //   }
  // }


  // Loop through all points to process one after the other
  for(int i=0; i<nrows; i++){
    pb.increment();

    // Initialize variables to store wthe mean coordinates of all neigbors with the
    // actual coordinates of the focal point from where the kernel starts moving
    double meanx = (double) X[i];
    double meany = (double) Y[i];
    double meanz = (double) pc(i,2);


    // Initialize variables to store the old coodinates with unrealistic values of -100
    double oldx = (double) meanx;
    double oldy = (double) meany;
    double oldz = (double) meanz;

    // Keep iterating as long as the centroid (or the maximum number of iterations) is not reached
    int IterCounter = 0;
    do {

      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      double sump = 0.0;

      // Increase the iteration counter
      IterCounter = IterCounter + 1;

      // Calculate cylinder dimensions based on point height
      const double r = H2CW_fac * oldz * 0.5;
      const double d = H2CW_fac * oldz;
      const double h = H2CL_fac * oldz;

      //
      int n_points = 0;
      // Rcout << "Starting getAround...\n";
      int* nn_idx = getAround(ht, oldx, oldy, r, mult, rngX[1], rngY[1], &n_points);
      // Rcout << "i: " << i << ", n_points: " << n_points << endl;


      // Loop through all points to identify the neighbors of the focal point
      for(int j=0; j < n_points; j++) {
        // for(int j=0; j<nrows; j++){

        const int idx = nn_idx[j];
        const double jx = (double) X[idx];
        const double jy = (double) Y[idx];
        const double jz = (double) pc(idx, 2);
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
            sumx = sumx + weight * jx;
            sumy = sumy + weight * jy;
            sumz = sumz + weight * jz;
            sump = sump + weight;
          }
          // If the option of a uniform kernel is set to true calculate the centroid
          // by summing up all coodinates and dividing by the number of points
          else
          {
            sumx = sumx + jx;
            sumy = sumy + jy;
            sumz = sumz + jz;
            sump = sump + 1.0;
          }
        }
      }
      meanx = sumx / sump;
      meany = sumy / sump;
      meanz = sumz / sump;

      if (meanx == oldx && meany == oldy && meanz == oldz)
        break;

      // Remember the coordinate means (kernel position) of previos iteration
      oldx = meanx;
      oldy = meany;
      oldz = meanz;

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

