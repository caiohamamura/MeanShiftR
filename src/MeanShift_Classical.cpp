#include <cmath>
#include "LittleFunctionsCollection.h"
#include <ANN/ANN.h>
#include "Progress.h"
#include <Rcpp.h>
using namespace Rcpp;


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
  const int n_dimensions = 2;  // X, Y
  ANNkd_tree	*the_tree;	// Search structure
  ANNpointArray data_pts 	= annAllocPts(nrows, n_dimensions);
  ANNidxArray nn_idx 		= new ANNidx[nrows];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[nrows];		// Allocate near neighbor dists
  int n_points;

  // now construct the points
  for(int i = 0; i < nrows; i++)
  {
    data_pts[i][0]=pc(i, 0);
    data_pts[i][1]=pc(i, 1);
  }
  the_tree = new ANNkd_tree( data_pts, nrows, n_dimensions);

  // Loop through all points to process one after the other
  for(int i=0; i<nrows; i++){
    pb.increment();

    // Initialize variables to store the mean coordinates of all neigbors with the
    // actual coordinates of the focal point from where the kernel starts moving
    double meanx = (double) pc(i,0);
    double meany = (double) pc(i,1);
    double meanz = (double) pc(i,2);


    // Initialize variables to store the old coodinates with unrealistic values of -100
    double oldx = (double) pc(i,0);
    double oldy = (double) pc(i,1);
    double oldz = (double) pc(i,2);

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

      // Query which points are within radius
      ANNpoint pt_query = annAllocPt(n_dimensions);
      pt_query[0] = oldx;
      pt_query[1] = oldy;
      //
      n_points = the_tree->annkFRSearch(pt_query, pow(r, 2), nrows, nn_idx, dists);

      // Loop through all points to identify the neighbors of the focal point
      for(int j=0; j < n_points; j++) {
        // for(int j=0; j<nrows; j++){

        const int idx = nn_idx[j];
        const double jx = (double) pc(idx, 0);
        const double jy = (double) pc(idx, 1);
        const double jz = (double) pc(idx, 2);
        // double jx = (double) pc(j, 0);
        // double jy = (double) pc(j, 1);
        // double jz = (double) pc(j, 2);

        // if(InCylinder(jx, jy, jz, r, h, oldx, oldy, oldz)){
        if(abs(jz - meanz) <= (h/2.0)){

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
    centroidx[i] = meanx;
    centroidy[i] = meany;
    centroidz[i] = meanz;

    if ((i % 1000) == 999) {
      try
      {
        Rcpp::checkUserInterrupt();
      }
      catch(Rcpp::internal::InterruptedException e)
      {
        // Dealloc stuff
        annDeallocPts(data_pts);
        delete the_tree;
        delete [] nn_idx;
        delete [] dists;
        return DataFrame::create();
      }
    }
  }

  // Dealloc stuff
  annDeallocPts(data_pts);
  delete the_tree;
  delete [] nn_idx;
  delete [] dists;

  // Return the result as a data.frame with XYZ-coordinates of all points and their corresponding centroids
  return DataFrame::create(_["X"]= pc(_,0),_["Y"]= pc(_,1),_["Z"]= pc(_,2),_["CtrX"]= centroidx,_["CtrY"]= centroidy,_["CtrZ"]= centroidz);
}

