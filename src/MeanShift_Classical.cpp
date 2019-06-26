#include <Rcpp.h>
#include <cmath>
#include "LittleFunctionsCollection.h"
#include <Progress.h>
#include <unordered_map>
#include <hashtable.h>
#include <CL/cl.hpp>
#include <CL/kernelMeanShift.h>

using namespace Rcpp;

void configureOpenCL() {



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
  // NumericVector Z = pc( _, 3 );
  int64_t maxX = ((NumericVector)Rcpp::range(X))[0];
  int64_t maxY = ((NumericVector)Rcpp::range(Y))[0];
  const int mult = maxY + 0.5;
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


  ////////////////////////////
  // OPENCL
  ////////////////////////////
  VECTOR_CLASS<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  VECTOR_CLASS<cl::Device> devices;
  platforms[0].getDevices(CL_DEVICE_TYPE_ALL, &devices);
  cl::Context context(devices);

  cl_int err = -1;


  cl::Program vectorWrapper(
      context,
      KERNELSTRING,
      false
  );
  err = vectorWrapper.build("-I");


  if (err == CL_BUILD_PROGRAM_FAILURE) {
    std::string buildlog = vectorWrapper.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    std::cerr << "Build log:"
              << buildlog << std::endl;
    return DataFrame::create();
  }

  cl::Buffer buf_ht(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (size_t)(sizeof(keyvalue) * ht->size), ht->table);
  cl::Buffer buf_ids(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int) * ht->size, ht->ids);
  cl::Buffer buf_X(context, X.begin(), X.end(), true);
  cl::Buffer buf_Y(context, Y.begin(), Y.end(), true);
  cl::Buffer buf_Z(context, Z.begin(), Z.end(), true);
  cl::Buffer buf_H2CW_fac(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (size_t)(sizeof(cl_double)), (void*)&H2CW_fac);
  cl::Buffer buf_H2CL_fac(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_double), (void*)&H2CL_fac);
  cl::Buffer buf_MaxIter(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int), (void*)&MaxIter);
  cl::Buffer buf_mult(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int), (void*)&mult);
  cl::Buffer buf_maxX(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_long), (void*)&maxX);
  cl::Buffer buf_maxY(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_long), (void*)&maxY);
  cl::Buffer buf_UniformKernel(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_bool), (void*)&UniformKernel);
  cl::Buffer buf_centroidX(context, CL_MEM_WRITE_ONLY, sizeof(cl_double)*nrows, NULL);
  cl::Buffer buf_centroidY(context, CL_MEM_WRITE_ONLY, sizeof(cl_double)*nrows, NULL);
  cl::Buffer buf_centroidZ(context, CL_MEM_WRITE_ONLY, sizeof(cl_double)*nrows, NULL);




  auto vectorAddKernel =
    cl::make_kernel<
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&,
      cl::Buffer&
    >( vectorWrapper, "cl_mean_shift" );

    vectorAddKernel(cl::EnqueueArgs(
        nrows),
        buf_ht,
        buf_ids,
        buf_X,
        buf_Y,
        buf_Z,
        buf_H2CW_fac,
        buf_H2CL_fac,
        buf_MaxIter,
        buf_mult,
        buf_maxX,
        buf_maxY,
        buf_UniformKernel,
        buf_centroidX,
        buf_centroidY,
        buf_centroidZ);

    cl::copy(buf_centroidX, std::begin(centroidx), std::end(centroidx));
    cl::copy(buf_centroidY, std::begin(centroidy), std::end(centroidy));
    cl::copy(buf_centroidZ, std::begin(centroidz), std::end(centroidz));

    clReleaseMemObject(buf_ht());
    clReleaseMemObject(buf_ids());
    clReleaseMemObject(buf_X());
    clReleaseMemObject(buf_Y());
    clReleaseMemObject(buf_Z());
    clReleaseMemObject(buf_H2CW_fac());
    clReleaseMemObject(buf_H2CL_fac());
    clReleaseMemObject(buf_MaxIter());
    clReleaseMemObject(buf_mult());
    clReleaseMemObject(buf_maxX());
    clReleaseMemObject(buf_maxY());
    clReleaseMemObject(buf_UniformKernel());
    clReleaseMemObject(buf_centroidX());
    clReleaseMemObject(buf_centroidY());
    clReleaseMemObject(buf_centroidZ());
  // Return the result as a data.frame with XYZ-coordinates of all points and their corresponding centroids
  return DataFrame::create(_["X"]= pc(_,0),_["Y"]= pc(_,1),_["Z"]= pc(_,2),_["CtrX"]= centroidx,_["CtrY"]= centroidy,_["CtrZ"]= centroidz);
}

