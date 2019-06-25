#include <hashtable.h>
#include <CL/LittleFunctionsCollection.h>


kernel void vectorAdd(global const keyvalue *table, global const int *ids, global const double *X, global const double *Y, global const double *Z, const global double *H2CW_fac, const global double *H2CL_fac, const global int *mult, const global double *maxX, const global double *maxY, const global int *size)
{
  hashtable *ht = (hashtable*)malloc(sizeof(hashtable));
  ht->table = table;
  ht->ids = ids;
  ht->size = *size;


  const int i = get_global_id(0);

  // Initialize variables to store wthe mean coordinates of all neigbors with the
  // actual coordinates of the focal point from where the kernel starts moving
  double meanx = (double) X[i];
  double meany = (double) Y[i];
  double meanz = (double) Z(i);


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
    const double r = *H2CW_fac * oldz * 0.5;
    const double d = *H2CW_fac * oldz;
    const double h = *H2CL_fac * oldz;

    //
    int n_points = 0;
    // Rcout << "Starting getAround...\n";
    int* nn_idx = getAround(ht, oldx, oldy, r, *mult, *maxX, *maxY, &n_points);
    // Rcout << "i: " << i << ", n_points: " << n_points << endl;


    // Loop through all points to identify the neighbors of the focal point
    for(int j=0; j < n_points; j++) {

        const int idx = nn_idx[j];
        const double jx = (double) X[idx];
        const double jy = (double) Y[idx];
        const double jz = (double) Z[idx];

        if(InCylinder(jx, jy, jz, r, h, oldx, oldy, oldz)) {
            if(UniformKernel == false) {
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
    }
}
