const std::string KERNELSTRING = R"====(
typedef	unsigned int uint32_t;
typedef	unsigned long long uint64_t;
typedef	long long int64_t;
// struct ptr_size {
//   int *ptr;
//   int count;
// };
//
// typedef struct ptr_size ptr_size_t;

struct KeyValue {
  uint64_t key;
  uint64_t pos;
  int count;
};

typedef struct KeyValue keyvalue;

struct HashTable {
  int size;
  uint64_t last_pos;
  __global int *ids;
  __global keyvalue *table;
};

typedef struct HashTable hashtable;


uint64_t hashCode(uint64_t x, const int size) {
  x = (x ^ (x >> 30)) * (uint64_t)(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * (uint64_t)(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return (x % size);
}


__global int * *ht_get(hashtable* ht, const uint64_t key, int *count) {
  //get the hash
  uint64_t hashIndex = hashCode(key, ht->size);

  //move in array until an empty
  while(ht->table[hashIndex].count != -1) {

    if(ht->table[hashIndex].key == key) {
      *count = ht->table[hashIndex].count;
      return &(ht->ids[ht->table[hashIndex].pos]);
    }

    //go to next cell
    ++hashIndex;

    //wrap around the table
    hashIndex %= ht->size;
  }

  return NULL;
}


int *getAround(hashtable *ht, const double x, const double y, const double dist, const int mult, const int64_t maxX, const int64_t maxY, int *oSize) {
  int n = 0;
  int k = 0;
  int64_t xMin = (x - dist);
  int64_t xMax = (x + dist);
  int64_t yMin = (y - dist);
  int64_t yMax = (y + dist);
  if (xMin < 0) xMin = 0;
  if (yMin < 0) yMin = 0;
  if (xMax > maxX) xMax = maxX;
  if (yMax > maxY) yMax = maxY;
  int out[500];
  for (int64_t curY = yMin; curY <= yMax; curY++) {
    for (int64_t curX = xMin; curX <= xMax; curX++) {
      uint64_t idx = (curX * mult) + curY;
      // Rcpp::Rcout << "\ridx: " << idx << "    ";
      int size = -1;
      // Rcpp::Rcout << "Getting from ht...\n";
      int* ids = ht_get(ht, idx, &size);
      n += size;
      // Rcpp::Rcout << "Reallocating to size:" << size << "...\n";
      // Rcpp::Rcout << "Putting values...\n";
      for (int it = 0; it < size;) {
        out[k++] = ids[it++];
      }
    }
  }
  *oSize = n;
  return out;
}

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


kernel void cl_mean_shift(global const keyvalue *table, global const int *ids, global const double *X, global const double *Y, global const double *Z, const global double *H2CW_fac, const global double *H2CL_fac, const global int *MaxIter, const global int *mult, const global double *maxX, const global double *maxY, const global int *size, const global bool *UniformKernel, global double* centroidx, global double* centroidy, global double* centroidz)
{
  hashtable ht[1];
  ht->table = table;
  ht->ids = ids;
  ht->size = *size;

  const int i = get_global_id(0);

  // Initialize variables to store wthe mean coordinates of all neigbors with the
  // actual coordinates of the focal point from where the kernel starts moving
  double meanx = (double) X[i];
  double meany = (double) Y[i];
  double meanz = (double) Z[i];


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
      }  while(IterCounter < *MaxIter);

      // Store the found position as the centroid position for the focal point
      centroidx[i] = meanx;
      centroidy[i] = meany;
      centroidz[i] = meanz;
}
)====";
