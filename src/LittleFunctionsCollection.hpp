#ifndef MEANSHIFTR_LITTLEFUNCTIONSCOLLECTION_HPP
#define MEANSHIFTR_LITTLEFUNCTIONSCOLLECTION_HPP
#include <Rcpp.h>
#include <cmath>

// Declarations of all the little functions used by the main functions

bool InCylinder(double PointX, double PointY, double PointZ, double Radius, double Height, double CtrX, double CtrY, double CtrZ);
double VerticalDistance(double Height, double CtrZ, double PointZ);
double VerticalMask(double Height, double CtrZ, double PointZ);
double EpanechnikovFunction(double Height, double CtrZ, double PointZ);
double GaussFunction(double Width, double CtrX, double CtrY, double PointX, double PointY);
double IntensityGaussFunction(double Width, double CtrI, double PointI);

namespace tools
{
    struct point{
        public:
            double x, y, z;
            int id;

            double distance(const point *pt) const {
                return (
                    std::sqrt(
                        std::pow(x - pt->x, 2) +
                        std::pow(y - pt->y, 2) +
                        std::pow(z - pt->z, 2)
                    )
                );
            }
    };
}

#endif /* MEANSHIFTR_LITTLEFUNCTIONSCOLLECTION_HPP */
