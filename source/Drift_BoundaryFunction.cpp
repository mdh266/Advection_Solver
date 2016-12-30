#include "../include/Drift_BoundaryFunction.hpp"

namespace Advection
{
using namespace dealii;

template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
                                     std::vector<double> &values,
                                     const unsigned int) const
{
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<values.size(); ++i)
    {
        if (points[i](0)<0.5 && points[i](1) < 0.5)
            values[i]=1.;
        else
            values[i]=0.;
    }
}
}
