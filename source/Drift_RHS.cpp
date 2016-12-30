#include "../include/Drift_RHS.hpp"

namespace Advection
{
using namespace dealii;

template <int dim>
void RHS<dim>::value_list(const std::vector<Point<dim> > &points,
                          std::vector<double> &values,
                          const unsigned int) const
{
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<values.size(); ++i)
        values[i]=0;
}

}
