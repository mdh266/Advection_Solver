#ifndef _DRIFT_BOUNDARYFUNCTION_H__
#define _DRIFT_BOUNDARYFUNCTION_H__

#include <deal.II/base/function.h>

namespace Advection
{
using namespace dealii;

template <int dim>
class BoundaryValues:  public Function<dim>
{
public:
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const;
};
}

#endif
