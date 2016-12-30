#ifndef _DRIFT_RHS_H__
#define _DRIFT_RHS_H__

#include <deal.II/base/function.h>

namespace Advection
{
using namespace dealii;

template <int dim>
class RHS:  public Function<dim>
{
public:
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double> &values,
                             const unsigned int component=0) const;
};

}



#endif
