#ifndef _FLOW_FIELD_H__
#define _FLOW_FIELD_H__

#include <deal.II/base/function.h>

namespace Advection
{
using namespace dealii;

template <int dim>
class Beta
{
public:
    Beta () {}
    void value_list (const std::vector<Point<dim> > &points,
                     std::vector<Point<dim> > &values) const;
};

}

#endif
