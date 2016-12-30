#include "../include/FlowField.hpp"


namespace Advection
{
using namespace dealii;

template <int dim>
void Beta<dim>::value_list(const std::vector<Point<dim> > &points,
                           std::vector<Point<dim> > &values) const
{
    Assert(values.size()==points.size(),
           ExcDimensionMismatch(values.size(),points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
    {
        //   if (points[i](0) > 0)
        //  {
        values[i](0) = +1.0;//-points[i](1);
        values[i](1) = +1.0;//points[i](0);
        /*        }
              else
                {
                  values[i] = Point<dim>();
                  values[i](0) = -points[i](1);
                }*/
    }
}

}
