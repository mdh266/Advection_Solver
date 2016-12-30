#ifndef _DG_TRANSPORT_H__
#define _DG_TRANSPORT_H__

#include <deal.II/base/function.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/vector.h>
#include "../source/Drift_RHS.cpp"
#include "../source/Drift_BoundaryFunction.cpp"
#include "../source/FlowField.cpp"

namespace Advection
{

using namespace dealii;

template <int dim>
class DGTransportEquation
{
public:
    DGTransportEquation();

    void assemble_cell_term(const FEValues<dim> &fe_v,
                            FullMatrix<double> &ui_vi_matrix,
                            Vector<double> &cell_vector) const;

    void assemble_boundary_term(const FEFaceValues<dim> &fe_v,
                                FullMatrix<double> &ui_vi_matrix,
                                Vector<double> &cell_vector) const;

    void assemble_interior_face_term(const FEFaceValuesBase<dim> &fe_v,
                                     const FEFaceValuesBase<dim> &fe_v_neighbor,
                                     FullMatrix<double> &ui_vi_matrix,
                                     FullMatrix<double> &ue_vi_matrix,
                                     FullMatrix<double> &ui_ve_matrix,
                                     FullMatrix<double> &ue_ve_matrix) const;
private:
    const Beta<dim> beta_function;
    const RHS<dim> rhs_function;
    const BoundaryValues<dim> boundary_function;
};

}

#endif
