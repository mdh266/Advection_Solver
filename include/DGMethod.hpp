#ifndef _DG_METHOD_H__
#define _DG_METHOD_H__


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>
#include <iostream>
#include <fstream>
#include "../source/DG_Transport.cpp"

namespace Advection
{
using namespace dealii;

template <int dim>
class DGMethod
{
public:
    DGMethod (const unsigned int degree);
    ~DGMethod ();

    void run ();

private:
    void setup_system ();
    void assemble_system ();
    void solve (Vector<double> &solution);
    void output_results () const;
    void output_system();

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    const unsigned int   degree;
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    const QGauss<dim>   quadrature;
    const QGauss<dim-1> face_quadrature;

    Vector<double>       solution2;
    Vector<double>       right_hand_side;

    const DGTransportEquation<dim> dg;
};

}


#endif
