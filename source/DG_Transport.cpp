#include "../include/DG_Transport.hpp"



namespace Advection
{
using namespace dealii;

template <int dim>
DGTransportEquation<dim>::DGTransportEquation ()
    :
    beta_function (),
    rhs_function (),
    boundary_function ()
{}


template <int dim>
void DGTransportEquation<dim>::assemble_cell_term(
    const FEValues<dim> &fe_v,
    FullMatrix<double> &ui_vi_matrix,
    Vector<double> &cell_vector) const
{
    const std::vector<double> &JxW = fe_v.get_JxW_values ();

    std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
    std::vector<double> rhs (fe_v.n_quadrature_points);

    beta_function.value_list (fe_v.get_quadrature_points(), beta);
    rhs_function.value_list (fe_v.get_quadrature_points(), rhs);

    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
        {
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
            {
                ui_vi_matrix(i,j) -= beta[point]*fe_v.shape_grad(i,point)*
                                     fe_v.shape_value(j,point) *
                                     JxW[point];
            }
            cell_vector(i) += rhs[point] * fe_v.shape_value(i,point) * JxW[point];
        }
}

template <int dim>
void DGTransportEquation<dim>::assemble_boundary_term(
    const FEFaceValues<dim> &fe_v,
    FullMatrix<double> &ui_vi_matrix,
    Vector<double> &cell_vector) const
{
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

    std::vector<Point<dim> > beta (fe_v.n_quadrature_points);
    std::vector<double> g(fe_v.n_quadrature_points);

    beta_function.value_list (fe_v.get_quadrature_points(), beta);
    boundary_function.value_list (fe_v.get_quadrature_points(), g);

    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
        const double beta_n=beta[point] * normals[point];

        // outflow boundary
        if (beta_n>0)
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
                for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
                    ui_vi_matrix(i,j) += beta_n *
                                         fe_v.shape_value(j,point) *
                                         fe_v.shape_value(i,point) *
                                         JxW[point];
        // inflow boundary
        else
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
                cell_vector(i) -= beta_n *
                                  g[point] *
                                  fe_v.shape_value(i,point) *
                                  JxW[point];
    }
}


template <int dim>
void DGTransportEquation<dim>::assemble_interior_face_term(
    const FEFaceValuesBase<dim> &fe_v,
    const FEFaceValuesBase<dim> &fe_v_neighbor,
    FullMatrix<double> &ui_vi_matrix,
    FullMatrix<double> &ue_vi_matrix,
    FullMatrix<double> &ui_ve_matrix,
    FullMatrix<double> &ue_ve_matrix) const
{
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

    std::vector<Point<dim> > beta (fe_v.n_quadrature_points);

    beta_function.value_list (fe_v.get_quadrature_points(), beta);

    for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
        const double beta_n=beta[point] * normals[point];

        // ouflow across interior cells
        if (beta_n>0)
        {
            //		std::cout << "outflow" << std::endl;
            // int_{beta * n^{-} > 0 } beta * v^{-} u^{-} dx
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
                for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
                    ui_vi_matrix(i,j) += beta_n *
                                         fe_v.shape_value(i,point) * //test
                                         fe_v.shape_value(j,point) * //trial
                                         JxW[point];

            // int_{beta * n^{-} > 0 } beta * v^{+} u^{-} dx
            for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
                for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
                    ui_ve_matrix(k,j) -= beta_n *
                                         fe_v_neighbor.shape_value(k,point) * // test
                                         fe_v.shape_value(j,point) * // trial
                                         JxW[point];
        }
        else
        {
            //		std::cout << "inflow" << std::endl;
            // int_{beta * n^{-} < 0 } beta * v^{-} u^{+} dx
            for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
                for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
                    ue_vi_matrix(i,l) += beta_n *
                                         fe_v.shape_value(i,point) * // test
                                         fe_v_neighbor.shape_value(l,point) * //trial
                                         JxW[point];

            // int_{beta * n^{-} < 0 } beta * v^{+} u^{+} dx
            for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
                for (unsigned int l=0; l<fe_v_neighbor.dofs_per_cell; ++l)
                    ue_ve_matrix(k,l) -= beta_n *
                                         fe_v_neighbor.shape_value(k,point) * //test
                                         fe_v_neighbor.shape_value(l,point) * //trial
                                         JxW[point];
        }
    } // for point
} // assemble_interior_face_term

} // namespace Advection



