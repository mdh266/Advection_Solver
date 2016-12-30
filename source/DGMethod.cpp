#include "../include/DGMethod.hpp"

namespace Advection
{
using namespace dealii;


template <int dim>
DGMethod<dim>::DGMethod (const unsigned int degree)
    :
    mapping (),
    degree(degree),
    fe (degree),
    dof_handler (triangulation),
    quadrature (degree+1),
    face_quadrature (degree+1),
    dg ()
{}


template <int dim>
DGMethod<dim>::~DGMethod ()
{
    dof_handler.clear ();
}


template <int dim>
void DGMethod<dim>::setup_system ()
{
    dof_handler.distribute_dofs (fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());

    DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit (sparsity_pattern);

    solution2.reinit (dof_handler.n_dofs());
    right_hand_side.reinit (dof_handler.n_dofs());
}



template <int dim>
void DGMethod<dim>::assemble_system()
{
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dofs (dofs_per_cell);
    std::vector<types::global_dof_index> dofs_neighbor (dofs_per_cell);

    const UpdateFlags update_flags = update_values
                                     | update_gradients
                                     | update_quadrature_points
                                     | update_JxW_values;

    const UpdateFlags face_update_flags = update_values
                                          | update_quadrature_points
                                          | update_JxW_values
                                          | update_normal_vectors;

    const UpdateFlags neighbor_face_update_flags = update_values;

    FEValues<dim> fe_v (
        mapping, fe, quadrature, update_flags);

    FEFaceValues<dim> fe_v_face (
        mapping, fe, face_quadrature, face_update_flags);

    FEFaceValues<dim> fe_v_face_neighbor (
        mapping, fe, face_quadrature, neighbor_face_update_flags);


    FullMatrix<double> ui_vi_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> ue_vi_matrix (dofs_per_cell, dofs_per_cell);

    FullMatrix<double> ui_ve_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> ue_ve_matrix (dofs_per_cell, dofs_per_cell);

    Vector<double>  cell_vector (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // LOOP OVER ALL THE CELLS
    for (; cell!=endc; ++cell)
    {
        ui_vi_matrix = 0;
        cell_vector = 0;

        fe_v.reinit (cell);

        // ASSEMBLE THE BODY INTEGRALS:
        //  \int_{omega_{e}} - \nabla v * beta  u dx
        //  \int_{omega_{e}} v f dx
        dg.assemble_cell_term(fe_v,
                              ui_vi_matrix,
                              cell_vector);

        //		ui_vi_matrix.print(std::cout);

        //			std::cout << std::endl;

        cell->get_dof_indices (dofs);
//				std::cout << "cell: " << cell << std::endl;
        // Now assemble the terms corresponding to the floxes
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            typename DoFHandler<dim>::face_iterator face=
                cell->face(face_no);

            //					std::cout << "face_no: " << face_no;
            // face is on the boundary
            if (face->at_boundary())
            {
                //							std::cout << " is on boundary " << std::endl;
                fe_v_face.reinit (cell, face_no);
                /*
                                dg.assemble_boundary_term(fe_v_face,
                                                          ui_vi_matrix,
                                                          cell_vector);
                  */
            }
            else // on interior face
            {
                Assert (cell->neighbor(face_no).state() == IteratorState::valid,
                        ExcInternalError());

                typename DoFHandler<dim>::cell_iterator neighbor=
                    cell->neighbor(face_no);

                //Thus we have the additional condition, that the
                // cell with the lower index does the work. In the rare
                // case that both cells have the same index, the cell with
                // lower level is selected.
                if ( (neighbor->index() > cell->index() )/* ||
                        (neighbor->level() < cell->level() &&
                        neighbor->index() == cell->index()))*/ )
                {
                    //	std::cout << " is on interior face" << std::endl;

                    // Here we know, that the neighbor is not coarser so we
                    // can use the usual @p neighbor_of_neighbor
                    // function. However, we could also use the more
                    // general @p neighbor_face_no function.
                    const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

                    ue_vi_matrix = 0;
                    ui_ve_matrix = 0;
                    ue_ve_matrix = 0;

                    fe_v_face.reinit (cell, face_no);
                    fe_v_face_neighbor.reinit (neighbor, neighbor2);

                    dg.assemble_interior_face_term(fe_v_face,
                                                   fe_v_face_neighbor,
                                                   ui_vi_matrix,
                                                   ue_vi_matrix,
                                                   ui_ve_matrix,
                                                   ue_ve_matrix);

                    neighbor->get_dof_indices (dofs_neighbor);

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            system_matrix.add(dofs[i], dofs_neighbor[j],
                                              ue_vi_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs[j],
                                              ui_ve_matrix(i,j));
                            system_matrix.add(dofs_neighbor[i], dofs_neighbor[j],
                                              ue_ve_matrix(i,j));
                        } // for j
                    } // for i
                } // if
                //:		else
                //	std::cout << " is already done" << std::endl;
            } // interior faces

        } // face_no

        // Now add on the into the system matrix the body integrals and flux inside cell
        for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
                system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i,j));

        // add on the boundary conditions and the rhs function
        for (unsigned int i=0; i<dofs_per_cell; ++i)
            right_hand_side(dofs[i]) += cell_vector(i);

    } // cell
} // assemble_system2



template <int dim>
void DGMethod<dim>::solve (Vector<double> &solution)
{
    SolverControl           solver_control (1000, 1e-12, false, false);
    SolverRichardson<>      solver (solver_control);

    PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;

    preconditioner.initialize(system_matrix, fe.dofs_per_cell);

    solver.solve (system_matrix, solution, right_hand_side,
                  preconditioner);
}


template <int dim>
void DGMethod<dim>::output_results () const
{
    std::ofstream eps_output ("grid.eps");

    GridOut grid_out;
    grid_out.write_eps (triangulation, eps_output);


    std::ofstream solution_output("solution.vtk");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution2, "u");

    data_out.build_patches (degree);

    data_out.write_vtk(solution_output);
}

template <int dim>
void DGMethod<dim>::output_system()
{
    std::ofstream output_system("Stiffness_Matrix.mtx");
    system_matrix.print_formatted(output_system);
    output_system.close();
    output_system.open("Rhs_Vector.vec");
    right_hand_side.print(output_system);
    output_system.close();
}

template <int dim>
void DGMethod<dim>::run ()
{
    GridGenerator::hyper_cube(triangulation,0,1);
    triangulation.refine_global (4);

    std::cout << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

    setup_system ();

    std::cout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

    Timer assemble_timer;
    assemble_system ();
    std::cout << "Time of assemble_system2: "
              << assemble_timer()
              << std::endl;

    //output_system();

    solve (solution2);

    output_results();
}
}



