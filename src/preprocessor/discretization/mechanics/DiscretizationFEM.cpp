#include "DiscretizationFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <stdexcept>
#ifdef WITH_EIGEN
#include "PolyhedralElementDirect.hpp"
#include "PolyhedralElementMSRSB.hpp"
#endif
#include "StandardFiniteElement.hpp"
#include "MeshStatsComputer.hpp"
#include "logger/ProgressBar.hpp"  // provides ProgressBar
#include "VTKWriter.hpp"                                       // debugging, provides io::VTKWriter


namespace discretization
{
using Point = angem::Point<3,double>;

DiscretizationFEM::DiscretizationFEM(const mesh::Mesh & grid, const FiniteElementConfig & config,
                                     const std::vector<int> & fracture_face_markers,
                                     const std::vector<int> & neumann_face_markers)
    : _grid(grid), _config( config )
{
  #ifndef WITH_EIGEN
  throw std::runtime_error("Cannot use DFEM method without linking to Eigen");
  #endif

  for (const int marker : fracture_face_markers)
    _fracture_face_orientation[marker] = {angem::Basis<3, double>(), false};
  for (const int marker : neumann_face_markers)
    _neumann_face_orientation[marker] = {angem::Basis<3, double>(), false};

  // analyze_cell_(_grid.cell(0));
  // if (_config.method != strong_discontinuity)
  // {
  //   auto cell = _grid.begin_active_cells();
  //   PolyhedralElementDirect de(*cell, _config);
  //   mesh::MeshStatsComputer stats(de.get_grid());
  //   std::cout << "Average edge length = " << stats.get_average_edge_length() << std::endl;
  // }

  _face_data.resize( _grid.n_faces() );
  _cell_data.resize( _grid.n_cells() );
  _frac_data.resize( _grid.n_faces() );
  logging::ProgressBar progress("Build Finite Elements", _grid.n_active_cells());
  angem::Basis<3, double> face_basis;
  // std::cout << std::endl;
  size_t item = 0;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    progress.set_progress(item++);

    const std::unique_ptr<FiniteElementBase> p_discr = build_element(*cell);

    FiniteElementData cell_fem_data = p_discr->get_cell_data();
    cell_fem_data.element_index = cell->index();
    _cell_data[cell->index()] = std::move(cell_fem_data);

    size_t iface = 0;
    for ( const mesh::Face * face : cell->faces() )
    {
      const size_t face_index = face->index();
      const bool is_fracture = _fracture_face_orientation.find( face->marker() ) !=
                               _fracture_face_orientation.end();
      const bool is_neumann = _neumann_face_orientation.find( face->marker() ) !=
                              _neumann_face_orientation.end();
      if (is_fracture)
        face_basis = get_basis_(*face, _fracture_face_orientation.find(face->marker())->second);
      if (is_neumann)
        face_basis = get_basis_(*face, _neumann_face_orientation.find(face->marker())->second);

      if (is_neumann || is_fracture)
        if ( _face_data[face->index()].points.empty() )
        {
          _face_data[face_index] = p_discr->get_face_data(iface, face_basis[2]);
          _face_data[face_index].element_index = face_index;
        }

      if (is_fracture)
        {
          _frac_data[face_index].push_back(p_discr->get_fracture_data(iface, face_basis[2]));
          _frac_data[face_index].back().element_index = cell->index();
        }

      iface++;
    }
  }

  progress.finalize();
}

void DiscretizationFEM::analyze_cell_(const mesh::Cell & cell)
{
  IO::VTKWriter::write_geometry(_grid, cell, "output/geometry-" + std::to_string(cell.index()) + ".vtk");

#ifdef WITH_EIGEN
  PolyhedralElementDirect de(cell, _grid, _config);
  // PolyhedralElementMSRSB de(cell, _config);
  de.save_shape_functions("output/shape_functions-" + std::to_string(cell.index())+ ".vtk");
  exit(0);

  StandardFiniteElement fe(cell);
  FiniteElementData an_data =  fe.get_cell_data();
  auto verts = cell.vertices();
  std::cout << "+======COORDINATES" << std::endl;
  // std::cout << "analytic" << std::endl;
  // for (size_t q=0; q<an_data.points.size(); ++q)
  // {
  //   Point p (0,00,0);
  //   for (size_t v=0; v<verts.size(); ++v)
  //   {
  //     p += _grid.vertex(verts[v]) * an_data.points[q].values[v];
  //   }
  //   std::cout << "p = " << p << std::endl;
  // }

  auto data =  de.get_cell_data();
  auto qpoints = de.get_cell_integration_points();
  std::cout << "============================" << std::endl;
  std::cout << "Recovery integration points:" << std::endl;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    Point p (0,00,0);
    for (size_t v=0; v<verts.size(); ++v)
    {
      assert( data.points[q].values.size() == verts.size() );
      p += _grid.vertex(verts[v]) * data.points[q].values[v];
    }

    std::cout << "diff = " << (p - qpoints[q]).norm()<< std::endl;
  }

  std::cout << "============================" << std::endl;
  std::cout << "Weights:" << std::endl;
  double wsum = 0;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    std::cout << data.points[q].weight << " ";
    wsum += data.points[q].weight;
  }
  std::cout << std::endl;

  std::cout << "weights sum = " << wsum << "  (should be " << data.center.weight << ")" << std::endl;


  // std::cout << "============================" << std::endl;
  // std::cout << "Values:" << std::endl;
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << "\n";
  //   for (size_t v=0; v<verts.size(); ++v)
  //   {
  //     std::cout << data.points[q].values[v] << " "
  //               << an_data.points[q].values[v] << " "
  //               << std::fabs( (data.points[q].values[v] - an_data.points[q].values[v]) / an_data.points[q].values[v]  )
  //               << std::endl;
  //   }
  // }

  std::cout << "============================" << std::endl;
  std::cout << "Gradients:" << std::endl;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    std::cout << q << "\n";
    for (size_t v=0; v<verts.size(); ++v)
    {
      std::cout <<  data.points[q].grads[v] << "\t";
    }
    std::cout << std::endl;
  }


  // std::cout << "======= shape functions ==========" << std::endl;
  // std::cout << "analytic" << std::endl;
  // data =  fe.get_cell_data();
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << ": ";
  //   for (size_t v=0; v<verts.size(); ++v)
  //     std::cout << data.points[q].values[v] << " ";
  //   std::cout << std::endl;
  // }
  // std::cout << "numeric" << std::endl;
  // data =  de.get_cell_data();
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << ": ";
  //   for (size_t v=0; v<verts.size(); ++v)
  //     std::cout << data.points[q].values[v] << " ";
  //   std::cout << std::endl;
  // }

  // This proves partition of unity
  std::cout << "======= partition of unity ===========" << std::endl;
  data =  de.get_cell_data();
  {
    double maxsum = 0;
    for (auto & qp : data.points)
    {
      double sum = 0;
      for (size_t v=0; v<qp.values.size(); ++v)
      {
        sum += qp.values[v];
      }
      maxsum = std::max(sum, maxsum);
    }
    std::cout << "numeric sum = " << std::scientific << std::fabs(maxsum - 1.0) << std::endl << std::defaultfloat;
  }

  std::cout << "======= patch test ========" << std::endl;
  {
    Point maxsum = {0,0,0};
    for (auto & qp : data.points)
    {
      for (size_t i=0; i<3; ++i)
      {
        double sum = 0;
        for (size_t v=0; v<qp.values.size(); ++v)
          sum += qp.grads[v][i];
        maxsum[i] = std::max(maxsum[i], std::fabs(sum));
      }

    }
    std::cout << "grad maxsum = " << maxsum << std::endl;
  }

  exit(0);
#endif
}

std::unique_ptr<FiniteElementBase> DiscretizationFEM::build_element(const mesh::Cell & cell)
{
  const bool need_fracture_values = false;
  const bool need_face_values = false;
  std::unique_ptr<FiniteElementBase> p_discr;
  if (_config.method == FEMMethod::polyhedral_finite_element)
  {
#ifdef WITH_EIGEN
    if (_config.solver == direct || _config.solver == cg)
      p_discr = std::make_unique<PolyhedralElementDirect>(cell, _grid, _config);
    else if (_config.solver == msrsb)
      p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _grid, _config);
#else
    throw std::runtime_error("Cannot use PDFEM method without linking to Eigen");
#endif

  }
  else if (_config.method == strong_discontinuity)
    p_discr = std::make_unique<StandardFiniteElement>(cell);
  else if (_config.method == mixed)
  {
    if (cell.vtk_id() == angem::GeneralPolyhedronID)
    {
#ifdef WITH_EIGEN
      if (_config.solver == direct || _config.solver == cg)
        p_discr = std::make_unique<PolyhedralElementDirect>(cell, _grid, _config);
      else if (_config.solver == msrsb)
        p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _grid, _config);
#else
    throw std::runtime_error("Cannot use PFEM method without linking to Eigen");
#endif
    }
    else
      p_discr = std::make_unique<StandardFiniteElement>(cell);
  }

  return p_discr;
}

const angem::Basis<3, double> &
DiscretizationFEM::get_basis_(const mesh::Face & face,
                              FaceOrientation &orientation) noexcept
{
  if (!orientation.assigned)
  {
    orientation.basis = face.polygon().plane().get_basis();
    orientation.assigned = true;
  }

  return orientation.basis;
}


}  // end namepsace discretization
