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
                                     const std::vector<int> & fracture_face_markers)
    : _grid(grid), _config( config ),
      _fracture_face_markers(fracture_face_markers.begin(), fracture_face_markers.end())
{
  #ifndef WITH_EIGEN
  throw std::runtime_error("Cannot use DFEM method without linking to Eigen");
  #endif

  // analyze_cell_(_grid.cell(546));
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
  // std::cout << std::endl;
  size_t item = 0;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    // std::cout << item++ << " (" << _grid.n_active_cells() << ")"<< std::endl;
    progress.set_progress(item++);

    const std::unique_ptr<FiniteElementBase> p_discr = build_element(*cell);

    FiniteElementData cell_fem_data = p_discr->get_cell_data();
    cell_fem_data.element_index = cell->index();
    _cell_data[cell->index()] = std::move(cell_fem_data);

    std::vector<FiniteElementData> face_data = p_discr->get_face_data();
    std::vector<FiniteElementData> frac_data = p_discr->get_fracture_data();
    size_t iface = 0;
    for ( const mesh::Face * face : cell->faces() )
    {
      const size_t face_index = face->index();
      if ( _face_data[face->index()].points.empty() )
      {
        face_data[iface].element_index = face_index;
       _face_data[face_index] = face_data[iface];
      }

      if ( _fracture_face_markers.find( face->marker() ) != _fracture_face_markers.end() )
        {
          frac_data[iface].element_index = cell->index();
          _frac_data[face_index].push_back(frac_data[iface]);
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
  PolyhedralElementDirect de(cell, _config);
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
  const bool need_face_values = true;
  const bool need_fracture_values = has_fracture_face_(cell);
  std::unique_ptr<FiniteElementBase> p_discr;
  if (_config.method == FEMMethod::polyhedral_finite_element)
  {
#ifdef WITH_EIGEN
    if (_config.solver == direct || _config.solver == cg)
      p_discr = std::make_unique<PolyhedralElementDirect>(cell, _config, need_face_values,
                                                          need_fracture_values);
    else if (_config.solver == msrsb)
      p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _config);
#else
    throw std::runtime_error("Cannot use PDFEM method without linking to Eigen");
#endif

  }
  else if (_config.method == strong_discontinuity)
    p_discr = std::make_unique<StandardFiniteElement>(cell, need_face_values, need_fracture_values);
  else if (_config.method == mixed)
  {
    if (cell.vtk_id() == angem::GeneralPolyhedronID)
    {
#ifdef WITH_EIGEN
      if (_config.solver == direct || _config.solver == cg)
        p_discr = std::make_unique<PolyhedralElementDirect>(cell, _config, need_face_values,
                                                            need_fracture_values);
      else if (_config.solver == msrsb)
        p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _config);
#else
    throw std::runtime_error("Cannot use PFEM method without linking to Eigen");
#endif
    }
    else
      p_discr = std::make_unique<StandardFiniteElement>(cell, need_face_values, need_fracture_values);
  }

  return p_discr;
}

bool DiscretizationFEM::has_fracture_face_(const mesh::Cell & cell)
{
  for (const mesh::Face * face : cell.faces() )
    if ( _fracture_face_markers.find( face->marker() ) != _fracture_face_markers.end()  )
      return true;
  return false;
}

}  // end namepsace discretization
