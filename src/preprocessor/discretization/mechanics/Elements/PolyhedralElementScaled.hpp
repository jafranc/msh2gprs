#pragma once

#include "PolyhedralElementBase.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

/** This class implements the scaling of PFEM shape functions
 * for the current polyhedron given a master element.
 * There is no polymorformism check (should be performed prior),
 * and there is no renumbering.
 * The idea of scaling is the mapping between the current element and
 * the master element, same as in stanfard FEM.
 */
class PolyhedralElementScaled : public PolyhedralElementBase {
 public:
  /**
   * @brief Constructor. Does not do anything
   * Input:
   * \param[in] cell : grid cell to be discretized
   * \param[in] master : polyhedral element built on a combinatorially-
   *                     identical cell. vertex and face numbering must
   *                     be the same.
   */
  PolyhedralElementScaled(const mesh::Cell & cell,
                          const mesh::Mesh & parent_grid,
                          PolyhedralElementBase & master,
                          const FiniteElementConfig & config);

  // Destructor
  virtual ~PolyhedralElementScaled() = default;

  // get face data
  FiniteElementData get_face_data(size_t iface,
                                  const angem::Basis<3,double> basis) override;


 protected:
  void build_fe_cell_data_() override;
  void build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                            FEPointData const & master,
                            FEPointData & current,
                            angem::Tensor2<3, double> & du_dx) const;

  void compute_detJ_and_invert_cell_jacobian_(const std::vector<angem::Point<3,double>> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ,
                                              std::vector<angem::Point<3,double>> const & vertex_coord) const;
  void update_shape_grads_(std::vector<angem::Point<3,double>> const & ref_grads,
                           angem::Tensor2<3, double> const & du_dx,
                           std::vector<angem::Point<3,double>> &shape_grads) const;



  PolyhedralElementBase & _master;

};
}  // end namespace discretization
