#pragma once

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <map>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <set>

#include "element.hpp"
#include "renum.hpp"
#include "transes.hpp"
#include "GElement.hpp"

#include "angem/Point.hpp"
#include "angem/PolyGroup.hpp"
#include "angem/Collisions.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Mesh.hpp"
#include <SimdataConfig.hpp>

struct GmConstraint
{
  double dPenalty;
  int nVertices;
  vector<int> vIndexes;
  vector<int> vVertices;
  vector<vector<double> > vvCoefficient;
};


struct RockProps
{
  int model;
  std::vector<double> v_props;
  double poro;
  double perm, perm_x, perm_y, perm_z;
  double thc, thc_x, thc_y, thc_z;
  double temp;
  double heat_capacity;

  double biot_plas;
  double biot_flow;
  double biot;
  double young;
  double poisson;
  double density;
  double poron;

  double temperature;
  double pressure;
  double volmult;

  double ref_pres;
  double ref_temp;
  double ref_stres;
  double ref_strain;

  double cohesion;
  double friction;
  double dilation;
  double thermal_expansion;
  double pore_thermal_expansion;

  vector<double> zmf;
  vector<double> stress;
};


struct PhysicalFace
{
  int ntype;
  int nface;
  int nmark;
  int nfluid;
  angem::Point<3,double> condition;
  double aperture;
  double conductivity;
};


struct SimpleWell
{
  vector<double> vRadiusPoisk;
  vector<double> vWellCoordinate;
  vector<int> vID;
  vector<int> vWi;
  double datum;
  string Type;
  double radius_poisk;
};


struct EmbeddedFracture
{
  std::vector<std::size_t>            cells;
  std::vector<angem::Point<3,double>> points;
  std::vector<double>                 dip;
  std::vector<double>                 strike;
  double                              cohesion;
  double                              friction_angle;
  double                              dilation_angle;
  // cells -> points
  // these two entries represent mesh within the frac
  mesh::SurfaceMesh<double>            mesh;
  std::vector<angem::Point<3,double>>  vVertices;
  // cells -> vertex indiced (fracture polygons)
  std::vector<std::vector<std::size_t>> vIndices;
  double                               aperture;     // m
  double                               conductivity;  // md-m
};




class SimData
{
public:
  // SimData(const string & inputstream, const SimdataConfig & config);
  SimData(mesh::Mesh & grid, const SimdataConfig & config);
  ~SimData();
  // void readSetupValues();
  void readTotalData();
  void readTotalTemp();

  void readGmshFile();
  void convertGmsh2Sim();

  void initilizeBoundaryConditions();

  void defineRockProperties();
  void defineEmbeddedFractureProperties();
  void computeCellClipping();
  // void mergeSmallFracCells();
  void definePhysicalFacets();
  void defineStressAndDispOnBoundary();

  void splitInternalFaces();

  void handleConnections();
  void computeReservoirTransmissibilities();
  void computeFracFracTran(const std::size_t                 frac,
                           const EmbeddedFracture          & efrac,
                           const mesh::SurfaceMesh<double> & mesh,
                           FlowData                        & frac_flow_data);
  void computeEDFMTransmissibilities(const std::vector<angem::PolyGroup<double>> & splits,
                                     const int   frac_ind);
  void computeInterEDFMTransmissibilities();
  // void computeTransBetweenDifferentEfracs();


  void createSimpleWells();

  void extractInternalFaces();
  std::size_t n_default_vars() const;

  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;

  angem::Point<3,double> get_permeability(const std::size_t cell) const;
  double get_volume_factor(const std::size_t cell) const;
  void meshFractures();

protected:
   void methodElementCenter(int nelem, vector<Gelement> &vsElement);
   void methodFaceNormalVector(int nelem, vector<Gelement> &vsElement);
   void methodChangeFacesNormalVector();

   void methodRandomRockProperties();
   double createLognormalDistribution(double E, double S);

  void handle_edfm_face_intersection(const std::size_t ifrac,
                                     const std::size_t jfrac,
                                     const std::vector<std::size_t> & icells,
                                     const std::vector<std::size_t> & jcells);

  void compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
                                              const std::vector<std::vector<std::size_t>> & polys,
                                              const std::vector<int>                      & markers,
                                              FlowData                                    & flow_data) const;
  // std::size_t get_flow_element_index(const std::size_t ifrac,
  //                                    const std::size_t ielement) const;


   int checkReservedBoundaryName(int nmarker)
   {
     if( nmarker > 1000000 ) return(-1);
     return(1);
   }
   renum * pRenum;

public:
  mesh::Mesh & grid;
  // int corner_cell;
  // vector<double> grade_total, temp_total;
  // int nNodes;
  // std::size_t nNodes;

  // double dNotNumber;

  // int nBndNodes;
  // vector<vector<double> > vvBndFaceNodesCoor;
  // vector<vector<int> >    vvBndFaceNodes;
  // vector<int>    vBndFaceCode;

  // vector<vector<double> > vvInputCoorNodes;
  // vector<vector<int> >    vvElementNodes;
  // vector<int> vElementCode;

  string outstream;

  // Internal Data
  // std::size_t nVertices;
  // vector< angem::Point<3,double> > vvVrtxCoords;
  // vector<vector<double> > vvVrtxDisp;
  // vector<int> vConstraintVertex;

  // int nCells;
  // vector<Gelement> vsCellCustom;
  std::vector<RockProps> vsCellRockProps;
  std::vector<std::string> rockPropNames;

  vector<EmbeddedFracture> vEfrac;
  FlowData flow_data;
  FlowData new_flow_data;

  std::size_t nDirichletFaces;
  std::size_t nDirichletNodes;
  std::size_t nNeumannFaces;

  std::unordered_map<uint256_t, PhysicalFace> boundary_faces;
  std::unordered_map<uint256_t, PhysicalFace> dfm_faces;
  vector<PhysicalFace> vsPhysicalBoundary;

  //wells
  vector<SimpleWell> vsWell;

  SimdataConfig config;

protected:
  StandardElements * pStdElement;
};
