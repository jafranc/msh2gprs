/************************************************************************/
/* DISCRETE FEATURE MODEL                                               */
/* by Mohammad Karimi-Fard (karimi@stanford.edu)                        */
/* March 2007, Stanford, CA.                                            */
/*                                                                      */
/* December 2012. Modifications by Timur Garipov                        */
/* Geomechanical interface. Convert Tetgen data to Karimi data           */
/************************************************************************/
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <algorithm>
#include <math.h>
#include <iterator>
#include <vector>
#include <set>
#include <time.h>

// #include "simdata.hpp"


struct FlowData
{
  std::vector<double> volumes, poro, depth;
  // regular transmissibilities
  std::vector<std::size_t> ielement, jelement;
  std::vector<double>      trans_ij, conduct_ij;
  // connections
  std::vector<int> connection_type;
  std::vector<int> connection_n;
};


class CalcTranses
{
public:
   CalcTranses();
  ~CalcTranses();
  void createKarimiData();
  void createKarimiApproximation();
  // this guy writes text output in a series of files
  void writeOutputFiles(const std::string & output_path) const;
  void extractData(FlowData & data) const;
  void init();


public:
  int NbNodes;
  int NbPolygons;
  int NbPolyhedra;
  int NbZones;
  int NbFracs;
  int NbOptions;
  //coordinates
  std::vector<double> vCoordinatesX;
  std::vector<double> vCoordinatesY;
  std::vector<double> vCoordinatesZ;
  //faces
  std::vector<int> vNbVFaces;
  std::vector<std::vector<int> > vvVFaces;
  std::vector<int> vCodePolygon;
  //elements
  std::vector<int> vNbFNodes;
  std::vector<std::vector<std::size_t> > vvFNodes;
  std::vector<int> vCodePolyhedron;
  //properties
  std::vector<int> vZoneCode;
  std::vector<double> vZVolumeFactor;
  std::vector<double> vTimurConnectionFactor;
  std::vector<double> vZPorosity;
  std::vector<int> vZPermCode;
  std::vector<double> vZPermeability;
  std::vector<double> vZConduction;
protected:
  double ABS(double v);
  void CheckIt(int TEST);
  void ProjectionA(double mx,double my,double mz,
			double px,double py,double pz,
			double ix,double iy,double iz,
			double nx,double ny,double nz,
			double *hx,double *hy,double *hz);
  void ProjectionB(double mx,double my,double mz,
			double ix,double iy,double iz,
			double ux,double uy,double uz,
			double *hx,double *hy,double *hz);
  void ComputeBasicGeometry();
  void ComputeControlVolumeList();
  void PrepareConnectionList();
  void ConstructConnectionList();
  void VolumeCorrection();
  void ComputeContinuityNode();
  void ComputeDirectionalPermeability();
  void ComputeTransmissibilityPart();
  void ComputeTransmissibilityList();

protected:
// SimData * pSim;

int		NbCVs,NbVolumes,NbInterfaces,NbEquations,NbFeatures,NbCF,NbConnections,NbIntersections;
int		NbTransmissibility,NbMetric;
int		OptionVC,OptionGO,OptionMC;

///// Grid information /////
double		*X,*Y,*Z;
int		*NbFNodes;
int		**FNodes;
int		*NbVFaces;
int		**VFaces;

int		*CodePolygon;
int		*CodePolyhedron;

///// Control volume numbers /////

int		*EQF;
int		*EQV;

///// Additional polygon information /////

double		*FArea;
double		*FXG,*FYG,*FZG;
double		*Fnx,*Fny,*Fnz;		// unit vector

///// Additional polyhedron information /////

double		*VVolume;
double		*VXG,*VYG,*VZG;

///// Definition of the control volumes /////

double		Tolerance;

int		*CVType;
int		*CVZone;
double		*CVx,*CVy,*CVz;
double		*CVVolume;

int		*ZoneCode;
double		*ZVolumeFactor;
double		*ZPorosity;
int		*ZPermCode;
double		**ZPermeability;
double		**ZConduction;
double		K1,K2,K3,K4,K5,K6;
double		Kx,Ky,Kz;

///// Definition of the connections /////

int		*ConType;
int		*ConN;
int		**ConCV;
double		**ConTr;

double		**ConGeom;
double		**ConMult;

double		**ConArea;
double		**ConPerm;
double		*ConP1x,*ConP1y,*ConP1z;
double		*ConP2x,*ConP2y,*ConP2z;
double		*ConIx,*ConIy,*ConIz;
double		*ConVx,*ConVy,*ConVz;
double		*Conhx,*Conhy,*Conhz;

///// Transmissibility List /////

int		*iTr,*jTr;
double		*Tij,SumTr, *TConductionIJ;

int		NbEdges;	// Actives
int		*ListV1,*ListV2;
int		*ListE1,*ListE2,*ListF;

double		TotalVolume,FaceArea,delta;
double		LocalDistance,t,xt,yt,zt;

clock_t		t1,t2,t3,t4,t5,t6,t7,t8,t9;
clock_t		Deb_Computing,Fin_Computing;

 public:
  double fracporo;
  static constexpr double transmissibility_conversion_factor =
      0.0085267146719160104986876640419948;
};
