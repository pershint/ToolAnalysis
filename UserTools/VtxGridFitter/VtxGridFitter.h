#ifndef VtxGridFitter_H
#define VtxGridFitter_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMath.h"
#include "Direction.h"
#include "Position.h"
#include "RecoVertex.h"

class VtxGridFitter: public Tool {


 public:

  VtxGridFitter();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  /// Clear true vertex and vertex seed list
  void Reset();


  /// Function that generates an array of directions approximately
  /// evenly distributed in all directions from the origin.  This is
  /// done by populating a unit sphere using Vogel's method.
  /// source: https://stackoverflow.com/questions/9600801/evenly-
  /// distributing-n-points-on-a-sphere
  void GenerateDirectionSeeds(int numDirs);
  
  /// Generate a list of time seeds to fit to in a time range defined
  double MinimumTime=-10.;
  double MaximumTime=40.;
  double TimeResolution=0.25;
  void GenerateRoughTimeSeeds();

  /// Generate a list of time seeds sampled from a gaussian at the given
  /// mean with TimeResolution standard deviation
  void GenerateFineTimeSeeds(double mean, double timeres, int nFineSeeds);

  /// Scales and translates the input vertex position with the provided parameters
  RecoVertex* TransformPosition(RecoVertex* InputVertex, double scale, Position translation);

  /// Given a RecoVertex with position parameters, return the direction
  /// and time with the highest FOM
  RecoVertex* FindBestVtxAtPos(RecoVertex* SeedPosition, std::vector<double> timeList,
          double coneweight, double vtxweight);

  void PushExtendedVertex(RecoVertex* vtx, bool savetodisk);

 private:
  TRandom3 frand;  ///< Random number generator
  int NumDirSeeds;

  RecoVertex* BestExtendedVertex = nullptr;
  RecoVertex* scaledSeedPos = nullptr;

  /// Vertex Geometry shared by Fitter tools
  VertexGeometry* fVtxGeo;

  std::vector<RecoDigit>* fDigitList = 0;
  std::vector<double> vSeedRTimeList;
  std::vector<double> vSeedFTimeList;
  std::vector<Direction> vSeedDirList;
  std::vector<RecoVertex>* vSeedVtxList = nullptr;

  /// verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity=-1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;

};


#endif
