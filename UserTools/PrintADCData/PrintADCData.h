#ifndef PrintADCData_H
#define PrintADCData_H

#include <string>
#include <iostream>
#include <thread>

#include "Tool.h"

#include "TApplication.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"


/**
 * \class PrintADCData
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class PrintADCData: public Tool {


 public:

  PrintADCData(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:

  bool visualize;
  std::string outputfile;
  TFile *file_out = nullptr;

  long totalentries=0;

  uint32_t RunNum;
  uint32_t SubrunNum;
  long EntryNum;
  std::map<unsigned long, std::vector<Waveform<uint16_t>> > RawADCData;
  
  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;

  // internal members:
  int SampleLength;
  std::vector<int> numberline;
  std::vector<int> upcastdata;
  int maxwfrmamp=0;
  int WaveformSource;
  
  // ROOT TApplication variables
  // ---------------------------
  TApplication* rootTApp=nullptr;
  double canvwidth, canvheight;
  TCanvas* mb_canv=nullptr;
  TGraph* mb_graph=nullptr;
  

};


#endif