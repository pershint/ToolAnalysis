#ifndef MrdPaddlePlot_H
#define MrdPaddlePlot_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TClonesArray.h"
#include "MRDSubEventClass.hh"      // a class for defining subevents
#include "MRDTrackClass.hh"         // a class for defining MRD tracks
// for drawing
#include "TApplication.h"
#include "TSystem.h"
//#include "TEveLine.h"
class TEveLine;

class MrdPaddlePlot: public Tool {
	
	public:
	MrdPaddlePlot();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	// Variables stored in Config file
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int verbosity=1;
	std::string gdmlpath;
	bool saveimages;
	const char* plotDirectory;  // where to save images and plots
	
	// Variables from ANNIEEVENT
	// ~~~~~~~~~~~~~~~~~~~~~~~~~
	std::string MCFile;     //-> currentfilestring
	uint32_t RunNumber;     // -> runnum     (keep as the type is different: TODO align types)
	uint32_t SubrunNumber;  // -> subrunnum  ( " " )
	uint32_t EventNumber;   // -> eventnum   ( " " )
	uint16_t MCTriggernum;  // -> triggernum ( " " )
	uint64_t MCEventNum;    // not yet in MRDTrackClass 
	
	// MRD TRACK DRAWING
	// ~~~~~~~~~~~~~~~~~~~~~~~~
	int numsubevs;
	int numtracksinev;
	bool highlight_true_paddles;
	std::map<unsigned long,int> channelkey_to_mrdpmtid;
	
	// we need to subtract the an offset to draw TEveLines over gdml
	Position buildingoffset;
	
	bool printTracks;
	bool printTClonesTracks;
	bool drawPaddlePlot;
	bool drawGdmlOverlay;
	bool drawStatistics;
	TApplication* rootTApp=nullptr;
	
	TCanvas* gdmlcanv=nullptr;
	TCanvas* mrdTrackCanv=nullptr;
	Double_t canvwidth = 700;
	Double_t canvheight = 600;
	TClonesArray* thesubeventarray;  // retrieve from track finder
	std::vector<TEveLine*> thiseventstracks;
	
	// Summary histograms on tracks found
	TH1D* hnumhclusters=nullptr;
	TH1D* hnumvclusters=nullptr;
	TH1D* hnumhcells=nullptr;
	TH1D* hnumvcells=nullptr;
	TH1D* hpaddleids=nullptr;
	TH1D* hpaddleinlayeridsh=nullptr;
	TH1D* hpaddleinlayeridsv=nullptr;
	TH1D* hdigittimes=nullptr;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
	
};


#endif
