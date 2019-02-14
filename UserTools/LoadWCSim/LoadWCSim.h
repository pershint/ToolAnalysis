/* vim:set noexpandtab tabstop=4 wrap */
#ifndef LoadWCSim_H
#define LoadWCSim_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TFile.h"
#include "TTree.h"
#include "wcsimT.h"
#include "Particle.h"
#include "Hit.h"
#include "Waveform.h"
#include "TriggerClass.h"
#include "Geometry.h"
#include "MRDspecs.hh"
#include "ChannelKey.h"
#include "Detector.h"
#include "BeamStatus.h"

namespace{
	//PMTs
	constexpr int ADC_CHANNELS_PER_CARD=4;
	constexpr int ADC_CARDS_PER_CRATE=20;
	constexpr int MT_CHANNELS_PER_CARD=4;
	constexpr int MT_CARDS_PER_CRATE=20;
	//LAPPDs
	constexpr int ACDC_CHANNELS_PER_CARD=30;
	constexpr int ACDC_CARDS_PER_CRATE=20;
	constexpr int ACC_CHANNELS_PER_CARD=8;
	constexpr int ACC_CARDS_PER_CRATE=20;
	//TDCs
	constexpr int TDC_CHANNELS_PER_CARD=32;
	constexpr int TDC_CARDS_PER_CRATE=6;
	//HV
	constexpr int CAEN_HV_CHANNELS_PER_CARD=16;
	constexpr int CAEN_HV_CARDS_PER_CRATE=10;
	constexpr int LECROY_HV_CHANNELS_PER_CARD=16;
	constexpr int LECROY_HV_CARDS_PER_CRATE=16;
	constexpr int LAPPD_HV_CHANNELS_PER_CARD=4; // XXX ??? XXX
	constexpr int LAPPD_HV_CARDS_PER_CRATE=10;  // XXX ??? XXX
}

class LoadWCSim: public Tool {
	
	public:
	
	LoadWCSim();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	// variables from config file
	/////////////////////////////
	int verbose=1;
	int HistoricTriggeroffset;
	int LappdNumStrips;           // number of Channels per LAPPD
	double LappdStripLength;      // [mm] for calculating relative x position for dual-ended readout
	double LappdStripSeparation;  // [mm] for calculating relative y position of each stripline
	
	// WCSim variables
	//////////////////
	TFile* file;
	TTree* wcsimtree;
	wcsimT* WCSimEntry; // from makeclass
	WCSimRootTrigger* atrigt, *atrigm, *atrigv;
	WCSimRootGeom* wcsimrootgeom;
	WCSimRootOptions* wcsimrootopts;
	
	// Misc Others
	//////////////
	long NumEvents;
	int numtankpmts;
	int numlappds;
	int nummrdpmts;
	int numvetopmts;
	
	// For constructing ToolChain Geometry
	//////////////////////////////////////
	void ConstructToolChainGeometry();
	std::map<int,unsigned long> lappd_tubeid_to_detectorkey;
	std::map<int,unsigned long> pmt_tubeid_to_channelkey;
	std::map<int,unsigned long> mrd_tubeid_to_channelkey;
	std::map<int,unsigned long> facc_tubeid_to_channelkey;
	// inverse
	std::map<unsigned long,int> detectorkey_to_lappdid;
	std::map<unsigned long,int> channelkey_to_pmtdit;
	std::map<unsigned long,int> channelkey_to_mrdpmtid;
	std::map<unsigned long,int> channelkey_to_faccpmtid;
	
	////////////////
	// things that will be filled into the store from this WCSim file.
	// note: filling everything in the right format is a complex process;
	// just do a simple filling here: this will be properly handled by the conversion
	// from WCSim to Raw and the proper RawReader Tools
	// bool MCFlag=true; 
	std::string MCFile;
	uint64_t MCEventNum;
	uint16_t MCTriggernum;
	uint32_t RunNumber;
	uint32_t SubrunNumber;
	uint32_t EventNumber; // will need to be tracked separately, since we flatten triggers
	TimeClass* EventTime;
	uint64_t EventTimeNs;
	std::vector<MCParticle>* MCParticles;
	std::map<unsigned long,std::vector<Hit>>* TDCData;
	std::map<unsigned long,std::vector<Hit>>* MCHits;
	std::vector<TriggerClass>* TriggerData;
	BeamStatusClass* BeamStatus;
	
};


#endif
