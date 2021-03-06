/* vim:set noexpandtab tabstop=4 wrap */
#ifndef GEOMETRYCLASS_H
#define GEOMETRYCLASS_H

#include<SerialisableObject.h>
#include "ChannelKey.h"
#include "Detector.h"
#include "Particle.h"
#include "Channel.h"
using namespace std;

enum class geostatus : uint8_t { FULLY_OPERATIONAL, TANK_ONLY, MRD_ONLY, };

class Geometry : public SerialisableObject{
	
	friend class boost::serialization::access;
	
	public:
	// Do we care to have the overloaded empty constructor?
	Geometry() : NextFreeChannelKey(0), NextFreeDetectorKey(0), Version(0.), tank_centre(Position(0,0,0)), tank_radius(0.), tank_halfheight(0.), pmt_enclosed_radius(0.), pmt_enclosed_halfheight(0.), mrd_width(0.), mrd_height(0.), mrd_depth(0.), mrd_start(0.), numtankpmts(0), nummrdpmts(0), numvetopmts(0), numlappds(0), Status(geostatus::FULLY_OPERATIONAL), Detectors(std::vector<std::map<unsigned long,Detector>* >{}) {
		serialise=true;
		RealDetectors.reserve(10);
	}
	
	Geometry(double ver, Position tankc, double tankr, double tankhh, double pmtencr, double pmtenchh, double mrdw, double mrdh, double mrdd, double mrds, int ntankpmts, int nmrdpmts, int nvetopmts, int nlappds, geostatus statin, std::vector<std::map<unsigned long,Detector>* >dets=std::vector<std::map<unsigned long,Detector>* >{});
	
	inline std::vector<std::map<unsigned long,Detector>* >* GetDetectors(){return &Detectors;}
	inline double GetVersion(){return Version;}
	inline geostatus GetStatus(){return Status;}
	inline Position GetTankCentre(){return tank_centre;}
	inline double GetTankRadius(){return tank_radius;}
	inline double GetTankHalfheight(){return tank_halfheight;}
	inline double GetPMTEnclosedRadius(){return pmt_enclosed_radius;}
	inline double GetPMTEnclosedHalfheight(){return pmt_enclosed_halfheight;}
	inline double GetFiducialCutRadius(){return fiducialradius;}
	inline double GetFiducialCutY(){return fiducialcuty;}
	inline double GetFiducialCutZ(){return fiducialcutz;}
	inline double GetMrdWidth(){return mrd_width;}
	inline double GetMrdHeight(){return mrd_height;}
	inline double GetMrdDepth(){return mrd_depth;}
	inline double GetMrdStart(){return mrd_start;}
	inline double GetMrdEnd(){return mrd_start+mrd_depth;}
	
	inline void SetVersion(double VersionIn){Version = VersionIn;}
	inline void SetStatus(geostatus StatusIn){Status=StatusIn;}
	inline void SetTankCentre(Position tank_centrein){tank_centre = tank_centrein;}
	inline void SetTankRadius(double tank_radiusIn){tank_radius = tank_radiusIn;}
	inline void SetTankHalfheight(double tank_halfheightIn){tank_halfheight = tank_halfheightIn;}
	inline void SetPMTEnclosedRadius(double pmt_enclosed_radiusIn){pmt_enclosed_radius = pmt_enclosed_radiusIn;}
	inline void SetPMTEnclosedHalfheight(double pmt_enclosed_halfheightIn){pmt_enclosed_halfheight = pmt_enclosed_halfheightIn;}
	inline void SetFiducialCutRadius(double fidcutradiusin){fiducialradius = fidcutradiusin;}
	inline void SetFiducialCutZ(double fidcutzin){fiducialcutz = fidcutzin;}
	inline void SetFiducialCutY(double fidcutyin){fiducialcuty = fidcutyin;}
	inline void SetMrdWidth(double mrd_widthIn){mrd_width = mrd_widthIn;}
	inline void SetMrdHeight(double mrd_heightIn){mrd_height = mrd_heightIn;}
	inline void SetMrdDepth(double mrd_depthIn){mrd_depth = mrd_depthIn;}
	inline void SetMrdStart(double mrd_startIn){mrd_start = mrd_startIn;}
	void SetDetectors(std::vector<std::map<unsigned long,Detector>* >DetectorsIn){
		Detectors = DetectorsIn;
		detectorcounts.clear();
	}
	
	unsigned long ConsumeNextFreeChannelKey(){
		unsigned long thefreechannelkey = NextFreeChannelKey;
		NextFreeChannelKey++;
		return thefreechannelkey;
	}
	unsigned long ConsumeNextFreeDetectorKey(){
		unsigned long thefreedetectorkey = NextFreeDetectorKey;
		NextFreeDetectorKey++;
		return thefreedetectorkey;
	}
	
	bool AddDetector(Detector detin){
		// Pass a pointer to it's owning geometry to this Detector
		// we need to do this before calling `emplace` as that must do a copy-construction
		detin.SetGeometryPtr(this);
		
		std::string thedetel = detin.GetDetectorElement();
		int detectorsetindex=-1;
		if(DetectorElements.count(thedetel)==0){
			// this is a new detector element - create a new entry in the Detectors vector
			RealDetectors.resize(RealDetectors.size()+1);
			Detectors.push_back(&RealDetectors.back());
			DetectorElements.emplace(thedetel,RealDetectors.size()-1); // maps det. element to vector index
			detectorsetindex = RealDetectors.size()-1;
		} else {
			// we already have a detector set for this detector element: add this detector to that set
			detectorsetindex = DetectorElements.at(thedetel);
		}
		if(Detectors.at(detectorsetindex)->count(detin.GetDetectorID())!=0){  // search it for this key
				std::cerr<<"Geometry error! AddDetector called with non-unique DetectorKey "
						 <<detin.GetDetectorID()<<std::endl;
				return false;
		} else {
			Detectors.at(detectorsetindex)->emplace(detin.GetDetectorID(), detin);
		}
		// increment counters of the number of detectors in each set
		if(detectorcounts.count(thedetel)==0){
			detectorcounts.emplace(thedetel,1);
		} else {
			detectorcounts.at(thedetel)++;
		}
		return true;
	}
	
	inline int GetNumDetectors(){return Detectors.size();}  // FIXME this is the num detector SETS
	Detector* GetDetector(unsigned long DetectorKey);
	Detector* ChannelToDetector(unsigned long ChannelKey);
	Channel* GetChannel(unsigned long ChannelKey);
	void InitChannelMap();
	
	int GetNumDetectorsInSet(std::string SetName){
		// if we've calculated it before, such as during filling with Geometry::AddDetector
		if(detectorcounts.count(SetName)!=0) return detectorcounts.at(SetName);
		// we may not have an existing count if Geometry::SetDetectors was used
		// in which case, scan the Detectors for this set and count its members
		int numdetectorsinthisset=-1;
		for(std::vector<std::map<unsigned long,Detector>>::iterator detsetit  = RealDetectors.begin();
																	detsetit != RealDetectors.end();
																	++detsetit){
			int setsize = detsetit->size();
			if(setsize==0){
				cerr<<"Geometry::GetNumTankPMTs ERROR: Detector set "
					<<std::distance(RealDetectors.begin(),detsetit)
					<<" had no Detectors!"<<endl;
					continue;
			}
			std::map<unsigned long,Detector>::iterator detectorsit = detsetit->begin();
			std::string detectorsetelement = detectorsit->second.GetDetectorElement();
			if(detectorsetelement==SetName){
				numdetectorsinthisset = setsize;
				break; // found the Tank detector set
			}
		}
		if(numdetectorsinthisset<0){
			cerr<<"Geometry::GetNumTankPMTs ERROR: Could not find a detector set with detectorElement \"Tank\"!"<<endl;
			return 0;
		}
		detectorcounts.emplace(SetName,numdetectorsinthisset);
		return numdetectorsinthisset;
	}
	
	/*  FIXME to make work with the new detector styling
	int GetNumTankPMTs(){
		if(numtankpmts!=0) return numtankpmts;
		for(auto adet : Detectors){
			unsigned long  chankey = adet.first; // subdetector 1 is ADC
			if(chankey.GetSubDetectorType()==subdetector::ADC) numtankpmts++;
		}
		return numtankpmts;
	}
	int GetNumMrdPMTs(){
		if(nummrdpmts!=0) return nummrdpmts;
		for(auto&& adet : Detectors){
			ChannelKey chankey = adet.first; // subdetector 0 is TDC...
			if(chankey.GetSubDetectorType()==subdetector::TDC){
				if( true FIXME ) nummrdpmts++;
			}
		}
		nummrdpmts -= 26; // FIXME
		return nummrdpmts;
	}
	// XXX how do we distinguish MRD vs Veto channels?
	int GetNumVetoPMTs(){
		if(numvetopmts!=0) return numvetopmts;
//		for(auto&& adet : Detectors){
//			ChannelKey chankey = adet.first; // subdetector 0 is TDC...
//			if(chankey.GetSubDetectorType()==subdetector::TDC){
//				if( false FIXME ) numvetopmts++;
//			}
//		}
		numvetopmts=26;  // FIXME
		return numvetopmts;
	}
	int GetNumLAPPDs(){
		// FIXME this will be inefficient if num lappds are actually 0!
		if(numlappds!=0) return numlappds;
		for(auto&& adet : Detectors){
			ChannelKey chankey = adet.first; // subdetector 2 is LAPPDs...
			if(chankey.GetSubDetectorType()==subdetector::LAPPD) numlappds++;
		}
		return numlappds;
	}
	*/
	
	bool GetTankContained(Particle part, int startstop=0){
		Position aVertex = (startstop==0) ? part.GetStopVertex() : part.GetStartVertex();
		bool tankcontained = (sqrt(pow(aVertex.X(),2.)+pow(aVertex.Z()-tank_centre.Z(),2.)) < tank_radius) &&
							 (abs(aVertex.Y()-tank_centre.Y())<tank_halfheight);
		return tankcontained;
	}
	bool GetTankContained(Position aVertex){
		bool tankcontained = (sqrt(pow(aVertex.X(),2.)+pow(aVertex.Z()-tank_centre.Z(),2.)) < tank_radius) &&
							 (abs(aVertex.Y()-tank_centre.Y())<tank_halfheight);
		return tankcontained;
	}
	
	bool GetMrdContained(Particle part, int startstop=0){
		Position aVertex = (startstop==0) ? part.GetStopVertex() : part.GetStartVertex();
		bool MrdContained = (abs(aVertex.X())<mrd_width/2.) &&
							(aVertex.Z()>mrd_start) && (aVertex.Z()<(mrd_start+mrd_depth)) &&
							(abs(aVertex.Y())<mrd_height/2.);
		return MrdContained;
	}
	bool GetMrdContained(Position aVertex){
		bool MrdContained = (abs(aVertex.X())<mrd_width/2.) &&
							(aVertex.Z()>mrd_start) && (aVertex.Z()<(mrd_start+mrd_depth)) &&
							(abs(aVertex.Y())<mrd_height/2.);
		return MrdContained;
	}
	
	bool Print(){
		int verbose=0;
		cout<<"Num Detectors : "<<Detectors.size()<<endl;
//		if(verbose){    // FIXME
//			cout<<"Detectors : {"<<endl;
//			for(auto&& adet : Detectors){
//				unsigned long tmp = adet->first;
//				cout<<"ChannelKey : "<<tmp.Print();
//				cout<<"Detector : "<<adet->second.Print();
//			}
//			cout<<"}"<<endl;
//		}
		cout<<"Version : "<<Version<<endl;
		cout<<"Status : "; PrintStatus(Status);
		cout<<"tank_centre : "; tank_centre.Print();
		cout<<"tank_radius : "<<tank_radius<<endl;
		cout<<"tank_halfheight : "<<tank_halfheight<<endl;
		cout<<"pmt_enclosed_radius : "<<pmt_enclosed_radius<<endl;
		cout<<"pmt_enclosed_halfheight : "<<pmt_enclosed_halfheight<<endl;
		cout<<"tank fiducial radius: "<<fiducialradius<<endl;
		cout<<"tank fiducial z cut: "<<fiducialcutz<<endl;
		cout<<"tank fiducial y cut: "<<fiducialcuty<<endl;
		cout<<"mrd_width : "<<mrd_width<<endl;
		cout<<"mrd_height : "<<mrd_height<<endl;
		cout<<"mrd_depth : "<<mrd_depth<<endl;
		cout<<"mrd_start : "<<mrd_start<<endl;
		cout<<"Number of tank PMTs : " << numtankpmts << endl;
		cout<<"Number of MRD PMTs : " << nummrdpmts << endl;
		cout<<"Number of veto PMTs : " << numvetopmts << endl;
		
		return true;
	}
	bool PrintStatus(geostatus status){
		switch(status){
			case (geostatus::FULLY_OPERATIONAL): cout<<"FULLY OPERATIONAL"<<endl; break;
			case (geostatus::TANK_ONLY): cout<<"TANK ONLY"<<endl; break;
			case (geostatus::MRD_ONLY) : cout<<"MRD ONLY"<<endl; break;
		}
		return true;
	}
	
	void PrintChannels();
	
	// helper functions
	Position GlobalToTankCentered(Position posin){ return (posin - tank_centre); }
	void CartesianToPolar(Position posin, double& R, double& Phi, double& Theta, bool tankcentered=false);
	
	private:
	unsigned long NextFreeDetectorKey;
	unsigned long NextFreeChannelKey;
	std::map<unsigned long,Detector*> ChannelMap;
	std::vector<std::map<unsigned long,Detector> > RealDetectors;
	std::vector<std::map<unsigned long,Detector>*> Detectors;
	std::map<std::string, int> detectorcounts;  // not stored
	std::map<std::string, int> DetectorElements;
	double Version;
	geostatus Status;
	Position tank_centre;
	double tank_radius;
	double tank_halfheight;
	double pmt_enclosed_radius;
	double pmt_enclosed_halfheight;
	double mrd_width;
	double mrd_height;
	double mrd_depth;
	double mrd_start;
	int numtankpmts;
	int nummrdpmts;
	int numvetopmts;
	int numlappds;
	double fiducialradius;
	double fiducialcutz;
	double fiducialcuty;
	
	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & RealDetectors;
			ar & Version;
			ar & Status;
			ar & tank_centre;
			ar & tank_radius;
			ar & tank_halfheight;
			ar & pmt_enclosed_radius;
			ar & pmt_enclosed_halfheight;
			ar & mrd_width;
			ar & mrd_height;
			ar & mrd_depth;
			ar & mrd_start;
			ar & fiducialradius;
			ar & fiducialcutz;
			ar & fiducialcuty;
		}
	}
};

#endif
