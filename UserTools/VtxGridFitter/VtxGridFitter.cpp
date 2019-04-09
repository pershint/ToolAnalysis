#include "VtxGridFitter.h"

VtxGridFitter::VtxGridFitter():Tool(){}


bool VtxGridFitter::Initialise(std::string configfile, DataModel &data){
  if(verbosity) cout<<"VtxGridFitter Tool: Initializing Tool VtxGridFitter"<<endl;
  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  coneweight = 1.0;
  vtxweight = 1.0;
  NumDirSeeds = 1000;
  numFineSeeds = 100;
  FineScaleReduction = 0.1;

  /// Get the Tool configuration variables
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("NumberDirSeeds", NumDirSeeds);
	m_variables.Get("ConeFitWeight", coneweight);
  m_variables.Get("ExtResidualWeight", vtxweight);
  m_variables.Get("FirstPassMinimumTime", MinimumTime);
  m_variables.Get("FirstPassMaximumTime", MaximumTime);
  m_variables.Get("FirstPassTimeResolution", TimeResolution);
  m_variables.Get("NumFineTimeSeeds", numFineSeeds);
  m_variables.Get("FineScalePosReduction", FineScaleReduction);

  return true;
}


bool VtxGridFitter::Execute(){

  // Retrive digits from RecoEvent
  auto get_ok = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(not get_ok){
    Log("VtxGridFitter Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
    return false;
  }
    
  // Get vertex seed candidates from the store
  auto get_seedlist = m_data->Stores.at("RecoEvent")->Get("vSeedVtxList", vSeedVtxList);
  if(!get_seedlist){ 
    Log("VtxGridFitter Tool: Error retrieving vertex seeds from RecoEvent!",v_error,verbosity);
    Log("VtxGridFitter Tool: Needs to run VtxSeedGenerator first!",v_error,verbosity);
    return false;
  }
 
  // Load digits to VertexGeometry
  fVtxGeo = VertexGeometry::Instance();
  fVtxGeo->LoadDigits(fDigitList);

  // ***** BEGIN FIRST ROUGH PASS OF ENTIRE TANK GRID *****
  //Generate the direction and time seeds
  Log("VtxGridFitter Tool: First pass of grid fitting...",v_message,verbosity); 
  double bestFOM = -1;
  this->GenerateRoughTimeSeeds();
  this->GenerateRoughDirectionSeeds(NumDirSeeds);

  //For each position in vSeedVtxList, run the best FOM search
  //On the output from this, check if it has the best FOM of all
  //positions
  RecoVertex* BestVtxCandidate; 
  for(int i=0;i<vSeedVtxList->size();i++){
    BestVtxCandidate = this->FindBestVtxAtPos(&(vSeedVtxList->at(i)),
            &vSeedRTimeList, &vSeedRoughDirList, coneweight, vtxweight);
    if(BestVtxCandidate->GetFOM() > bestFOM) BestExtendedVertex = BestVtxCandidate;
  }

  // **** BEGIN SECOND FINE PASS OF REDUCED TANK GRID *****
  // Generate the fine time seeds
  // TODO: Add these into the config
  Log("VtxGridFitter Tool: Second pass of grid fitting...",v_message,verbosity); 
  this->GenerateFineTimeSeeds(BestVtxCandidate->GetTime(), TimeResolution,
          numFineSeeds);
  this->GenerateFineDirectionSeeds(NumDirSeeds, BestExtendedVertex->GetDirection());
  // Re-run grid search with position grid shrunk and centered about the
  // Best fit position from the first pass
  for(int i=0;i<vSeedVtxList->size();i++){
    RecoVertex* thisSeedPos = &(vSeedVtxList->at(i));
    scaledSeedPos = this->TransformPosition(thisSeedPos,FineScaleReduction,
            BestVtxCandidate->GetPosition());
    RecoVertex* BestVtxCandidate = this->FindBestVtxAtPos(scaledSeedPos,
            &vSeedFTimeList, &vSeedFineDirList, coneweight, vtxweight);
    if(BestVtxCandidate->GetFOM() > bestFOM) BestExtendedVertex = BestVtxCandidate;
  }
  
  // Push highest FOM vertex into the BoostStore
  this->PushExtendedVertex(BestExtendedVertex,true);
  this->Reset();

  return true;
}


bool VtxGridFitter::Finalise(){
  delete scaledSeedPos;
  delete BestExtendedVertex;
  return true;
}

void VtxGridFitter::Reset(){
  scaledSeedPos->Reset();
  BestExtendedVertex->Reset();
}

// Add extended vertex to RecoEvent store
void VtxGridFitter::PushExtendedVertex(RecoVertex* vtx, bool savetodisk) {  
  // push vertex to RecoEvent store
  Log("VtxGridFitter Tool: Push extended vertex to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("ExtendedVertex", vtx, savetodisk);
}

void VtxGridFitter::GenerateRoughTimeSeeds(){
  Log("VtxGridFitter Tool: Generating Rough Time Seed Array",v_debug,verbosity); 
  for(double timeStep=MinimumTime;timeStep<=MaximumTime;timeStep+=TimeResolution){
    vSeedRTimeList.push_back(timeStep);
  }
}

void VtxGridFitter::GenerateFineTimeSeeds(double mean, double timeres, int numFineSeeds){
  Log("VtxGridFitter Tool: Generating Fine Time Seed Array",v_debug,verbosity); 
  for(int i=0;i<numFineSeeds;i++){
	double timeStep = frand.Gaus(mean, timeres); // time is smeared with 100 ps time resolution. Harded-coded for now.
    vSeedFTimeList.push_back(timeStep);
  }
}

void VtxGridFitter::GenerateRoughDirectionSeeds(int numDirs){
  numDirs=numDirs*2;
  vSeedRoughDirList.clear();
  double offset = 2./numDirs;
  double increment = TMath::Pi()*(3. - sqrt(5.));
  for(int i=0;i<numDirs;i++){
    double y = ((i*offset)-1.) + (offset/2.);
    double r = sqrt(1 - pow(y,2));
    double phi = ((i+1) % numDirs) * increment;
    double z = sin(phi) * r;
    if(z<0) continue; //Only search forward pointing direction
    double x = cos(phi) * r;
    Direction thisDir(x,y,z);
    //Now, rotate this direction into the best rough direction candidate
    vSeedRoughDirList.push_back(thisDir);
  };
}

void VtxGridFitter::GenerateFineDirectionSeeds(int numDirs,Direction seedDir){
  vSeedFineDirList.clear();
  logmessage = " Direction from first pass is: " + std::to_string(seedDir.X()) +
      ","+std::to_string(seedDir.Y())+","+std::to_string(seedDir.Z());
  Log(logmessage,v_debug,verbosity); 
  //Generate the matrix that rotates our smears to the seedDir reference frame
  TVector3 seedDirT3(seedDir.X(), seedDir.Y(), seedDir.Z());
  TVector3 local(0.,0.,1.);
  TMatrixD seedRotator = this->Rotateatob(local,seedDirT3);
  for(int i=0;i<numDirs;i++){
    //First, random sample a phi angle
    double phi = 2.*TMath::Pi()*gRandom->Uniform();
	  //Sample a radius, assuming a standard deviation with
    //the radius given by the approx. surface area covered by each direction vector
    double coneradius = atan(2*TMath::Pi()/numDirs);
    double r = abs(frand.Gaus(0,coneradius));
    double x = r*cos(phi);
    double y = r*sin(phi);
    double z = sqrt(1 - pow(x,2) - pow(y,2));
    TVector3 seedSmear(x,y,z);
    seedSmear = seedSmear.Unit();
    TVector3 seedSmearT3 = seedRotator*seedSmear;
    Direction smearedSeedDir(seedSmearT3.X(), seedSmearT3.Y(), seedSmearT3.Z());
    vSeedFineDirList.push_back(smearedSeedDir);
  };
}

RecoVertex* VtxGridFitter::TransformPosition(RecoVertex* InputVertex, 
        double scale, Position translation){
  Position thisPos = InputVertex->GetPosition();
  Position scaledPos;
  scaledPos.SetX(thisPos.X()*scale + translation.X());
  scaledPos.SetY(thisPos.Y()*scale + translation.Y());
  scaledPos.SetZ(thisPos.Z()*scale + translation.Z());
  RecoVertex* ScaledVtx = new RecoVertex();
  ScaledVtx->SetVertex(scaledPos,InputVertex->GetTime());
  ScaledVtx->SetDirection(InputVertex->GetDirection());
  return ScaledVtx;
}


RecoVertex* VtxGridFitter::FindBestVtxAtPos(RecoVertex* fSeedVertex, 
        std::vector<double>* timeList, std::vector<Direction>* dirList, double coneweight, double vtxweight){
  //Initialize variables holding best fit values
  RecoVertex* BestVertex = new RecoVertex();
  double bestCombTime = -999.;
  Direction bestCombDirection;  
  double bestCombFOM = -1.;
  //Load up our Optimizer, which has access to the FOM calculations
  MinuitOptimizer* myOptimizer = new MinuitOptimizer();

  //Loop over all directions, finding fit FOMs for each
  for(int i=0;i<dirList->size();i++){
    double coneFOM = -1;
    Direction SeedDir = dirList->at(i);
    fSeedVertex->SetDirection(SeedDir);
    myOptimizer->LoadVertex(fSeedVertex);
    // calculate residuals
    Position SeedPos = fSeedVertex->GetPosition();
    this->fVtxGeo->CalcExtendedResiduals(SeedPos.X(), SeedPos.Y(),
            SeedPos.Z(), 0.0, SeedDir.X(), SeedDir.Y(), SeedDir.Z());
    myOptimizer->LoadVertexGeometry(fVtxGeo); 
    //Get this direction's cone FOM
    myOptimizer->ConePropertiesFoM(coneFOM);
    if(coneFOM<0.0) continue;
    //For this direction, find the best vertex time and vtxFOM
    double vtxParam=0;
    double bestVtxFOM = -1.;
    double bestTimeForDir = -999.;
    for(int i=0;i<timeList->size();i++){
      double vtxTime = timeList->at(i);
      double thisVtxFOM = -1.;
      myOptimizer->TimePropertiesLnL(vtxTime,vtxParam,thisVtxFOM);
      if(thisVtxFOM>bestVtxFOM){
        bestTimeForDir = vtxTime;
        bestVtxFOM = thisVtxFOM;
      }
    }
    if(bestVtxFOM<0.0) continue;
    //Finally, see if the two FOMs give the best total FOM
    double totalFOM = (coneweight*coneFOM + vtxweight*bestVtxFOM)/
        (coneweight + vtxweight);
    if(totalFOM > bestCombFOM){
      bestCombDirection = SeedDir;
      bestCombFOM = totalFOM;
      bestCombTime = bestTimeForDir;
    }
  }
  delete myOptimizer; myOptimizer = 0;
  logmessage = " BestCombinedFOM is: " + std::to_string(bestCombFOM);
  Log(logmessage,v_debug,verbosity); 
  logmessage = " BestFittedTime is: " + std::to_string(bestCombTime);
  Log(logmessage,v_debug,verbosity); 
  //Set the best Best total FOM results as the best vertex
  BestVertex->SetVertex(fSeedVertex->GetPosition(),bestCombTime);
  BestVertex->SetDirection(bestCombDirection);
  BestVertex->SetFOM(bestCombFOM,1,1);
  return BestVertex;
}

TMatrixD VtxGridFitter::Rotateatob(TVector3 a, TVector3 b){
  //Get rotation matrix that rotates the z-axis of one coordinate system
  // to vector b.
  //b must be normalized.
  TVector3 au = a.Unit();
  TVector3 bu = b.Unit();
  TVector3 c = au.Cross(bu);
  TVector3 cu = c.Unit();
  double cosa = au.Dot(bu);
  double s = sqrt(1- (cosa*cosa));
  double C = 1-cosa;
  TMatrixD rotn(3,3);
  rotn[0][0]= (cu.X()*cu.X()*C + cosa);
  rotn[0][1]= (cu.X()*cu.Y()*C)-(cu.Z()*s);
  rotn[0][2]= (cu.X()*cu.Z()*C)+(cu.Y()*s);
  rotn[1][0]= (cu.Y()*cu.X()*C)+(cu.Z()*s);
  rotn[1][1]= (cu.Y()*cu.Y()*C)+cosa;
  rotn[1][2]= (cu.Y()*cu.Z()*C)-(cu.X()*s);
  rotn[2][0]= (cu.Z()*cu.X()*C)-(cu.Y()*s);
  rotn[2][1]= (cu.Z()*cu.Y()*C)+(cu.X()*s);
  rotn[2][2]= (cu.Z()*cu.Z()*C)+cosa;
  return rotn;
} 
