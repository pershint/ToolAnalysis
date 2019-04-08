#include "VtxGridFitter.h"

VtxGridFitter::VtxGridFitter():Tool(){}


bool VtxGridFitter::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool VtxGridFitter::Execute(){
  this->Reset();

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
  //TODO: Add these inputs to config file
  NumDirSeeds = 1000;
  double coneweight = 1.0;
  double vtxweight = 1.0;

  double bestFOM = -1;
  this->GenerateRoughTimeSeeds();
  this->GenerateDirectionSeeds(NumDirSeeds);

  //For each position in vSeedVtxList, run the best FOM search
  //On the output from this, check if it has the best FOM of all
  //positions
  RecoVertex* BestVtxCandidate; 
  for(int i=0;i<vSeedVtxList->size();i++){
    BestVtxCandidate = this->FindBestVtxAtPos(&(vSeedVtxList->at(i)),
            vSeedRTimeList, coneweight, vtxweight);
    if(BestVtxCandidate->GetFOM() > bestFOM) BestExtendedVertex = BestVtxCandidate;
  }

  // **** BEGIN SECOND FINE PASS OF REDUCED TANK GRID *****
  // Generate the fine time seeds
  // TODO: Add these into the config
  int numFineSeeds = 100;
  double FineScaleReduction = 0.1;
  this->GenerateFineTimeSeeds(BestVtxCandidate->GetTime(), TimeResolution,
          numFineSeeds);
  // Re-run grid search with position grid shrunk and centered about the
  // Best fit position from the first pass
  for(int i=0;i<vSeedVtxList->size();i++){
    RecoVertex* thisSeedPos = &(vSeedVtxList->at(i));
    scaledSeedPos = this->TransformPosition(thisSeedPos,FineScaleReduction,
            BestVtxCandidate->GetPosition());
    RecoVertex* BestVtxCandidate = this->FindBestVtxAtPos(scaledSeedPos,
            vSeedFTimeList, coneweight, vtxweight);
    if(BestVtxCandidate->GetFOM() > bestFOM) BestExtendedVertex = BestVtxCandidate;
  }
  
  // Push highest FOM vertex into the BoostStore
  this->PushExtendedVertex(BestExtendedVertex,true);

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
  for(double timeStep=MinimumTime;timeStep<=MaximumTime;timeStep+=TimeResolution){
    vSeedRTimeList.push_back(timeStep);
  }
}

void VtxGridFitter::GenerateFineTimeSeeds(double mean, double timeres, int numFineSeeds){
  for(int i=0;i<numFineSeeds;i++){
	double timeStep = frand.Gaus(mean, timeres); // time is smeared with 100 ps time resolution. Harded-coded for now.
    vSeedFTimeList.push_back(timeStep);
  }
}

void VtxGridFitter::GenerateDirectionSeeds(int numDirs){
  vSeedDirList.clear();
  double offset = 2./numDirs;
  double increment = TMath::Pi()*(3. - sqrt(5.));
  for(int i=0;i<numDirs;i++){
    double y = ((i*offset)-1.) + (offset/2.);
    double r = sqrt(1 - pow(y,2));
    double phi = ((i+1) % numDirs) * increment;

    double x = cos(phi) * r;
    double z = sin(phi) * r;

    Direction thisDir(x,y,z);
    vSeedDirList.push_back(thisDir);
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
        std::vector<double> timeList, double coneweight, double vtxweight){
  //Initialize variables holding best fit values
  RecoVertex* BestVertex = new RecoVertex();
  double bestCombTime = -999.;
  Direction bestCombDirection;  
  double bestCombFOM = -1.;

  //Load up our Optimizer, which has access to the FOM calculations
  MinuitOptimizer* myOptimizer = new MinuitOptimizer();
  myOptimizer->LoadVertexGeometry(fVtxGeo); 

  //Loop over all directions, finding fit FOMs for each
  for(int i=0;i<vSeedDirList.size();i++){
    double coneFOM = -1;
    Direction SeedDir = vSeedDirList.at(i);
    fSeedVertex->SetDirection(SeedDir);
    myOptimizer->LoadVertex(fSeedVertex);
    // calculate residuals
    Position SeedPos = fSeedVertex->GetPosition();
    this->fVtxGeo->CalcExtendedResiduals(SeedPos.X(), SeedPos.Y(),
            SeedPos.Z(), 0.0, SeedDir.X(), SeedDir.Y(), SeedDir.Z());
    
    //Get this direction's cone FOM
    myOptimizer->ConePropertiesFoM(coneFOM);
    
    //For this direction, find the best vertex time and vtxFOM
    double vtxParam=0;
    double bestVtxFOM = -1.;
    double bestTimeForDir;
    for(int i=0;i<timeList.size();i++){
      double vtxTime = timeList.at(i);
      double thisVtxFOM = -1.;
      myOptimizer->TimePropertiesLnL(vtxTime,vtxParam,thisVtxFOM);
      if(thisVtxFOM>bestVtxFOM){
        bestTimeForDir = vtxTime;
        bestVtxFOM = thisVtxFOM;
      }
    }
    //Finally, see if the two FOMs give the best total FOM
    double totalFOM = (coneweight*coneFOM + vtxweight*bestVtxFOM)/
        (coneweight + vtxweight);
    if(totalFOM > bestCombFOM){
      bestCombDirection = SeedDir;
      bestCombFOM = totalFOM;
      bestCombTime = bestTimeForDir;
    }
  }
  delete myOptimizer;
  //Set the best Best total FOM results as the best vertex
  BestVertex->SetVertex(fSeedVertex->GetPosition(),bestCombTime);
  BestVertex->SetDirection(bestCombDirection);
  BestVertex->SetFOM(bestCombFOM,1,1);
  return BestVertex;
}
