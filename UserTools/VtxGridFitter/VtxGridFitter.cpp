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

  return true;
}


bool VtxGridFitter::Finalise(){

  return true;
}
