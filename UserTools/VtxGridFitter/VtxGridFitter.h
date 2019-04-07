#ifndef VtxGridFitter_H
#define VtxGridFitter_H

#include <string>
#include <iostream>

#include "Tool.h"

class VtxGridFitter: public Tool {


 public:

  VtxGridFitter();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:





};


#endif
