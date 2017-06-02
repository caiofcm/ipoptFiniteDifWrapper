/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     cpp_example.cpp
 Revision: $$
 Contents: example for class myADOLC_NPL for interfacing with Ipopt
 
 Copyright (c) Andrea Walther
   
 This file is part of ADOL-C. This software is provided as open source.
 Any use, reproduction, or distribution of the software constitutes 
 recipient's acceptance of the terms of the accompanying license file.
 
 This code is based on the file corresponding file cpp_example.cpp contained 
 in the Ipopt package with the authors:  Carl Laird, Andreas Waechter   
----------------------------------------------------------------------------*/

#include "optexample1.h"
#include "cfcmINLP.h"

using namespace Ipopt;
using namespace std;

int main(int argv, char* argc[]){
  myExample1 aow;
  OptInterface *optI = &aow;
  SmartPtr<TNLP> nlp = new cfcmINLP(optI);
  createAppAndRun(nlp);
}