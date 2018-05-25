//File: DebugStyle.cxx
//Brief: Configuration for drawing basic plots.  Should make the axis labels larger, 
//disable the statistics box, and possibly set a default color palette.   
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOT_STYLES_DEBUGSTYLE_CXX
#define PLOT_STYLES_DEBUGSTYLE_CXX

#include "TStyle.h"
#include "TROOT.h"

namespace util
{
  void SetDebugStyle()
  {
    gStyle->SetOptTitle(1); //Turns on drawing the canvas title?
    gStyle->SetOptStat(111111); //Draw underflows and overflows in histograms
    gStyle->SetOptTitle(1); //Turns on drawing the canvas title?
    gStyle->SetTitleSize(0.04); //Make the canvas title bigger
    gStyle->SetLabelSize(0.04, "X"); //Make x axis labels bigger
    gStyle->SetLabelSize(0.04, "Y"); //Make y axis labels bigger
    gStyle->SetTitleSize(0.04, "X"); //Make x axis titles bigger
    gStyle->SetTitleSize(0.04, "Y"); //Make y axis titles bigger
    gStyle->SetHistLineWidth(2); //Default in ROOT is 1

    gROOT->ForceStyle();
  }
}

#endif //PLOT_STYLES_DEBUGSTYLE_CXX
