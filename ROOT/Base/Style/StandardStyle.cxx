//File: StandardStyle.cxx
//Brief: Configuration for drawing basic plots.  Should make the axis labels larger, 
//disable the statistics box, and possibly set a default color palette.   
//Author: Andrew Olivier aolivier@ur.rochester.edu

#include "TStyle.h"
#include "TROOT.h"

#ifndef PLOT_STYLES_STANDARDSTYLE_CXX
#define PLOT_STYLES_STANDARDSTYLE_CXX

namespace util
{
  void SetStandardStyle()
  {
    //Histograms
    gStyle->SetOptTitle(1); //Turns on drawing the canvas title?
    gStyle->SetTitleSize(0.04); //Make the canvas title bigger
    gStyle->SetOptStat(0); //Disable the statistics box on histograms
    gStyle->SetLabelSize(0.04, "X"); //Make x axis labels bigger
    gStyle->SetLabelSize(0.04, "Y"); //Make y axis labels bigger
    gStyle->SetTitleSize(0.04, "X"); //Make x axis titles bigger
    gStyle->SetTitleSize(0.04, "Y"); //Make y axis titles bigger
    gStyle->SetTitleOffset(1.4, "Y"); //Move the y axis titles farther away from the axis than the default distance
    gStyle->SetHistLineWidth(2); //Default in ROOT is 1

    //Make sure this style is applied
    gROOT->ForceStyle();
  }
}

#endif //PLOT_STYLES_STANDARDSTYLE_CXX
