#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TMath.h"
void run_gainscan(const char* list="input_lists/test_gain_list.txt",
                  const char* out_csv="gains_long.csv",
                  const char* outdir="plots") {
  gSystem->AddIncludePath("-I/opt/homebrew/opt/cfitsio/include");
  gSystem->AddLinkedLibs("-L/opt/homebrew/opt/cfitsio/lib -lcfitsio");
  gROOT->ProcessLine(Form(".x src/gainscan.cc++(\"%s\",\"%s\",\"%s\")", list, out_csv, outdir));
}