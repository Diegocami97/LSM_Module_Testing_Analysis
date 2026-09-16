// Diego Venegas-Vargas
// DAMIC-M Collaboration
// PlotClusterChargeSpectrum.cc
//
// Diagnostic macro: histograms the total charge (in electrons) of every
// cluster in a filelist, with NO energy-window cut. Used to find out which
// low-multiplicity charge values (e.g. ~1, ~4, ~15 e-) actually exist in a
// given dataset before picking targets for MakeExampleClusterFigure.cc.

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"

double E_PER_e_CCS = 3.74; // eV -> electrons conversion

static inline void EnsureDirCCS(const TString& path){
  gSystem->mkdir(path, kTRUE);
}

static inline void SaveCanvasCCS(TCanvas* c, const TString& basePath){
  if (!c) return;
  c->Modified();
  c->Update();
  c->SaveAs(basePath + ".png");
}

void PlotClusterChargeSpectrum(const char* filelist = "test_list.txt",
                                const char* moduleUsed = "DMXX",
                                const char* outDir = "out",
                                double qMax_e = 50.0,
                                double binWidth_e = 0.25)
{
  std::ifstream fin(filelist);
  if (!fin){
    std::cerr << "[ERROR] Cannot open file list: " << filelist << "\n";
    return;
  }
  TChain chain("clustersRec");
  std::string line;
  int nFilesAdded = 0;
  while (std::getline(fin, line)){
    if (line.empty()) continue;
    chain.Add(line.c_str());
    ++nFilesAdded;
  }
  fin.close();
  if (nFilesAdded == 0){
    std::cerr << "[ERROR] No files found in " << filelist << "\n";
    return;
  }

  std::vector<float>* Energy = nullptr;
  chain.SetBranchAddress("Energy", &Energy);

  const int nBins = std::max(1, int(qMax_e / binWidth_e));
  TH1F* hCharge = new TH1F("hClusterCharge_e", "Cluster Charge Spectrum",
                            nBins, 0.0, qMax_e);

  Long64_t nClusters = 0;
  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i){
    chain.GetEntry(i);
    for (size_t n = 0; n < Energy->size(); ++n){
      const double q_e = Energy->at(n) * 1000.0 / E_PER_e_CCS;
      hCharge->Fill(q_e);
      ++nClusters;
    }
  }

  std::cout << "[INFO] Processed " << nEntries << " entries, "
            << nClusters << " clusters total (all energies, no window cut).\n";

  // Print low-charge bin contents so the user can read off populated values.
  std::cout << "[INFO] Charge spectrum below 20 e- (bin center : counts):\n";
  for (int ib = 1; ib <= hCharge->GetNbinsX(); ++ib){
    const double x = hCharge->GetXaxis()->GetBinCenter(ib);
    if (x > 20.0) break;
    const double c = hCharge->GetBinContent(ib);
    if (c > 0) std::cout << "  " << x << " e- : " << c << "\n";
  }

  TString baseTag = Form("%s/Example_Clusters/%s", outDir, moduleUsed);
  EnsureDirCCS(baseTag);

  gStyle->SetOptStat(1);
  TCanvas* c = new TCanvas("cChargeSpectrum", "Cluster Charge Spectrum", 900, 600);
  c->cd();
  gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
  hCharge->SetTitle(";Cluster charge [e-];Counts");
  hCharge->SetLineColor(kBlue+2);
  hCharge->SetFillColorAlpha(kBlue+1, 0.25);
  hCharge->Draw("HIST");
  { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.045);
    lab.DrawLatex(0.50, 0.965, Form("Cluster Charge Spectrum - %s", moduleUsed)); }
  SaveCanvasCCS(c, baseTag + "/charge_spectrum");

  std::cout << "[DONE] Saved: " << baseTag << "/charge_spectrum.png\n";
}
