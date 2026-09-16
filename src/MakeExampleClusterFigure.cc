// Diego Venegas-Vargas
// DAMIC-M Collaboration
// Builds a Fig-6-style figure: real single-cluster pixel maps side by side
// at decreasing electron multiplicity (e.g. ~15, ~4, 1 e-), illustrating how
// spatial cluster shape collapses as multiplicity drops.
//
// Given a list of target charges (in electrons), this scans every cluster in
// a filelist (no energy-window cut) and, for each target, keeps the single
// cluster whose total charge is closest to that target. It then draws each
// winning cluster's real pixel map (cropped to its own bounding box, padded
// by a couple of pixels for context) in a shared 1xN canvas.
//
// Run PlotClusterChargeSpectrum.cc first to see which charge values are
// actually populated in your data before choosing targetsCSV.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"

#ifdef __CLING__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#endif

double E_PER_e_MEC = 3.74; // eV -> electrons conversion

static inline void EnsureDirMEC(const TString& path){
  gSystem->mkdir(path, kTRUE);
}

static inline void SaveCanvasMEC(TCanvas* c, const TString& basePath){
  if (!c) return;
  c->Modified();
  c->Update();
  c->SaveAs(basePath + ".png");
}

static std::vector<double> ParseTargetsCSV(const char* csv){
  std::vector<double> out;
  std::stringstream ss(csv);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) continue;
    out.push_back(std::stod(tok));
  }
  return out;
}

// One selected cluster's data, deep-copied out of the TChain buffers
// (buffers get overwritten on the next GetEntry, so this must be a copy).
struct ClusterMatch {
  bool found = false;
  double target_e = 0.0;
  double q_e = 0.0;
  double diff_e = std::numeric_limits<double>::infinity();
  Long64_t entry = -1;
  size_t clusterIdx = 0;
  double posX = 0.0, posY = 0.0;
  int npix = 0;
  std::vector<int> px, py;
  std::vector<float> pE;
};

// Build a TH2F of one cluster's pixel map, cropped to its own bounding box
// padded by padPix pixels on each side, with bin content in electrons/pixel.
static TH2F* ClusterToTH2F_e(const ClusterMatch& m, int padPix, const char* name){
  if (m.px.empty()){
    return new TH2F(name, name, 1, -0.5, 0.5, 1, -0.5, 0.5);
  }
  int min_x = *std::min_element(m.px.begin(), m.px.end()) - padPix;
  int max_x = *std::max_element(m.px.begin(), m.px.end()) + padPix;
  int min_y = *std::min_element(m.py.begin(), m.py.end()) - padPix;
  int max_y = *std::max_element(m.py.begin(), m.py.end()) + padPix;

  TH2F* h2 = new TH2F(name, name,
                       max_x - min_x + 1, min_x - 0.5, max_x + 0.5,
                       max_y - min_y + 1, min_y - 0.5, max_y + 0.5);
  for (size_t i = 0; i < m.px.size(); ++i){
    const double pE_e = m.pE[i] * 1000.0 / E_PER_e_MEC;
    h2->Fill(m.px[i], m.py[i], pE_e);
  }
  return h2;
}

void MakeExampleClusterFigure(const char* filelist = "test_list.txt",
                               const char* targetsCSV = "15,4,1",
                               const char* moduleUsed = "DMXX",
                               const char* outDir = "out",
                               int padPix = 2,
                               double zMax_e = 5.0)
{
  std::vector<double> targets = ParseTargetsCSV(targetsCSV);
  if (targets.empty()){
    std::cerr << "[ERROR] No targets parsed from targetsCSV='" << targetsCSV << "'\n";
    return;
  }

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

  std::vector<std::vector<int>>* pixels_x = nullptr;
  std::vector<std::vector<int>>* pixels_y = nullptr;
  std::vector<std::vector<float>>* pixels_E = nullptr;
  std::vector<float>* PosX = nullptr;
  std::vector<float>* PosY = nullptr;
  std::vector<float>* Energy = nullptr;
  std::vector<int>* Npix = nullptr;

  chain.SetBranchAddress("pixels_x", &pixels_x);
  chain.SetBranchAddress("pixels_y", &pixels_y);
  chain.SetBranchAddress("pixels_E", &pixels_E);
  chain.SetBranchAddress("PosX", &PosX);
  chain.SetBranchAddress("PosY", &PosY);
  chain.SetBranchAddress("Energy", &Energy);
  chain.SetBranchAddress("Npix", &Npix);

  std::vector<ClusterMatch> matches(targets.size());
  for (size_t t = 0; t < targets.size(); ++t) matches[t].target_e = targets[t];

  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i){
    chain.GetEntry(i);
    const size_t nC = Energy->size();
    for (size_t n = 0; n < nC; ++n){
      const double q_e = Energy->at(n) * 1000.0 / E_PER_e_MEC;
      for (size_t t = 0; t < targets.size(); ++t){
        const double diff = std::fabs(q_e - targets[t]);
        if (diff < matches[t].diff_e){
          ClusterMatch& m = matches[t];
          m.found = true;
          m.q_e = q_e;
          m.diff_e = diff;
          m.entry = i;
          m.clusterIdx = n;
          m.posX = PosX->at(n);
          m.posY = PosY->at(n);
          m.npix = Npix->at(n);
          m.px = pixels_x->at(n);
          m.py = pixels_y->at(n);
          m.pE = pixels_E->at(n);
        }
      }
    }
  }

  for (size_t t = 0; t < targets.size(); ++t){
    const ClusterMatch& m = matches[t];
    if (!m.found){
      std::cerr << "[WARN] No cluster found at all for target " << targets[t] << " e-\n";
      continue;
    }
    std::cout << "[MATCH] target=" << targets[t] << " e- -> found q=" << m.q_e
               << " e- (diff=" << m.diff_e << ") at entry=" << m.entry
               << " cluster=" << m.clusterIdx
               << " PosX=" << m.posX << " PosY=" << m.posY
               << " Npix=" << m.npix << "\n";
    if (m.diff_e > 5.0){
      std::cerr << "[WARN] Closest match for target " << targets[t]
                << " e- is " << m.diff_e << " e- away - may not be representative.\n";
    }
  }

  TString baseTag = Form("%s/Example_Clusters/%s", outDir, moduleUsed);
  EnsureDirMEC(baseTag);

  gStyle->SetOptStat(0);
  const int N = (int)targets.size();
  TCanvas* c = new TCanvas("cExampleClusters", "Example Cluster Multiplicity", 500 * N, 550);
  c->Divide(N, 1);

  for (int t = 0; t < N; ++t){
    c->cd(t + 1);
    gPad->SetRightMargin(0.18); gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.18); gPad->SetBottomMargin(0.12);

    const ClusterMatch& m = matches[t];
    TH2F* h2 = ClusterToTH2F_e(m, padPix, Form("hExampleCluster_%d", t));
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("x [pix]");
    h2->GetYaxis()->SetTitle("y [pix]");
    h2->SetMinimum(0.0);
    h2->SetMaximum(zMax_e);
    h2->GetZaxis()->SetTitle("e-/pixel");
    h2->Draw("COLZ");

    TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.08);
    lab.DrawLatex(0.5, 0.94, m.found ? Form("~%.1f e-", m.q_e) : "no match");
  }

  SaveCanvasMEC(c, baseTag + "/example_cluster_multiplicity");
  std::cout << "[DONE] Saved: " << baseTag << "/example_cluster_multiplicity.png\n";
}
