// ParameterScanEnergyPeaks.C
// Reads a file list where each line is a ROOT file named ..._<image_id>_ext<ext>.root.
// Files with the same grid point (same prefix before _<id>_ext<ext>.root) and same ext
// are combined, then front/back energy spectra are filled and fitted (double Gaussian)
// like MakeComposite_Fe55_Clusters_Extensions.cc. Outputs CSV and comparison plots.
//
// Usage:
//   root -l 'ParameterScanEnergyPeaks.C("input_lists/Image_5_optimal.txt", 5, "Parameter_Scan_Fe")'
//
// Arguments: filelist, image_number, output_dir
// If the list contains paths from multiple modules (e.g. .../DM01_VR/... and .../DM02_VR/...),
// module is inferred from the parent directory of each file; results are grouped and output per module.

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#ifdef __CLING__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<int>>+;
#pragma link C++ class std::vector<std::vector<float>>+;
#endif

// --------------- Config (match MakeComposite_Fe55_Clusters_Extensions.cc) ---------------
const double ENE_MIN       = 4.5;
const double ENE_MAX       = 7.0;
const double SIGMA_FRONT   = 1.0;
const double G1_LO = 5.3, G1_HI = 6.1;

// --------------- Plot Y-axis ranges (manual; edit to fix each panel) ---------------
// Energy plot: 2 panels per canvas — [0]=top (Kβ), [1]=bottom (Kα)
const double energy_ymin[2] = { 6.2,  5.6 };
const double energy_ymax[2] = { 6.8,  6.2 };
// Calibration plot: 3 panels — [0]=Gain, [1]=Single-e σ, [2]=Dark current
const double cal_ymin[3]    = { 1.2,   0.12,   0.4 };
const double cal_ymax[3]    = { 1.5,  0.25,  6 };

// --------------- Parse basename: ..._<grid_point>_ext<1-4>.root -> (grid_point, ext_num) ---------------
// Grid point = full prefix before _ext (e.g. ..._test_1, ..._test_2 are different grid points).
bool ParseRootName(const TString& fullpath, TString& grid_point, int& ext_num) {
  TString base = gSystem->BaseName(fullpath.Data());
  int iext = base.Index("_ext");
  if (iext < 0 || iext + 5 > base.Length()) return false;
  char c = base[iext+4];
  if (c < '1' || c > '4') return false;
  ext_num = c - '0';
  grid_point = base(0, iext);  // full prefix before _ext: test_1, test_2, test_3, test_4 -> different grid points
  return true;
}

// --------------- Infer module from path: parent directory of file (e.g. DM02_VR) ---------------
TString GetModuleFromPath(const TString& fullpath) {
  TString dir = gSystem->DirName(fullpath.Data());
  return gSystem->BaseName(dir);
}

// --------------- Group file list: (module, grid_point, ext) -> list of paths ---------------
// Key format: "module|grid_point_eN" so files from different modules are not mixed.
std::map<std::string, std::vector<TString>> GroupFileList(const char* filelist) {
  std::ifstream f(filelist);
  std::map<std::string, std::vector<TString>> groups;
  std::string line_str;
  while (std::getline(f, line_str)) {
    TString line = line_str.c_str();
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0 || line[0] == '#') continue;
    if (!line.EndsWith(".root")) continue;
    TString grid;
    int ext = 0;
    if (!ParseRootName(line, grid, ext)) continue;
    TString module = GetModuleFromPath(line);
    if (module.Length() == 0) module = "unknown";
    std::string key = std::string(module.Data()) + "|" + std::string(grid.Data()) + "_e" + std::to_string(ext);
    groups[key].push_back(line);
  }
  return groups;
}

// --------------- Fill front/back energy histos from chain (same logic as .cc) ---------------
void FillFrontBackHists(TChain* chain, int image_number,
                       TH1F* hFront, TH1F* hBack) {
  std::vector<float> *PosX = nullptr, *PosY = nullptr;
  std::vector<float> *rmsx = nullptr, *rmsy = nullptr;
  std::vector<float> *Energy = nullptr;

  chain->SetBranchAddress("Energy", &Energy);
  chain->SetBranchAddress("fSTD_X", &rmsx);
  chain->SetBranchAddress("fSTD_Y", &rmsy);
  chain->SetBranchAddress("PosX",   &PosX);

  const Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    chain->GetEntry(i);
    const size_t nC = (PosX ? PosX->size() : (Energy ? Energy->size() : 0));
    if (nC == 0) continue;
    for (size_t n = 0; n < nC; n++) {
      double E = Energy->at(n);
      if (E <= ENE_MIN || E >= ENE_MAX) continue;
      double sx = rmsx->at(n), sy = rmsy->at(n);
      double sxy = std::sqrt(0.5 * (sx*sx + sy*sy));
      double sCut = sxy;
      if (image_number == 6 || image_number == 9) sCut = sx;
      if (image_number == 7 || image_number == 10) sCut = sy;
      if (sCut < SIGMA_FRONT)
        hFront->Fill(E);
      else
        hBack->Fill(E);
    }
  }
}

// --------------- Fit double Gaussian (first peak [G1_LO,G1_HI], second dynamic); return mu1, mu2 ---------------
void FitDoubleGauss(TH1F* h, double& mu1, double& mu1_err, double& mu2, double& mu2_err) {
  mu1 = 5.9; mu1_err = 0; mu2 = 6.4; mu2_err = 0;
  if (!h || h->GetEntries() < 10) return;
  TF1* g1 = new TF1("g1", "gaus", G1_LO, G1_HI);
  g1->SetParameters(h->GetMaximum(), 5.9, 0.12);
  g1->SetParLimits(2, 0.03, 0.5);
  int st = h->Fit(g1, "RQ0S");
  if (st != 0) { delete g1; return; }
  double m1 = g1->GetParameter(1), s1 = g1->GetParameter(2);
  mu1 = m1; mu1_err = g1->GetParError(1);
  double g2_left = TMath::Max(ENE_MIN, m1 + 2.*s1);
  double g2_right = TMath::Min(ENE_MAX, g2_left + 1.0);
  if (g2_right <= g2_left) g2_right = TMath::Min(ENE_MAX, g2_left + 0.5);
  TF1* g2 = new TF1("g2", "gaus", g2_left, g2_right);
  g2->SetParameters(h->GetMaximum()/3., (g2_left+g2_right)/2., 0.12);
  g2->SetParLimits(2, 0.03, 0.5);
  st = h->Fit(g2, "RQ0S");
  if (st == 0) { mu2 = g2->GetParameter(1); mu2_err = g2->GetParError(1); }
  delete g2;
  delete g1;
}

// --------------- One result row ---------------
struct PeakRow {
  TString module;
  TString grid_point;
  int ext;
  double mu_KA_f, mu_KB_f, mu_KA_b, mu_KB_b;
  double mu_KA_f_err, mu_KB_f_err, mu_KA_b_err, mu_KB_b_err;
  // Info tree (same as MakeComposite_Fe55_Clusters_Extensions.cc)
  double gainext, gainext_err, sigext, sigext_err, dcext, dcext_err;
};

// --------------- Read gain, single-e resolution, DC from info tree (first file): average over all entries ---------------
// Gain: mean over entries; error on mean = sqrt(sum(ge_i^2))/N.
// Sigma: convert to electrons (sig / gain) per entry, then average; output in e-.
// DC: convert to e/pix/day via DC * (1/gain) * (86400/11208) per entry, then average.
const double DC_SEC_TO_DAY = 86400.0 / 11208.0;  // factor for e/pix/s -> e/pix/day
void ReadInfoFromFirstFile(const char* first_file, int ext_num,
                           double& gain, double& gain_err, double& sig, double& sig_err,
                           double& dc, double& dc_err) {
  gain = 0; gain_err = 0; sig = 0; sig_err = 0; dc = 0; dc_err = 0;
  TFile* finfo = TFile::Open(first_file);
  if (!finfo || finfo->IsZombie()) { if (finfo) delete finfo; return; }
  TTree* tinfo = (TTree*)finfo->Get("info");
  if (!tinfo || tinfo->GetEntries() == 0) { finfo->Close(); delete finfo; return; }
  double gv = 0, ge = 0, sv = 0, se = 0, dv = 0, de = 0;
  const char* gain_br     = Form("ext%d_gainext",     ext_num);
  const char* gainerr_br  = Form("ext%d_gainext_err", ext_num);
  const char* sig_br      = Form("ext%d_sigext",      ext_num);
  const char* sigerr_br   = Form("ext%d_sigext_err", ext_num);
  const char* dc_br       = Form("ext%d_dcext",      ext_num);
  const char* dcerr_br    = Form("ext%d_dcext_err",  ext_num);
  if (tinfo->FindBranch(gain_br))    tinfo->SetBranchAddress(gain_br,    &gv);
  if (tinfo->FindBranch(gainerr_br)) tinfo->SetBranchAddress(gainerr_br, &ge);
  if (tinfo->FindBranch(sig_br))     tinfo->SetBranchAddress(sig_br,     &sv);
  if (tinfo->FindBranch(sigerr_br))  tinfo->SetBranchAddress(sigerr_br,  &se);
  if (tinfo->FindBranch(dc_br))     tinfo->SetBranchAddress(dc_br,      &dv);
  if (tinfo->FindBranch(dcerr_br))  tinfo->SetBranchAddress(dcerr_br,   &de);
  const Long64_t nEntries = tinfo->GetEntries();
  double sum_g = 0, sum_ge2 = 0;
  double sum_sige = 0, sum_sige_err2 = 0;   // sigma in electrons (sig * gain) per entry, then average
  double sum_dc_day = 0, sum_dc_day_err2 = 0;  // DC in e/pix/day per entry, then average
  for (Long64_t i = 0; i < nEntries; i++) {
    tinfo->GetEntry(i);
    sum_g += gv;
    sum_ge2 += ge * ge;
    // Sigma in electrons: sig_e = sig / gain; error propagation (relative): (err_sige/sig_e)^2 = (se/sv)^2 + (ge/gv)^2
    if (gv > 0 && sv >= 0) {
      double sig_e = sv / gv;
      sum_sige += sig_e;
      double rel2 = (se > 0 && sv > 0 ? (se / sv) * (se / sv) : 0) + (ge / gv) * (ge / gv);
      sum_sige_err2 += rel2 * sig_e * sig_e;
    }
    // DC in e/pix/day: dc_day = dc * (1/gain) * (86400/11208); error propagation (relative): (de/dv)^2 + (ge/gv)^2
    if (gv > 0 && dv >= 0) {
      double dc_day = dv * DC_SEC_TO_DAY / gv;
      sum_dc_day += dc_day;
      double rel2 = (de > 0 && dv > 0 ? (de / dv) * (de / dv) : 0) + (ge / gv) * (ge / gv);
      sum_dc_day_err2 += rel2 * dc_day * dc_day;
    }
  }
  gain = sum_g / nEntries;
  gain_err = (nEntries > 0 && sum_ge2 >= 0) ? std::sqrt(sum_ge2) / nEntries : 0;
  sig  = (nEntries > 0) ? sum_sige / nEntries : 0;
  sig_err = (nEntries > 0 && sum_sige_err2 >= 0) ? std::sqrt(sum_sige_err2) / nEntries : 0;
  dc   = (nEntries > 0) ? sum_dc_day / nEntries : 0;
  dc_err = (nEntries > 0 && sum_dc_day_err2 >= 0) ? std::sqrt(sum_dc_day_err2) / nEntries : 0;
  finfo->Close();
  delete finfo;
}

// --------------- Main macro ---------------
void ParameterScanEnergyPeaks(const char* filelist = "input_lists/Image_5_optimal.txt",
                              const int image_number = 5,
                              const char* output_dir = "Parameter_Scan_Fe") {
  gSystem->mkdir(output_dir, true);
  auto groups = GroupFileList(filelist);
  if (groups.empty()) {
    std::cerr << "No valid files in " << filelist << std::endl;
    return;
  }

  // Per-module sorted grid points (key is "module|grid_point_eN")
  std::map<std::string, std::set<std::string>> module_gp_set;
  for (const auto& kv : groups) {
    const std::string& key = kv.first;
    size_t pipe = key.find('|');
    if (pipe == std::string::npos) continue;
    size_t ep = key.find("_e", pipe);
    if (ep == std::string::npos) continue;
    std::string module = key.substr(0, pipe);
    std::string grid_point = key.substr(pipe+1, ep - (pipe+1));
    module_gp_set[module].insert(grid_point);
  }
  std::map<std::string, std::vector<std::string>> module_grid_points;
  for (auto& kv : module_gp_set) {
    module_grid_points[kv.first].assign(kv.second.begin(), kv.second.end());
    std::sort(module_grid_points[kv.first].begin(), module_grid_points[kv.first].end());
  }

  std::vector<PeakRow> rows;
  const int nBins = 50;

  for (const auto& kv : groups) {
    const std::string& key = kv.first;
    const std::vector<TString>& paths = kv.second;
    size_t pipe = key.find('|');
    if (pipe == std::string::npos) continue;
    size_t ep = key.find("_e", pipe);
    if (ep == std::string::npos) continue;
    std::string module = key.substr(0, pipe);
    std::string grid_point = key.substr(pipe+1, ep - (pipe+1));
    int ext_num = std::stoi(key.substr(ep+2));
    TChain* chain = new TChain("clustersRec");
    for (const auto& p : paths) chain->Add(p.Data());
    if (chain->GetEntries() == 0) { delete chain; continue; }
    TH1F* hFront = new TH1F("hF", "", nBins, ENE_MIN, ENE_MAX);
    TH1F* hBack  = new TH1F("hB", "", nBins, ENE_MIN, ENE_MAX);
    FillFrontBackHists(chain, image_number, hFront, hBack);
    delete chain;
    PeakRow r;
    r.module = module.c_str();
    r.grid_point = grid_point.c_str();
    r.ext = ext_num;
    r.gainext = 0; r.gainext_err = 0; r.sigext = 0; r.sigext_err = 0; r.dcext = 0; r.dcext_err = 0;
    ReadInfoFromFirstFile(paths[0].Data(), ext_num, r.gainext, r.gainext_err, r.sigext, r.sigext_err, r.dcext, r.dcext_err);
    FitDoubleGauss(hFront, r.mu_KA_f, r.mu_KA_f_err, r.mu_KB_f, r.mu_KB_f_err);
    FitDoubleGauss(hBack,  r.mu_KA_b, r.mu_KA_b_err, r.mu_KB_b, r.mu_KB_b_err);
    rows.push_back(r);
    delete hFront;
    delete hBack;
    std::cout << "  " << module << " " << grid_point << " ext" << ext_num
              << " Kα_f=" << r.mu_KA_f << " Kβ_f=" << r.mu_KB_f
              << " Kα_b=" << r.mu_KA_b << " Kβ_b=" << r.mu_KB_b << std::endl;
  }

  // CSV: one row per (scan_id, module, extension); scan_id is per-module; include calibration from info tree
  TString csv_path = TString(output_dir) + "/energy_peaks_scan.csv";
  std::ofstream csv(csv_path.Data());
  csv << "scan_id,module,extension,mu_KA_front,mu_KB_front,mu_KA_back,mu_KB_back,gain,gain_err,sigext,sigext_err,dcext,dcext_err\n";
  for (const auto& mg : module_grid_points) {
    const std::string& module = mg.first;
    const std::vector<std::string>& grid_points = mg.second;
    for (size_t gi = 0; gi < grid_points.size(); gi++) {
      TString sid = TString::Format("S%d", (int)gi+1);
      for (int e = 1; e <= 4; e++) {
        TString ext_name = TString::Format("Ext%d", e);
        double ka_f = 0, kb_f = 0, ka_b = 0, kb_b = 0;
        double gain = 0, gain_err = 0, sig = 0, sig_err = 0, dc = 0, dc_err = 0;
        for (const auto& r : rows) {
          if (r.module != module.c_str() || r.grid_point != grid_points[gi].c_str() || r.ext != e) continue;
          ka_f = r.mu_KA_f; kb_f = r.mu_KB_f; ka_b = r.mu_KA_b; kb_b = r.mu_KB_b;
          gain = r.gainext; gain_err = r.gainext_err; sig = r.sigext; sig_err = r.sigext_err;
          dc = r.dcext; dc_err = r.dcext_err;
          break;
        }
        csv << sid.Data() << "," << module.c_str() << "," << ext_name.Data()
            << "," << ka_f << "," << kb_f << "," << ka_b << "," << kb_b
            << "," << gain << "," << gain_err << "," << sig << "," << sig_err << "," << dc << "," << dc_err << "\n";
      }
    }
  }
  csv.close();
  std::cout << "Wrote " << csv_path.Data() << std::endl;

  // Comparison plots like parameter_scan_comparison.py: one figure (Front, Back) with
  // x = (Module, Extension), two panels (K#beta top, K#alpha bottom), one series per scan with jitter
  TString plot_base = TString(output_dir) + "/scan_plots";
  gSystem->mkdir(plot_base.Data(), true);
  const double E_KA = 5.898, E_KB = 6.490;
  const int nExt = 4;
  const int gap = 1;
  int colors[] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+1, kYellow+2};
  int markers[] = {20, 21, 22, 23, 24, 25, 26, 27};

  // Ordered modules and minimum scan count (so we have same number of points per module)
  std::vector<std::string> modules;
  int nScansMin = 999;
  for (const auto& mg : module_grid_points) {
    modules.push_back(mg.first);
    int n = (int)mg.second.size();
    if (n < nScansMin) nScansMin = n;
  }
  std::sort(modules.begin(), modules.end());
  if (nScansMin == 999) nScansMin = 0;
  int nModules = (int)modules.size();
  int nX = nModules * nExt;
  if (nX == 0 || nScansMin == 0) { std::cout << "Nothing to plot." << std::endl; return; }

  // X positions: for (module_i, ext_j), x = mi*(nExt+gap)+j  -> 0,1,2,3, 5,6,7,8 for 2 modules
  double xmax_plot = (nModules - 1) * (nExt + gap) + nExt - 0.5;
  std::vector<double> x_edges(nX + 1);
  x_edges[0] = -0.5;
  for (int b = 1; b <= nX; b++) {
    int mi = (b - 1) / nExt;
    int ei = (b - 1) % nExt;
    x_edges[b] = mi * (nExt + gap) + ei + 0.5;
  }

  auto get_value = [&](const std::string& mod, const std::string& gp, int ext, bool use_front, bool is_ka) -> double {
    for (const auto& r : rows) {
      if (r.module != mod.c_str() || r.grid_point != gp.c_str() || r.ext != ext) continue;
      if (use_front) return is_ka ? r.mu_KA_f : r.mu_KB_f;
      return is_ka ? r.mu_KA_b : r.mu_KB_b;
    }
    return 0.;
  };

  // Module ID for x-axis: use part before "_" (e.g. DM01_VR -> DM01)
  auto module_id = [](const std::string& m) -> TString {
    TString s = m.c_str();
    int idx = s.Index("_");
    return (idx > 0) ? s(0, idx) : s;
  };

  auto draw_scan_summary = [&](const char* title, bool use_front, const char* fname) {
    TCanvas* c = new TCanvas("c", title, 1200, 720);
    c->Divide(1, 2, 0, 0);
    // Top panel: K#beta (no x-axis labels; tight margin below)
    c->cd(1);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.05);
    gPad->SetGrid();
    TH1F* frame_top = new TH1F("frame_top", ";;#mu [keV]", nX, &x_edges[0]);
    frame_top->SetStats(0);
    frame_top->SetMinimum(energy_ymin[0]);
    frame_top->SetMaximum(energy_ymax[0]);
    frame_top->GetYaxis()->CenterTitle();
    frame_top->GetYaxis()->SetTitle("Mean Peak Energy [keV]");
    frame_top->GetYaxis()->SetTitleSize(0.055);
    frame_top->GetYaxis()->SetLabelSize(0.045);
    frame_top->GetXaxis()->SetLabelSize(0);
    frame_top->GetXaxis()->SetTitleSize(0);
    frame_top->GetXaxis()->SetNdivisions(0);
    frame_top->Draw("AXIS");
    TLine* l_kb = new TLine(-0.5, E_KB, xmax_plot, E_KB); l_kb->SetLineStyle(2); l_kb->SetLineColor(kMagenta+1); l_kb->Draw();
    TLatex* lab_kb = new TLatex(0.01, 0.92, "K#beta"); lab_kb->SetNDC(); lab_kb->SetTextSize(0.065); lab_kb->Draw();
    for (int mi = 1; mi < nModules; mi++) {
      double sep = mi * (nExt + gap) - gap / 2.0 - 0.5;
      TLine* v = new TLine(sep, energy_ymin[0], sep, energy_ymax[0]); v->SetLineStyle(3); v->SetLineWidth(1); v->Draw();
    }
    TLegend* leg = new TLegend(0.48, 0.68, 0.96, 0.92);
    leg->SetNColumns(TMath::Min(nScansMin, 8));
    leg->SetTextSize(0.042);
    leg->SetTextFont(42);
    leg->SetFillStyle(1001);
    leg->SetBorderSize(1);
    for (int si = 0; si < nScansMin; si++) {
      double jitter = (nScansMin > 1) ? (-0.15 + 0.30 * (double)si / (nScansMin - 1)) : 0.;
      std::vector<double> vx, vy;
      for (int mi = 0; mi < nModules; mi++) {
        const std::string& mod = modules[mi];
        const std::vector<std::string>& gps = module_grid_points.at(mod);
        if (si >= (int)gps.size()) continue;
        const std::string& gp = gps[si];
        for (int e = 0; e < nExt; e++) {
          double x = mi * (nExt + gap) + e + jitter;
          double y = get_value(mod, gp, e + 1, use_front, false);
          vx.push_back(x); vy.push_back(y);
        }
      }
      if (vx.empty()) continue;
      TGraph* g = new TGraph((int)vx.size(), &vx[0], &vy[0]);
      g->SetMarkerStyle(markers[si % 8]);
      g->SetMarkerColor(colors[si % 8]);
      g->SetMarkerSize(1.1);
      g->Draw("P SAME");
      leg->AddEntry(g, TString::Format("S%d", si + 1), "p");
    }
    // leg->Draw();
    c->cd(2);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.14);
    gPad->SetGrid();
    TH1F* frame_bot = new TH1F("frame_bot", ";Module, Extension; ", nX, &x_edges[0]);
    frame_bot->SetStats(0);
    frame_bot->SetMinimum(energy_ymin[1]);
    frame_bot->SetMaximum(energy_ymax[1]);
    frame_bot->GetYaxis()->CenterTitle();
    frame_bot->GetYaxis()->SetTitle("Mean Peak Energy [keV]");
    frame_bot->GetYaxis()->SetTitleSize(0.055);
    frame_bot->GetYaxis()->SetLabelSize(0.045);
    frame_bot->GetXaxis()->SetTitleSize(0.055);
    frame_bot->GetXaxis()->SetLabelSize(0.045);
    for (int i = 1; i <= nX; i++) {
      int mi = (i - 1) / nExt;
      int e = (i - 1) % nExt;
      frame_bot->GetXaxis()->SetBinLabel(i, TString::Format("%s-EXT%d", module_id(modules[mi]).Data(), e + 1));
    }
    frame_bot->Draw("AXIS");
    TLine* l_ka = new TLine(-0.5, E_KA, xmax_plot, E_KA); l_ka->SetLineStyle(2); l_ka->SetLineColor(kMagenta+1); l_ka->Draw();
    TLatex* lab_ka = new TLatex(0.01, 0.92, "K#alpha"); lab_ka->SetNDC(); lab_ka->SetTextSize(0.065); lab_ka->Draw();
    for (int mi = 1; mi < nModules; mi++) {
      double sep = mi * (nExt + gap) - gap / 2.0 - 0.5;
      TLine* v = new TLine(sep, energy_ymin[1], sep, energy_ymax[1]); v->SetLineStyle(3); v->SetLineWidth(1); v->Draw();
    }
    for (int si = 0; si < nScansMin; si++) {
      double jitter = (nScansMin > 1) ? (-0.15 + 0.30 * (double)si / (nScansMin - 1)) : 0.;
      std::vector<double> vx, vy;
      for (int mi = 0; mi < nModules; mi++) {
        const std::string& mod = modules[mi];
        const std::vector<std::string>& gps = module_grid_points.at(mod);
        if (si >= (int)gps.size()) continue;
        const std::string& gp = gps[si];
        for (int e = 0; e < nExt; e++) {
          double x = mi * (nExt + gap) + e + jitter;
          double y = get_value(mod, gp, e + 1, use_front, true);
          vx.push_back(x); vy.push_back(y);
        }
      }
      if (vx.empty()) continue;
      TGraph* g = new TGraph((int)vx.size(), &vx[0], &vy[0]);
      g->SetMarkerStyle(markers[si % 8]);
      g->SetMarkerColor(colors[si % 8]);
      g->SetMarkerSize(1.1);
      g->Draw("P SAME");
    }
    c->cd(0);
    TLatex* tit = new TLatex(0.5, 0.98, title); tit->SetNDC(); tit->SetTextAlign(22); tit->SetTextSize(0.038); tit->Draw();
    c->SaveAs(plot_base + "/" + fname);
    delete c;
  };

  draw_scan_summary("Fe-55  Front Events", true,  "scan_energy_peaks_front.png");
  draw_scan_summary("Fe-55 Back Events",  false, "scan_energy_peaks_back.png");

  // Calibration (info tree): Gain, single-e resolution, DC — same comparison format (Module × Extension, one series per scan)
  auto get_cal = [&](const std::string& mod, const std::string& gp, int ext, int which) -> double {
    for (const auto& r : rows) {
      if (r.module != mod.c_str() || r.grid_point != gp.c_str() || r.ext != ext) continue;
      if (which == 0) return r.gainext;
      if (which == 1) return r.sigext;
      return r.dcext;
    }
    return 0.;
  };
  const char* cal_titles[] = { "Gain [ADU/e^{-}]", "Single-e Resolution [e^{-}]", "Dark current [e^{-}/pix/Day]" };
  TCanvas* cCal = new TCanvas("cCal", "Calibration (info tree)", 1200, 960);
  cCal->Divide(1, 3, 0, 0);
  for (int pad = 0; pad < 3; pad++) {
    cCal->cd(pad + 1);
    gPad->SetTopMargin(pad == 0 ? 0.05 : 0.03);
    gPad->SetBottomMargin(pad == 2 ? 0.14 : 0.02);
    gPad->SetGrid();
    TH1F* frame_cal = new TH1F(Form("frame_cal_%d", pad), ";Module, Extension;", nX, &x_edges[0]);
    frame_cal->SetStats(0);
    frame_cal->SetTitle(TString::Format(";Module, Extension;%s", cal_titles[pad]));
    frame_cal->SetMinimum(cal_ymin[pad]);
    frame_cal->SetMaximum(cal_ymax[pad] > cal_ymin[pad] ? cal_ymax[pad] : cal_ymin[pad] + 1.);
    frame_cal->GetYaxis()->SetTitleSize(0.055);
    frame_cal->GetYaxis()->SetLabelSize(0.045);
    frame_cal->GetXaxis()->SetTitleSize(0.055);
    frame_cal->GetXaxis()->SetLabelSize(0.045);
    frame_cal->Draw("AXIS");
    for (int i = 1; i <= nX; i++) {
      int mi = (i - 1) / nExt;
      int e = (i - 1) % nExt;
      frame_cal->GetXaxis()->SetBinLabel(i, TString::Format("%s-EXT%d", module_id(modules[mi]).Data(), e + 1));
    }
    for (int mi = 1; mi < nModules; mi++) {
      double sep = mi * (nExt + gap) - gap / 2.0 - 0.5;
      TLine* v = new TLine(sep, frame_cal->GetMinimum(), sep, frame_cal->GetMaximum());
      v->SetLineStyle(3); v->SetLineWidth(1); v->Draw();
    }
    TLegend* leg_cal = (pad == 0) ? new TLegend(0.48, 0.68, 0.96, 0.92) : nullptr;
    if (leg_cal) {
      leg_cal->SetNColumns(TMath::Min(nScansMin, 8));
      leg_cal->SetTextSize(0.042);
      leg_cal->SetTextFont(42);
      leg_cal->SetFillStyle(1001);
      leg_cal->SetBorderSize(1);
    }
    for (int si = 0; si < nScansMin; si++) {
      double jitter = (nScansMin > 1) ? (-0.15 + 0.30 * (double)si / (nScansMin - 1)) : 0.;
      std::vector<double> vx, vy;
      for (int mi = 0; mi < nModules; mi++) {
        const std::string& mod = modules[mi];
        const std::vector<std::string>& gps = module_grid_points.at(mod);
        if (si >= (int)gps.size()) continue;
        for (int e = 0; e < nExt; e++) {
          double x = mi * (nExt + gap) + e + jitter;
          double y = get_cal(mod, gps[si], e + 1, pad);
          vx.push_back(x); vy.push_back(y);
        }
      }
      if (vx.empty()) continue;
      TGraph* g = new TGraph((int)vx.size(), &vx[0], &vy[0]);
      g->SetMarkerStyle(markers[si % 8]);
      g->SetMarkerColor(colors[si % 8]);
      g->SetMarkerSize(1.1);
      g->Draw("P SAME");
      if (leg_cal) leg_cal->AddEntry(g, TString::Format("S%d", si + 1), "p");
    }
    // if (leg_cal) leg_cal->Draw();
  }
  cCal->cd(0);
  TLatex* titCal = new TLatex(0.5, 0.98, ""); titCal->SetNDC(); titCal->SetTextAlign(22); titCal->SetTextSize(0.038); titCal->Draw();
  cCal->SaveAs(plot_base + "/compare_calibration.png");
  delete cCal;
  std::cout << "Plots saved in " << plot_base.Data() << std::endl;
}
