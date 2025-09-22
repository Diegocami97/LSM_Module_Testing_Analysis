// MakeComposite_Fe55_Clusters.C
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TPad.h"
#include "TImage.h"
#include "TSpectrum.h"


// --- For cling vector dictionaries ---
#ifdef __CLING__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#endif

// -------------------------- Config ---------------------------------
int    NBINS_COMP_X   = 401;     // composite grid bins per axis
int    NBINS_COMP_Y   = 401;     // composite grid bins per axis
static const double RANGE_COMP   = 4.0;     // +/- range (pixels)
static const int    BINS_PER_PIX = 25;      // sampling density for pixel maps
static const double E_PER_ADU    = 3.8;     // eV->ADU scaling used in printouts
static const double ENE_MIN      = 4.5;     // energy window [keV]
static const double ENE_MAX      = 7.0;
static const double SIGMA_FRONT_CUT = 1;  // front/back cut
// --------------------- Small helpers --------------------------------


// // clamp utility
// template <typename T>
// static inline T clamp(T v, T lo, T hi){ return std::max(lo, std::min(v, hi)); }

// estimate local sigma from histogram around mu (±window keV)
static double estimate_local_sigma(TH1* h, double mu, double window_keV){
  const int ib = h->GetXaxis()->FindBin(mu);
  const int halfw = std::max(1, int(window_keV / h->GetBinWidth(1)));
  double sw=0, sx=0, sxx=0;
  for (int k=ib-halfw; k<=ib+halfw; ++k){
    if (k<1 || k>h->GetNbinsX()) continue;
    const double x = h->GetXaxis()->GetBinCenter(k);
    const double w = h->GetBinContent(k);
    sw += w; sx += w*x; sxx += w*x*x;
  }
  if (sw<=0) return 0.12; // fallback (keV)
  const double mean = sx/sw;
  const double var  = std::max(1e-6, sxx/sw - mean*mean);
  return std::sqrt(var);
}

static inline void EnsureDir(const TString& path){
  gSystem->mkdir(path, kTRUE);
}

void SaveCanvasHiRes(TCanvas* c, const TString& base, int W=4800, int H=3200) {
  gSystem->mkdir(gSystem->DirName(base), kTRUE);
  std::unique_ptr<TImage> img(TImage::Create());
  img->FromPad(c, W, H);                // re-render at W×H
  img->WriteImage((base + ".png").Data());
  c->SaveAs((base + ".png").Data());    // crisp vector copy
}

static inline TString BaseLeaf(const char* path){
  TString s(path);
  int slash = s.Last('/');
  if (slash >= 0) s = s(slash+1, s.Length()-slash-1);
  int dot = s.Last('.');
  if (dot >= 0) s = s(0, dot);
  return s;
}

static inline int ParseExtNum(const TString& p){
  // find "_ext<digit>"
  int idx = p.Index("_ext");
  if (idx < 0 || idx+4 >= p.Length()) return -1;
  char c = p[idx+4];
  if (c<'1' || c>'4') return -1;
  return c - '0';
}

// Save canvas as PNG+PDF quickly
static inline void SaveCanvas(TCanvas* c, const TString& basePath){
  if (!c) return;
  c->SaveAs(basePath + ".png");
//   c->SaveAs(basePath + ".pdf");
}

// Build a TH2F from pixel lists, with half-integer edges (bin centers at ints)
static TH2F* VectorsToTH2F(const std::vector<int>& x,
                           const std::vector<int>& y,
                           const std::vector<float>& val,
                           const char* name="h2")
{
  if (x.size()!=y.size() || y.size()!=val.size() || x.empty()){
    // Return a tiny empty hist to avoid null checks
    return new TH2F(name, name, 1, -0.5, 0.5, 1, -0.5, 0.5);
  }
  int min_x = *std::min_element(x.begin(), x.end());
  int max_x = *std::max_element(x.begin(), x.end());
  int min_y = *std::min_element(y.begin(), y.end());
  int max_y = *std::max_element(y.begin(), y.end());

  TH2F* h2 = new TH2F(name, name,
                      max_x - min_x + 1, min_x - 0.5, max_x + 0.5,
                      max_y - min_y + 1, min_y - 0.5, max_y + 0.5);
  for (size_t i=0;i<x.size();++i) h2->Fill(x[i], y[i], val[i]);
  return h2;
}

// Rectangle integral without changing axis ranges (robust)
static double SumRect(TH2D* h, double xlo, double xhi, double ylo, double yhi){
  auto ax = h->GetXaxis();
  auto ay = h->GetYaxis();
  const double eps = 1e-9;
  int x1 = ax->FindBin(xlo + eps), x2 = ax->FindBin(xhi - eps);
  int y1 = ay->FindBin(ylo + eps), y2 = ay->FindBin(yhi - eps);
  return h->Integral(x1, x2, y1, y2);
}

// ProjectX over a Y slab [ylo,yhi]
static TH1D* ProjXRect(TH2D* h, const char* name, double ylo, double yhi){
  auto ay = h->GetYaxis();
  int y1 = ay->FindBin(ylo + 1e-9);
  int y2 = ay->FindBin(yhi - 1e-9);
  return h->ProjectionX(name, y1, y2);
}

// ProjectY over an X slab [xlo,xhi]
static TH1D* ProjYRect(TH2D* h, const char* name, double xlo, double xhi){
  auto ax = h->GetXaxis();
  int x1 = ax->FindBin(xlo + 1e-9);
  int x2 = ax->FindBin(xhi - 1e-9);
  return h->ProjectionY(name, x1, x2);
}

// -------------------- results container per extension ----------------
struct ExtResults {
    // bookkeeping
    int    nFiles = 0;
    Long64_t totalClusters = 0;
    int    nAcceptedEnergy = 0;  // number passing energy window (for bookkeeping)
    int    nAcceptedFront  = 0;
    int    nAcceptedBack   = 0;

    // composites & spectra
    TH2D*  hFront  = nullptr;
    TH2D*  hBack   = nullptr;
    TH1F*  hEnergy = nullptr;
    TH1F*  hEnergy_front = nullptr;
    TH1F*  hEnergy_back = nullptr;
    TH2F*  henergy_sigmaxy = nullptr;
    TH2F*  henergy_sigmaxy2 = nullptr;

    // projections (normalized copies for compare tiles)
    TH1D*  projX   = nullptr;    // back: ProjectionX over y ∈ [-4,-3]
    TH1D*  projY   = nullptr;    // back: ProjectionY over x ∈ [-4,-3]
    double projY_int_e = 0.0;    // "hleft" integral [e-]
    double projX_int_e = 0.0;    // "hbelow" integral [e-]
    double projY_mean=0, projY_rms=0, projY_skew=0;
    double projX_mean=0, projX_rms=0, projX_skew=0;

    // integrals (e-)
    double frontIntegral_e = 0.0;
    double backIntegral_e  = 0.0;

    // CTI fractions (front)
    double centerCharge = 0.0;
    double fAbove=0.0, fBelow=0.0, fRight=0.0, fLeft=0.0;

    // energy fits (separate gaussians in windows)
    double A1=0, A1_err=0;   // peak heights for Kα
    double A2=0, A2_err=0;   // peak heights for Kβ
    double mu1=0, mu1_err=0, s1=0, s1_err=0;
    double mu2=0, mu2_err=0, s2=0, s2_err=0;
    double chi2=0; int ndf=0; double p=0;

    // energy fits (separate gaussians in windows)
    double A3=0, A3_err=0;   // peak heights for Kα
    double A4=0, A4_err=0;   // peak heights for Kβ
    double mu3=0, mu3_err=0, s3=0, s3_err=0;
    double mu4=0, mu4_err=0, s4=0, s4_err=0;
    double chi2_2=0; int ndf_2=0; double p_2=0;
};

// --------------------- per-extension processing ----------------------
static ExtResults ProcessOneExtension(const std::vector<TString>& files,
                                      int extNum,
                                      int image_number,
                                      const TString& outDir,
                                      const char* moduleUsed)
{
  ExtResults R;
  R.nFiles = (int)files.size();
  if (R.nFiles==0) return R;

  // Prepare output folder for this extension
  TString extDir = Form("%s/ext%d", outDir.Data(), extNum);
  EnsureDir(extDir);

  // Build a chain with only those files
  TChain chain("clustersRec");
  for (auto& f : files) chain.Add(f);

  // Branches
  std::vector<std::vector<int>>   *pixels_x=nullptr, *pixels_y=nullptr;
  std::vector<std::vector<float>> *pixels_E=nullptr;
  std::vector<float> *PosX=nullptr, *PosY=nullptr;
  std::vector<float> *rmsxy=nullptr, *rmsx=nullptr, *rmsy=nullptr;
  std::vector<float> *Energy=nullptr;

  chain.SetBranchAddress("pixels_x", &pixels_x);
  chain.SetBranchAddress("pixels_y", &pixels_y);
  chain.SetBranchAddress("pixels_E", &pixels_E);
  chain.SetBranchAddress("PosX",     &PosX);
  chain.SetBranchAddress("PosY",     &PosY);
  chain.SetBranchAddress("wSTD_XY",  &rmsxy);
  chain.SetBranchAddress("fSTD_X",   &rmsx);
  chain.SetBranchAddress("fSTD_Y",   &rmsy);
  chain.SetBranchAddress("Energy",   &Energy);

  // Image sizes based on image_number (as per your prior convention)
  int nColumns=0, nRows=0;
  if (image_number==4){ nColumns=6400; nRows=1600; }
  else if (image_number==5){ nColumns=1280; nRows=480; }
  else if (image_number==6){ nColumns=640;  nRows=1600; NBINS_COMP_Y=41;}
  else if (image_number==7){ nColumns=640;  nRows=1600; NBINS_COMP_X=41;}
  else { nColumns=640; nRows=640; }

  // Hists
  R.hFront  = new TH2D(Form("hFront_ext%d",extNum),  "Composite Cluster Front Events",
  NBINS_COMP_X, -RANGE_COMP, RANGE_COMP,
  NBINS_COMP_Y, -RANGE_COMP, RANGE_COMP);
  R.hBack   = new TH2D(Form("hBack_ext%d",extNum),   "Composite Cluster Back Events",
  NBINS_COMP_X, -RANGE_COMP, RANGE_COMP,
  NBINS_COMP_Y, -RANGE_COMP, RANGE_COMP);
  R.hEnergy = new TH1F(Form("hEnergy_ext%d",extNum), "hClusterE",
                       50, ENE_MIN, ENE_MAX);
  R.hEnergy_front = new TH1F(Form("hEnergy_ext%d",extNum), "hClusterE",
                       50, ENE_MIN, ENE_MAX);
  R.hEnergy_back = new TH1F(Form("hEnergy_ext%d",extNum), "hClusterE",
                       50, ENE_MIN, ENE_MAX);

  // Additional distributions (saved later)
  TH1F* hPosX = new TH1F(Form("hPosX_ext%d",extNum), "PosX", 20, 0, 1);   // range updated below
  TH1F* hPosY = new TH1F(Form("hPosY_ext%d",extNum), "PosY", 20, 0, 1);
  TH2F* hPosXY= new TH2F(Form("hPosXY_ext%d",extNum),"PosXY", 10, 0, 1, 10, 0, 1);
  TH1F* hSigmaxy = new TH1F(Form("hSigmaxy_ext%d",extNum), "Sigmaxy", 50, 0, 2.0);
  TH1F* hSigmax  = new TH1F(Form("hSigmax_ext%d", extNum),  "Sigmax",  50, 0, 2.0);
  TH1F* hSigmay  = new TH1F(Form("hSigmay_ext%d", extNum),  "Sigmay",  50, 0, 2.0);
  TH2F* hposx_sigmaxy = nullptr;
  R.henergy_sigmaxy= new TH2F(Form("hE_sigmaxy_ext%d",extNum),
                                  "Energy vs #sigma_{xy}", 50, ENE_MIN, ENE_MAX, 50, 0, 2.0);

  

  // finalize axis ranges for the helper hists
  hPosX->GetXaxis()->SetLimits(0, nColumns);
  hPosY->GetXaxis()->SetLimits(0, nRows);
  hPosXY = new TH2F(Form("hPosXY_ext%d",extNum), "PosXY", std::max(1,nColumns/100), 0, nColumns,
                                                   std::max(1,nRows/100),   0, nRows);
  hposx_sigmaxy = new TH2F(Form("hposx_sigmaxy_ext%d",extNum),
                           "posx_sigmaxy", nColumns, 0, nColumns, 50, 0, 2.0);

  // Loop
  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t i=0;i<nEntries;++i){
    chain.GetEntry(i);
    const size_t nC = PosX->size();
    R.totalClusters += nC;

    for (size_t n=0;n<nC;++n){
      const double E = Energy->at(n);
      const double x = PosX->at(n);
      const double y = PosY->at(n);
      const double sx= rmsx->at(n);
      const double sy= rmsy->at(n);
      const double sxy= std::sqrt(0.5*(sx*sx + sy*sy));

      const auto& px = pixels_x->at(n);
      const auto& py = pixels_y->at(n);
      const auto& pv = pixels_E->at(n);

      // book-keeping distributions for events in energy window
      if (E>ENE_MIN && E<ENE_MAX){
        hSigmaxy->Fill(sxy); hSigmax->Fill(sx); hSigmay->Fill(sy);
        hPosX->Fill(x); hPosY->Fill(y);
        hPosXY->Fill(x,y,E);
        hposx_sigmaxy->Fill(x, sxy, E);
        R.henergy_sigmaxy->Fill(E, sxy);
        R.nAcceptedEnergy++;
      }

      // choose variable for front/back split
      double sCutVar = sxy;
      if (image_number==6) sCutVar = sx;
      if (image_number==7) sCutVar = sy;

      // FRONT selection
      if (E>ENE_MIN && E<ENE_MAX && sCutVar < SIGMA_FRONT_CUT){
        R.hEnergy->Fill(E);
        R.hEnergy_front->Fill(E);
        TH2F* hview = VectorsToTH2F(px, py, pv, "h2_front_tmp");
        // Accumulate into composite centered on (x,y)
        for (int ix=1; ix<=NBINS_COMP_X; ++ix){
          for (int iy=1; iy<=NBINS_COMP_Y; ++iy){
            const double gx = x + R.hFront->GetXaxis()->GetBinLowEdge(ix);
            const double gy = y + R.hFront->GetYaxis()->GetBinLowEdge(iy);
            const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (BINS_PER_PIX*BINS_PER_PIX);
            R.hFront->SetBinContent(ix, iy, R.hFront->GetBinContent(ix,iy) + bc);
          }
        }
        delete hview;
        R.nAcceptedFront++;
      }

      // BACK selection
      if (E>ENE_MIN && E<ENE_MAX && sCutVar > SIGMA_FRONT_CUT){
        R.hEnergy_back->Fill(E);

        TH2F* hview = VectorsToTH2F(px, py, pv, "h2_back_tmp");
        for (int ix=1; ix<=NBINS_COMP_X; ++ix){
          for (int iy=1; iy<=NBINS_COMP_Y; ++iy){
            const double gx = x + R.hBack->GetXaxis()->GetBinLowEdge(ix);
            const double gy = y + R.hBack->GetYaxis()->GetBinLowEdge(iy);
            const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (BINS_PER_PIX*BINS_PER_PIX);
            R.hBack->SetBinContent(ix, iy, R.hBack->GetBinContent(ix,iy) + bc);
          }
        }
        delete hview;
        R.nAcceptedBack++;
      }
    } // clusters
  } // entries

  // Normalize composites
  if (R.nAcceptedFront>0) R.hFront->Scale(1.0/R.nAcceptedFront);
  if (R.nAcceptedBack >0) R.hBack ->Scale(1.0/R.nAcceptedBack);

  // Integrals in electrons
  R.frontIntegral_e = R.hFront->Integral()*1000.0/E_PER_ADU;
  R.backIntegral_e  = R.hBack ->Integral()*1000.0/E_PER_ADU;

  // --------- FRONT: CTI fractions w.r.t center ----------
  // regions relative to pixel boxes:
  // center: x,y in [-0.5,0.5]
  // above : x in [-0.5,0.5], y in [0.5,1.5]
  // below : x in [-0.5,0.5], y in [-1.5,-0.5]
  // right : x in [0.5,1.5], y in [-0.5,0.5]
  // left  : x in [-1.5,-0.5],y in [-0.5,0.5]
  if (R.hFront){
    const double c = SumRect(R.hFront, -0.5, 0.5,  -0.5,  0.5);
    const double a = SumRect(R.hFront, -0.5, 1.5,   0.5,  2.5);
    const double b = SumRect(R.hFront, -0.5, 0.5,  -1.5, -0.5);
    const double r = SumRect(R.hFront,  1.5, -0.5,  2.5,  0.5);
    const double l = SumRect(R.hFront, -1.5,-0.5,  -0.5,  0.5);
    R.centerCharge = c;
    if (c>0){
      R.fAbove = a/c; R.fBelow = b/c; R.fRight = r/c; R.fLeft = l/c;
    }
  }

  // ---------------- BACK: projections & stats ----------------
  if (R.hBack){
    // Projections on the slab [-4,-3] pixels on the other axis
    TH1D* hbelow = ProjXRect(R.hBack, Form("hbelow_ext%d",extNum), -4.0, -3.0);
    TH1D* hleft  = ProjYRect(R.hBack, Form("hleft_ext%d", extNum), -4.0, -3.0);

    // Store stats BEFORE normalization
    R.projY_int_e = hleft ->Integral()*1000.0/E_PER_ADU;
    R.projX_int_e = hbelow->Integral()*1000.0/E_PER_ADU;

    R.projY_mean = hleft->GetMean();  R.projY_rms = hleft->GetRMS();  R.projY_skew = hleft->GetSkewness();
    R.projX_mean = hbelow->GetMean(); R.projX_rms = hbelow->GetRMS(); R.projX_skew = hbelow->GetSkewness();

    // Normalize for overlays in compare canvases
    if (hbelow->Integral()>0) hbelow->Scale(1.0/hbelow->Integral());
    if (hleft ->Integral()>0) hleft ->Scale(1.0/hleft ->Integral());
    R.projX = (TH1D*)hbelow->Clone(Form("projX_ext%d",extNum));
    R.projY = (TH1D*)hleft ->Clone(Form("projY_ext%d",extNum));

    // Per-extension canvas: composite + projections
    gStyle->SetOptStat(0);
    TCanvas* cBack = new TCanvas(Form("cBack_ext%d",extNum),"Back composite + projections",1200,550);
    cBack->Divide(2,1);
    // Left: heatmap
    cBack->cd(1);
    gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetLogz();
    TH2D* hb = (TH2D*)R.hBack->Clone();
    hb->GetXaxis()->SetRangeUser(-4,4); hb->GetYaxis()->SetRangeUser(-4,4);
    hb->SetTitle("");
    hb->GetXaxis()->SetTitle("#Delta x [pix]");
    hb->GetYaxis()->SetTitle("#Delta y [pix]");
    hb->Draw("COLZ");
    // Title
    { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Back-side Composite Cluster - %s - ext %d", moduleUsed, extNum)); }

    // Right: projections with legend
    cBack->cd(2);
    gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16);
    R.projX->SetTitle(";Pixel offset [pix];Normalized counts");
    R.projX->SetLineColor(kBlack); 
    R.projX->SetLineWidth(3);
    R.projX->SetFillColorAlpha(kBlack, 0.08);
    R.projX->SetMarkerStyle(20);
    R.projX->SetMarkerColor(kBlack);

    R.projY->SetLineColor(kBlue+1); 
    R.projY->SetLineWidth(3);
    R.projX->SetFillColorAlpha(kBlue+1, 0.08);
    R.projX->SetMarkerStyle(20);
    R.projX->SetMarkerColor(kBlue+1);

    double ymax = std::max(R.projX->GetMaximum(), R.projY->GetMaximum())*1.15;
    R.projX->SetMaximum(ymax); R.projX->SetMinimum(0.0);
    R.projX->Draw("HIST"); R.projY->Draw("HIST SAME");
    auto legB = new TLegend(0.57,0.72,0.89,0.88);
    legB->AddEntry(R.projX,"ProjX (Y: -3 #rightarrow -2 px)","l");
    legB->AddEntry(R.projY,"ProjY (X: -3 #rightarrow -2 px)","l");
    legB->Draw();
    // Stats box (integrals & moments)
    // auto paveB = new TPaveText(0.16,0.60,0.54,0.88,"NDC");
    // paveB->SetFillStyle(0); paveB->SetBorderSize(0); paveB->SetTextFont(42); paveB->SetTextSize(0.035);
    // paveB->AddText("BACK-SIDE EVENTS CTI METRICS");
    // paveB->AddText(Form("ProjY int [e^-]: %.1f | mean=%.3f RMS=%.3f skew=%.3f",
    //                     R.projY_int_e, R.projY_mean, R.projY_rms, R.projY_skew));
    // paveB->AddText(Form("ProjX int [e^-]: %.1f | mean=%.3f RMS=%.3f skew=%.3f",
    //                     R.projX_int_e, R.projX_mean, R.projX_rms, R.projX_skew));
    // paveB->Draw();
    SaveCanvas(cBack, Form("%s/backside", extDir.Data()));
    delete cBack;
  }

  // -------------- FRONT: composite + region bars ----------------
  {
    gStyle->SetOptStat(0);
    // gStyle->SetNumberContours(255);

    TCanvas* cFront = new TCanvas(Form("cFront_ext%d",extNum),"Front composite + CTI bars",1200,550);
    cFront->Divide(2,1);

    // Left: heatmap with squares
    cFront->cd(1);
    gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetLogz();
    TH2D* hf=(TH2D*)R.hFront->Clone();
    hf->GetXaxis()->SetRangeUser(-4,4); hf->GetYaxis()->SetRangeUser(-4,4);
    hf->SetTitle("");
    hf->GetXaxis()->SetTitle("#Delta x [pix]");
    hf->GetYaxis()->SetTitle("#Delta y [pix]");
    hf->Draw("COLZ");
    // color boxes
    auto drawBoxOnly=[&](double x1,double y1,double x2,double y2, Color_t lc){
      TBox* b=new TBox(x1,y1,x2,y2); b->SetLineColor(lc); b->SetLineWidth(3); b->SetFillStyle(0); b->Draw("same"); return b; };
    TBox* boxCenter = drawBoxOnly(-0.5,-0.5, 0.5, 0.5, kBlack);
    TBox* boxAbove  = drawBoxOnly(-0.5, 1.5, 0.5, 2.5, kBlue+1);
    TBox* boxBelow  = drawBoxOnly(-0.5,-1.5, 0.5,-0.5, kRed+1);
    TBox* boxRight  = drawBoxOnly( 1.5,-0.5, 2.5, 0.5, kGreen+2);
    TBox* boxLeft   = drawBoxOnly(-1.5,-0.5,-0.5, 0.5, kMagenta);

    { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Front-Side Composite Cluster - %s - Ext %d", moduleUsed, extNum)); }

    // Right: region bars + legend with numbers
    cFront->cd(2);
    gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.18);
    TH1F* frame=new TH1F(Form("frameF_ext%d",extNum),";Region;Fraction (w.r.t. center)",5,0.5,5.5);
    frame->GetXaxis()->SetBinLabel(1,"Above");
    frame->GetXaxis()->SetBinLabel(2,"Below");
    frame->GetXaxis()->SetBinLabel(3,"Right");
    frame->GetXaxis()->SetBinLabel(4,"Left");
    frame->GetXaxis()->SetBinLabel(5,"Center=1.0");
    double yMaxFrac = 1.25*std::max(1.0, std::max({R.fAbove,R.fBelow,R.fRight,R.fLeft}));
    frame->SetMinimum(0.0); frame->SetMaximum(yMaxFrac); frame->Draw("AXIS");
    TLine* ref=new TLine(0.5,1.0,5.5,1.0); ref->SetLineStyle(2); ref->SetLineColor(kGray+2); ref->Draw();
    auto drawBar=[&](int bin,double y, Color_t col){ double x1=frame->GetXaxis()->GetBinLowEdge(bin), x2=frame->GetXaxis()->GetBinUpEdge(bin);
      TBox* b=new TBox(x1,0.0,x2,y); b->SetFillColorAlpha(col,0.35); b->SetLineColor(col); b->SetLineWidth(3); b->Draw("same"); return b; };
    TBox* bAbove  = drawBar(1,R.fAbove, kBlue+1);
    TBox* bBelow  = drawBar(2,R.fBelow, kRed+1);
    TBox* bRight  = drawBar(3,R.fRight, kGreen+2);
    TBox* bLeft   = drawBar(4,R.fLeft,  kMagenta);
    TBox* bCenter = drawBar(5,1.0,      kGray+1);

    auto legF = new TLegend(0.18, 0.58, 0.55, 0.88);
    legF->SetBorderSize(0); legF->SetFillStyle(0); legF->SetTextFont(42); legF->SetTextSize(0.035);
    legF->AddEntry(bCenter, Form("Center charge = %.0f", R.centerCharge), "f");
    legF->AddEntry(bAbove,  Form("Above  = %.3f", R.fAbove),  "f");
    legF->AddEntry(bBelow,  Form("Below  = %.3f", R.fBelow),  "f");
    legF->AddEntry(bRight,  Form("Right  = %.3f", R.fRight),  "f");
    legF->AddEntry(bLeft,   Form("Left   = %.3f", R.fLeft ),  "f");
    legF->Draw();
    SaveCanvas(cFront, Form("%s/frontside", extDir.Data()));
    delete cFront;
  }

  // -------------- ENERGY: double-Gaussian style plot - Front----------------
  {
    TCanvas* cE = new TCanvas(Form("cE_ext%d",extNum),"Cluster Energy (double Gaussian)",900,600);
    cE->cd(); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

    // matplotlib-like bars
    R.hEnergy_front->SetTitle(";Energy [keV];Counts");
    R.hEnergy_front->SetFillColorAlpha(kGray+1, 0.6);
    R.hEnergy_front->SetLineColor(kGray+2);
    R.hEnergy_front->SetBarWidth(0.95);
    R.hEnergy_front->SetBarOffset(0.025);
    R.hEnergy_front->Draw("BAR");
    auto outline=(TH1F*)R.hEnergy_front->Clone(Form("hE_outline_ext%d",extNum));
    outline->SetFillStyle(0); outline->SetLineColor(kBlack); outline->SetLineWidth(2); outline->Draw("BAR SAME");

    // // fit two per-peak windows
    // TF1* g1 = new TF1(Form("g1_ext%d",extNum),"gaus",5.0,6.4);
    // TF1* g2 = new TF1(Form("g2_ext%d",extNum),"gaus",6.0,6.8);
    // g1->SetParameters(R.hEnergy_front->GetMaximum(), 6.4, 0.10);
    // g2->SetParameters(R.hEnergy_front->GetMaximum()/3.0, 6.40, 0.12);
    // g1->SetParLimits(2,0.03,0.40); 

    // g2->SetParLimits(2,0.03,0.40);
    // R.hEnergy_front->Fit(g1,"RQ0S"); 
    
    // R.hEnergy_front->Fit(g2,"RQ0S");

    // g1->SetLineColor(kRed+1); g1->SetLineWidth(3);
    // g2->SetLineColor(kBlue+1); g2->SetLineWidth(3);
    // g1->Draw("SAME"); g2->Draw("SAME");

    // fit two per-peak windows
    TF1* g1 = new TF1(Form("g1_ext%d",extNum),"gaus",4.5,6.3); 
    g1->SetParameters(R.hEnergy_front->GetMaximum(), 6.4, 0.10);
    g1->SetParLimits(2,0.03,0.40);
    R.hEnergy_front->Fit(g1,"RQ0S");

    // Get the mean and 2*sigma from first fit
    double mu1 = g1->GetParameter(1);
    double sigma1 = g1->GetParameter(2);
    double two_sig = sigma1 * 2;  

    // Set second peak window based on first peak
    double g2_left = mu1 + two_sig;
    double g2_right = g2_left + 1;

    TF1* g2 = new TF1(Form("g2_ext%d",extNum),"gaus",g2_left,g2_right);
    g2->SetParameters(R.hEnergy_front->GetMaximum()/3.0, (g2_left+g2_right)/2.0, 0.12);  // center initial guess
    // g2->SetParameters(R.hEnergy_front->GetMaximum()/3.0, 5.25, 0.12);  // center initial guess
    g2->SetParLimits(2,0.03,0.40);
    R.hEnergy_front->Fit(g2,"RQ0S");

    double mu2 = g2->GetParameter(1);
    double sigma2 = g2->GetParameter(2);

    g1->SetLineColor(kRed+1); g1->SetLineWidth(3);
    g2->SetLineColor(kBlue+1); g2->SetLineWidth(3);
    g1->Draw("SAME"); g2->Draw("SAME");

    // extract params & manual chi2 over union of windows
    R.mu1=g1->GetParameter(1); R.mu1_err=g1->GetParError(1);
    R.s1 =g1->GetParameter(2); R.s1_err =g1->GetParError(2);
    R.mu2=g2->GetParameter(1); R.mu2_err=g2->GetParError(1);
    R.s2 =g2->GetParameter(2); R.s2_err =g2->GetParError(2);
    R.A1 = g1->GetParameter(0); R.A1_err = g1->GetParError(0);
    R.A2 = g2->GetParameter(0); R.A2_err = g2->GetParError(0);


    auto inG1=[&](double x){ return x>=5.3 && x<=6.1; };
    auto inG2=[&](double x){ return x>=6.0 && x<=6.8; };
    R.chi2=0; R.ndf=0;
    for (int ib=1; ib<=R.hEnergy_front->GetNbinsX(); ++ib){
      const double x = R.hEnergy_front->GetXaxis()->GetBinCenter(ib);
      if (!(inG1(x) || inG2(x))) continue;
      const double y = R.hEnergy_front->GetBinContent(ib);
      const double yfit = g1->Eval(x) + g2->Eval(x);
      const double var = y + 1e-6;
      R.chi2 += (y-yfit)*(y-yfit)/var;
      ++R.ndf;
    }
    R.ndf = std::max(1, R.ndf - 6);
    R.p = TMath::Prob(R.chi2, R.ndf);

    // Legend with splitline labels
    auto leg = new TLegend(0.58,0.63,0.89,0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(g1,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu1,R.mu1_err,R.s1,R.s1_err),"l");
    leg->AddEntry(g2,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu2,R.mu2_err,R.s2,R.s2_err),"l");
    leg->AddEntry(R.hEnergy_front,"Data","f");
    leg->Draw();

    // per-tile title (consistency with compare)
    { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Front-Side Energy Spectrum - %s - Ext %d", moduleUsed, extNum)); }

    SaveCanvas(cE, Form("%s/cluster_energy", extDir.Data()));
    delete cE;
  }

  // -------------- ENERGY: double-Gaussian style plot - Back ----------------
  {
    TCanvas* cE_back = new TCanvas(Form("cE_ext%d",extNum),"Cluster Energy (double Gaussian)",900,600);
    cE_back->cd(); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

    // matplotlib-like bars
    R.hEnergy_back->SetTitle(";Energy [keV];Counts");
    R.hEnergy_back->SetFillColorAlpha(kGray+1, 0.6);
    R.hEnergy_back->SetLineColor(kGray+2);
    R.hEnergy_back->SetBarWidth(0.95);
    R.hEnergy_back->SetBarOffset(0.025);
    R.hEnergy_back->Draw("BAR");
    auto outline=(TH1F*)R.hEnergy_back->Clone(Form("hE_outline_ext%d",extNum));
    outline->SetFillStyle(0); outline->SetLineColor(kBlack); outline->SetLineWidth(2); outline->Draw("BAR SAME");

    // fit two per-peak windows
    TF1* g3 = new TF1(Form("g3_ext%d",extNum),"gaus",4.5,6.3);  // broader range to 6.3
    g3->SetParameters(R.hEnergy_back->GetMaximum(), 6.4, 0.10);
    g3->SetParLimits(2,0.03,0.40);
    R.hEnergy_back->Fit(g3,"RQ0S");

    // Get the mean and 2*sigma from first fit
    double mu3 = g3->GetParameter(1);
    double sigma3 = g3->GetParameter(2);
    double two_sig = sigma3 * 2;  

    // Set second peak window based on first peak
    double g4_left = mu3 + two_sig;
    double g4_right = g4_left + 1;

    TF1* g4 = new TF1(Form("g4_ext%d",extNum),"gaus",g4_left,g4_right);
    g4->SetParameters(R.hEnergy_back->GetMaximum()/3.0, (g4_left+g4_right)/2.0, 0.12);  // center initial guess
    g4->SetParLimits(2,0.03,0.40);
    R.hEnergy_back->Fit(g4,"RQ0S");

    double mu2 = g4->GetParameter(1);
    double sigma2 = g4->GetParameter(2);

    g3->SetLineColor(kRed+1); g3->SetLineWidth(3);
    g4->SetLineColor(kBlue+1); g4->SetLineWidth(3);
    g3->Draw("SAME"); g4->Draw("SAME");

    // extract params & manual chi2 over union of windows
    R.mu3=g3->GetParameter(1); R.mu3_err=g3->GetParError(1);
    R.s3 =g3->GetParameter(2); R.s3_err =g3->GetParError(2);
    R.mu4=g4->GetParameter(1); R.mu4_err=g4->GetParError(1);
    R.s4 =g4->GetParameter(2); R.s4_err =g4->GetParError(2);
    R.A3 = g3->GetParameter(0); R.A3_err = g3->GetParError(0);
    R.A4 = g4->GetParameter(0); R.A4_err = g4->GetParError(0);

    double ymaxPanel = 1.2 * std::max({
      R.hEnergy_back->GetMaximum(),
      R.A3,
      R.A4
    });
    R.hEnergy_back->SetMaximum(ymaxPanel);
    outline->SetMaximum(ymaxPanel);

    // NEED TO LOOK AT THIS MORE
    auto inG3=[&](double x){ return x>=5.3 && x<=6.1; };
    auto inG4=[&](double x){ return x>=6.0 && x<=6.8; };
    R.chi2_2=0; R.ndf_2=0;
    for (int ib=1; ib<=R.hEnergy_back->GetNbinsX(); ++ib){
      const double x = R.hEnergy_back->GetXaxis()->GetBinCenter(ib);
      if (!(inG3(x) || inG4(x))) continue;
      const double y = R.hEnergy_back->GetBinContent(ib);
      const double yfit = g3->Eval(x) + g4->Eval(x);
      const double var = y + 1e-6;
      R.chi2_2 += (y-yfit)*(y-yfit)/var;
      ++R.ndf_2;
    }
    R.ndf_2 = std::max(1, R.ndf_2 - 6);
    R.p_2 = TMath::Prob(R.chi2_2, R.ndf_2);

    // Legend with splitline labels
    auto leg = new TLegend(0.58,0.63,0.89,0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(g3,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu3,R.mu3_err,R.s3,R.s3_err),"l");
    leg->AddEntry(g4,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu4,R.mu4_err,R.s4,R.s4_err),"l");
    leg->AddEntry(R.hEnergy_back,"Data","f");
    leg->Draw();

    // per-tile title (consistency with compare)
    { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Back-Side Energy Spectrum - %s - Ext %d", moduleUsed, extNum)); }

    SaveCanvas(cE_back, Form("%s/cluster_energy_back", extDir.Data()));
    delete cE_back;
  }

  // -------------- ENERGY: double-Gaussian style plot - Front & Back ----------------
  {
    TCanvas* cE_fb = new TCanvas(Form("cE_fb_ext%d",extNum),
                              "Cluster Energy: Front & Back", 900, 600);
    cE_fb->cd(); gPad->SetGridx(); gPad->SetGridy();
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);


    // ----------- FRONT HIST ---------
    R.hEnergy_front->SetTitle(";Energy [keV];Counts");
    R.hEnergy_front->SetFillColorAlpha(kRed,0.6);
    R.hEnergy_front->SetLineColor(kRed+2);
    R.hEnergy_front->SetBarWidth(0.95);
    R.hEnergy_front->SetBarOffset(0);
    R.hEnergy_front->Draw("BAR");

    auto frontOutline = (TH1F*)R.hEnergy_front->Clone("hE_outline_front");
    frontOutline->SetFillStyle(0); frontOutline->SetLineColor(kRed);
    frontOutline->SetLineWidth(1);
    frontOutline->Draw("BAR SAME");

    // ----------- BACK HIST ---------
    R.hEnergy_back->SetLineColor(kBlue+1);
    R.hEnergy_back->SetFillColorAlpha(kBlue+1, 0.6);
    R.hEnergy_back->SetBarWidth(0.95);
    R.hEnergy_back->SetBarOffset(0);  
    R.hEnergy_back->Draw("BAR SAME");

    auto backOutline = (TH1F*)R.hEnergy_back->Clone("hE_outline_back");
    backOutline->SetFillStyle(0); backOutline->SetLineColor(kBlue);
    backOutline->SetLineWidth(1);
    backOutline->Draw("BAR SAME");


    // draw fitted components using stored params
    TF1* g1=new TF1(Form("cEg1_%d",extNum),"gaus", R.mu1-4*R.s1, R.mu1+4*R.s1);
    TF1* g2=new TF1(Form("cEg2_%d",extNum),"gaus", R.mu2-4*R.s2, R.mu2+4*R.s2);
    TF1* g3=new TF1(Form("cEg3_%d",extNum),"gaus", R.mu3-4*R.s3, R.mu3+4*R.s3);
    TF1* g4=new TF1(Form("cEg4_%d",extNum),"gaus", R.mu4-4*R.s4, R.mu4+4*R.s4);

    // fallback if a fit failed
    const double A1 = (R.A1>0 ? R.A1 : backOutline->GetMaximum());
    const double A2 = (R.A2>0 ? R.A2 : 0.3*backOutline->GetMaximum());
    const double A3 = (R.A3>0 ? R.A3 : backOutline->GetMaximum());
    const double A4 = (R.A4>0 ? R.A4 : 0.3*backOutline->GetMaximum());

    g1->SetParameters(R.A1, R.mu1, std::max(0.03, R.s1));
    g2->SetParameters(R.A2, R.mu2, std::max(0.03, R.s2));
    g1->SetNpx(600); g2->SetNpx(600);
    g3->SetParameters(R.A3, R.mu3, std::max(0.03, R.s3));
    g4->SetParameters(R.A4, R.mu4, std::max(0.03, R.s4));
    g3->SetNpx(600); g4->SetNpx(600);

    // ---- DRAW FITS (already existing objects: g1/g2/g3/g4) ----
    g1->SetLineColor(kRed+1);   g1->SetLineWidth(4); g1->SetLineStyle(1); g1->Draw("SAME");
    g2->SetLineColor(kRed-7);   g2->SetLineWidth(4); g2->SetLineStyle(1); g2->Draw("SAME");
    g3->SetLineColor(kBlue+1);  g3->SetLineWidth(2); g3->SetLineStyle(1); g3->Draw("SAME"); 
    g4->SetLineColor(kBlue-7);   g4->SetLineWidth(2); g4->SetLineStyle(1); g4->Draw("SAME");



    /// Left legend for "Front data"
    auto leg_front = new TLegend(0.145, 0.65, 0.385, 0.89);
    leg_front->SetBorderSize(1);
    leg_front->SetFillStyle(1001);
    leg_front->SetTextFont(42);
    leg_front->SetTextSize(0.030);

    leg_front->AddEntry(R.hEnergy_front, "Front", "f");
    leg_front->AddEntry(g1, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu1, R.mu1_err, R.s1, R.s1_err), "l");
    leg_front->AddEntry(g2, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu2, R.mu2_err, R.s2, R.s2_err), "l");

    // Right legend for "Back data"
    auto leg_back = new TLegend(0.64, 0.65, 0.88, 0.89);
    leg_back->SetBorderSize(1);
    leg_back->SetFillStyle(1001);
    leg_back->SetTextFont(42);
    leg_back->SetTextSize(0.030);

    leg_back->AddEntry(R.hEnergy_back, "Back", "f");
    leg_back->AddEntry(g3, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu3, R.mu3_err, R.s3, R.s3_err), "l");
    leg_back->AddEntry(g4, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R.mu4, R.mu4_err, R.s4, R.s4_err), "l");

    // Draw both legends
    leg_front->Draw();
    leg_back->Draw();

    // TITLE
    { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Front & Back Cluster Energy - %s - Ext %d", moduleUsed, extNum));
    }

    SaveCanvas(cE_fb, Form("%s/cluster_energy_fb", extDir.Data()));
    delete cE_fb;
  }

  // ----------------- Extra plots per extension ------------------
  {
    // posX vs sigmaXY
    
    TCanvas* c1=new TCanvas(Form("c_posx_sigmaxy_ext%d",extNum),"posx vs sigmaxy",800,600);
    c1->cd(); gPad->SetLogz(); hposx_sigmaxy->SetTitle("Mean x vs #sigma_{xy}");
    hposx_sigmaxy->GetXaxis()->SetTitle("Mean x [col]");
    hposx_sigmaxy->GetYaxis()->SetTitle("#sigma_{xy} [pix]");
    hposx_sigmaxy->Draw("COLZ");
    SaveCanvas(c1, Form("%s/posx_vs_sigmaXY", extDir.Data())); delete c1;

    // energy vs sigmaXY
    TCanvas* c2=new TCanvas(Form("c_E_sigmaxy_ext%d",extNum),"Energy vs sigmaxy",800,600);
    c2->cd(); R.henergy_sigmaxy->SetTitle("Cluster Energy vs #sigma_{xy}");
    R.henergy_sigmaxy->GetXaxis()->SetTitle("Energy [keV]");
    R.henergy_sigmaxy->GetYaxis()->SetTitle("#sigma_{xy} [pix]");
    R.henergy_sigmaxy->Draw("COLZ");
    SaveCanvas(c2, Form("%s/energy_vs_sigmaXY", extDir.Data())); delete c2;

    // sigma hists
    TCanvas* c3=new TCanvas(Form("c_sigmaxy_ext%d",extNum),"sigmaXY",700,500);
    hSigmaxy->SetTitle("#sigma_{xy} Distribution;#sigma_{xy} [pix];Counts");
    hSigmaxy->SetLineColor(kRed); hSigmaxy->SetFillColor(kRed); hSigmaxy->SetFillStyle(3001);
    hSigmaxy->Draw("hist"); SaveCanvas(c3, Form("%s/sigmaXY", extDir.Data())); delete c3;

    TCanvas* c4=new TCanvas(Form("c_sigmax_ext%d",extNum),"sigmaX",700,500);
    hSigmax->SetTitle("#sigma_{x} Distribution;#sigma_{x} [pix];Counts");
    hSigmax->SetLineColor(kRed); hSigmax->SetFillColor(kRed); hSigmax->SetFillStyle(3001);
    hSigmax->Draw("hist"); SaveCanvas(c4, Form("%s/sigmaX", extDir.Data())); delete c4;

    TCanvas* c5=new TCanvas(Form("c_sigmay_ext%d",extNum),"sigmaY",700,500);
    hSigmay->SetTitle("#sigma_{y} Distribution;#sigma_{y} [pix];Counts");
    hSigmay->SetLineColor(kRed); hSigmay->SetFillColor(kRed); hSigmay->SetFillStyle(3001);
    hSigmay->Draw("hist"); SaveCanvas(c5, Form("%s/sigmaY", extDir.Data())); delete c5;

    // positions
    TCanvas* c6=new TCanvas(Form("c_posXY_ext%d",extNum),"PosXY",800,700);
    hPosXY->SetTitle("Cluster Position Distribution;X [col];Y [row]");
    hPosXY->Draw("BOXCOLZ"); SaveCanvas(c6, Form("%s/posXY", extDir.Data())); delete c6;

    TCanvas* c7=new TCanvas(Form("c_posX_ext%d",extNum),"PosX",700,500);
    hPosX->SetTitle("Cluster X Position;X [col];Counts");
    hPosX->SetLineColor(kBlue); hPosX->SetFillColor(kBlue); hPosX->SetFillStyle(3001);
    hPosX->Draw("hist"); SaveCanvas(c7, Form("%s/posX", extDir.Data())); delete c7;

    TCanvas* c8=new TCanvas(Form("c_posY_ext%d",extNum),"PosY",700,500);
    hPosY->SetTitle("Cluster Y Position;Y [row];Counts");
    hPosY->SetLineColor(kBlue); hPosY->SetFillColor(kBlue); hPosY->SetFillStyle(3001);
    hPosY->Draw("hist"); SaveCanvas(c8, Form("%s/posY", extDir.Data())); delete c8;
  }

  return R;
}

// ---------------------- main entry point ----------------------------
void MakeComposite_Fe55_Clusters_Extensions(const char* filelist="test_list.txt",
                                 const int   image_number=4,
                                 const char* moduleUsed="DMXX",
                                 const char* outDir="out")
{
  // Read list and bucket by extension
  std::ifstream fin(filelist);
  if (!fin){
    std::cerr << "[ERROR] Cannot open file list: " << filelist << "\n";
    return;
  }
  std::map<int, std::vector<TString>> byExt; // 1..4
  std::string line;
  while (std::getline(fin, line)){
    if (line.empty()) continue;
    TString p(line.c_str());
    int ext = ParseExtNum(p);
    if (ext>=1 && ext<=4) byExt[ext].push_back(p);
  }
  fin.close();

  // Output root folder: out/<Module>/<leaf of listfile>
  TString leaf = BaseLeaf(filelist);
  TString baseTag = Form("%s/Fe_55_Analysis_Plots/%s/%s", outDir, moduleUsed, leaf.Data());
  EnsureDir(baseTag);
  EnsureDir(baseTag + "/compare");

  // Process each extension, collect results
  ExtResults R[5]; // 1..4
  for (int e=1;e<=4;++e){
    if (byExt.count(e)==0 || byExt[e].empty()){
      std::cout << "[INFO] No files for ext" << e << "\n";
      continue;
    }
    std::cout << "[RUN] Processing ext" << e << " with " << byExt[e].size() << " files\n";
    R[e] = ProcessOneExtension(byExt[e], e, image_number, baseTag, moduleUsed);
  }

  // ---------------- Comparison 2x2: common scales -------------------
  // Back heatmap z-range
  double zminB=1e-6, zmaxB=0.0;
  double ymaxProj=0.0;
  for (int e=1;e<=4;++e){
    if (!R[e].hBack) continue;
    zmaxB = std::max(zmaxB, R[e].hBack->GetMaximum());
    if (R[e].projX) ymaxProj = std::max({ymaxProj, R[e].projX->GetMaximum(), R[e].projY->GetMaximum()});
  }
  if (zmaxB<=0) zmaxB=1.0;

  // Front heatmap z-range
  double zminF=1e-6, zmaxF=0.0;
  for (int e=1;e<=4;++e){
    if (!R[e].hFront) continue;
    zmaxF = std::max(zmaxF, R[e].hFront->GetMaximum());
  }
  if (zmaxF<=0) zmaxF=1.0;

  // Energy y-range
  double ymaxEnergy=0.0;
  for (int e=1;e<=4;++e){
    if (!R[e].hEnergy) continue;
    ymaxEnergy = std::max(ymaxEnergy, R[e].hEnergy->GetMaximum());
  }
  if (ymaxEnergy<=0) ymaxEnergy=1.0;

  // ----------------- (1) Back-side comparison 2×2 -------------------
  {
    TCanvas* c = new TCanvas("compare_back","Back 2x2", 2800, 2200);
    c->Divide(2,2);
    int pad=1;
    for (int e=1; e<=4; ++e){
      if (!R[e].hBack) { ++pad; continue; }
      c->cd(pad++);
      // two subpads per tile
      TPad* left = new TPad(Form("bL_%d",e),"",0.0,0.0,0.57,1.0); left->Draw();
      TPad* right= new TPad(Form("bR_%d",e),"",0.58,0.0,1.0,1.0); right->Draw();

      left->cd(); left->SetRightMargin(0.16); left->SetLeftMargin(0.12); left->SetBottomMargin(0.12); left->SetLogz();
      TH2D* h = (TH2D*)R[e].hBack->Clone();
      h->GetXaxis()->SetRangeUser(-4,4); h->GetYaxis()->SetRangeUser(-4,4);
      h->SetMinimum(zminB); h->SetMaximum(zmaxB);
      h->GetXaxis()->SetTitle("#Delta x [pix]");
      h->GetYaxis()->SetTitle("#Delta y [pix]");
      h->SetTitle("");
      h->Draw("COLZ");
      { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
        lab.DrawLatex(0.50,0.965,Form("Back-side Composite Cluster - %s - Ext %d", moduleUsed, e)); }

      right->cd(); right->SetGridx(); right->SetGridy(); right->SetLeftMargin(0.16); right->SetBottomMargin(0.16);
      if (R[e].projX && R[e].projY){
        TH1D* px=(TH1D*)R[e].projX->Clone(); TH1D* py=(TH1D*)R[e].projY->Clone();
        px->GetXaxis()->SetTitle("Pixel offset [pix]");
        px->GetYaxis()->SetTitle("Normalized counts");
        px->GetYaxis()->SetTitleOffset(2.3);
        px->SetMaximum(1.15*ymaxProj); px->SetMinimum(0.0);
        px->SetLineColor(kBlack);
        px->SetLineWidth(3);
        px->SetFillColorAlpha(kBlack, 0.08);
        px->SetMarkerStyle(20);
        px->SetMarkerColor(kBlack);
        py->SetLineColor(kBlue+1);
        py->SetLineWidth(3);
        py->SetFillColorAlpha(kBlue+1, 0.08);
        py->SetMarkerStyle(20);
        py->SetMarkerColor(kBlue+1);

        px->Draw("HIST"); py->Draw("HIST SAME");
        auto leg = new TLegend(0.45,0.8,0.89,0.88);
        leg->SetBorderSize(1);
        leg->SetFillStyle(1001);
        leg->SetTextFont(42);
        leg->SetTextSize(0.038);
        leg->AddEntry(px,"ProjX (Y: -3 #rightarrow -2 px)","l");
        leg->AddEntry(py,"ProjY (X: -3 #rightarrow -2 px)","l");
        leg->Draw();
      }
    }
    SaveCanvas(c, Form("%s/compare/compare_back_2x2", baseTag.Data()));
    delete c;
  }


  // ----------------- (2) Front-side comparison 2×2 ------------------
  {
    TCanvas* c = new TCanvas("compare_front","Front 2x2", 2800, 2200);
    c->Divide(2,2);
    int pad=1;
    for (int e=1; e<=4; ++e){
      if (!R[e].hFront) { ++pad; continue; }
      c->cd(pad++);
      TPad* left = new TPad(Form("fL_%d",e),"",0.0,0.0,0.57,1.0); left->Draw();
      TPad* right= new TPad(Form("fR_%d",e),"",0.57,0.0,1.0,1.0); right->Draw();

      // LEFT: heatmap + squares + title
      left->cd(); left->SetRightMargin(0.16); left->SetLeftMargin(0.12); left->SetBottomMargin(0.12); left->SetLogz();
      TH2D* h=(TH2D*)R[e].hFront->Clone();
      h->GetXaxis()->SetRangeUser(-4,4); h->GetYaxis()->SetRangeUser(-4,4);
      h->SetMinimum(zminF); h->SetMaximum(zmaxF);
      h->GetXaxis()->SetTitle("#Delta x [pix]");
      h->GetYaxis()->SetTitle("#Delta y [pix]");
      h->SetTitle("");
      h->Draw("COLZ");
      auto box=[&](double x1,double y1,double x2,double y2,Color_t c){ TBox*b=new TBox(x1,y1,x2,y2); b->SetLineColor(c); b->SetLineWidth(6); b->SetFillStyle(0); b->Draw("same"); return b; };
      TBox* bC = box(-0.5,-0.5,0.5,0.5,kBlack);
      TBox* bA = box(-0.5,1.5,0.5,2.5,kBlue+1);
      TBox* bB = box(-0.5,-1.5,0.5,-0.5,kRed+1);
      TBox* bR = box(1.5,-0.5,2.5,0.5,kGreen+2);
      TBox* bL = box(-1.5,-0.5,-0.5,0.5,kMagenta);
      { TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
        lab.DrawLatex(0.50,0.965,Form("Front-side Composite Cluster - %s - ext %d", moduleUsed, e)); }

      // RIGHT: bars + legend with colored boxes + ratios
      right->cd(); right->SetGridx(); right->SetGridy(); right->SetLeftMargin(0.16); right->SetBottomMargin(0.18);
      TH1F* frame=new TH1F(Form("frameF_cmp_%d",e),";Region;Fraction (w.r.t. center)",5,0.5,5.5);
      frame->GetXaxis()->SetBinLabel(1,"Above");
      frame->GetXaxis()->SetBinLabel(2,"Below");
      frame->GetXaxis()->SetBinLabel(3,"Right");
      frame->GetXaxis()->SetBinLabel(4,"Left");
      frame->GetXaxis()->SetBinLabel(5,"Center");
      double ymax = 1.25*std::max(1.0, std::max({R[e].fAbove,R[e].fBelow,R[e].fRight,R[e].fLeft}));
      frame->SetMinimum(0.0); frame->SetMaximum(ymax); frame->Draw("AXIS");
      TLine* ref=new TLine(0.5,1.0,5.5,1.0); ref->SetLineStyle(2); ref->SetLineColor(kGray+2); ref->Draw();
      auto drawBar=[&](int bin,double y,Color_t col){ double x1=frame->GetXaxis()->GetBinLowEdge(bin), x2=frame->GetXaxis()->GetBinUpEdge(bin);
        TBox* b=new TBox(x1,0.0,x2,y); b->SetFillColorAlpha(col,0.35); b->SetLineColor(col); b->SetLineWidth(3); b->Draw("same"); return b; };
      TBox* BA=drawBar(1,R[e].fAbove,kBlue+1);
      TBox* BB=drawBar(2,R[e].fBelow,kRed+1);
      TBox* BR=drawBar(3,R[e].fRight,kGreen+2);
      TBox* BL=drawBar(4,R[e].fLeft, kMagenta);
      TBox* BC=drawBar(5,1.0,       kGray+1);

      auto leg = new TLegend(0.18, 0.58, 0.55, 0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
      leg->AddEntry(BC, Form("Center charge = %.0f", R[e].centerCharge), "f");
      leg->AddEntry(BA, Form("Above  = %.3f", R[e].fAbove),  "f");
      leg->AddEntry(BB, Form("Below  = %.3f", R[e].fBelow),  "f");
      leg->AddEntry(BR, Form("Right  = %.3f", R[e].fRight),  "f");
      leg->AddEntry(BL, Form("Left   = %.3f", R[e].fLeft ),  "f");
      leg->Draw();
    }
    SaveCanvas(c, Form("%s/compare/compare_front_2x2", baseTag.Data()));
    delete c;
  }

  // ----------------- (3) Energy comparison 2×2 ----------------------
  {
        TCanvas* c = new TCanvas("compare_energy","Energy 2x2", 1400, 1000);
        c->Divide(2,2);
        int pad=1;
        for (int e=1;e<=4;++e){
        if (!R[e].hEnergy_front){ ++pad; continue; }
        c->cd(pad++); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

        TH1F* h=(TH1F*)R[e].hEnergy_front->Clone();
        h->SetTitle(";Energy [keV];Counts");
        h->SetFillColorAlpha(kGray+1,0.60); 
        h->SetLineColor(kGray+2); 
        h->SetBarWidth(0.95); 
        h->SetBarOffset(0.025);
        h->SetMaximum(1.15 * ymaxEnergy);
        h->GetYaxis()->SetNoExponent(true);   // avoid "×10^3" label
        h->GetYaxis()->SetMaxDigits(3);
        h->Draw("BAR");
        auto outline=(TH1F*)h->Clone(Form("hEout_cmp_%d",e)); outline->SetFillStyle(0); outline->SetLineColor(kBlack); outline->SetLineWidth(2); 
        outline->Draw("BAR SAME");

        // draw fitted components using stored params
        TF1* g1=new TF1(Form("cEg1_%d",e),"gaus",R[e].mu1-4*R[e].s1, R[e].mu1+4*R[e].s1);
        TF1* g2=new TF1(Form("cEg2_%d",e),"gaus",R[e].mu2-4*R[e].s2, R[e].mu2+4*R[e].s2);

        // fallback if a fit failed
        const double A1 = (R[e].A1>0 ? R[e].A1 : h->GetMaximum());
        const double A2 = (R[e].A2>0 ? R[e].A2 : 0.3*h->GetMaximum());

        g1->SetParameters(A1, R[e].mu1, std::max(0.03, R[e].s1));
        g2->SetParameters(A2, R[e].mu2, std::max(0.03, R[e].s2));
        g1->SetNpx(600); g2->SetNpx(600);
        g1->SetLineColor(kRed+1);  g1->SetLineWidth(3);
        g2->SetLineColor(kBlue+1); g2->SetLineWidth(3);

        g1->Draw("SAME");
        g2->Draw("SAME");

        auto leg=new TLegend(0.58,0.63,0.89,0.88);
        leg->AddEntry(g1,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu1,R[e].mu1_err,R[e].s1,R[e].s1_err),"l");
        leg->AddEntry(g2,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu2,R[e].mu2_err,R[e].s2,R[e].s2_err),"l");
        leg->AddEntry(h,"Data","f"); leg->Draw();

        TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
        lab.DrawLatex(0.50,0.965,Form("Front-side Energy Spectrum - %s - Ext %d", moduleUsed, e));
    }
    SaveCanvas(c, Form("%s/compare/compare_energy_2x2", baseTag.Data()));
    delete c;
  }

// ----------------- (4) Energy comparison Front and Back 2×2 ----------------------
  {
        TCanvas* c = new TCanvas("compare_energy","Energy 2x2", 1400, 1000);
        c->Divide(2,2);
        int pad=1;
        for (int e=1;e<=4;++e){
        if (!R[e].hEnergy_front){ ++pad; continue; }
        c->cd(pad++); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
        
        // ----------- FRONT HIST ---------
        TH1F* h_front=(TH1F*)R[e].hEnergy_front->Clone();
        h_front->SetTitle(";Energy [keV];Counts");
        h_front->SetFillColorAlpha(kRed+1,0.5);
        h_front->SetLineColor(kRed+2);
        h_front->SetBarWidth(0.95);
        h_front->SetBarOffset(0.05);
        h_front->SetMaximum(1.15 * ymaxEnergy);
        h_front->GetYaxis()->SetNoExponent(true);   // avoid "×10^3" label
        h_front->GetYaxis()->SetMaxDigits(3);
        h_front->Draw("BAR");

        auto frontOutline = (TH1F*)R[e].hEnergy_front->Clone("hE_outline_front");
        frontOutline->SetFillStyle(0); frontOutline->SetLineColor(kBlack);
        frontOutline->SetLineWidth(1);
        frontOutline->Draw("BAR SAME");

        // ----------- BACK HIST ---------
        TH1F* h_back=(TH1F*)R[e].hEnergy_back->Clone();
        h_back->SetTitle(";Energy [keV];Counts");
        h_back->SetFillColorAlpha(kBlue+1,0.5);
        h_back->SetLineColor(kBlue+2);
        h_back->SetBarWidth(0.95);
        h_back->SetBarOffset(0.05);
        h_back->SetMaximum(1.15 * ymaxEnergy);
        h_back->GetYaxis()->SetNoExponent(true);   // avoid "×10^3" label
        h_back->GetYaxis()->SetMaxDigits(3);
        h_back->Draw("BAR SAME");

        auto backOutline = (TH1F*)R[e].hEnergy_back->Clone("hE_outline_back");
        backOutline->SetFillStyle(0); backOutline->SetLineColor(kBlack);
        backOutline->SetLineWidth(1);
        backOutline->Draw("BAR SAME");


        // draw fitted components using stored params
        TF1* g1=new TF1(Form("cEg1_%d",e),"gaus", R[e].mu1-4*R[e].s1, R[e].mu1+4*R[e].s1);
        TF1* g2=new TF1(Form("cEg2_%d",e),"gaus", R[e].mu2-4*R[e].s2, R[e].mu2+4*R[e].s2);
        TF1* g3=new TF1(Form("cEg3_%d",e),"gaus", R[e].mu3-4*R[e].s3, R[e].mu3+4*R[e].s3);
        TF1* g4=new TF1(Form("cEg4_%d",e),"gaus", R[e].mu4-4*R[e].s4, R[e].mu4+4*R[e].s4);

        // fallback if a fit failed
        const double A1 = (R[e].A1>0 ? R[e].A1 : backOutline->GetMaximum());
        const double A2 = (R[e].A2>0 ? R[e].A2 : 0.3*backOutline->GetMaximum());
        const double A3 = (R[e].A3>0 ? R[e].A3 : backOutline->GetMaximum());
        const double A4 = (R[e].A4>0 ? R[e].A4 : 0.3*backOutline->GetMaximum());

        g1->SetParameters(R[e].A1, R[e].mu1, std::max(0.03, R[e].s1));
        g2->SetParameters(R[e].A2, R[e].mu2, std::max(0.03, R[e].s2));
        g1->SetNpx(600); g2->SetNpx(600);
        g3->SetParameters(R[e].A3, R[e].mu3, std::max(0.03, R[e].s3));
        g4->SetParameters(R[e].A4, R[e].mu4, std::max(0.03, R[e].s4));
        g3->SetNpx(600); g4->SetNpx(600);
   

        // ---- DRAW FITS (already existing objects: g1/g2/g3/g4) ----
        g1->SetLineColor(kRed+1);   g1->SetLineWidth(4); g1->SetLineStyle(1); g1->Draw("SAME");
        g2->SetLineColor(kRed+2);   g2->SetLineWidth(2); g2->SetLineStyle(1); g2->Draw("SAME");
        g3->SetLineColor(kBlue+1);  g3->SetLineWidth(4); g3->SetLineStyle(1); g3->Draw("SAME"); 
        g4->SetLineColor(kBlue+2);  g4->SetLineWidth(2); g4->SetLineStyle(1); g4->Draw("SAME");

        // ensure headroom for curves too
        // double ymaxPanel = 1.2 * std::max( { frontOutline->GetMaximum(), A1, A2, backOutline->GetMaximum(), A3, A4 } );
        // frontOutline->SetMaximum(ymaxPanel);

        g1->Draw("SAME");
        g2->Draw("SAME");
        g3->Draw("SAME");
        g4->Draw("SAME");

        auto leg_front = new TLegend(0.145, 0.65, 0.405, 0.89); // Left top position
        leg_front->SetBorderSize(1);
        leg_front->SetFillStyle(1001);
        leg_front->SetTextFont(42);
        leg_front->SetTextSize(0.030);
        leg_front->AddEntry(R[e].hEnergy_front, "#sigma_{xy}<1 (front)", "f");
        leg_front->AddEntry(g1, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu1, R[e].mu1_err, R[e].s1, R[e].s1_err), "l");
        leg_front->AddEntry(g2, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu2, R[e].mu2_err, R[e].s2, R[e].s2_err), "l");


        auto leg_back = new TLegend(0.63, 0.65, 0.89, 0.89); // Right top position
        leg_back->SetBorderSize(1);
        leg_back->SetFillStyle(1001);
        leg_back->SetTextFont(42);
        leg_back->SetTextSize(0.030);
        leg_back->AddEntry(R[e].hEnergy_back, "#sigma_{xy}>1 (back)", "f");
        leg_back->AddEntry(g3, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu3, R[e].mu3_err, R[e].s3, R[e].s3_err), "l");
        leg_back->AddEntry(g4, Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu4, R[e].mu4_err, R[e].s4, R[e].s4_err), "l");

      
        leg_front->Draw();
        leg_back->Draw();

        TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
        lab.DrawLatex(0.50,0.965,Form("Front and Back Energy Spectrum - %s - Ext %d", moduleUsed, e));
    }
    SaveCanvas(c, Form("%s/compare/compare_energy_front_back_2x2", baseTag.Data()));
    delete c;
  }
// ----------------- (5) Energy vs. Sigmaxy comparison 2×2 ------------------
  {
    TCanvas* c = new TCanvas("compare_front","Front 2x2", 2800, 2200);
    c->Divide(2,2);
    int pad=1;
    for (int e=1; e<=4; ++e){
      if (!R[e].henergy_sigmaxy) { ++pad; continue; }
      c->cd(pad++);  
      R[e].henergy_sigmaxy->SetTitle("");
      R[e].henergy_sigmaxy->GetXaxis()->SetTitle("Energy [keV]");
      R[e].henergy_sigmaxy->GetYaxis()->SetTitle("#sigma_{xy} [pix]");
      R[e].henergy_sigmaxy->Draw("COLZ");
      {TLatex lab; lab.SetNDC(); lab.SetTextAlign(22); lab.SetTextSize(0.05);
      lab.DrawLatex(0.50,0.965,Form("Cluster Energy vs #sigma_{xy} - %s - Ext %d", moduleUsed, e));}
    }
    SaveCanvas(c, Form("%s/compare/compare_E_Sxy_2x2", baseTag.Data()));
    delete c;
  }

  // -------------------- Summary log Energy --------------------
  std::ofstream energy_log(baseTag + "/energy_summary.csv");
  energy_log << "Extension,location,energypeak1,energypeak2\n";
  for (int e=1; e<=4; ++e){
    if (R[e].nFiles==0) continue;
    energy_log 
      << e << ","
      << 0 << ","
      << R[e].mu1 << " +/- " << R[e].mu1_err << ","
      << R[e].mu2 << " +/- " << R[e].mu2_err << "\n";
    energy_log 
      << e << ","
      << 1 << ","
      << R[e].mu3 << " +/- " << R[e].mu3_err << ","
      << R[e].mu4 << " +/- " << R[e].mu4_err << "\n";
  }
  energy_log.close();

  std::ofstream CTI_log(baseTag + "/CTI_summary.csv");
  CTI_log << "Extension,location,fraction,Mean,RMS,Skewness,Integral\n";
  for (int e=1; e<=4; ++e){
    if (R[e].nFiles==0) continue;
    CTI_log 
      << e << ","
      << 0 << ","
      << R[e].fRight << ","
      << R[e].projX_mean << ","
      << R[e].projX_rms << ","
      << R[e].projX_skew << ","
      << R[e].projX_int_e << "\n";
    CTI_log 
      << e << ","
      << 1 << ","
      << R[e].fAbove << ","
      << R[e].projY_mean << ","
      << R[e].projY_rms << ","
      << R[e].projY_skew << ","
      << R[e].projY_int_e << "\n";
  }
  CTI_log.close();

  // -------------------- main log --------------------
  std::ofstream main_log(baseTag + "/summary.log");
  main_log << "Summary for " << baseTag << "  (Module: " << moduleUsed << ")\n";
  main_log << "================================================================\n";
  for (int e=1; e<=4; ++e){
    if (R[e].nFiles==0) continue;
    main_log << "----------------------------------------------------------------\n";
    main_log << "Extension: ext" << e << "\n";
    main_log << "Files combined: " << R[e].nFiles << "\n";
    main_log << "Total clusters: " << R[e].totalClusters << "\n";
    main_log << "Accepted (energy window): " << R[e].nAcceptedEnergy << "\n";
    main_log << "Accepted (front-side):    " << R[e].nAcceptedFront  << "\n";
    main_log << "Accepted (back-side):     " << R[e].nAcceptedBack   << "\n";
    main_log << "Front integral [e-]: " << R[e].frontIntegral_e << "\n";
    main_log << "Back  integral [e-]: " << R[e].backIntegral_e  << "\n";
    main_log << "CTI fractions for front events (w.r.t. center)\n";
    main_log << "  center charge: " << R[e].centerCharge << "\n";
    main_log << "  above=" << R[e].fAbove
        << " below=" << R[e].fBelow
        << " right=" << R[e].fRight
        << " left="  << R[e].fLeft  << "\n";
    main_log << "Back-side projections (slabs at [-4,-3] px)\n";
    main_log << "  ProjY (X slab) integral [e-]: " << R[e].projY_int_e
        << " | mean=" << R[e].projY_mean
        << " RMS="    << R[e].projY_rms
        << " skew="   << R[e].projY_skew << "\n";
    main_log << "  ProjX (Y slab) integral [e-]: " << R[e].projX_int_e
        << " | mean=" << R[e].projX_mean
        << " RMS="    << R[e].projX_rms
        << " skew="   << R[e].projX_skew << "\n";
    main_log << "Energy fit (two separate Gaussians in per-peak windows)\n";
    main_log << "  mu1=" << R[e].mu1 << " +/- " << R[e].mu1_err
        << "  sigma1=" << R[e].s1 << " +/- " << R[e].s1_err << "\n";
    main_log << "  mu2=" << R[e].mu2 << " +/- " << R[e].mu2_err
        << "  sigma2=" << R[e].s2 << " +/- " << R[e].s2_err << "\n";
    main_log << "  chi2/ndf=" << R[e].chi2 << "/" << R[e].ndf
        << "  p=" << R[e].p << "\n";
  }
  main_log.close();

  std::cout << "[DONE] Outputs under: " << baseTag << "\n";
}
