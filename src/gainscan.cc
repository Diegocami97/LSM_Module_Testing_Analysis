// gainscan.C  — run as a ROOT macro (interpreted or ACLiC-compiled).
//
// Usage examples:
//   root -l -b -q 'gainscan.C(" /path/to/campaign_root ")'
//   root -l -b -q 'gainscan.C+(" /path/to/campaign_root ", "gains_long.csv", "plots")'
//
// Notes:
//  - Requires CFITSIO. The macro tries to load it dynamically.
//  - Reads per-extension keys: MEGAIN1 (gain), EMEGAIN1 (err), fallback GAIN_1.
//  - Parses module (e.g., DM07) and image label (Image_3_500, Image_3_1000, 5, 6, 7)
//    from the directory/file path via regex.
//  - Makes per-module and cross-module “plot” plots.
//
// If CFITSIO doesn't autoload, ensure libcfitsio is discoverable, or set:
//   gSystem->Load("/full/path/to/libcfitsio.dylib");  (macOS)
//   gSystem->Load("/full/path/to/libcfitsio.so");     (Linux)

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

#include <vector>
#include <string>
#include <set>
#include <map>
#include <regex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <limits>
#include <cmath>

// gSystem->AddIncludePath("-I/opt/homebrew/opt/cfitsio/include");
// gSystem->AddLinkedLibs("-L/opt/homebrew/opt/cfitsio/lib -lcfitsio");
#if __has_include(<fitsio.h>)
  #include <fitsio.h>
#elif __has_include(<cfitsio/fitsio.h>)
  #include <cfitsio/fitsio.h>
#else
  #error "fitsio.h not found. Provide CFITSIO include path at runtime (see run command)."
#endif


// #include <fitsio.h>



// // ---- CFITSIO setup (place this ABOVE any fitsio include) ----
// #ifdef __CLING__
// // Tell ACLiC/Cling where to find headers and which lib to load at runtime.
// // Change /opt/homebrew to your actual CFITSIO prefix if different.
// #pragma cling add_include_path("/opt/homebrew/include")
// // Safer to load the exact dylib path:
// #pragma cling load("/opt/homebrew/lib/libcfitsio.dylib")
// #endif

// // Include the header once (some installs put it under cfitsio/)
// #if __has_include(</opt/homebrew/opt/cfitsio/fitsio.h>)
//   #include <cfitsio/fitsio.h>
// #elif __has_include(<fitsio.h>)
//   #include </opt/homebrew/opt/cfitsio/fitsio.h>
// #else
//   #error "fitsio.h not found. Point to your CFITSIO include dir (brew --prefix cfitsio)."
// #endif





// ---------- CFITSIO helpers ----------
static void CFITSIO_CHECK(int status, const char* where) {
  if (status) {
    char msg[FLEN_STATUS]; msg[0] = '\0';
    fits_get_errstatus(status, msg);
    std::cerr << "CFITSIO error in " << where << ": " << msg << std::endl;
    gSystem->Exit(status);
  }
}

// forward declaration so collectFromList can call it
static void collectFits(const std::string& dir, std::vector<std::string>& out);


static bool read_key_double(fitsfile* f, const char* key, double& out) {
  int status = 0;
  char cmt[FLEN_COMMENT] = {0};
  fits_read_key(f, TDOUBLE, key, &out, cmt, &status);
  if (status) {
    status = 0;
    long tmp = 0;
    fits_read_key(f, TLONG, key, &tmp, cmt, &status);
    if (status) return false;
    out = (double)tmp;
  }
  return true;
}

// Read INT-like (keep your fixed read_key_int)

static bool read_indexed_sigma(fitsfile* f, int k, double& sig, double& err, std::string& src) {
  const std::string keyS = "MESIG" + std::to_string(k);
  const std::string keyE = "EMESIG" + std::to_string(k);
  if (read_key_double(f, keyS.c_str(), sig)) {
    src = keyS;
    double e; if (read_key_double(f, keyE.c_str(), e)) err = e; else err = 0.0;
    return true;
  }
  return false; // we DON'T synthesize sigma from anything else
}

// Read INT-like header keys (e.g., EXTVER)
static bool read_key_int(fitsfile* f, const char* key, int& out) {
  int status = 0;
  char cmt[FLEN_COMMENT] = {0};

  // Try as TINT
  if (fits_read_key(f, TINT, key, &out, cmt, &status) == 0) return true;

  // Fallback: TLONG
  status = 0;
  long v = 0;
  if (fits_read_key(f, TLONG, key, &v, cmt, &status) == 0) {
    out = static_cast<int>(v);
    return true;
  }
  return false;
}


// Read an indexed gain pair: MEGAINk / EMEGAINk (fallback to GAIN_k)
static bool read_indexed_gain(fitsfile* f, int k, double& gain, double& err, std::string& src) {
  const std::string keyG = "MEGAIN" + std::to_string(k);
  const std::string keyE = "EMEGAIN" + std::to_string(k);
  if (read_key_double(f, keyG.c_str(), gain)) {
    src = keyG;
    double e;
    if (read_key_double(f, keyE.c_str(), e)) err = e; else err = 0.0;
    return true;
  }
  const std::string keyG2 = "GAIN_" + std::to_string(k);
  if (read_key_double(f, keyG2.c_str(), gain)) {
    src = keyG2;
    err = 0.0;
    return true;
  }
  return false;
}


// // Read integer-like header keys (EXTVER, etc.)
// static bool read_key_int(fitsfile* f, const char* key, int& out) {
//   int status = 0;
//   char cmt[FLEN_COMMENT] = {0};
//   fits_read_key(f, TINT, key, &out, cmt, &status);
//   if (status) {
//     status = 0;
//     long v = 0;
//     fits_read_key(f, TLONG, key, &v, cmt, &status);
//     if (status) return false;
//     out = int(v);
//   }
//   return true;
// }



static std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b-1])) --b;
  return s.substr(a, b-a);
}

static bool read_key_string(fitsfile* f, const char* key, std::string& out) {
  int status = 0;
  char val[FLEN_VALUE] = {0};
  char cmt[FLEN_COMMENT] = {0};
  fits_read_key(f, TSTRING, key, val, cmt, &status);
  if (status) return false;
  std::string s = val;
  s.erase(std::remove(s.begin(), s.end(), '\''), s.end());
  out = trim(s);
  return true;
}

// Read a .txt list of paths (one per line). Lines may be:
//  - full path to a .fits/.fz file
//  - a directory (we’ll recurse inside it)
//  - blank or starting with '#' (ignored)
static void collectFromList(const std::string& listfile, std::vector<std::string>& out) {
  std::ifstream in(listfile);
  if (!in) {
    std::cerr << "Cannot open list file: " << listfile << "\n";
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    // strip inline comments and whitespace
    auto hash = line.find('#'); if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;

    TString p = line.c_str();
    if (gSystem->AccessPathName(p)) {
      std::cerr << "Path not found (from list): " << line << "\n";
      continue;
    }
    FileStat_t st; gSystem->GetPathInfo(p, st);
    if ((st.fMode & S_IFDIR) != 0) {
      collectFits(line, out);            // recurse into directory
    } else {
      TString low = p; low.ToLower();
      if (low.EndsWith(".fits") || low.EndsWith(".fz")) {
        out.emplace_back(line);          // single FITS/FZ file
      } else {
        std::cerr << "Skipping non-FITS path in list: " << line << "\n";
      }
    }
  }
}


// ---------- Records ----------
struct Rec {
  std::string module_id;   // e.g., DM07
  std::string image_label; // e.g., Image_3_500 / Image_3_1000 / Image_5 / Image_6 / Image_7
  int         skip_count = -1;
  int         ccd_id = -1; // 1..4
  std::string extname;
  double      gain = std::numeric_limits<double>::quiet_NaN();
  double      gain_err = 0.0;
  std::string gain_source; // MEGAIN1 or GAIN_1
  std::string file_base;
  std::string file_path;
  double      sigma_adu = std::numeric_limits<double>::quiet_NaN();
  double      sigma_err_adu = 0.0;
  std::string sigma_source; // MESIGk
};

// ---------- Parsing (module/image from path) ----------
struct ParsedPath {
  std::string module_id;
  std::string image_label;
  int skip_count = -1;
};

static ParsedPath parse_from_path(const std::string& path) {
  ParsedPath out;

  // Find "DMxx" in path components
  std::regex re_mod(R"(\/(DM\d{2})\/)");
  std::smatch mmod;
  if (std::regex_search(path, mmod, re_mod)) {
    out.module_id = mmod.str(1);
  }

  // Image_3_500 / Image_3_1000
  std::regex re_img_skips(R"(Image_(\d+)[^0-9]?(\d{3,4}))");
  std::regex re_img_only(R"(Image_(\d+))");
  std::smatch mimg;
  if (std::regex_search(path, mimg, re_img_skips)) {
    out.image_label = "Image_" + mimg.str(1) + "_" + mimg.str(2);
    out.skip_count  = std::stoi(mimg.str(2));
  } else if (std::regex_search(path, mimg, re_img_only)) {
    out.image_label = "Image_" + mimg.str(1);
  }
  return out;
}

// Try infer CCD id from FILEEXT, EXTNAME=EXT#, or filename "EXT#"
static int infer_ccd_id(const std::string& fileext_str,
                        const std::string& extname,
                        const std::string& fname,
                        int fallback_seq) {
  // FILEEXT (often "1", "2", "3", "4")
  if (!fileext_str.empty()) {
    try {
      int v = std::stoi(trim(fileext_str));
      if (v >= 1 && v <= 4) return v;
    } catch (...) {}
  }
  // EXTNAME like "EXT1"
  std::regex re_extn(R"(^EXT\s*([1-4]))", std::regex::icase);
  std::smatch m;
  if (std::regex_search(extname, m, re_extn)) {
    return std::stoi(m.str(1));
  }
  // Filename contains "EXT#"
  std::regex re_fn(R"(EXT([1-4]))");
  if (std::regex_search(fname, m, re_fn)) {
    return std::stoi(m.str(1));
  }
  // Fallback sequence within file 1..4
  return std::max(1, std::min(4, fallback_seq));
}


// Canonical image order (for jitter + legend order)
static const std::vector<std::string>& image_order() {
  static const std::vector<std::string> ord = {
    "Image_3_500", "Image_3_1000", "Image_5", "Image_6", "Image_7"
  };
  return ord;
}
static int image_rank(const std::string& label) {
  const auto& ord = image_order();
  auto it = std::find(ord.begin(), ord.end(), label);
  if (it == ord.end()) return (int)ord.size(); // unknown → last
  return int(it - ord.begin());
}
static int image_marker(int rank) {
  static int styles[] = {20, 21, 22, 23, 33, 34, 29, 25};
  if (rank < 0) rank = 0;
  if (rank >= (int)(sizeof(styles)/sizeof(int))) rank = (int)(sizeof(styles)/sizeof(int)) - 1;
  return styles[rank];
}
static double image_jitter(int rank) {
  static double off[] = {-0.20, -0.10, 0.0, +0.10, +0.20, +0.26, -0.26, 0.30};
  if (rank < 0) rank = 0;
  if (rank >= (int)(sizeof(off)/sizeof(double))) rank = 0;
  return off[rank];
}
static int module_color(const std::string& mod) {
  static std::vector<std::string> order;
  auto it = std::find(order.begin(), order.end(), mod);
  if (it == order.end()) { order.push_back(mod); it = std::prev(order.end()); }
  int idx = int(it - order.begin());
  static int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kGray+2};
  return colors[idx % (int)(sizeof(colors)/sizeof(int))];
}

// ---------- Directory recursion (collect .fits/.fz) ----------
static void collectFits(const std::string& dir, std::vector<std::string>& out) {
  TSystemDirectory sd("sd", dir.c_str());
  TList* list = sd.GetListOfFiles();
  if (!list) return;
  TIter next(list);
  while (TObject* obj = next()) {
    auto* f = dynamic_cast<TSystemFile*>(obj);
    if (!f) continue;
    TString name = f->GetName();
    if (name=="." || name=="..") continue;
    TString full = dir; if (!full.EndsWith("/")) full += "/"; full += name;
    if (f->IsDirectory()) {
      collectFits(full.Data(), out);
    } else {
      TString low = name; low.ToLower();
      if (low.EndsWith(".fits") || low.EndsWith(".fz")) {
        out.emplace_back(full.Data());
      }
    }
  }
}

// ---------- FITS scan per file ----------
// static void scan_file(const std::string& fpath, std::vector<Rec>& out) {
//   fitsfile* f = nullptr;
//   int status = 0;
//   fits_open_file(&f, fpath.c_str(), READONLY, &status);
//   if (status) { std::cerr << "Skip (open fail): " << fpath << "\n"; status = 0; return; }

//   int n_hdu = 0;
//   fits_get_num_hdus(f, &n_hdu, &status); CFITSIO_CHECK(status, "fits_get_num_hdus");

//   ParsedPath pp = parse_from_path(fpath);
//   int seq_ccd = 0; // fallback order (1..4) within a file for per-CCD case

//   for (int i = 1; i <= n_hdu; ++i) {
//     int hdutype = 0;
//     fits_movabs_hdu(f, i, &hdutype, &status); CFITSIO_CHECK(status, "fits_movabs_hdu");
//     if (hdutype != IMAGE_HDU) continue;

//     // Identify context
//     std::string extname, fileext; read_key_string(f, "EXTNAME", extname); read_key_string(f, "FILEEXT", fileext);
//     int extver = 0; (void)read_key_int(f, "EXTVER", extver);

//     // Probe which indexed gains exist in this header
//     int present_mg = 0;
//     bool present_idx[5] = {false,false,false,false,false};
//     for (int k=1;k<=4;++k) {
//       double tmp;
//       if (read_key_double(f, ("MEGAIN"+std::to_string(k)).c_str(), tmp) ||
//           read_key_double(f, ("GAIN_"+std::to_string(k)).c_str(),  tmp)) {
//         present_idx[k] = true; present_mg++;
//       }
//     }
//     bool fileext_has_list = (fileext.find(',') != std::string::npos);

//     // Branch: aggregated (Image_3) vs per-CCD (Image_5/6/7)
//     if (present_mg >= 2 || fileext_has_list) {
//       // --- Aggregated: emit 4 rows (k=1..4) using MEGAINk/EMEGAINk (or GAIN_k)
//       for (int k=1;k<=4;++k) {
//         double gain = std::numeric_limits<double>::quiet_NaN(), err = 0.0;
//         std::string src;
//         if (!read_indexed_gain(f, k, gain, err, src)) continue; // skip missing k
//         Rec r;
//         r.module_id   = pp.module_id.empty() ? "UNKNOWN" : pp.module_id;
//         r.image_label = pp.image_label.empty() ? "Image_?" : pp.image_label;
//         r.skip_count  = pp.skip_count;
//         r.ccd_id      = k;
//         r.extname     = extname;        // keep as-is (e.g., "CALIBRATED")
//         r.gain        = gain;
//         r.gain_err    = err;
//         r.gain_source = src;

//         size_t slash = fpath.find_last_of("/\\");
//         r.file_base = (slash == std::string::npos) ? fpath : fpath.substr(slash+1);
//         r.file_path = fpath;

//         out.push_back(std::move(r));
//       }
//       // Avoid duplicating if multiple IMAGE_HDUs appear (unlikely for your Image_3)
//       break;
//     } else {
//       // --- Per-CCD: pick the single indexed gain present (or MEGAIN1/GAIN_1)
//       int k_guess = 0;
//       for (int k=1;k<=4;++k) { if (present_idx[k]) { k_guess = k; break; } }

//       // Infer CCD id: prefer FILEEXT, then EXTVER/EXTNAME, then k_guess, then sequential
//       int ccd = -1;
//       // from FILEEXT if numeric
//       if (!fileext.empty()) {
//         try { int v = std::stoi(trim(fileext)); if (v>=1 && v<=4) ccd = v; } catch(...) {}
//       }
//       if (ccd < 0 && extver >= 1 && extver <= 4) ccd = extver;
//       if (ccd < 0) {
//         std::regex re_extn(R"(^EXT\s*([1-4]))", std::regex::icase);
//         std::smatch m; if (std::regex_search(extname, m, re_extn)) ccd = std::stoi(m.str(1));
//       }
//       if (ccd < 0 && k_guess >= 1 && k_guess <= 4) ccd = k_guess;
//       if (ccd < 0) ccd = std::max(1, std::min(4, ++seq_ccd));

//       // Read gain for that CCD index
//       double gain = std::numeric_limits<double>::quiet_NaN(), err = 0.0;
//       std::string src;
//       int k_to_read = (k_guess>=1 && k_guess<=4) ? k_guess : ccd;
//       if (!read_indexed_gain(f, k_to_read, gain, err, src)) {
//         // Conservative fallback
//         if (read_key_double(f, "MEGAIN1", gain)) { src = "MEGAIN1"; (void)read_key_double(f, "EMEGAIN1", err); }
//         else if (read_key_double(f, "GAIN_1", gain)) { src = "GAIN_1"; err = 0.0; }
//         else continue;
//       }

//       Rec r;
//       r.module_id   = pp.module_id.empty() ? "UNKNOWN" : pp.module_id;
//       r.image_label = pp.image_label.empty() ? "Image_?" : pp.image_label;
//       r.skip_count  = pp.skip_count;
//       r.ccd_id      = ccd;
//       r.extname     = extname;
//       r.gain        = gain;
//       r.gain_err    = err;
//       r.gain_source = src;

//       size_t slash = fpath.find_last_of("/\\");
//       r.file_base = (slash == std::string::npos) ? fpath : fpath.substr(slash+1);
//       r.file_path = fpath;

//       out.push_back(std::move(r));
//     }
//   }

//   fits_close_file(f, &status); CFITSIO_CHECK(status, "fits_close_file");
// }

static void scan_file(const std::string& fpath, std::vector<Rec>& out) {
  fitsfile* f = nullptr;
  int status = 0;
  fits_open_file(&f, fpath.c_str(), READONLY, &status);
  if (status) { std::cerr << "Skip (open fail): " << fpath << "\n"; status = 0; return; }

  int n_hdu = 0;
  fits_get_num_hdus(f, &n_hdu, &status); CFITSIO_CHECK(status, "fits_get_num_hdus");

  ParsedPath pp = parse_from_path(fpath);
  int seq_ccd = 0;

  for (int i = 1; i <= n_hdu; ++i) {
    int hdutype = 0;
    fits_movabs_hdu(f, i, &hdutype, &status); CFITSIO_CHECK(status, "fits_movabs_hdu");
    if (hdutype != IMAGE_HDU) continue;

    std::string extname, fileext; read_key_string(f, "EXTNAME", extname); read_key_string(f, "FILEEXT", fileext);
    int extver = 0; (void)read_key_int(f, "EXTVER", extver);

    // Probe for indexed keys present in this header
    int present_mg = 0, present_ms = 0;
    bool present_idx_gain[5] = {false,false,false,false,false};
    bool present_idx_sig [5]  = {false,false,false,false,false};
    for (int k=1;k<=4;++k) {
      double tmp;
      if (read_key_double(f, ("MEGAIN"+std::to_string(k)).c_str(), tmp) ||
          read_key_double(f, ("GAIN_"+std::to_string(k)).c_str(),  tmp)) {
        present_idx_gain[k] = true; present_mg++;
      }
      if (read_key_double(f, ("MESIG"+std::to_string(k)).c_str(), tmp)) {
        present_idx_sig[k] = true; present_ms++;
      }
    }
    bool fileext_has_list = (fileext.find(',') != std::string::npos);

    // -------- Aggregated header (Image_3_*): emit up to 4 rows
    if (present_mg >= 2 || present_ms >= 2 || fileext_has_list) {
      for (int k=1;k<=4;++k) {
        // Gain
        double gain = std::numeric_limits<double>::quiet_NaN(), gerr = 0.0; std::string gsrc;
        if (!(read_indexed_gain(f, k, gain, gerr, gsrc))) continue; // skip CCDs with no gain

        // Sigma (optional)
        double sig = std::numeric_limits<double>::quiet_NaN(), serr = 0.0; std::string ssrc;
        (void)read_indexed_sigma(f, k, sig, serr, ssrc); // ok if missing

        Rec r;
        r.module_id   = pp.module_id.empty() ? "UNKNOWN" : pp.module_id;
        r.image_label = pp.image_label.empty() ? "Image_?" : pp.image_label;
        r.skip_count  = pp.skip_count;
        r.ccd_id      = k;
        r.extname     = extname;
        r.gain        = gain;
        r.gain_err    = gerr;
        r.gain_source = gsrc;
        r.sigma_adu   = sig;
        r.sigma_err_adu = serr;
        r.sigma_source  = ssrc;

        size_t slash = fpath.find_last_of("/\\");
        r.file_base = (slash == std::string::npos) ? fpath : fpath.substr(slash+1);
        r.file_path = fpath;

        out.push_back(std::move(r));
      }
      break; // avoid duplicate emission if strange extra HDUs exist
    }

    // -------- Per-CCD header (Images 5/6/7): one row
    int k_guess = 0;
    for (int k=1;k<=4;++k) if (present_idx_gain[k] || present_idx_sig[k]) { k_guess = k; break; }

    // Infer CCD
    int ccd = -1;
    if (!fileext.empty()) { try { int v = std::stoi(trim(fileext)); if (v>=1 && v<=4) ccd = v; } catch(...) {} }
    if (ccd < 0 && extver >= 1 && extver <= 4) ccd = extver;
    if (ccd < 0) {
      std::regex re_extn(R"(^EXT\s*([1-4]))", std::regex::icase);
      std::smatch m; if (std::regex_search(extname, m, re_extn)) ccd = std::stoi(m.str(1));
    }
    if (ccd < 0 && k_guess >= 1 && k_guess <= 4) ccd = k_guess;
    if (ccd < 0) ccd = std::max(1, std::min(4, ++seq_ccd));

    // Read gain for k
    int k_to_read = (k_guess>=1 && k_guess<=4) ? k_guess : ccd;
    double gain = std::numeric_limits<double>::quiet_NaN(), gerr = 0.0; std::string gsrc;
    if (!read_indexed_gain(f, k_to_read, gain, gerr, gsrc)) {
      if (read_key_double(f, "MEGAIN1", gain)) { gsrc = "MEGAIN1"; (void)read_key_double(f, "EMEGAIN1", gerr); }
      else if (read_key_double(f, "GAIN_1", gain)) { gsrc = "GAIN_1"; gerr = 0.0; }
      else continue; // no usable gain
    }

    // Read sigma for k (optional)
    double sig = std::numeric_limits<double>::quiet_NaN(), serr = 0.0; std::string ssrc;
    (void)read_indexed_sigma(f, k_to_read, sig, serr, ssrc);

    Rec r;
    r.module_id   = pp.module_id.empty() ? "UNKNOWN" : pp.module_id;
    r.image_label = pp.image_label.empty() ? "Image_?" : pp.image_label;
    r.skip_count  = pp.skip_count;
    r.ccd_id      = ccd;
    r.extname     = extname;
    r.gain        = gain;
    r.gain_err    = gerr;
    r.gain_source = gsrc;
    r.sigma_adu   = sig;
    r.sigma_err_adu = serr;
    r.sigma_source  = ssrc;

    size_t slash = fpath.find_last_of("/\\");
    r.file_base = (slash == std::string::npos) ? fpath : fpath.substr(slash+1);
    r.file_path = fpath;
    out.push_back(std::move(r));
  }

  fits_close_file(f, &status); CFITSIO_CHECK(status, "fits_close_file");
}

// Simpler scan: only MEGAIN1/GAIN_1 from each IMAGE_HDU, no aggregation

static void scan_file_basic(const std::string& fpath, std::vector<Rec>& out) {
  fitsfile* f = nullptr;
  int status = 0;
  fits_open_file(&f, fpath.c_str(), READONLY, &status);
  if (status) {
    std::cerr << "Skip (open fail): " << fpath << "\n";
    status = 0; return;
  }
  int n_hdu = 0;
  fits_get_num_hdus(f, &n_hdu, &status); CFITSIO_CHECK(status, "fits_get_num_hdus");

  ParsedPath pp = parse_from_path(fpath);
  int seq_ccd = 0; // fallback order

  for (int i = 1; i <= n_hdu; ++i) {
    int hdutype = 0;
    fits_movabs_hdu(f, i, &hdutype, &status); CFITSIO_CHECK(status, "fits_movabs_hdu");
    if (hdutype != IMAGE_HDU) continue;

    double gain = 0.0, gain_err = 0.0;
    std::string source;
    if (read_key_double(f, "MEGAIN1", gain)) {
      source = "MEGAIN1";
      if (!read_key_double(f, "EMEGAIN1", gain_err)) gain_err = 0.0;
    } else if (read_key_double(f, "GAIN_1", gain)) {
      source = "GAIN_1";
      gain_err = 0.0;
    } else {
      continue; // no usable gain in this HDU
    }

    std::string extname, fileext;
    read_key_string(f, "EXTNAME", extname);
    read_key_string(f, "FILEEXT", fileext);

    int ccd = infer_ccd_id(fileext, extname, fpath, ++seq_ccd);

    Rec r;
    r.module_id   = pp.module_id.empty() ? "UNKNOWN" : pp.module_id;
    r.image_label = pp.image_label.empty() ? "Image_?" : pp.image_label;
    r.skip_count  = pp.skip_count;
    r.ccd_id      = ccd;
    r.extname     = extname;
    r.gain        = gain;
    r.gain_err    = gain_err;
    r.gain_source = source;

    // Basename
    size_t slash = fpath.find_last_of("/\\");
    r.file_base = (slash == std::string::npos) ? fpath : fpath.substr(slash+1);
    r.file_path = fpath;

    out.push_back(std::move(r));
  }

  fits_close_file(f, &status); CFITSIO_CHECK(status, "fits_close_file");
}

// ---------- CSV ----------
static void write_csv(const std::string& csv_path, const std::vector<Rec>& rows) {
  std::ofstream os(csv_path);
//   os << "module_id,image_label,skip_count,ccd_id,extname,gain_adupere,gain_err_adupere,gain_source,file_basename,file_path\n";
  os << "module_id,image_label,skip_count,ccd_id,extname,"
      "gain_adupere,gain_err_adupere,gain_source,"
      "file_basename,file_path,"
      "sigma_adu,sigma_err_adu,sigma_source\n";
  os << std::setprecision(16) << std::fixed;
  for (const auto& r : rows) {
    os << r.module_id << ','
       << r.image_label << ','
       << ((r.skip_count>=0)? std::to_string(r.skip_count) : "") << ','
       << r.ccd_id << ','
       << '"' << r.extname << '"' << ','
       << r.gain << ','
       << r.gain_err << ','
       << r.gain_source << ','
       << '"' << r.file_base << '"' << ','
       << '"' << r.file_path << '"' << ','
       << r.sigma_adu << ','
       << r.sigma_err_adu << ','
       << r.sigma_source << '\n';

    
  }
  std::cerr << "Wrote CSV: " << csv_path << "  (" << rows.size() << " rows)\n";
}

// ---------- Utilities ----------
static std::vector<std::string> unique_modules(const std::vector<Rec>& rows) {
  std::set<std::string> s; for (auto& r: rows) s.insert(r.module_id);
  return std::vector<std::string>(s.begin(), s.end());
}
static std::vector<std::string> unique_images(const std::vector<Rec>& rows) {
  std::set<std::string> s; for (auto& r: rows) s.insert(r.image_label);
  std::vector<std::string> v(s.begin(), s.end());
  std::stable_sort(v.begin(), v.end(), [](const std::string& a, const std::string& b){
    int ra = image_rank(a), rb = image_rank(b);
    if (ra != rb) return ra < rb;
    return a < b;
  });
  return v;
}
static TH1F* make_frame(double ymin, double ymax, const char* title) {
  static int counter = 0;
  auto* h = new TH1F(Form("frame_%d", counter++), "", 100, 0.5, 4.5);
  h->SetStats(0);
  h->SetTitle(title);

  // Axis titles
  h->GetXaxis()->SetTitle("CCD / Extension");
  h->GetYaxis()->SetTitle("Gain [ADU/e^{-}]");

  // Axis ranges
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);

  // Ticks & labels
  h->GetXaxis()->SetNdivisions(404, false);
  h->LabelsOption("h", "X");                 // 1 = X axis (2=Y, 3=Z)
  h->GetXaxis()->SetBit(TAxis::kLabelsHori);
//   h->GetXaxis()->ClearBit(TAxis::kLabelsVert);         // force horizontal tick labels
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.40);       // avoid overlap with ticks
  h->GetYaxis()->SetTitleOffset(1.40);

  // Bin labels
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(1.0), "EXT1");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(2.0), "EXT2");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(3.0), "EXT3");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(4.0), "EXT4");
  return h;
}

// Frame for sigma (same style; different Y title)
static TH1F* make_frame_sigma(double ymin, double ymax, const char* title) {
  static int counter = 0;
  auto* h = new TH1F(Form("frame_sigma_%d", counter++), "", 100, 0.5, 4.5);
  h->SetStats(0);
  h->SetTitle(title);
  h->GetXaxis()->SetTitle("CCD / Extension");
  h->GetYaxis()->SetTitle("Single-e Resolution [ADU]");
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetXaxis()->SetNdivisions(404, false);
  h->LabelsOption("h", "X");
  h->GetXaxis()->SetBit(TAxis::kLabelsHori);
//   h->GetXaxis()->ClearBit(TAxis::kLabelsVert);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.30);
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(1.0), "EXT1");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(2.0), "EXT2");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(3.0), "EXT3");
  h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(4.0), "EXT4");
  return h;
}

static TGraphErrors* make_graph_sigma(const std::vector<Rec>& rows,
                                      const std::string& module_id,
                                      const std::string& image_label,
                                      int marker_style, int color) {
  std::vector<double> xs, ys, exs, eys;
  int rank = image_rank(image_label);
  for (const auto& r : rows) {
    if (r.module_id != module_id) continue;
    if (r.image_label != image_label) continue;
    if (!(r.ccd_id >= 1 && r.ccd_id <= 4)) continue;
    if (!std::isfinite(r.sigma_adu)) continue;
    xs.push_back(r.ccd_id + image_jitter(rank));
    ys.push_back(r.sigma_adu);
    exs.push_back(0.0);
    eys.push_back(r.sigma_err_adu);
  }
  if (xs.empty()) return nullptr;
  auto* gr = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
  gr->SetMarkerStyle(marker_style);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetMarkerSize(1.1);
  return gr;
}

static std::pair<double,double> y_range_sigma(const std::vector<Rec>& rows) {
  double lo =  1e30, hi = -1e30;
  for (auto& r: rows) if (std::isfinite(r.sigma_adu)) { lo = std::min(lo, r.sigma_adu); hi = std::max(hi, r.sigma_adu); }
  if (lo > hi) { lo = 0.9; hi = 1.1; } // harmless default
  double pad = 0.04 * std::max(1.0, hi - lo);
  return {lo - pad, hi + pad};
}

static std::pair<double,double> y_range(const std::vector<Rec>& rows) {
  double lo =  1e30, hi = -1e30;
  for (auto& r: rows) if (std::isfinite(r.gain)) { lo = std::min(lo, r.gain); hi = std::max(hi, r.gain); }
  if (lo > hi) { lo = 0.9; hi = 1.1; }
  double pad = 0.02 * std::max(1.0, hi - lo);
  return {lo - pad, hi + pad};
}
static TGraphErrors* make_graph(const std::vector<Rec>& rows,
                                const std::string& module_id,
                                const std::string& image_label,
                                int marker_style, int color) {
  std::vector<double> xs, ys, exs, eys;
  int rank = image_rank(image_label);
  for (const auto& r : rows) {
    if (r.module_id != module_id) continue;
    if (r.image_label != image_label) continue;
    if (!(r.ccd_id >= 1 && r.ccd_id <= 4)) continue;
    xs.push_back(r.ccd_id + image_jitter(rank));
    ys.push_back(r.gain);
    exs.push_back(0.0);
    eys.push_back(r.gain_err);
  }
  if (xs.empty()) return nullptr;
  auto* gr = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
  gr->SetMarkerStyle(marker_style);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetMarkerSize(1.1);
  return gr;
}

// ---------- MAIN ENTRY ----------
void gainscan(const char* in_path = ".", const char* out_csv = "gains_long.csv", const char* outdir = "plots") {
  // Try to load CFITSIO
  if (gSystem->Load("cfitsio") < 0) {
    if (gSystem->Load("libcfitsio") < 0) {
      std::cerr << "Warning: could not load CFITSIO by name. If this fails later, "
                   "load your lib explicitly (gSystem->Load(\"/full/path/libcfitsio.so\"))\n";
    }
  }

  gStyle->SetOptStat(0);
  gROOT->SetBatch(kFALSE);
  gSystem->mkdir(outdir, kTRUE);

  // Gather files
    std::vector<std::string> files;
    {
    TString p = in_path;
    if (gSystem->AccessPathName(p)) {
        std::cerr << "Input path not found: " << p << "\n";
    } else {
        FileStat_t st; gSystem->GetPathInfo(p, st);
        if ((st.fMode & S_IFDIR) != 0) {
        collectFits(p.Data(), files);                           // recurse directory
        } else {
        TString low = p; low.ToLower();
        if (low.EndsWith(".txt")) {
            collectFromList(p.Data(), files);                     // read list of paths
        } else if (low.EndsWith(".fits") || low.EndsWith(".fz")) {
            files.emplace_back(p.Data());                         // single FITS/FZ file
        } else {
            std::cerr << "Unsupported input: " << p
                    << " (expect a directory, .txt list, or .fits/.fz)\n";
        }
        }
    }
    }
    if (files.empty()) {
    std::cerr << "No .fits/.fz files found from input: " << in_path << "\n";
    return;
    }


  // Scan
  std::vector<Rec> rows;
  rows.reserve(files.size()*4);
  for (const auto& f : files) scan_file(f, rows);

  if (rows.empty()) {
    std::cerr << "No usable HDUs with MEGAIN1/GAIN_1 found.\n";
    return;
  }

  // Sort (module, image_rank, image_label, ccd)
  std::sort(rows.begin(), rows.end(), [](const Rec& a, const Rec& b){
    if (a.module_id != b.module_id) return a.module_id < b.module_id;
    int ra = image_rank(a.image_label), rb = image_rank(b.image_label);
    if (ra != rb) return ra < rb;
    if (a.image_label != b.image_label) return a.image_label < b.image_label;
    return a.ccd_id < b.ccd_id;
  });

  write_csv(out_csv, rows);

  // Build module/image lists + Y-range
  auto modules = unique_modules(rows);
  auto images  = unique_images(rows);
  auto yr = y_range(rows);
  double ymin = yr.first, ymax = yr.second;

  // Per-module plots
  for (const auto& mod : modules) {
    TCanvas c(("c_"+mod).c_str(), ("plot "+mod).c_str(), 900, 650);
    auto* frame = make_frame(ymin, ymax, (mod + " : Gain vs CCD").c_str());
    frame->Draw();

   // Comfortable margins so labels/titles don't collide
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.08);

    // Legend INSIDE the frame (top-right)
    TLegend leg(0.70, 0.70, 0.93, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.033);

    for (size_t i = 0; i < images.size(); ++i) {
      const auto& label = images[i];
      auto* gr = make_graph(rows, mod, label, image_marker((int)i), kBlack);
      if (!gr) continue;
      gr->Draw("P SAME");
      leg.AddEntry(gr, label.c_str(), "p");
    }
    leg.Draw();

    TString pdir = outdir; TString png = pdir + "/" + mod.c_str() + "_plot.png";
    TString pdf  = pdir + "/" + mod.c_str() + "_plot.pdf";
    c.SaveAs(png); c.SaveAs(pdf);
    }

  // All-modules grand plot (color by module, marker by image)
{
    TCanvas c("c_all", "All modules", 1000, 700);
    auto* frame = make_frame(ymin, ymax, "All modules : Gain vs CCD");
    frame->Draw();

    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.08);

    // Module colors legend (top-left, inside)
    TLegend leg_mod(0.14, 0.70, 0.31, 0.88);
    leg_mod.SetBorderSize(0);
    leg_mod.SetFillStyle(0);
    leg_mod.SetTextSize(0.033);

    // Image markers legend (top-right, inside)
    TLegend leg_img(0.70, 0.70, 0.93, 0.88);
    leg_img.SetBorderSize(0);
    leg_img.SetFillStyle(0);
    leg_img.SetTextSize(0.033);

    // TLegend leg_mod(0.10, 0.72, 0.28, 0.92); leg_mod.SetBorderSize(0); leg_mod.SetFillStyle(0);
    // TLegend leg_img(0.72, 0.72, 0.92, 0.92); leg_img.SetBorderSize(0); leg_img.SetFillStyle(0);

    // Build once so legends don't duplicate
    bool added_mod_legend = false;
    bool added_img_legend = false;

    // Draw module × image graphs
    for (size_t im = 0; im < images.size(); ++im) {
      const auto& label = images[im];
      int mstyle = image_marker((int)im);
      for (const auto& mod : modules) {
        int col = module_color(mod);
        auto* gr = make_graph(rows, mod, label, mstyle, col);
        if (!gr) continue;
        gr->Draw("P SAME");

        if (!added_mod_legend) {
          // Add module legend entries once
          auto* dm = new TGraphErrors(); dm->SetMarkerStyle(20); dm->SetMarkerColor(col); dm->SetLineColor(col);
          leg_mod.AddEntry(dm, mod.c_str(), "p");
        }
      }
      // Add image legend (marker shapes) once
      if (!added_img_legend) {
        auto* di = new TGraphErrors(); di->SetMarkerStyle(mstyle); di->SetMarkerColor(kBlack);
        leg_img.AddEntry(di, label.c_str(), "p");
      }
      added_mod_legend = true;
    }

    leg_mod.Draw(); leg_img.Draw();

    TString pdir = outdir; TString png = pdir + "/ALL_modules_plot.png";
    TString pdf  = pdir + "/ALL_modules_plot.pdf";
    c.SaveAs(png); c.SaveAs(pdf);
}

    auto yrS  = y_range_sigma(rows);
    double yminS = yrS.first;
    double ymaxS = yrS.second;

    // Per-module sigma plots
    for (const auto& mod : modules) {
        TCanvas cS(("c_sigma_"+mod).c_str(), ("Sigma "+mod).c_str(), 900, 650);
        auto* frameS = make_frame_sigma(yminS, ymaxS, (mod + " : Single-e Resolution vs CCD").c_str());
        frameS->Draw();

        // margins
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.08);

        // Legend INSIDE the frame (top-right) — same box as gain
        TLegend legS(0.70, 0.70, 0.93, 0.88);
        legS.SetBorderSize(0);
        legS.SetFillStyle(0);
        legS.SetTextSize(0.033);

        for (size_t i = 0; i < images.size(); ++i) {
            const auto& label = images[i];
            auto* gr = make_graph_sigma(rows, mod, label, image_marker((int)i), kBlack);
            if (!gr) continue;
            gr->Draw("P SAME");
            legS.AddEntry(gr, label.c_str(), "p");
        }
        legS.Draw();

        TString pdir = outdir; TString png = pdir + "/" + mod.c_str() + "_sigma_plot.png";
        TString pdf  = pdir + "/" + mod.c_str() + "_sigma_plot.pdf";
        cS.SaveAs(png); cS.SaveAs(pdf);
    }

    // All-modules sigma plot (color by module, marker by image)
    {
    TCanvas cS("c_all_sigma", "All modules (sigma)", 1000, 700);
    auto* frameS = make_frame_sigma(yminS, ymaxS, "All modules : Single-e Resolution vs CCD");
    frameS->Draw();

    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.08);

    // Module colors legend (top-left, inside)
    TLegend leg_modS(0.14, 0.70, 0.31, 0.88);
    leg_modS.SetBorderSize(0);
    leg_modS.SetFillStyle(0);
    leg_modS.SetTextSize(0.033);

    // Image markers legend (top-right, inside)
    TLegend leg_imgS(0.70, 0.70, 0.93, 0.88);
    leg_imgS.SetBorderSize(0);
    leg_imgS.SetFillStyle(0);
    leg_imgS.SetTextSize(0.033);

    // Draw module × image graphs
    for (size_t im = 0; im < images.size(); ++im) {
        const auto& label = images[im];
        int mstyle = image_marker((int)im);
        for (const auto& mod : modules) {
        int col = module_color(mod);
        if (auto* gr = make_graph_sigma(rows, mod, label, mstyle, col)) gr->Draw("P SAME");
        }
    }

    // Build legends AFTER drawing (one entry per module / per image)
    for (const auto& mod : modules) {
        int col = module_color(mod);
        auto* dm = new TGraphErrors(); dm->SetMarkerStyle(20); dm->SetMarkerColor(col); dm->SetLineColor(col);
        leg_modS.AddEntry(dm, mod.c_str(), "p");
    }
    for (size_t im = 0; im < images.size(); ++im) {
        int mstyle = image_marker((int)im);
        auto* di = new TGraphErrors(); di->SetMarkerStyle(mstyle); di->SetMarkerColor(kBlack); di->SetLineColor(kBlack);
        leg_imgS.AddEntry(di, images[im].c_str(), "p");
    }

    leg_modS.Draw(); leg_imgS.Draw();

    TString pdir = outdir; TString png = pdir + "/ALL_modules_sigma_plot.png";
    TString pdf  = pdir + "/ALL_modules_sigma_plot.pdf";
    cS.SaveAs(png); cS.SaveAs(pdf);
    }




  std::cerr << "Done. CSV: " << out_csv << "  Plots dir: " << outdir << "\n";
}
