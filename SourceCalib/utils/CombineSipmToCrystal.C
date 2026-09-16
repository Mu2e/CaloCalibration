// Combine per-SiPM ADC histograms into per-crystal ADC histograms.
//
// Each crystal is read out by two SiPMs (an even/odd pair): sipm_2*i and
// sipm_2*i+1 are summed bin-by-bin into crystal i. Reads
// <top>/sipm_ADC/sipm_<n> from the input file and writes
// <top>/crystals_ADC/cry_<i> into the output file, preserving the on-disk
// directory layout used elsewhere in this package (see
// MakeAnalysisTree_main.cc's "SourceAna/crystals_ADC/cry_" convention).
//
// Usage:
//   root -l -b -q 'CombineSipmToCrystal.C("/exp/mu2e/data/users/hjafree/new_seed.root")'
//   root -l -b -q 'CombineSipmToCrystal.C("in.root", "out.root", "sourceana")'

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <regex>
#include <map>

// Find a subdirectory of `parent` whose name matches `name` case-insensitively.
static TDirectory* FindDirCI(TDirectory* parent, const char* name) {
    TString target(name);
    target.ToLower();
    for (TObject* obj : *parent->GetListOfKeys()) {
        TKey* key = (TKey*)obj;
        TString cname(key->GetClassName());
        if (!cname.Contains("TDirectory")) continue;
        TString kname(key->GetName());
        TString lname(kname);
        lname.ToLower();
        if (lname == target) {
            return (TDirectory*)parent->Get(kname);
        }
    }
    return nullptr;
}

void CombineSipmToCrystal(const char* inputPath  = "/exp/mu2e/data/users/hjafree/new_seed.root",
                           const char* outputPath = "",
                           const char* topDirName  = "sourceana") {
    TString outPath(outputPath);
    if (outPath.IsNull()) {
        TString in(inputPath);
        Ssiz_t dot = in.Last('.');
        outPath = (dot == kNPOS) ? in + "_crystals" : in(0, dot) + "_crystals" + in(dot, in.Length() - dot);
    }

    TFile* fin = TFile::Open(inputPath, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: could not open input file " << inputPath << std::endl;
        return;
    }

    TDirectory* topIn = FindDirCI(fin, topDirName);
    if (!topIn) {
        std::cerr << "ERROR: could not find top-level directory \"" << topDirName
                   << "\" in " << inputPath << std::endl;
        return;
    }
    TDirectory* sipmDir = FindDirCI(topIn, "sipm_ADC");
    if (!sipmDir) {
        std::cerr << "ERROR: could not find \"sipm_ADC\" directory under \""
                   << topDirName << "\"" << std::endl;
        return;
    }

    // Discover all sipm_<n> histograms and their indices.
    std::map<int, TString> sipmNames;
    std::regex sipmRegex("^sipm_([0-9]+)$");
    for (TObject* obj : *sipmDir->GetListOfKeys()) {
        TKey* key = (TKey*)obj;
        std::string kname(key->GetName());
        std::smatch m;
        if (std::regex_match(kname, m, sipmRegex)) {
            sipmNames[std::stoi(m[1].str())] = kname;
        }
    }
    if (sipmNames.empty()) {
        std::cerr << "ERROR: no sipm_<n> histograms found under \"" << topDirName
                   << "/sipm_ADC\"" << std::endl;
        return;
    }

    int nSipms = sipmNames.rbegin()->first + 1; // assumes contiguous 0..N-1
    int nCrystals = nSipms / 2;
    std::cout << "Found " << sipmNames.size() << " sipm histograms (max index "
              << sipmNames.rbegin()->first << ") -> building " << nCrystals
              << " crystal histograms." << std::endl;
    if (nSipms % 2 != 0) {
        std::cerr << "WARNING: number of sipm indices (0.." << sipmNames.rbegin()->first
                   << ") is odd; the last sipm (" << sipmNames.rbegin()->first
                   << ") will be dropped." << std::endl;
    }

    TFile* fout = new TFile(outPath, "RECREATE");
    TDirectory* topOut = fout->mkdir(topDirName);
    TDirectory* cryDir = topOut->mkdir("crystals_ADC");

    int nWritten = 0;
    for (int c = 0; c < nCrystals; ++c) {
        int i0 = 2 * c;
        int i1 = 2 * c + 1;
        if (!sipmNames.count(i0) || !sipmNames.count(i1)) {
            std::cerr << "WARNING: missing sipm_" << i0 << " or sipm_" << i1
                       << "; skipping crystal " << c << std::endl;
            continue;
        }
        TH1F* h0 = (TH1F*)sipmDir->Get(sipmNames[i0]);
        TH1F* h1 = (TH1F*)sipmDir->Get(sipmNames[i1]);
        if (!h0 || !h1) {
            std::cerr << "WARNING: could not read sipm_" << i0 << " or sipm_" << i1
                       << "; skipping crystal " << c << std::endl;
            continue;
        }

        TString cryName = Form("cry_%d", c);
        TH1F* hCry = (TH1F*)h0->Clone(cryName);
        hCry->SetTitle(cryName);
        hCry->SetDirectory(nullptr);
        hCry->Add(h1);

        cryDir->cd();
        hCry->Write(cryName);
        delete hCry;
        ++nWritten;
    }

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Wrote " << nWritten << " crystal histograms to " << outPath
               << " under \"" << topDirName << "/crystals_ADC\"" << std::endl;
}