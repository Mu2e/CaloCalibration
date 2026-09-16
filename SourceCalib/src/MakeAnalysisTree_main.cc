#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include "CaloCalibration/SourceCalib/inc/2dcontour.hh"
#include "CaloCalibration/SourceCalib/inc/mcinfo.hh"// extracting mc information -- need to be removed for real data
#include <algorithm>

#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/scratch/users/hjafree/fixedgeom_both_disks.root";//change file here


std::pair<TH1F*, TFile*> get_data_histogram(int cryNum, bool isSiPMRun) {
    TString histPath = isSiPMRun ? "SourceAna/sipm_ADC/sipm_" : "SourceAna/crystals_ADC/cry_";
    //uncomment line below and comment the one above if cut and count needs to be applied to the mc truth histograms
    //histPath = isSiPMRun ? "SourceAna/crystals_edep_truth/cry_" : "SourceAna/crystals_ADC/cry_";

    TFile *f = new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get(histPath + crystalNumber);
		hist->SetDirectory(0);
    return std::make_pair(hist, f);
}
std::pair<double, double> ComputeHistogramStats(TH1F* hist) {
    double sum = 0;
    double weightedSum = 0;
    double weightedSumSq = 0;

    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double content = hist->GetBinContent(i);
        double center = hist->GetBinCenter(i);

        sum += content;
        weightedSum += content * center;
        weightedSumSq += content * center * center;
    }

    double mean = (sum > 0) ? weightedSum / sum : 0;
    double variance = (sum > 0) ? (weightedSumSq / sum) - (mean * mean) : 0;
    double stddev = (variance > 0) ? std::sqrt(variance) : 0;

    return std::make_pair(mean, stddev);
}


int main(int argc, char* argv[]) {
    std::cout << "========== Welcome to the Mu2e Source Calibration Analysis ==========" << std::endl;
    std::vector<TString> allowedParams = {
    	"peak" , "alpha", "n_full" , "n_1st", "n_2nd","beta"
    }; 
    auto isInvalid = [&](TString input){
    	input.ToLower();
    	for (const auto& p : allowedParams) {
    		if (input == p) return false;
    	}
    	return true;
    };

    if (argc < 6) {
        std::cerr << "[ERROR] Missing arguments.\n";
        std::cerr << "Usage: " << argv[0] << " <start_cry> <end_cry> <alg> <disk> <readout: sipm|crystal> [flags: overlay, contour, mc, singlesided]\n";
        return 1;
    }

    // 2. Parse Fixed Positional Arguments
    int anacrys_start = std::atoi(argv[1]);
    int anacrys_end   = std::atoi(argv[2]);
    TString alg       = argv[3];
    int disk          = std::atoi(argv[4]);

    TString readout = argv[5];
    readout.ToLower();
    if (readout != "sipm" && readout != "crystal") {
        std::cerr << "[ERROR] Invalid readout '" << argv[5] << "' - must be 'sipm' or 'crystal'\n";
        return 1;
    }
    bool isSiPMRun = (readout == "sipm");

    // 3. Initialize Flags
    bool doOverlay = false;
    bool contour   = false;
    TString xSelect = "Peak"; 
    TString ySelect = "Width"; 
    bool isMC      = false;
    for (int i = 6; i < argc; i++) {
        TString arg = argv[i];
        arg.ToLower(); 
        
        if (arg.Contains("overlay")) doOverlay = true;
        if (arg.Contains("contour")){
        	contour   = true;
        	if (i+2 < argc){
        		xSelect = argv[i+1];
        		ySelect = argv[i+2];
        			if (isInvalid(xSelect)||isInvalid(ySelect)){
        				std::cerr<<"\n[ERROR] Invalid contour variables "<<std::endl;
        				std::cerr << "Allowed options are: \n";
        				for (auto &p : allowedParams) std::cerr<< " - " << p << "\n";
        				std::cerr <<std::endl;
        				return 1;
        			}
        		i += 2;
        	} else {
        		std::cerr<< "[ERROR] contour flag requires two variables (eg Peak, Width, Alpha, N_Full, N_1st, N_2nd)"<<std::endl;
        		}
        }
        
        if (arg.Contains("mc"))      isMC      = true;  
        if(arg.Contains("singlesided")) SourceFitter::singlesided = true;
        if (arg.Contains("help")) {
            std::cout << "Usage: " << argv[0] << " <start_cry> <end_cry> <alg> <disk> <readout: sipm|crystal> [flags: overlay, contour, mc, singlesided]\n";
            std::cout << "Example: " << argv[0] << " 0 10 fit 0 sipm overlay\n";
            return 0;
        }
         	
    }
//------------------
//Caphri crystals
//------------------
std::vector<int> lysoSiPMs = {1164,1165,1220,1221,1218,1219,1274,1275};
std::vector<int> lysoCrystals = {582,610,609,637};
std::cout << "\n[INFO] Starting loop from ID" << anacrys_start << " to " << anacrys_end <<std::endl;
    
// --------------------------------------------------------
// MC TRUTH BLOCK (Cut-and-Count + Fitting Comparison)
// --------------------------------------------------------
if (isMC) {
    std::cout << "[INFO] Running in MC mode, using MC file paths.\n";
    std::cout << "=======================================\n";
    std::cout << "  Running MC Truth Analysis (Cut-Count + Fit)\n";
    std::cout << "=======================================\n";

    TFile *truthtable = new TFile("truthinfo.root", "RECREATE");

    // --- Cut-and-Count Tree ---
    TTree *trueinfo = new TTree("trueinfo", "MC Truth Information");
    Int_t cryNumparam, tot_evts, mainpeak, first_espeak, second_espeak, background;
    Int_t compton1, compton2, compton3;
    Float_t frmainpeak, frfirst_espeak, frsecond_espeak, frbackground;
    Float_t frcompton1, frcompton2, frcompton3;
    trueinfo->Branch("cryNumparam",     &cryNumparam,     "cryNumparam/I");
    trueinfo->Branch("tot_evts",        &tot_evts,        "tot_evts/I");
    trueinfo->Branch("mainpeak",        &mainpeak,        "mainpeak/I");
    trueinfo->Branch("first_espeak",    &first_espeak,    "first_espeak/I");
    trueinfo->Branch("second_espeak",   &second_espeak,   "second_espeak/I");
    trueinfo->Branch("background",      &background,      "background/I");
    trueinfo->Branch("compton1",        &compton1,        "compton1/I");
    trueinfo->Branch("compton2",        &compton2,        "compton2/I");
    trueinfo->Branch("compton3",        &compton3,        "compton3/I");
    trueinfo->Branch("frmainpeak",      &frmainpeak,      "frmainpeak/F");
    trueinfo->Branch("frfirst_espeak",  &frfirst_espeak,  "frfirst_espeak/F");
    trueinfo->Branch("frsecond_espeak", &frsecond_espeak, "frsecond_espeak/F");
    trueinfo->Branch("frbackground",    &frbackground,    "frbackground/F");
    trueinfo->Branch("frcompton1",      &frcompton1,      "frcompton1/F");
    trueinfo->Branch("frcompton2",      &frcompton2,      "frcompton2/F");
    trueinfo->Branch("frcompton3",      &frcompton3,      "frcompton3/F");

    // --- RooFit Tree for MC Truth Spectrum Fits ---
    TTree *truthfit = new TTree("truthfit", "MC Truth Spectrum Fits (RooFit)");
    Int_t truth_cryNum, truth_ndof;
    Float_t truth_peak_mev, truth_fr_full, truth_fr_1st, truth_fr_2nd;
    Float_t truth_fr_c1, truth_fr_c2, truth_fr_c3;
    Float_t truth_width, truth_alpha, truth_chi2;
    Float_t truth_fr_ebk;
    truthfit->Branch("crystalNo",  &truth_cryNum,   "truth_cryNum/I");
    truthfit->Branch("peak_mev",   &truth_peak_mev, "truth_peak_mev/F");
    truthfit->Branch("fr_full",    &truth_fr_full,  "truth_fr_full/F");
    truthfit->Branch("fr_1st",     &truth_fr_1st,   "truth_fr_1st/F");
    truthfit->Branch("fr_2nd",     &truth_fr_2nd,   "truth_fr_2nd/F");
    truthfit->Branch("fr_c1",      &truth_fr_c1,    "truth_fr_c1/F");
    truthfit->Branch("fr_c2",      &truth_fr_c2,    "truth_fr_c2/F");
    truthfit->Branch("fr_c3",      &truth_fr_c3,    "truth_fr_c3/F");
    truthfit->Branch("fr_ebk",      &truth_fr_ebk,    "truth_fr_ebk/F");
    truthfit->Branch("width",      &truth_width,    "truth_width/F");
    truthfit->Branch("alpha",      &truth_alpha,    "truth_alpha/F");
    truthfit->Branch("chi2",       &truth_chi2,     "truth_chi2/F");
    truthfit->Branch("ndof",       &truth_ndof,     "truth_ndof/I");

    for (int cryNum = anacrys_start; cryNum < anacrys_end; cryNum++) {
        auto [hist, file] = get_data_histogram(cryNum, isSiPMRun);
        if (!hist) {
            if (file) { file->Close(); delete file; }
            continue;
        }
        hist->SetDirectory(0);

        // --- 1. Cut-and-count MC truth (existing) ---
        cryNumparam = cryNum;
        mcinfo truth;
        truth.RunMCTruth(
            hist, cryNum, disk, trueinfo,
            cryNumparam, tot_evts, mainpeak, first_espeak, second_espeak, background,
            frmainpeak, frfirst_espeak, frsecond_espeak, frbackground,
            compton1, compton2, compton3, frcompton1, frcompton2, frcompton3
        );

        // --- 2. Fit MC truth spectrum (NEW) ---
        truth_cryNum = cryNum;
        truth.FitMCTruthSpectrum(hist, cryNum, truthfit,
                                 truth_peak_mev, truth_fr_full, truth_fr_1st, truth_fr_2nd,
                                 truth_fr_c1, truth_fr_c2, truth_fr_c3,truth_fr_ebk,
                                 truth_width, truth_alpha, truth_chi2, truth_ndof,false);


        file->Close();
        delete file;
        delete hist;
    }

    mcinfo summary;
    summary.FinalizeMCSummary(trueinfo, truthfit, false);

    truthtable->cd();
    trueinfo->Write();
    truthfit->Write();  
    truthtable->Write();
    truthtable->Close();

    std::cout << "[INFO] MC Truth analysis complete. Saved to truthinfo.root" << std::endl;
    std::cout << "[INFO] MC Truth fits available in truthfit tree for comparison" << std::endl;

     return 0;
}
// --------------------------------------------------------
// FITTING BLOCK (with MC Truth Comparison if available)
// --------------------------------------------------------
  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents, convergencestatus;
  Float_t fpeak, peakerrorhigh,peakerrorlo,redpeak,fsigma, chiSq, fstpeak,
  scdpeak,fcbalphaparam,fcbndegparam,fcbalphaRparam,fcbndegRparam,comCnstparam,combetaparam,frFullparam,frFrstparam,
  frScndparam,crystalNoparam,frBKGparam,frcomptonparam1,frcomptonparam2,frcomptonparam3,
  pval,h_means,h_stddevs,unreducedchi2,fval,mparam,etaparam,widtherrorhigh,
  widtherrorlo,errbarhigh,errbarlo,evtfullerrorhigh,evtfullerrorlo,
  Esparam,Aparam,Aperr,Bparam,Berr,Cparam,Cerr;
  Int_t ndof;

  // MC Truth reference branches (for comparison diagnostics)
  Float_t truth_peak_mev_ref, truth_fr_full_ref, truth_fr_1st_ref, truth_fr_2nd_ref;
  Float_t recon_peak_mev;  // Reconstructed peak converted to MeV
  Int_t has_truth_match;   // Flag: was MC truth data found for this crystal

  TTree *covar = new TTree("covar","Covariance Plot");
  covar->Branch("nEvents", &nEvents,"nEvents/I");
  covar->Branch("Peak", &fpeak,"fpeak/F");
  covar->Branch("PeakErrHigh", &peakerrorhigh,"peakerrorhigh/F");
  covar->Branch("PeakErrLo", &peakerrorlo,"peakerrorlo/F");
  covar->Branch("CompositePeak", &redpeak,"redpeak/F");
  covar->Branch("Width", &fsigma,"fsigma/F");
  covar->Branch("WidthErrHigh", &widtherrorhigh,"widtherrorhigh/F");
  covar->Branch("WidthErrLo", &widtherrorlo,"widtherrorlo/F");
  covar->Branch("ChiSq", &chiSq,"chiSq/F");
  covar->Branch("1stPeak", &fstpeak,"fstpeak/F");
  covar->Branch("2ndPeak", &scdpeak,"scdpeak/F");
  covar->Branch("Alpha", &fcbalphaparam,"fcbalphaparam/F");
  covar->Branch("AlphaR", &fcbalphaRparam,"fcbalphaRparam/F");
  covar->Branch("NdegR",  &fcbndegRparam, "fcbndegRparam/F");
  covar->Branch("Ndeg", &fcbndegparam,"fcbndegparam/F");
  covar->Branch("comCnst", &comCnstparam,"comCnstparam/F");
  covar->Branch("combeta", &combetaparam,"combetaparam/F");
  covar->Branch("frFull", &frFullparam,"frFullparam/F");
  covar->Branch("frFrst", &frFrstparam,"frFrstparam/F");
  covar->Branch("frScnd", &frScndparam,"frScndparam/F");
  covar->Branch("frBKG", &frBKGparam,"frBKGparam/F");
  covar->Branch("frcompton1", &frcomptonparam1,"frcomptonparam1/F");
  covar->Branch("frcompton2", &frcomptonparam2,"frcomptonparam2/F");
  covar->Branch("frcompton3", &frcomptonparam3,"frcomptonparam3/F");
  covar->Branch("crystalNo", &crystalNoparam,"crystalNoparam/F");
  covar->Branch("convgstatus", &convergencestatus,"convergencestatus/I");
  covar->Branch("pval", &pval,"pval/F");
  covar->Branch("h_means", &h_means,"h_means/F");
  covar->Branch("h_stddevs", &h_stddevs,"h_stddevs/F");
  covar->Branch("unreducedchi2", &unreducedchi2,"unreducedchi2/F");
  covar->Branch("fval", &fval,"fval/F");
  covar->Branch("m", &mparam,"mparam/F");
  covar->Branch("eta", &etaparam,"etaparam/F");
  covar->Branch("ndof", &ndof,"ndof/I");
  covar->Branch("errbarhigh", &errbarhigh,"errbarhigh/F");
  covar->Branch("errbarlo", &errbarlo,"errbarlo/F");
  covar->Branch("evtfullerrorhigh", &evtfullerrorhigh,"evtfullerrorhigh/F");
  covar->Branch("evtfullerrorlo", &evtfullerrorlo,"evtfullerrorlo/F");
  covar->Branch("Esparam", &Esparam,"Esparam/I");
  covar->Branch("A",    &Aparam, "Aparam/F");
  covar->Branch("AErr", &Aperr,  "Aperr/F");
  covar->Branch("B",    &Bparam, "Bparam/F");
  covar->Branch("BErr", &Berr,   "Berr/F");
  covar->Branch("C",    &Cparam, "Cparam/F");
  covar->Branch("CErr", &Cerr,   "Cerr/F");

  // MC Truth reference branches (diagnostics for peak misattribution)
  covar->Branch("truth_peak_mev",     &truth_peak_mev_ref,    "truth_peak_mev_ref/F");
  covar->Branch("truth_fr_full",      &truth_fr_full_ref,     "truth_fr_full_ref/F");
  covar->Branch("truth_fr_1st",       &truth_fr_1st_ref,      "truth_fr_1st_ref/F");
  covar->Branch("truth_fr_2nd",       &truth_fr_2nd_ref,      "truth_fr_2nd_ref/F");
  covar->Branch("recon_peak_mev",     &recon_peak_mev,        "recon_peak_mev/F");
  covar->Branch("has_truth_match",    &has_truth_match,       "has_truth_match/I");

  // Load MC truth fits for reference (if available)
  TFile *truthfile = nullptr;
  TTree *truthfit_ref = nullptr;
  Int_t truth_cryNum_ref;
  Float_t truth_peak_mev_load, truth_fr_full_load, truth_fr_1st_load, truth_fr_2nd_load;
  std::map<int, std::tuple<Float_t, Float_t, Float_t, Float_t>> truth_data_map;

  // Try to open MC truth file for comparison
  truthfile = TFile::Open("truthinfo.root", "READ");
  if (truthfile && truthfile->IsOpen()) {
      truthfit_ref = (TTree*)truthfile->Get("truthfit");
      if (truthfit_ref) {
          std::cout << "[INFO] Loaded MC truth fits from truthinfo.root for comparison\n";
          truthfit_ref->SetBranchAddress("crystalNo",   &truth_cryNum_ref);
          truthfit_ref->SetBranchAddress("peak_mev",    &truth_peak_mev_load);
          truthfit_ref->SetBranchAddress("fr_full",     &truth_fr_full_load);
          truthfit_ref->SetBranchAddress("fr_1st",      &truth_fr_1st_load);
          truthfit_ref->SetBranchAddress("fr_2nd",      &truth_fr_2nd_load);

          // Cache truth data by crystal number
          for (Long64_t i = 0; i < truthfit_ref->GetEntries(); i++) {
              truthfit_ref->GetEntry(i);
              truth_data_map[truth_cryNum_ref] = std::make_tuple(
                  truth_peak_mev_load, truth_fr_full_load, truth_fr_1st_load, truth_fr_2nd_load
              );
          }
      }
  }

  auto start_bin = high_resolution_clock::now();
  for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    auto [hSum, file] = get_data_histogram(cryNum, isSiPMRun);
    auto [mean, stddev] = ComputeHistogramStats(hSum);
		h_means   = mean;
		h_stddevs = stddev;
		
	bool islyso = false;
	if (isSiPMRun){
		if (std::find(lysoSiPMs.begin(),lysoSiPMs.end(),cryNum) != lysoSiPMs.end()) islyso = true;
	}
	else{
		if(std::find(lysoCrystals.begin(), lysoCrystals.end(), cryNum) != lysoCrystals.end()) islyso = true; 
	}
	
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(hSum, alg, cryNum, covar, nEvents,convergencestatus, fpeak, peakerrorhigh,peakerrorlo, redpeak,
                    fsigma,widtherrorhigh,widtherrorlo, chiSq, fstpeak, scdpeak,
                    fcbalphaparam, fcbndegparam, comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam, crystalNoparam, frBKGparam,
                    frcomptonparam1,frcomptonparam2,frcomptonparam3,pval,h_means,
                    h_stddevs,unreducedchi2,
                    fval,mparam,etaparam,ndof,contour,xSelect,ySelect,errbarhigh,
                    errbarlo, evtfullerrorhigh,evtfullerrorlo,Esparam,Aparam,Aperr,Bparam,Berr,Cparam,Cerr,fcbalphaRparam,fcbndegRparam,islyso);
    file->Close();
    delete file;

    // --- MC Truth Comparison (if available) ---
    has_truth_match = 0;
    truth_peak_mev_ref = -1.0;
    truth_fr_full_ref = -1.0;
    truth_fr_1st_ref = -1.0;
    truth_fr_2nd_ref = -1.0;
    recon_peak_mev = -1.0;

    if (truth_data_map.count(cryNum)) {
        has_truth_match = 1;
        auto [t_peak, t_ff, t_f1, t_f2] = truth_data_map[cryNum];
        truth_peak_mev_ref = t_peak;
        truth_fr_full_ref = t_ff;
        truth_fr_1st_ref = t_f1;
        truth_fr_2nd_ref = t_f2;

        // Convert reconstructed peak to MeV using fitted calibration
        recon_peak_mev = mparam * fpeak;  // mparam is MeV/ADC, fpeak is ADC

        // Diagnostic: flag potential peak misattribution
        if (frFrstparam < 0.05 * truth_fr_1st_ref && truth_fr_1st_ref > 0.1) {
            std::cout << "[WARNING] Crystal " << cryNum << " peak misattribution suspected!\n";
            std::cout << "          Truth 1st esc frac: " << truth_fr_1st_ref
                      << " | Recon: " << frFrstparam << "\n";
        }
    }

    delete fit;
    delete hSum;
}
// --------------------------------------------------------
// OVERLAY BLOCK (Optional)
// --------------------------------------------------------
if (doOverlay) {
	for (int cryNum = anacrys_start; cryNum + 1 < anacrys_end; cryNum += 2) {

    auto [hist_even, file_even] = get_data_histogram(cryNum, isSiPMRun);
    auto [hist_odd, file_odd]  = get_data_histogram(cryNum + 1, isSiPMRun);
     hist_even->SetDirectory(0);
     hist_odd->SetDirectory(0);
    hist_even->SetLineColor(kBlue);
    hist_even->SetLineWidth(2);
    hist_odd->SetLineColor(kRed);
    hist_odd->SetLineWidth(2);
    TH1F* residual = (TH1F*)hist_odd->Clone(Form("residual_%d_%d", cryNum, cryNum+1));
    residual->SetDirectory(0);
    residual->Reset();
    int nBins = hist_odd->GetNbinsX();
    	for (int i = 1; i <= nBins; ++i) {
        double odd   = hist_odd->GetBinContent(i);
        double even  = hist_even->GetBinContent(i);
        double denom = sqrt(odd + even);
        double value = (denom > 0) ? (odd - even) / denom : 0;
        residual->SetBinContent(i, value);

        double err_odd  = hist_odd->GetBinError(i);
        double err_even = hist_even->GetBinError(i);
        double err = denom > 0 ? sqrt(err_odd*err_odd + err_even*err_even) / denom : 0;
        residual->SetBinError(i, err);
    }
    TCanvas* cOverlay = new TCanvas(Form("cOverlay_%d_%d", cryNum, cryNum+1),
                                    Form("Even/Odd Overlay %d & %d", cryNum, cryNum+1),
                                    800, 800);
    cOverlay->Divide(1, 2, 0, 0);
    cOverlay->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    hist_even->Draw("hist");
    hist_even->GetYaxis()->SetRangeUser(0, 5000);
    hist_odd->Draw("hist same");
    hist_odd->GetYaxis()->SetRangeUser(0, 5000);

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(hist_even, Form("Even SiPM (%d)", cryNum), "l");
    leg->AddEntry(hist_odd, Form("Odd SiPM (%d)", cryNum+1), "l");
    leg->Draw();

    TPaveText* statsText = new TPaveText(0.15, 0.7, 0.4, 0.88, "NDC");
    statsText->SetFillColor(0);
    statsText->SetTextSize(0.03);
    statsText->AddText(Form("Even Entries: %.0f", hist_even->GetEntries()));
    statsText->AddText(Form("Odd Entries: %.0f", hist_odd->GetEntries()));
    statsText->Draw();

    // --- Bottom pad: residual ---
    cOverlay->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3);
    residual->SetTitle("Normalised Residual (Residual / sqrt(Bin Count))");
    residual->GetXaxis()->SetTitle("ADC");
    residual->GetXaxis()->SetLabelSize(0.06);
    residual->GetXaxis()->SetTitleSize(0.07);
    residual->GetYaxis()->SetTitle("Counts");
    residual->GetYaxis()->SetLabelSize(0.06);
    residual->GetYaxis()->SetTitleSize(0.07);
    residual->Draw("hist");

    gPad->Update();
    if (auto stats = (TPaveStats*)residual->FindObject("stats")) {
        stats->SetX1NDC(0.7);
        stats->SetX2NDC(0.95);
        stats->SetY1NDC(0.15);
        stats->SetY2NDC(0.45);
        stats->SetTextSize(0.07);
    }
    // --- Save to a unique ROOT file ---
	TString oName = "mu2e_simu_hist_" + std::to_string(cryNum) + "_" + std::to_string(cryNum+1) + ".root";
    TFile* outputFile = new TFile(oName, "RECREATE");
    cOverlay->Write();
    outputFile->Close();
    delete outputFile;
    delete cOverlay;
    delete residual;
		delete leg;
		delete statsText;
		delete hist_even;
		delete hist_odd;
		}
}      
  auto end_bin = high_resolution_clock::now();
// --------------------------------------------------------
// SUMMARY & PLOTTING
// --------------------------------------------------------  
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots, anacrys_start, anacrys_end, isSiPMRun);

  // Cleanup MC truth file if loaded
  if (truthfile && truthfile->IsOpen()) {
      truthfile->Close();
      std::cout << "[INFO] MC truth comparison file closed.\n";
  }

  table->cd();
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();

  std::cout << "[INFO] Reconstruction fits written to arXivTable.root\n";
  if (truth_data_map.size() > 0) {
      std::cout << "[INFO] " << truth_data_map.size() << " crystals have MC truth comparison data\n";
  }

auto printConv = [](std::ostream& out, const std::string& label, int count, const std::vector<int>& crystals, bool printDetails) {
    
        out << label << ": " << count;
        if (printDetails && !crystals.empty()) {
            out << " (";
            for (size_t i = 0; i < crystals.size(); ++i) {
                out << crystals[i] << (i == crystals.size() - 1 ? "" : ", ");
            }
            out << ")";
        }
        out << "\n";
    };

    std::ofstream logFile("Fit_Summary.txt");
    if (!logFile.is_open()) {
        std::cerr << "[ERROR] Could not open Fit_Summary.txt for writing!" << std::endl;
    }

    auto writeSummary = [&](std::ostream& out, bool detailed) {
        out << "\n================ Refit Summary =================\n";
        out << "First fit converged: " << SourceFitter::nFirstFitConverged << "\n";

        printConv(out, "Second fit converged", SourceFitter::nSecondFitConverged, SourceFitter::crystalsSecondFitConverged, detailed);
        printConv(out, "Third fit converged",  SourceFitter::nThirdFitConverged,  SourceFitter::crystalsThirdFitConverged, detailed);
        
        if (detailed && !SourceFitter::thirdFitRetryCount.empty()) {
            out << "\nRetry attempts per crystal:\n";
            for (const auto &p : SourceFitter::thirdFitRetryCount) {
                out << "  Crystal " << p.first << ": " << p.second << " attempts\n";
            }
            out << "\n";
        }
        
        printConv(out, "Asymmetric Errors Found", SourceFitter::nAsymErrors, SourceFitter::crystalsWithAsymErrors, detailed);
        out << "Bad chi2: " << SourceFitter::badchi2;
        if (detailed && !SourceFitter::crystalswithbadchi2.empty()) {
            out << " (";
            for (size_t i = 0; i < SourceFitter::crystalswithbadchi2.size(); ++i) {
                int cid = SourceFitter::crystalswithbadchi2[i];
                out << cid << ": chi2=" << SourceFitter::badchi2Values[cid];
                if (i != SourceFitter::crystalswithbadchi2.size()-1) out << ", ";
            }
            out << ")";
        }
        out << "\n";
        printConv(out, "Hesse errors used", SourceFitter::nHesseFallbacks, SourceFitter::crystalsHesseFallback, detailed);
        
        if (!SourceFitter::convFailures.empty()) {
            if (detailed) {
                out << "Crystals with non-zero convergence status:\n";
                for (const auto& [cryNo, status] : SourceFitter::convFailures) {
                    out << "  Crystal " << cryNo << " : status = " << status << "\n";
                }
            } else {
                out << "Total convergence failures: " << SourceFitter::convFailures.size() << "\n";
            }
        }
        out << "================================================\n";
    };
	//can switch to true if want a list of crystals for each counter displayed in terminal
    writeSummary(std::cout, false);

    if (logFile.is_open()) {
        std::cout << "(Detailed summary with crystal lists saved to Fit_Summary.txt)\n";
        writeSummary(logFile, true);
        logFile.close();
    }

  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
