/*#include "CaloCalibration/SourceCalib/inc/validationfitter.hh"

using namespace RooFit;

namespace CaloSourceCalib {

    void ValidationFitter::ValidateCrystal(TH1F* h_spec, int crystalNo, std::string txtFilePath) {
        
        // =========================================================
        // 1. READ ALL SAVED PARAMETERS FROM THE TXT FILE
        // =========================================================
        double loadedPeak = 0.0;
        double loadedWidth = 0.0; 
        double loadedEvtsFull = 0.0;
        double loadedEvtsFst = 0.0;
        double loadedEvtsScd = 0.0;
        double loadedEvtsBkg = 0.0;
        double loadedComCnst = 0.0;
        double loadedComBeta = 0.0;

        std::ifstream infile(txtFilePath);
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open file: " << txtFilePath << std::endl;
            return;
        }

        std::string line;
        while (std::getline(infile, line)) {
            // Skip the header line or any empty lines
            if (line.empty() || line[0] == '#') continue;

            std::stringstream ss(line);
            
            // Temporary variables (all doubles to handle formatting)
            double txtCrysID, txtPeak, txtPeakErrHi, txtPeakErrLo, txtUnredChi2;
            double txtNdof, txtNEvents, txtWidth, txtWidthErrHi, txtWidthErrLo;
            double txtChi2, txtEvtsFull, txtEvtsFst, txtEvtsScd;
            
            // New variables for the background
            double txtEvtsBkg, txtComCnst, txtComBeta, txtExtraEvents;

            // Extract all 18 columns to keep perfect alignment
            if (ss >> txtCrysID >> txtPeak >> txtPeakErrHi >> txtPeakErrLo 
                   >> txtUnredChi2 >> txtNdof >> txtNEvents 
                   >> txtWidth >> txtWidthErrHi >> txtWidthErrLo 
                   >> txtChi2 >> txtEvtsFull >> txtEvtsFst >> txtEvtsScd
                   >> txtEvtsBkg >> txtComCnst >> txtComBeta >> txtExtraEvents) {
                
                if (static_cast<int>(txtCrysID) == crystalNo) {
                    loadedPeak = txtPeak;
                    loadedWidth = txtWidth;
                    loadedEvtsFull = txtEvtsFull;
                    loadedEvtsFst = txtEvtsFst;
                    loadedEvtsScd = txtEvtsScd;
                    
                    // Save the new background variables
                    loadedEvtsBkg = txtEvtsBkg;
                    loadedComCnst = txtComCnst;
                    loadedComBeta = txtComBeta;
                    break; 
                }
            }
        }
        infile.close();

      

        // =========================================================
        // 2. SETUP ROOFIT VARIABLES (LOCKED TO SAVED VALUES)
        // =========================================================
        RooRealVar crysADC("crysADC", "ADC [counts]", 40, 120);
        RooRealVar m_e("m_e", "electron energy in MeV", 0.511);
        RooRealVar E0("E0", "energy offset [MeV]", 0.0); 
        
        // Lock the Peak and Width from the text file
        RooRealVar fullPeak("fullPeak", "Full peak [ADC]", loadedPeak);
        fullPeak.setConstant(kTRUE); 
        
        RooRealVar fullWidth("fullWidth", "Full width [MeV]", loadedWidth); 
        fullWidth.setConstant(kTRUE);

        // Rebuild linked formulas
        RooFormulaVar eta("eta", "ADC/MeV", "fullPeak / (6.13-E0)", RooArgSet(fullPeak, E0));
        RooFormulaVar m("m","Mev/ADC", "(6.13-E0) / fullPeak", RooArgSet(fullPeak, E0));
        RooFormulaVar Es("Es", "calibrated energy [MeV]", "m*crysADC + E0", RooArgSet(m, crysADC, E0));
        
        // Shape constants (Hardcoded to what you used in the main script)
        RooRealVar fcbalpha("fcbalpha", "alpha", 0.5); fcbalpha.setConstant(kTRUE);
        RooRealVar fcbndeg("fcbndeg", "n", 5); fcbndeg.setConstant(kTRUE);

        // Peak energy definitions
        RooRealVar Egamma("Egamma", "Full peak [MeV]", 6.13);
        RooRealVar fstesc("fstesc", "first peak [MeV]", 5.619);
        RooRealVar scdesc("scdesc", "second peak [MeV]", 5.108);

        // Gamma Shapes
        RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
        RooCBShape firsErg("firsErg", "Single escape", Es, fstesc, fullWidth, fcbalpha, fcbndeg);
        RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha, fcbndeg);

        // Background Shape (Logistic)
        RooRealVar comCnst("comCnst", "Compton Constant", loadedComCnst);
        RooRealVar combeta("combeta", "Compton Beta", loadedComBeta); 
        RooRealVar bkg_fixed("bkg_fixed", "Background offset", 0.0); 
        RooGenericPdf comPdf("comPdf", "logistic background",
            "pow(1.0 + exp((Es - comCnst)/combeta), -1.0) + bkg_fixed",
            RooArgSet(Es, comCnst, combeta, bkg_fixed));

        // =========================================================
        // 3. YIELDS & FINAL MODEL
        // =========================================================
        
        // Lock the gamma yields to your saved values!
        RooRealVar evtsFull("evtsFull", "Full peak yield", loadedEvtsFull);
        evtsFull.setConstant(kTRUE);
        
        RooRealVar evtsFrst("evtsFrst", "First escape yield", loadedEvtsFst);
        evtsFrst.setConstant(kTRUE);
        
        RooRealVar evtsScnd("evtsScnd", "Second escape yield", loadedEvtsScd);
        evtsScnd.setConstant(kTRUE);
        
        RooRealVar evtsbkg("evtsbkg", "Background Yield", loadedEvtsBkg);

        RooAddPdf fitFun("fitFun", "Total Validation Model",
            RooArgList(fullErg, firsErg, secdErg, comPdf), 
            RooArgList(evtsFull, evtsFrst, evtsScnd, evtsbkg) 
        );

        // =========================================================
        // 4. FIT (Only Background Changes) & PLOT
        // =========================================================
        RooDataHist chSpec("chSpec", "Data", crysADC, h_spec);
        
        // One very fast fit pass. Since all gammas are locked, this ONLY scales the background!
        fitFun.fitTo(chSpec, RooFit::PrintLevel(-1), RooFit::Range(40, 115.2)); 

        TCanvas *can = new TCanvas("can", "Validation Canvas", 800, 600);
        RooPlot *chFrame = crysADC.frame(RooFit::Title(Form("Validation: SiPM %d", crystalNo)));
        
        chSpec.plotOn(chFrame, RooFit::MarkerColor(kBlack), RooFit::MarkerSize(0.5));
        fitFun.plotOn(chFrame, RooFit::LineColor(kRed));
        fitFun.plotOn(chFrame, RooFit::Components(fullErg), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed));
        fitFun.plotOn(chFrame, RooFit::Components(firsErg), RooFit::LineColor(kViolet), RooFit::LineStyle(kDashed));
        fitFun.plotOn(chFrame, RooFit::Components(secdErg), RooFit::LineColor(kCyan),   RooFit::LineStyle(kDashed));
        fitFun.plotOn(chFrame, RooFit::Components(comPdf),  RooFit::LineColor(kBlue),   RooFit::LineStyle(kDashed));
        
        chFrame->Draw();
        
        // Save 
        can->SaveAs(Form("Validation_Crystal_%d.root", crystalNo));
        
        delete can;
        std::cout << "[ValidationFitter] Successfully overlaid data for Crystal " << crystalNo << std::endl;
    }
}*/
