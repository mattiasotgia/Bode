/**
 * @file Analysis.cpp
 * @author Mattia Sotgia (mattiasotgia01@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include<iostream>

#include<TBox.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TGaxis.h>
#include<TGraphErrors.h>
#include<TLatex.h>
#include<TLegend.h>
#include<TLine.h>

#include"Bode/Analysis.h"
#include"ErrorAnalysis.h"
#include"LabPlot.h" // set_atlas_style() called from here
#include"Logger.h"

Bode::Bode(){
    // set_atlas_style(tsize);
    fGain = new TGraphErrors();
    fPhase = new TGraphErrors();
    fGainFit = new TF1("gainfit", _gainfit, fmin, fmax);
    fPhaseFit = new TF1("phasefit", _phasefit, fmin, fmax);
}

Bode::Bode(System_t sys){
    fSystem = sys;
    // set_atlas_style(tsize);

    SetSystem(sys);

    fGain = new TGraphErrors();
    fPhase = new TGraphErrors();
    fGainFit = new TF1("gainfit", _gainfit, fmin, fmax);

    fPhaseFit = new TF1("phasefit", _phasefit, fmin, fmax);
}


Bode::Bode(System_t sys, const char *filename, Option_t *option){
    fSystem = sys;
    // set_atlas_style(tsize);

    SetSystem(sys);

    fGain = new TGraphErrors();
    fPhase = new TGraphErrors();
    fGainFit = new TF1("gainfit", _gainfit, fmin, fmax);
    fPhaseFit = new TF1("phasefit", _phasefit, fmin, fmax);
}

void Bode::Plot(bool plotphase, bool plotgain){
    fFigure = new TCanvas("fFigure", "", 800, 600);
    TLatex *text = new TLatex();
    TLegend *legend = new TLegend();
    
    fGainPad = new TPad("fGainPad", "", 0, 0, 1, 1);
    fPhasePad = new TPad("fPhasePad", "", 0, 0, 1, 1);
    fPhasePad->SetFillStyle(4000);

    if(plotgain){
        if(plotphase){
            set_atlas_style(tsize, false);
        } else {
            set_atlas_style(tsize);
        }
        _apply_LineColor();

        fGainPad->cd();
        fGain->Draw("ap");
        fGainFit->Draw("same");
        if(plotphase){
            fGainPad->SetTicky(0);
            
            Double_t xmin = fGainPad->GetUxmin();
            Double_t xmax = fGainPad->GetUxmax();
            Double_t dx = (xmax - xmin) / 0.68; // 10 percent margins left and right
            Double_t ymin = fPhase->GetHistogram()->GetMinimum();
            Double_t ymax = fPhase->GetHistogram()->GetMaximum();
            Double_t dy = (ymax - ymin) / 0.79; // 10 percent margins top and bottom
            fPhasePad->Range(xmin-0.16*dx, ymin-0.16*dy, xmax+0.16*dx, ymax+0.05*dy);

            fPhasePad->cd();

            fPhase->Draw("p");
            fPhaseFit->Draw("same");

            Style_t tfont = fGain->GetHistogram()->GetYaxis()->GetTitleFont();
            Float_t tsize = fGain->GetHistogram()->GetYaxis()->GetTitleSize();
            Style_t lfont = fGain->GetHistogram()->GetYaxis()->GetLabelFont();
            Float_t lsize = fGain->GetHistogram()->GetYaxis()->GetLabelSize();

            TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
            axis->SetTitle("Phase");
            axis->SetTitleOffset(1.5);
            axis->SetTitleFont(tfont);
            axis->SetTitleSize(tsize);
            axis->SetLabelFont(lfont);
            axis->SetLabelSize(lsize);
            axis->SetMaxDigits(1);
            axis->Draw();
            gPad->Update();
        }
    } else {
        // plotgain is false

        set_atlas_style(tsize);
        _apply_LineColor();

        fPhasePad->cd();
        fPhase->Draw("ap");
        fPhaseFit->Draw("same");
        gPad->Update();
    }
}

void Bode::PlotGain(){
    Plot(false);
}

void Bode::PlotPhase(){
    Plot(true, false);
}

Bool_t Bode::SetFreqVec(std::vector<Double_t> Freq, std::vector<Double_t> ErrFreq){
    fPointFreq = Freq.data(); 
    fPErrFreq = ErrFreq.data();

    if (Freq.size()!=ErrFreq.size()){
        printf("%s", Logger::error("array size of Freq and its Error do not match!"));
        return false;
    }

    Int_t kSizePG = sizeof(fPointFreq)/sizeof(Double_t);
    Int_t kSizeEPG = sizeof(fPErrFreq)/sizeof(Double_t);

    if(Freq.size()!=kSizePG || ErrFreq.size()!=kSizeEPG){
        printf("%s", Logger::error(Form("size of Freq input vec: %lu, size of Freq copy vec: %d", Freq.size(), kSizePG)));
        printf("%s", Logger::error(Form("size of ErrFreq input vec: %lu, size of ErrFreq copy vec: %d", Freq.size(), kSizePG)));
        return false;
    }
    return true;
}

Bool_t Bode::SetGainVec(std::vector<Double_t> Gain, std::vector<Double_t> ErrGain){
    fPointGain = Gain.data(); 
    fPErrGain = ErrGain.data();

    if (Gain.size()!=ErrGain.size()){
        printf("%s", Logger::error("array size of Gain and its Error do not match!"));
        return false;
    }

    Int_t kSizePG = sizeof(fPointGain)/sizeof(Double_t);
    Int_t kSizeEPG = sizeof(fPErrGain)/sizeof(Double_t);

    if(Gain.size()!=kSizePG || ErrGain.size()!=kSizeEPG){
        printf("%s", Logger::error(Form("size of Gain input vec: %lu, size of Gain copy vec: %d", Gain.size(), kSizePG)));
        printf("%s", Logger::error(Form("size of ErrGain input vec: %lu, size of ErrGain copy vec: %d", Gain.size(), kSizePG)));
        return false;
    }
    return true;
}

Bool_t Bode::SetPhaseVec(std::vector<Double_t> Phase, std::vector<Double_t> ErrPhase){
    fPointPhase = Phase.data(); 
    fPErrPhase = ErrPhase.data();

    if (Phase.size()!=ErrPhase.size()){
        printf("%s", Logger::error("array size of Phase and its Error do not match!"));
        return false;
    }

    Int_t kSizePG = sizeof(fPointPhase)/sizeof(Double_t);
    Int_t kSizeEPG = sizeof(fPErrPhase)/sizeof(Double_t);

    if(Phase.size()!=kSizePG || ErrPhase.size()!=kSizeEPG){
        printf("%s", Logger::error(Form("size of Phase input vec: %lu, size of Phase copy vec: %d", Phase.size(), kSizePG)));
        printf("%s", Logger::error(Form("size of ErrPhase input vec: %lu, size of ErrPhase copy vec: %d", Phase.size(), kSizePG)));
        return false;
    }
    return true;
}

void Bode::SetSystem(System_t sys){
    fSystem = sys;
    switch(sys.Hash()){
        case op_amp_lp:
        case OP_AMP_LP:
            _gainfit = "[0]/sqrt(1+pow(x/[1], 2))";
            _phasefit = "";
            break;
        case op_amp_hp:
        case OP_AMP_HP:
            _gainfit = "";
            _phasefit = "";
            break;
        case rlc_lp:
        case RLC_LP:
            _gainfit = "";
            _phasefit = "";
            break;
        case rlc_hp:
        case RLC_HP:
            _gainfit = "";
            _phasefit = "";
            break;
        case rlc_bp:
        case RLC_BP:
            _gainfit = "";
            _phasefit = "";
            break;
        default:
            fprintf(stderr, "%s", Logger::warning(Form("System_t option '%s' not recognised", sys.Data())));
    }
}

Bode::~Bode(){
    fprintf(stderr, "%s\n", Logger::warning("Deleted obj Bode"));
}