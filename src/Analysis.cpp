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
#include<vector>

#include<TBox.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TGaxis.h>
#include<TGraphErrors.h>
#include<TLatex.h>
#include<TLegend.h>
#include<TLine.h>

#include"Analysis.h"
#include"ErrorAnalysis.h"
#include"LabPlot.h" // set_atlas_style() called from here
#include"Logger.h"

Bode::Bode(System_t sys){
    fSystem = sys;
    // set_atlas_style(tsize);

    SetSystem(sys);

    // We cant do much more, since its impossible to allocate memory for graphs, use functions!
}


Bode::Bode(System_t sys, const char *filename, Option_t *option){
    fSystem = sys;
    // set_atlas_style(tsize);

    SetSystem(sys);

    ReadInput(filename, option);

}

void Bode::Plot(bool plotphase, bool plotgain){
    fFigure = new TCanvas("fFigure", "", 800, 600);
    fFigure->cd();
    TLatex *text = new TLatex();
    TLegend *legend = new TLegend();
    
    fGainPad = new TPad("fGainPad", "", 0, 0, 1, 1);
    fGainPad->SetLogx();
    fGainPad->SetLogy();
    fPhasePad = new TPad("fPhasePad", "", 0, 0, 1, 1);
    fPhasePad->SetLogx();
    fPhasePad->SetFillStyle(4000);

    TGaxis *axis;

    if(plotgain){
        if(plotphase){
            set_atlas_style(tsize, false);
        } else {
            set_atlas_style(tsize);
        }
        _apply_LineColor();

        fGainPad->cd();
        fGain->Draw("ap");
        if(_hasfitted) { fGainFit->Draw("same"); }
        fGainPad->Modified();
        fGainPad->Draw();
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

            fPhase->Draw("ap");
            if(_hasfitted) { fPhaseFit->Draw("same"); }

            Style_t tfont = fGain->GetHistogram()->GetYaxis()->GetTitleFont();
            Float_t tsize = fGain->GetHistogram()->GetYaxis()->GetTitleSize();
            Style_t lfont = fGain->GetHistogram()->GetYaxis()->GetLabelFont();
            Float_t lsize = fGain->GetHistogram()->GetYaxis()->GetLabelSize();

            axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
            axis->SetTitle("Phase");
            axis->SetTitleOffset(1.5);
            axis->SetTitleFont(tfont);
            axis->SetTitleSize(tsize);
            axis->SetLabelFont(lfont);
            axis->SetLabelSize(lsize);
            axis->SetMaxDigits(1);
            axis->Draw();
            fPhasePad->Modified();
            fPhasePad->Draw();
            gPad->Update();
        }
    } else {
        // plotgain is false

        set_atlas_style(tsize);
        _apply_LineColor();

        fPhasePad->cd();
        fPhase->Draw("ap");
        if(_hasfitted) { fPhaseFit->Draw("same"); }
        fPhasePad->Modified();
        fPhasePad->Draw();
        gPad->Update();
    }
    fGainPad->Draw();
    fPhasePad->Draw();
    fFigure->Draw();
}

void Bode::PlotGain(){
    Plot(false);
}

void Bode::PlotPhase(){
    Plot(true, false);
}

double get_VRangeErr(double errPercent, int partitions, double range1){return errPercent * partitions *  range1;}
double get_TRangeErr(double range1, double errPercent = 0.0016, int partition = 10){return range1 * errPercent * partition;}
double get_HErr(double Vin, double Vout, double eVin, double eVout){ return sqrt(pow(eVout / Vin, 2) + pow(eVin * Vout / pow(Vin, 2), 2));}
double get_phi(double T, double dt){return 2 * M_PI * dt / T;}
double get_phiErr(double T, double dt, double eT, double edt){return 2 * M_PI * sqrt(pow(edt/T, 2) + pow(dt * eT/(pow(T, 2)), 2));}


Bool_t Bode::ReadInput(const char *filename, Option_t *option){

    std::ifstream data(filename);

    std::vector<double> tmpfreq, tmpfreqerr, tmpgain, tmpgainerr, tmpphase, tmpphaseerr;
    double Vin, fsVin, Vout, fsVout, T, fsT, dt, fsdt;


    for(int i=0; data >> Vin >> fsVin >> Vout >> fsVout >> T >> fsT >> dt >> fsdt ; i++){
        
        double eVin, eVout;
        if(fsVin<=0.01){
            eVin = get_VRangeErr(0.045, 8, fsVin)/sqrt(3);
        }else{
            eVin = get_VRangeErr(0.035, 8, fsVin)/sqrt(3);
        }
        if(fsVout<=0.01){
            eVout = get_VRangeErr(0.045, 8, fsVout)/sqrt(3);
        }else{
            eVout = get_VRangeErr(0.035, 8, fsVout)/sqrt(3);
        }
        double eT = get_TRangeErr(fsT)/sqrt(3);
        double edt = get_TRangeErr(fsdt)/sqrt(3);

        tmpfreq.push_back(1/T);
        tmpfreqerr.push_back(eT/pow(T, 2));
        tmpgain.push_back(Vout/Vin);
        tmpgainerr.push_back(get_HErr(Vin, Vout, eVin, eVout));
        tmpphase.push_back(get_phi(T, dt));
        tmpphaseerr.push_back(get_phiErr(T, dt, eT, edt));
    }

    SetGainVec(tmpgain, tmpgainerr);
    SetPhaseVec(tmpphase, tmpphaseerr);
    SetFreqVec(tmpfreq, tmpfreqerr);
    SetFunctions();

    return true;
}

Bool_t Bode::SetFreqVec(std::vector<Double_t> Freq, std::vector<Double_t> ErrFreq){
    fPointFreq = Freq.data(); 
    fPErrFreq = ErrFreq.data();

    if (fNpoints == -1){
        fNpoints = Freq.size();
    }else{
        if(Freq.size() != fNpoints){
            printf("%s", Logger::error("freq. array size does not match previous array size."));
            return false;
        }
    }
    

    if (Freq.size()!=ErrFreq.size()){
        printf("%s", Logger::error("array size of Freq and its Error do not match!"));
        return false;
    }

    return true;
}

Bool_t Bode::SetFunctions(){

    fGain = new TGraphErrors(fNpoints, fPointFreq, fPointGain, fPErrFreq, fPErrGain);
    fPhase = new TGraphErrors(fNpoints, fPointFreq, fPointPhase, fPErrFreq, fPErrPhase);
    fGainFit = new TF1("gain_fit", _gainfit);
    fPhaseFit = new TF1("phase_fit", _phasefit);

    return true;
}

Bool_t Bode::SetGainVec(std::vector<Double_t> Gain, std::vector<Double_t> ErrGain){
    fPointGain = Gain.data(); 
    fPErrGain = ErrGain.data();

    if (fNpoints == -1){
        fNpoints = Gain.size();
    }else{
        if(Gain.size() != fNpoints){
            printf("%s", Logger::error("gain array size does not match previous array size."));
            return false;
        }
    }

    if (Gain.size()!=ErrGain.size()){
        printf("%s", Logger::error("array size of Gain and its Error do not match!"));
        return false;
    }

    return true;
}

void Bode::SetParGain(Double_t gain, Double_t cutoff, Double_t Q){

    if(_islowhighpass){
        fGainFit->SetParameters(gain, cutoff);
    }else{
        fGainFit->SetParameters(gain, cutoff, Q);
    }

    return;
}

void Bode::SetParPhase(Double_t gain, Double_t cutoff, Double_t Q){

    if(_islowhighpass){
        fPhaseFit->SetParameters(gain, cutoff);
    }else{
        fPhaseFit->SetParameters(gain, cutoff, Q);
    }

    return;
}

Bool_t Bode::SetPhaseVec(std::vector<Double_t> Phase, std::vector<Double_t> ErrPhase){
    fPointPhase = Phase.data(); 
    fPErrPhase = ErrPhase.data();

    if (fNpoints == -1){
        fNpoints = Phase.size();
    }else{
        if(Phase.size() != fNpoints){
            printf("%s", Logger::error("phase array size does not match previous array size."));
            return false;
        }
    }

    if (Phase.size()!=ErrPhase.size()){
        printf("%s", Logger::error("array size of Phase and its Error do not match!"));
        return false;
    }

    return true;
}

void Bode::SetSystem(System_t sys){
    fSystem = sys;
    switch(sys.Hash()){
        case lowpass:
            _gainfit = "[0]/sqrt(1+pow(x/[1], 2))"; // [0] gain, [1] cutoff
            _phasefit = "";
            break;
        case highpass:
            _gainfit = "[0]/sqrt(1+pow([1]/x, 2))"; // [0] gain, [1] cutoff
            _phasefit = "";
            break;
        case bandpass:
            _gainfit = "[0]/sqrt(1+pow([2], 2)*pow(x/[1]-[1]/x, 2))"; // [0] gain, [1] cutoff/peak frequency, [2] Q factor
            _phasefit = "";
            _islowhighpass = false;
            break;
        default:
            fprintf(stderr, "%s", Logger::warning(Form("System_t option '%s' not recognised!\n"
            "Available options are: \n\t\"lowpass\"\n\t\"highpass\"\n\t\"bandpass\"\n ", sys.Data())));
    }
}

Bool_t Bode::FitGain(Option_t *option, Option_t *goption, Axis_t xmin, Axis_t xmax){

    fGain->Fit("gain_fit");
    _hasfitted = true;

    return true;
}

Bool_t Bode::FitPhase(Option_t *option, Option_t *goption, Axis_t xmin, Axis_t xmax){

    fPhase->Fit("phase_fit");
    _hasfitted = true;

    return true;
}

Bode::~Bode(){
    fprintf(stderr, "%s\n", Logger::warning("Deleted obj. Bode"));
}