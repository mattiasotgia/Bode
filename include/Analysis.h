/**
 * @file Analysis.h
 * @author Mattia Sotgia (mattiasotgia01@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-01
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef LABTOOLS_BODE_Analysis
#define LABTOOLS_BODE_Analysis

#include<iostream>
#include<vector>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraphErrors.h>
#include<TPad.h>
#include<TString.h>

// typedefs
typedef int NPar_t;
typedef TString System_t;


class Bode{
private:
    // miscellaneous
    enum {
        OP_AMP_HP = 1991226620,
        OP_AMP_LP = 4138710268,
        RLC_HP    = 3899251227,
        RLC_LP    = 3899252251,
        RLC_BP    = 3899248667,
        op_amp_hp = 1982541912,
        op_amp_lp = 4130025560,
        rlc_hp    = 3966631483,
        rlc_lp    = 3966630459,
        rlc_bp    = 3966628923
    };

    System_t            fSystem = "";

    bool                _isfunctioncalled = false; ///> check if user called for function different from base;
    bool                _residualOn = false;

    // font size for calling ATLASStyle
    Size_t              tsize = 29;

    // phase and gain fit formula
    const char         *_gainfit = "";  // todo: fill function
    const char         *_phasefit = ""; // todo: fill function
    NPar_t              _CutoffPar;
    NPar_t              _GainPar;

    // computational values;
    Double_t            gGWB;       ///> Gain Bandwidth coefficient
    Double_t            gErrGWB;    ///> Gain Bandwidth coefficient error
    Double_t            gCutoff;    ///> Cutoff value 
    Double_t            gErrCutoff; ///> Cutoff value error
    Double_t            gGain;      ///> Gain value
    Double_t            gErrGain;   ///> Gain value error

    /// graphical objects []
    TCanvas            *fFigure;
    TPad               *fGainPad;
    TPad               *fPhasePad;
    TGraphErrors       *fGain;
    TGraphErrors       *fPhase;
    TF1                *fGainFit; 
    TF1                *fPhaseFit;

    /// function variables declaration
    Double_t            fmin = (0.0);   ///> minimum for frequency range
    Double_t            fmax = (1.0);   ///> maximum for frequency range
    Double_t            fNpoints;        ///> points in graph
    Double_t           *fPointGain;     ///>[fNpoints] array for gain points
    Double_t           *fPErrGain;      ///>[fNpoints] array for error gain points
    Double_t           *fPointPhase;    ///>[fNpoints] array for phase points
    Double_t           *fPErrPhase;     ///>[fNpoints] array for error phase points
    Double_t           *fPointFreq;
    Double_t           *fPErrFreq;

    // miscellaneous
    void                _apply_LineColor();
public:
    Bode();
    Bode(System_t sys);                                             ///> sys: OP_AMP: 0x1 high pass, 0x2 low pass; RLC 0x101 high pass, 0x102 low pass, 0x103 band pass;
    Bode(System_t sys, const char *filename, Option_t *option="");  ///> sys: OP_AMP: 0x1 high pass, 0x2 low pass; RLC 0x101 high pass, 0x102 low pass, 0x103 band pass;
    ~Bode();
    Bool_t              FitCorrelated(Option_t *option="", Option_t *goption="", Axis_t xmin=0, Axis_t xmax=0);
    inline Double_t     GetCutoff()     const { return gCutoff; }
    inline Double_t     GetErrCutoff()  const { return gErrCutoff; }
    inline Double_t     GetErrGain()    const { return gErrGain; }
    inline Double_t     GetErrGWB()     const { return gErrGWB; }
    inline Double_t     GetGain()       const { return gGain; }
    inline Double_t     GetGWB()        const { return gGWB; }
    void                Plot(bool plotphase = true, bool plotgain = true);
    void                PlotGain();
    void                PlotPhase();
    Bool_t              ReadInput(const char *filename, Option_t *option="");       ///> read input for both phase and gain data 
    // bool                ReadInputGain(const char *filename, Option_t *option="");   ///> read input for gain data
    // bool                ReadInputPhase(const char *filename, Option_t *option="");  ///> read input for phase data
    // bool                ReadInputRDF()  // TO BE IMPLEMENTED
    inline void         SetCutoffNpar(NPar_t npar = 0)  { _CutoffPar = npar; }
    inline void         SetGainNpar(NPar_t npar = 1)    { _GainPar = npar; }
    Bool_t              SetFreqVec(std::vector<Double_t> Freq, std::vector<Double_t> ErrFreq);
    // void                SetGainFunction(const char *formula, Option_t *option="");
    Bool_t              SetGainVec(std::vector<Double_t> Gain, std::vector<Double_t> ErrGain);
    void                SetParGain(Double_t *params);
    void                SetParPhase(Double_t *params);
    // void                SetPhaseFunction(const char *formula, Option_t *option="");
    Bool_t              SetPhaseVec(std::vector<Double_t> Phase, std::vector<Double_t> ErrPhase);
    void                SetLogx();
    void                SetLogy();
    void         SetSystem(System_t sys);
    inline void         SetResidual(bool residual = true) { _residualOn = residual; }
};

void Bode::_apply_LineColor(){
    fGain->SetTitle(";Frequency;Gain");
    fPhase->SetTitle(";Frequency;Phase");
    fGainFit->SetLineColor(kBlack);
    fPhase->SetLineColor(kRed);
    fPhase->SetMarkerColor(kRed);
    fPhaseFit->SetMarkerColor(kRed);
    fPhaseFit->SetLineStyle(kDashed);
}

#endif