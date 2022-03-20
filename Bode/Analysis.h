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

#ifndef BODE_Analysis
#define BODE_Analysis

#include<iostream>
#include<vector>
#include<string>

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
        lowpass     = 244089597,    // "lowpass"
        highpass    = 2141088893,   // "highpass"
        bandpass    = 4273958204    // "bandpass"

    };

    System_t            fSystem = "";

    bool                _isfunctioncalled = false; ///> check if user called for function different from base;
    bool                _residualOn = false;
    bool                _islowhighpass = true;
    bool                _hasfittedgain = false;
    bool                _hasfittedphase = false;

    // font size for calling ATLASStyle
    Size_t              tsize = 25;

    // phase and gain fit formula
    const char         *_gainfit  = ""; 
    const char         *_phasefit = "";
    NPar_t              _CutoffPar = 1;
    NPar_t              _GainPar = 0;
    NPar_t              _QPar = 2;

    const char         *label = "Preliminary";

    // computational values;
    Double_t            gGBW;       ///> Gain Bandwidth coefficient
    Double_t            gErrGBW;    ///> Gain Bandwidth coefficient error
    Double_t            gCutoff;    ///> Cutoff value 
    Double_t            gErrCutoff; ///> Cutoff value error
    Double_t            gGain;      ///> Gain value
    Double_t            gErrGain;   ///> Gain value error

    /// graphical objects []
    TGraphErrors       *fGain;
    TGraphErrors       *fPhase;
    TF1                *fGainFit; 
    TF1                *fPhaseFit;

    Float_t             legendX1 = 0.2;
    Float_t             legendY1 = 0.2;
    Float_t             legendX2 = 0.5;
    Float_t             legendY2 = 0.35;

    /// function variables declaration
    Double_t            fmin = (0.0);   ///> minimum for frequency range
    Double_t            fmax = (1.0);   ///> maximum for frequency range
    Double_t            fNpoints = -1;       ///> points in graph
    std::vector<double> fPointGain;     ///>[fNpoints] array for gain points
    std::vector<double> fPErrGain;      ///>[fNpoints] array for error gain points
    std::vector<double> fPointPhase;    ///>[fNpoints] array for phase points
    std::vector<double> fPErrPhase;     ///>[fNpoints] array for error phase points
    std::vector<double> fPointFreq;     ///>[fNpoints] array for frequency points
    std::vector<double> fPErrFreq;      ///>[fNpoints] array for error frequency points

public:
    // Bode();
    Bode(System_t sys);                                             ///> sys: OP_AMP: 0x1 high pass, 0x2 low pass; RLC 0x101 high pass, 0x102 low pass, 0x103 band pass;
    Bode(System_t sys, const char *filename, Option_t *option="");  ///> sys: OP_AMP: 0x1 high pass, 0x2 low pass; RLC 0x101 high pass, 0x102 low pass, 0x103 band pass;
    ~Bode();
    Bool_t              FitGain(Option_t *option="", Option_t *goption="", Axis_t xmin=0, Axis_t xmax=0);
    Bool_t              FitPhase(Option_t *option="", Option_t *goption="", Axis_t xmin=0, Axis_t xmax=0);
    Bool_t              FitCorrelated(Option_t *option="", Option_t *goption="", Axis_t xmin=0, Axis_t xmax=0);
    inline Double_t     GetCutoff()     const { return gCutoff; }
    inline Double_t     GetErrCutoff()  const { return gErrCutoff; }
    inline Double_t     GetErrGain()    const { return gErrGain; }
    inline Double_t     GetErrGBW()     const { return gErrGBW; }
    inline Double_t     GetGain()       const { return gGain; }
    inline Double_t     GetGBW()        const { return gGBW; }
    void                Plot(const char *filename = "", bool plotphase = true, bool plotgain = true);
    void                PlotGain(const char *filename = "");
    void                PlotPhase(const char *filename = "");
    Bool_t              ReadInput(const char *filename = "", Option_t *option="");       ///> read input for both phase and gain data 
    // bool                ReadInputGain(const char *filename, Option_t *option="");   ///> read input for gain data
    // bool                ReadInputPhase(const char *filename, Option_t *option="");  ///> read input for phase data
    // bool                ReadInputRDF()  // TO BE IMPLEMENTED
    inline void         SetCutoffNpar(NPar_t npar = 1)  { _CutoffPar = npar; }
    inline void         SetGainNpar(NPar_t npar = 0)    { _GainPar = npar; }
    Bool_t              SetFreqVec(std::vector<Double_t> Freq, std::vector<Double_t> ErrFreq);
    Bool_t              SetFunctions();
    // void                SetGainFunction(const char *formula, Option_t *option="");
    Bool_t              SetGainVec(std::vector<Double_t> Gain, std::vector<Double_t> ErrGain);
    inline void         SetLabel(Option_t *fmt) { label = fmt; }
    void                SetParGain(Double_t gain, Double_t cutoff, Double_t Q = -1);
    void                SetParPhase(Double_t gain, Double_t cutoff, Double_t Q = -1);
    // void                SetPhaseFunction(const char *formula, Option_t *option="");
    Bool_t              SetPhaseVec(std::vector<Double_t> Phase, std::vector<Double_t> ErrPhase);
    void                SetSystem(System_t sys);
    inline void         SetResidual(bool residual = true) { _residualOn = residual; }
};

#endif
