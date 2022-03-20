/**
 * @file SimEngine.h
 * @author Mattia Sotgia (mattiasotgia01@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2022-03-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef BODEDATASIM_SimEngine
#define BODEDATASIM_SimEngine

#include<Rtypes.h>
#include<RtypesCore.h>

#include<TFormula.h>
#include<TRandom.h>

#include"Logger.h"
#include"Bode/Analysis.h"

typedef TString System_t;

class SimEngine{
private:
    ULong_t             seed = 0;
    const char         *printfname; ///> printing on file tmp
    System_t            fltype;

    Double_t            gCutoff = -99999;
    Double_t            gGain   = -99999;
    Double_t            gQ      = -99999;

    const char         *_gainformula  = ""; 
    const char         *_phaseformula = "";

    Bool_t              _islowhighpass = true;

    TRandom            *gen;    ///> random generator
    TFormula           *fFunctionFilter;

    enum {
        lowpass     = 244089597,    // "lowpass"
        highpass    = 2141088893,   // "highpass"
        bandpass    = 4273958204    // "bandpass"

    };

public:
    SimEngine();
    Bool_t              DataSim(const char *filename = "datasim.txt");
    void                GenLowNoise();
    void                GenHighNoise();
    void                SetCutoff(Double_t cutoff)  { gCutoff = cutoff; }
    void                SetGain(Double_t gain)      { gGain = gain; }
    void                SetQ(Double_t Q)             { gQ = Q; }
    void                SetFilterType(System_t filter = "lowpass");
    ~SimEngine();
};

SimEngine::~SimEngine(){
    fprintf(stderr, "%s\n", Logger::warning("Deleted obj. SimEngine"));
}


#endif