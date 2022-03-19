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

#include<TRandom.h>

#include"Logger.h"
#include"Bode/Analysis.h"

typedef TString System_t;

class SimEngine{
private:
    ULong_t             seed;
    const char         *printfname = "datasimulated.txt"; ///> printing on file tmp
    System_t            fltype;

    Double_t            gCutoff;
    Double_t            gGain;
    Double_t            gQ;

    TRandom            *gen;    ///> random generator

    enum {
        lowpass     = 244089597,    // "lowpass"
        highpass    = 2141088893,   // "highpass"
        bandpass    = 4273958204    // "bandpass"

    };

public:
    SimEngine(): gen{TRandom()} {};
    Bool_t              DataSim(const char *filename = printfname);
    void                GenLowNoise();
    void                GenHighNoise();
    Double_t            GetRandom();
    void                SetCutoff(Double_t cutoff)  { gCutoff = cutoff; }
    void                SetGain(Double_t gain)      { gGain = gain; }
    void                SetQ(Doule_t Q)             { gQ = Q; }
    void                SetFilterType(System_t filter = "lowpass");
    ~SimEngine();
};

SimEngine::~SimEngine(){
    fprintf(stderr, "%s\n", Logger::warning("Deleted obj. SimEngine"));
}


#endif