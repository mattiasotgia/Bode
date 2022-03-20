/**
 * @file Simulate.cpp
 * @author Mattia Sotgia (mattiasotgia01@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2022-03-19
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include<TFormula.h>
#include<TRandom.h>
#include<TMath.h>

#include"BodeDataSim/SimEngine.h"

SimEngine::SimEngine(){ gen = new TRandom(); }

void SimEngine::SetFilterType(System_t filter){
    fltype = filter;
    switch(fltype.Hash()){
        case lowpass:
            _gainformula = "[0]/sqrt(1+pow(x/[1], 2))"; // [0] gain, [1] cutoff
            _phaseformula = "";
            break;
        case highpass:
            _gainformula = "[0]/sqrt(1+pow([1]/x, 2))"; // [0] gain, [1] cutoff
            _phaseformula = "";
            break;
        case bandpass:
            _gainformula = "[0]/sqrt(1+pow([2], 2)*pow(x/[1]-[1]/x, 2))"; // [0] gain, [1] cutoff/peak frequency, [2] Q factor
            _phaseformula = "";
            _islowhighpass = false;
            break;
        default:
            fprintf(stderr, "%s", Logger::warning(Form("System_t option '%s' not recognised!\n"
            "Available options are: \n\t\"lowpass\"\n\t\"highpass\"\n\t\"bandpass\"\n ", filter.Data())));
    }
}