#pragma once

#include "parameters.h"

bool loadInputData(Files* files, InputData* inputData, Parameters* parameters);

//for debugging
void output2dArrayToFile(InputData* inputData, Parameters* parameters);
void outputArrayToFile(InputData* inputData, Parameters* parameters);
void outputLysesPerBlockToFile(InputData* inputData, Parameters* parameters);