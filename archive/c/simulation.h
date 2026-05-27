#pragma once

void distributeMolecules(NodeGrid* grid, RankGrid* rankGrid, int totalMolecules);
void simulationLoop(NodeGrid* grid, Parameters* parameters, InputData* inputData, RankGrid* rankGrid, Files* files);
void addMultipleUnboundToFiber(Fiber* fiber, Parameters* parameters, int time, int count);
void stop(int processNum);
void debug_checkAllMolecules(NodeGrid* grid, Parameters* parameters, char* file, int line);
//void debug_checkLocalMolecules(NodeGrid* grid, char* file, int line);