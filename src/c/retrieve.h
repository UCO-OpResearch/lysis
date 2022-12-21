#pragma once

#include "nodeGrid.h"

void printUnboundMoleculesToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName);
void printBoundMoleculesToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName);
void printIsDegradedToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName);
//void test_MarkForInspection(NodeGrid* grid);