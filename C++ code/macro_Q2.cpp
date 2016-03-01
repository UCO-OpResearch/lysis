/*******************************************************************************
 *******************************************************************************
 *** File:   macro_Q2.cpp
 *** Author: Dr. Brittany Bannish
 *** Converted to C++ by Dr. Bradley Paynter
 ***
 *** Created on December 18, 2015, 4:28 PM
 ***
 *** Runs the macroscale model in a clot
 *** with 72.7 nm diameter fibers and
 *** pore size. 1.0135 uM.
 *** FB conc. = 8.8 uM
 *******************************************************************************
 *******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>

extern "C" {
#include "kiss.h"
}

#define RIGHT 1
#define LEFT -1
#define UP 2
#define DOWN -2
#define OUT 3
#define IN -3
    
using namespace std;

const bool verbose = true;
/******************************************************************************
 ** Physical Parameters
 **
 ******************************************************************************/
// The tPA binding rate. Units of inverse (micromolar*sec)
const float bindingRate = 1.0e-2; // kon
// Pore size (distance betwen nodes), measured in centimeters
const float poreSize = 1.0135e-4; // delx
// Diffusion coefficient, measured in cm^2/s
const float diffusionCoeff = 5.0e-7; // Diff
// Concentration of binding sites in micromolar
const float bindingSites = 4.27e+2; // bs
/*
 * Distance from the start of one fiber to the next, in microns 
 * because distance between nodes is 1.0135 micron 
 * and diameter of 1 fiber is 0.0727 micron
 */
const float gridSymmetryDistance = 1.0862; // dist

/******************************************************************************
 ** Model Parameters
 **
 ******************************************************************************/
// The number of lattice nodes in each (horizontal) row
const int nodesInRow = 93; // N
// The number of lattice nodes in each (vertical) column
const int nodesInColumn = 121; // F
// Fibers in a full row
const int fullRow = 3*nodesInRow - 1;
// all right and out fibers in a row
const int xzRow = 2*nodesInRow - 1;
// The total number of fibers in the model
const int totalFibers = (2*nodesInRow - 1)*nodesInColumn
                                       + nodesInRow*(nodesInColumn - 1); // num
/*
 * 1st node in vertical direction containing fibers.
 * So if firstFiberRow=10, then rows 1-9 have no fibers,
 * there's one more row of fiber-free planar vertical edges,
 * and then the row starting with the firstFiberRow-th (e.g. 10th) vertical node
 * is a full row of fibers
 */
const int firstFiberRow = 29 - 1; // Ffree-1
// The last edge number without fibrin
const int lastGhostFiber = (3*nodesInRow - 1)*(firstFiberRow) - 1; // enoFB-1
/* The total number of tPA molecules:
 *      43074 is Colin's [tPA]=0.6 nM
 *      86148 is Colin's [tPA]=1.2 nM
 */
const int totalMolecules = 43074; // M
// The probability of moving. Make sure it is small enough that we've converged.
const float movingProbability = 0.2; // q

/******************************************************************************
 ** Experimental Parameters
 **
 ******************************************************************************/
// The number of independent trials to be run
const int totalTrials = 10; // stats
// Total running time for model in seconds
const int totalTime = 10*60; // tf
// The length of one timestep, in seconds
const double timeStep =
        movingProbability*pow(poreSize,2)/(12*diffusionCoeff); // tstep
// The total number of timesteps
const int totalTimeSteps = totalTime / timeStep;
// Seed for the random number generator
UINT_LEAST32_T seed = 912309035; // seed
UINT_LEAST32_T state[4] = { 129281, 362436069, 123456789, seed}; // state
/******************************************************************************
 ** Data Parameters
 **
 ******************************************************************************/
// Data file names
const string UnbindingTimeFile = "tsectPAPLG135_Q2.dat";
const string lysisTimeFile = "lysismat_PLG135_Q2.dat";
const string totalLysesFile = "lenlysisvect_PLG135_Q2.dat";
// Data size parameters
const unsigned int lysisBlocks = 100;
const unsigned int unbindsPerBlock = 500;
const unsigned int maxLysesPerBlock = 283;
/******************************************************************************
 ** The following variables will hold the microscale data
 ** which will be read in from the named data files
 ******************************************************************************/

double UnbindingTimeDistribution[101]; // tsec1;

/*
 * lysismat_PLG135_Q2.dat is a matrix with column corresponding to 
 * bin number (1-100) and with entries equal to the lysis times obtained 
 * in that bin. 
 * An entry of 6000 means lysis didn't happen.
 * i.e. the first column, lysismat(:,1), gives the lysis times 
 * for the first 100 tPA leaving times.
 */
double lysisTime[maxLysesPerBlock][lysisBlocks]; // lysismat

/*
 * lenlysisvect_PLG135_Q2.dat saves the first row entry in each column of 
 * lysismat_PLG135_Q2.dat where lysis did not occur, 
 * i.e., the first entry there's a 6000
 * i.e., out of the unbinds in the nth percentile (with respect to time),
 */
int lysesPerBlock[lysisBlocks]; // lenlysismat

/******************************************************************************
 ** Model State Variables
 **
 ******************************************************************************/
// The index of the fiber that molecule j is bound to
unsigned short location[totalMolecules]; // V(1,:)

// Whether or not molecule j is bound
bool bound[totalMolecules]; // V(2,:)

// Timestep when the molecule will bind/unbind
unsigned int unbindingTime[totalMolecules]; // bind/t_leave

// Degradation state of each edge.
bool degraded[totalFibers]; // degrade

// Timestep when the fiber will degrade
int degradeTime[totalFibers]; // t_degrade


/******************************************************************************
 ** Methods
 **
 ******************************************************************************/

/*
 * The following three methods read in the microscale data
 * from the appropriate files.
 */
bool readLysisTimeFromFile() {
    // Open data file for reading
    ifstream reader(lysisTimeFile);
    // Check if the file can be opened
    if (! reader) {
        // If not, return an error
        cerr << "Failed to open file " << lysisTimeFile << "." << endl;
        return false;
    } else {
        // If the file was successfully opened
        int i = 0;
        string line;
        while (getline(reader, line)) {
            if (i >= maxLysesPerBlock)
                cerr << "Extra rows in " << lysisTimeFile << " discarded." << endl;
            else {
                istringstream stringReader(line);
                int j = 0;
                double a;
                // Read in the data one element at a time
                while (stringReader >> a) {
                    if (j >= lysisBlocks)
                        cerr << "Extra entries in row " << i 
                                << " of " << lysisTimeFile << " discarded." << endl;
                    else
                        lysisTime[i][j++] = a;
                }
                if (j < lysisBlocks) {
                    cerr << "Insufficient entries in row " << i 
                            << " of " << lysisTimeFile << "." << endl;
                    return false;
                }
            }
            i++;
        }
        // Close the file
        reader.close();
        if (i < maxLysesPerBlock) {
            cerr << "Insufficient rows in " << lysisTimeFile << ":("
                << i << ", " << maxLysesPerBlock << ")" << endl;
            return false;
        }
    }
    return true;
}

/*
 * This method reads the unbinding time given in a file.
 * Extra entries are discarded, while an amount of files
 * too small will close the file.
 */
bool readUnbindingTimeFromFile() {
    int i = 0;
    // Open data file for reading
    ifstream reader(UnbindingTimeFile);
    // Check if the file can be opened
    if (! reader) {
        // If not, return an error
        cerr << "Failed to open file " << UnbindingTimeFile << "." << endl;
        return false;
    } else {
        // If the file was successfully opened
        double a;
        // Read in the data one T at a time
        while (reader >> a) {
            if (i >= 101)
                cerr << "Extra entries in " << UnbindingTimeFile
                        << " discarded." << endl;
            else
                UnbindingTimeDistribution[i++] = a;
        }
        // Close the file
        reader.close();
        if (i < 101) {
            cerr << "Insufficient entries in " << UnbindingTimeFile << "." << endl;
            return false;
        }
    }
    return true;
}

/*
 * This method reads the Lyses per block from
 * a file. Extra files will be removed, while
 * an amount of files too small will close the file.
 */
bool readLysesPerBlockFromFile() {
    int i = 0;
    // Open data file for reading
    ifstream reader(totalLysesFile);
    // Check if the file can be opened
    if (! reader) {
        // If not, return an error
        cerr << "Failed to open file " << totalLysesFile << "." << endl;
        return false;
    } else {
        // If the file was successfully opened
        int a;
        // Read in the data one T at a time
        while (reader >> a) {
            if (i >= lysisBlocks)
                cerr << "Extra entries in " << totalLysesFile
                << " discarded." << endl;
            else
                lysesPerBlock[i++] = a;
        }
        // Close the file
        reader.close();
        if (i < lysisBlocks) {
            cerr << "Insufficient entries in " << totalLysesFile << "." << endl;
            return false;
        }
    }
    return true;
}

/*
 * This method reads in the appropriate data from disk.
 * It returns "true" if all data read in correctly, 
 * "false" if any IO errors occur. 
 */
bool readData() {
    bool success;
    success = readLysisTimeFromFile();
    if (!success)
        return false;
    success = readUnbindingTimeFromFile();
    if (!success)
        return false;
    success = readLysesPerBlockFromFile();
    if (!success)
        return false;
    return true;
}

/*
 * Set up the random number generator 'kiss32' with the values from 'state'
 */
void initializeRandomGenerator() {
    set_kiss32_(state);
    get_kiss32_(state);
}

/*
 * Returns a random integer between 0 and N-1 inclusive.
 */
unsigned int randomInt(int N) {
    double t = urcw1_();
    // cout << t << " -> " << (unsigned int)(t*N) << endl;
    return (unsigned int)(t*N);
}

/*
 * Returns the (single-dimensional) index of the fiber 
 * leaving node (nodeX, nodeY) in direction 'fiber'.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberIndex(unsigned short nodeX, unsigned short nodeY, char direction) {
    if ((nodeX < 0) || (nodeY < 0))
        throw invalid_argument("Index must be non-negative.");
    else if (nodeX >= nodesInRow)
        throw invalid_argument("Index out of bounds.");
    else if (nodeY >= nodesInColumn)
        throw invalid_argument("Index out of bounds.");
    else
        switch (direction) {
            case DOWN:
                return (nodeY == 0) ? -1 : fiberIndex(nodeX, nodeY-1, UP);
            case LEFT:
                return (nodeX == 0) ? -1 : fiberIndex(nodeX-1, nodeY, RIGHT);
            case IN:
                return fiberIndex(nodeX, nodeY, OUT);
            case UP:
                return (nodeY == nodesInColumn - 1) ? -1 :
                    nodeY * fullRow + xzRow + nodeX;
            case RIGHT:
                return (nodeX == nodesInRow - 1) ? -1 :
                    nodeY * fullRow + nodeX * 2 + 1;
            case OUT:
                return nodeY * fullRow + nodeX * 2;
            default:
                return -1;
        }
}

/*
 * Returns the x-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeX(unsigned short fiberIndex) {
    if (fiberIndex < 0)
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;
    return (rowPosition >= xzRow) ? rowPosition - xzRow :
        (unsigned short)(rowPosition / 2);
}

/*
 * Returns the y-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeY(unsigned short fiberIndex) {
    if (fiberIndex < 0)
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)
        throw invalid_argument("Index out of bounds.");
    return (unsigned short)(fiberIndex / fullRow);
}

/*
 * Returns the direction of a fiber given the 
 * (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberDirection(unsigned short fiberIndex) {
    if (fiberIndex < 0)
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;
    return (rowPosition >= xzRow) ? UP : (rowPosition % 2 == 0 ? OUT : RIGHT);
}

/*
 * Gets the neighbors surrounding the selected fiber, given the
 * index of a fiber.
 */
void getNeighbors(unsigned short fiber, unsigned short neighbors[]) {
    unsigned short x = nodeX(fiber);
    unsigned short y = nodeY(fiber);
    unsigned short dir = fiberDirection(fiber);
    switch (dir) {
        case UP:
            neighbors[0] = (x != 0) ? fiberIndex(x,y,LEFT) : fiberIndex(x,y,RIGHT);
            neighbors[1] = fiberIndex(x, y, OUT);
            neighbors[2] = fiberIndex(x, y, IN);
            neighbors[3] = (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT);
            neighbors[4] = (x != 0) ? fiberIndex(x,y+1,LEFT) : fiberIndex(x,y+1,RIGHT);
            neighbors[5] = fiberIndex(x, y+1, OUT);
            neighbors[6] = fiberIndex(x, y+1, IN);
            neighbors[7] = (x != nodesInRow-1) ? fiberIndex(x, y+1, RIGHT) : fiberIndex(x, y+1, LEFT);
            break;
        case RIGHT:
            neighbors[0] = (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP);
            neighbors[1] = fiberIndex(x, y, OUT);
            neighbors[2] = fiberIndex(x, y, IN);
            neighbors[3] = (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN);
            neighbors[4] = (y != 0) ? fiberIndex(x+1, y, DOWN) : fiberIndex(x+1, y, UP);
            neighbors[5] = fiberIndex(x, y, OUT);
            neighbors[6] = fiberIndex(x, y, IN);
            neighbors[7] = (y != nodesInColumn-1) ? fiberIndex(x+1, y, UP) : fiberIndex(x+1, y, DOWN);
            break;
        case OUT:
            neighbors[0] = (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT); // z
            neighbors[1] = (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT); // z+1
            neighbors[2] = (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT); // z
            neighbors[3] = (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT); // z+1
            neighbors[4] = (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // z
            neighbors[5] = (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // z+1
            neighbors[6] = (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // z
            neighbors[7] = (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // z+1
            break;
        default:
            throw invalid_argument("Something went wrong while finding neighbors.");
    }
}

/*
 * Output all parameter (constant) values to stout.
 */
void printParameters() {
    cout << "Model Parameters:" << endl;
    cout << "   (N)         Nodes per row                   = " << nodesInRow << endl;
    cout << "   (F)         Nodes per column                = " << nodesInColumn << endl;
    cout << "   (Ffree-1)   First row of fibers             = " << firstFiberRow << endl;
    cout << "   (num)       Total fibers                    = " << totalFibers << endl;
    cout << "   (M)         Total molecules                 = " << totalMolecules << endl;
    cout << "   (enoFB-1)   Index of the last ghost fiber   = " << lastGhostFiber << endl;
    cout << "   (tstep)     Length of timestep              = " << timeStep << " sec" << endl;
    cout << "   (num_t)     Total timesteps                 = " << totalTimeSteps << endl;
    cout << "   (tf)        Total time                      = " << totalTime << " sec" << endl;
    cout << "   (seed)      Random number generator seed    = " << seed << endl;
    cout << "   (q)         Molecule moving probability     = " << movingProbability << endl;
    cout << "   (delx)      Distance between nodes          = " << poreSize << " cm" << endl;
    cout << "   (kon)       Binding rate                    = " << bindingRate << " micromolar*sec" << endl;
    cout << "   (bs)        Concentration of binding sites  = " << bindingSites << " micromolar" << endl;
}

/*
 * Sets the model state variables to their initial values, that is,
 * All fibers are undegraded and have no degrade time set
 * All molecules are unbound, have no timer set, and their locations
 * are uniformly distributed on the ghost fibers
 */
void initializeVariables() {
    for (int i = 0; i < totalFibers; i++) {
        degraded[i] = false;
        degradeTime[i] = totalTimeSteps + 1;
    }
    for (int i = 0; i < totalMolecules; i++) {
        bound[i] = false;
        location[i] = randomInt(lastGhostFiber + 1);
        unbindingTime[i] = totalTimeSteps + 1;
    }
}

/*
 * Degrades fibers for a given value.
 */
void degradeFibers(unsigned int t) {
    for (int i = lastGhostFiber+1; i < totalFibers; i++) {
        if (degradeTime[i] == t) {
            degraded[i] = true;
            for (int j = 0; j < totalMolecules; j++) {
                if (location[j] == i) {
                    bound[j] = false;
                    unbindingTime[j] = t;
                }
            }
        }
    }
}

/*
 * Runs the BloodClotting Model.
 */
void runModel() {
    for (unsigned int t = 0; t < totalTimeSteps; t++) {
        if ((t % 400000 == 0) && verbose)
            cout << "Time: " << t * timeStep << " sec" << endl;
        degradeFibers(t);
        for (unsigned int j = 0; j < totalMolecules; j++) {
            if (bound[j] && (unbindingTime[j] == t)) {
                bound[j] = false;
                unbindingTime[j] = t - log(urcw1_()) / (bindingRate * bindingSites);
            }
            if (!bound[j]) {
                if (urcw1_() <= 1-movingProbability) {
                    
                }
            }
        }
    }
}

/******************************************************************************
 ** Main
 **
 ******************************************************************************/
int main(int argc, char** argv) {
    if (verbose) cout << "Read in data.............................";
    bool success = readData();
    if (verbose) cout << (success ? "Done." : "Failed!") << endl;
    cout << "Output obtained using code macro_Q2.cpp" << endl;
    printParameters();
    if (verbose) cout << "Initializing random number generator.....";
    initializeRandomGenerator();
    if (verbose) cout << "Done." << endl;
    if (verbose) cout << "Initializing variables...................";
    initializeVariables();
    if (verbose) cout << "Done." << endl;
    return 0;
}


template <class T>
struct heap {
    
};

