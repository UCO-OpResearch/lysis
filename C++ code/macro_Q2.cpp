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
#include <cmath>
#include <fstream>
#include <sstream>

extern "C" {
#include "kiss.h"
}
    
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
// The total number of fibers in the model
const int numberOfFibers = (2*nodesInRow - 1)*nodesInColumn 
                                       + nodesInRow*(nodesInColumn - 1); // num
/*
 * 1st node in vertical direction containing fibers.
 * So if firstFiberRow=10, then rows 1-9 have no fibers,
 * there's one more row of fiber-free planar vertical edges,
 * and then the row starting with the firstFiberRow-th (e.g. 10th) vertical node
 * is a full row of fibers
 */
const int firstFiberRow = 29; // Ffree
// The last edge number without fibrin
const int lastGhostFiber = (3*nodesInRow - 1)*(firstFiberRow - 1); // enoFB
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
const int numberOfTrials = 10; // stats
// Total running time for model in seconds
const int totalTime = 10*60; // tf
// The length of one timestep, in seconds
const double timeStep =
        movingProbability*pow(poreSize,2)/(12*diffusionCoeff); // tstep
// The total number of timesteps
const int numberOfTimeSteps = totalTime / timeStep;
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

double UnbindingTimeDistribution[101];

/*
 * lysismat_PLG135_Q2.dat is a matrix with column corresponding to 
 * bin number (1-100) and with entries equal to the lysis times obtained 
 * in that bin. 
 * An entry of 6000 means lysis didn't happen.
 * i.e. the first column, lysismat(:,1), gives the lysis times 
 * for the first 100 tPA leaving times.
 */
double lysisTime[maxLysesPerBlock][lysisBlocks];

/*
 * lenlysisvect_PLG135_Q2.dat saves the first row entry in each column of 
 * lysismat_PLG135_Q2.dat where lysis did not occur, 
 * i.e., the first entry there's a 6000
 * i.e., out of the unbinds in the nth percentile (with respect to time),
 */
int lysesPerBlock[lysisBlocks];

/******************************************************************************
 ** Model State Variables
 **
 ******************************************************************************/
// The index of the fiber 


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
 * 
 */
int main(int argc, char** argv) {
    cout << "Read in data.....";
    bool success = readData();
    cout << (success ? "Done." : "Failed!") << endl;
    cout << "Parameters:" << endl;
    cout << "   N = " << nodesInRow << endl;
    cout << "   F = " << nodesInColumn << endl;
    cout << "   Ffree = " << firstFiberRow << endl;
    cout << "   num = " << numberOfFibers << endl;
    cout << "   M = " << totalMolecules << endl;
    cout << "   endFB = " << lastGhostFiber << endl;
    cout << "Obtained using code macro_Q2.cpp" << endl;
    cout << "Initializing random number generator.....";
    initializeRandomGenerator();
    cout << "Done." << endl;
    
    return 0;
}

