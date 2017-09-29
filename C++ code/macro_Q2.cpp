/*******************************************************************************
 *******************************************************************************
 *** File:   macro_Q2.cpp
 *** Author: Dr. Brittany Bannish
 *** Converted to C++ by Dr. Bradley Paynter and C. Beadle
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
#include <ncurses.h>
#include <stdexcept>
#include <fstream>
#include <sstream>

// This is the pseudo-random number generator code. It is written in C.
extern "C" {
#include "kiss.h"
}

// These constants allow the code to be more human-readable
#define RIGHT 1
#define LEFT -1
#define UP 2
#define DOWN -2
#define OUT 3
#define IN -3
    
using namespace std;

/******************************************************************************
 ******************************************************************************
 ****
 **** PARAMETERS
 ****
 ******************************************************************************
 ******************************************************************************/

// This flag determines how much information the code outputs
const bool verbose = true;
/******************************************************************************
 ** Physical Parameters
 ** i.e., those parameters that are measured from the real world
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
 ** Model Size Parameters
 ** i.e., those parameters selected to define 
 **       the physical properties of the model
 ******************************************************************************/
// The number of lattice nodes in each (horizontal) row
const int nodesInRow = 30; // N (93)
// The number of lattice nodes in each (vertical) column
const int nodesInColumn = 11; // F (121)
// Fibers in a full row
const int fullRow = 3*nodesInRow - 1;
// all right and out fibers in a row
const int xzRow = 2*nodesInRow - 1;
// The total number of fibers in the model
const int totalFibers = (2*nodesInRow - 1)*nodesInColumn
                                       + nodesInRow*(nodesInColumn - 1); // num
/*
 * Only part of the grid is occupied by the clot (fibers). The bottom part of the 
 * grid is occupied with 'ghost' fibers which serve as starting locations for the
 * tPA molecules.
 * The constant 'firstFiberRow' is the index of the 1st node in vertical 
 * direction containing fibers. So if firstFiberRow=10, then rows 1-9 have no 
 * fibers, there's one more row of fiber-free planar vertical edges, and then the 
 * row starting with the firstFiberRow-th (e.g. 10th) vertical node is a full row 
 * of fibers
 */
const int firstFiberRow = 4 - 1; // Ffree-1 (29-1)
// The last edge number without fibrin
const int lastGhostFiber = (3*nodesInRow - 1)*(firstFiberRow) - 1; // enoFB-1
/* The total number of tPA molecules:
 *      43074 is Colin's [tPA]=0.6 nM
 *      86148 is Colin's [tPA]=1.2 nM
 */
const int totalMolecules = 500; // M (43074)
// The probability of moving. Make sure it is small enough that we've converged.
const float movingProbability = 0.5; // q (0.2)

/******************************************************************************
 ** Model Time Parameters
 ** i.e., those parameters selected to define 
 **       the running time and granularity of the model
 ******************************************************************************/
// The number of independent trials to be run
const int totalTrials = 1; // stats (10)
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
 ** Output Parameters
 **
 ******************************************************************************/
const int gridSize = 1;

/******************************************************************************
 ******************************************************************************
 ****
 **** VARIABLES
 ****
 ******************************************************************************
 ******************************************************************************/


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
unsigned int bindUnbindTime[totalMolecules]; // bind/t_leave

// Degradation state of each edge.
bool degraded[totalFibers]; // degrade

// Timestep when the fiber will degrade
int degradeTime[totalFibers]; // t_degrade

// A linked list of molecules located at the current fiber
// Each fiber points to the first molecule located on it.
// Each molecule points to the next molecule located on the same fiber.
// nextMolecule is -1 if it is the last in the list.
unsigned int firstMolecule[totalFibers];
int nextMolecule[totalMolecules];


/******************************************************************************
 ** Output Data Storage Variables
 **
 ******************************************************************************/
 
//
int successfulBinds;

//
int unsuccessfulBinds;

//For the set of vectors for saveData.
//struct counterSet{
//	short locationCheck;
//	bool boundCheck;
//	int timeCheck;
//	bool degradedCheck;
//	int totalBindingsCheck;
//	int totalaAttemptsCheck;
//}
// A vector being used for processData.
//struct dvSet{
//	int x;
//	int y;
//	bool deg;
//}

//struct checkSet{
//	bool degCheck;
//}
//
//vector<counterSet> dataVec;
//
//vector<dvSet> degreeVector;
//
//vector<checkSet> degradeCheckVector;


/******************************************************************************
 ******************************************************************************
 ****
 **** METHODS
 ****
 ******************************************************************************
 ******************************************************************************/

/******************************************************************************
 ** Initialization methods 
 **
 ******************************************************************************/

/*
 * The following three methods read in the microscale data
 * from the appropriate files.
 * Each method returns "false" if the data file cannot be read, 
 * or if it does not contain the correct amount of data.
 * Each method returns "true" if the data is read in correctly.
 */
bool readLysisTimeFromFile() {
    ifstream reader(lysisTimeFile);													// Open data file for reading
    if (! reader) {																	// Check if the file can be opened
        cerr << "Failed to open file " << lysisTimeFile << "." << endl;
        return false;
    } else {																		// If the file was successfully opened
        int i = 0; 																		// Make a var to store the current line index
        string line; 																	// Make a var to hold the current line
        while (getline(reader, line)) { 												// Get the next line from the file
            if (i >= maxLysesPerBlock) 														// Check that we haven't exceeded the expected number of lines
                cerr << "Extra rows in " << lysisTimeFile << " discarded." << endl;
            else {
                istringstream stringReader(line); 											// Set up the line for reading
                int j = 0; 																	// Make a var to store the current entry's index
                double a; 																	// Make a var to store the current entry
                while (stringReader >> a) { 												// Get the next entry from the current line
                    if (j >= lysisBlocks) 														// Check that we haven't exceeded the expected number of entries in this line
                        cerr << "Extra entries in row " << i 
                                << " of " << lysisTimeFile << " discarded." << endl;
                    else
                        lysisTime[i][j++] = a; 													// Store the data in the appropriate position in the array
                }
                if (j < lysisBlocks) {														// Check if we've reached the end of the line with insufficient entries
                    cerr << "Insufficient entries in row " << i 
                            << " of " << lysisTimeFile << "." << endl;
                    return false;
                }
            }
            i++;																			// Increment the line counter
        }
        reader.close();																	// Close the file
        if (i < maxLysesPerBlock) {														// Check if we've reached the end of the file with insufficient rows
            cerr << "Insufficient rows in " << lysisTimeFile << ":("
                << i << ", " << maxLysesPerBlock << ")" << endl;
            return false;
        }
    }
    return true;																	// Everything was read in correctly
}

/*
 * This method reads the unbinding time given in a file.
 * Extra entries are discarded, while an amount of files
 * too small will close the file.
 */
bool readUnbindingTimeFromFile() {
    ifstream reader(UnbindingTimeFile);												// Open data file for reading
    if (! reader) {																	// Check if the file can be opened
        cerr << "Failed to open file " << UnbindingTimeFile << "." << endl;				// If not, return an error
        return false;
    } else {																		// If the file was successfully opened
    	int i = 0;																		// Make a var to store the current entry index
        double a;																		// Make a var to store the current entry
        while (reader >> a) {															// Read in the data one at a time
            if (i >= 101)																	// Check if we've exceeded the expected number of entries
                cerr << "Extra entries in " << UnbindingTimeFile
                        << " discarded." << endl;
            else
                UnbindingTimeDistribution[i++] = a;											// Store the data in the appropriate position in the array
        }
        reader.close();																	// Close the file
        if (i < 101) {
            cerr << "Insufficient entries in " << UnbindingTimeFile << "." << endl;		// Check if we've reached the end of the file with insufficient entries
            return false;
        }
    }
    return true;																	// Everything was read in correctly
}

/*
 * This method reads the Lyses per block from
 * a file. Extra files will be removed, while
 * an amount of files too small will close the file.
 */
bool readLysesPerBlockFromFile() {																		
    ifstream reader(totalLysesFile);												// Open data file for reading
    if (! reader) {																	// Check if the file can be opened
        cerr << "Failed to open file " << totalLysesFile << "." << endl;				// If not, return an error
        return false;
    } else {																		// If the file was successfully opened
    	int i = 0;																		// Make a var to store the current entry index
        int a;																			// Make a var to store the current entry
        while (reader >> a) {															// Read in the data one entry at a time
            if (i >= lysisBlocks)															// Check if we've exceeded the expected number of entries
                cerr << "Extra entries in " << totalLysesFile
                << " discarded." << endl;
            else
                lysesPerBlock[i++] = a;														// Store the current entry in the array
        }
        reader.close();																	// Close the file
        if (i < lysisBlocks) {															// Check if we've reached the end of the file with insufficient entries
            cerr << "Insufficient entries in " << totalLysesFile << "." << endl;
            return false;
        }
    }
    return true;																	// Everything was read in correctly
}

/*
 * This method reads in the appropriate data from disk.
 * It returns "true" if all data read in correctly, 
 * "false" if any IO errors occur. 
 */
bool readData() {
    bool success;
    success = readLysisTimeFromFile(); 												// Read the first data file
    if (!success)  																	// Check if it was read successfully
        return false;
    success = readUnbindingTimeFromFile();											// Read the next data file
    if (!success)  																	// Check if it was read successfully
        return false; 
    success = readLysesPerBlockFromFile();											// Read the next data file
    if (!success)  																	// Check if it was read successfully
        return false; 
    return true; 																	// If all are successful, then readData returns true
}

/*
 * Set up the random number generator 'kiss32' with the values from 'state'
 */
void initializeRandomGenerator() {
    set_kiss32_(state); 															// Initialize the random number generator at "state"
    get_kiss32_(state);																// Read the current state of the random number generator
}

/*
 * Returns a random integer between 0 and N-1 inclusive.
 */
unsigned int randomInt(int N) {
    double t = urcw1_();															// Call the method from 'kiss32' and store the resulting double
    return (unsigned int)(t*N);														// Cast the double to an unsigned int and return it
}



/******************************************************************************
 ** Fiber/Node indexing methods 
 **
 ******************************************************************************
 * In this model, nodes and fibers are arranged as follows
 *    /       /       /       /       /       /       /       /       /       /
 *   O-------O-------O-------O-------O-------O-------O-------O-------O-------O
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *   O-------O-------O-------O-------O-------O-------O-------O-------O-------O
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *   O-------O-------O-------O-------O-------O-------O-------O-------O-------O
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *   O-------O-------O-------O-------O-------O-------O-------O-------O-------O
 *  /       /       /       /       /       /       /       /       /       /
 *
 * Fibers can be referenced by their node O and their direction from that node
 *
 * The directions are
 *        UP  IN   
 *         | /
 *         |/   
 *  LEFT---O---RIGHT
 *        /|    
 *       / |
 *    OUT  DOWN
 *
 * Note that, in order to reduce the problem to just 2 dimensions,
 * the IN and OUT fibers are combined into a single fiber.
 *
 * Nodes are given x- and y-coordinates. 
 * A node's x-coordinate gives number of nodes to its left in its row
 * A node's y-coordinate gives the number of nodes below in its column
 * i.e.,
 *    /       /       /       /       /       /       /       /       /       /
 * (3,0)---(3,1)---(3,2)---(3,3)---(3,4)---(3,5)---(3,6)---(3,7)---(3,8)---(3,9)
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 * (2,0)---(2,1)---(2,2)---(2,3)---(2,4)---(2,5)---(2,6)---(2,7)---(2,8)---(2,9)
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 * (1,0)---(1,1)---(1,2)---(1,3)---(1,4)---(1,5)---(1,6)---(1,7)---(1,8)---(1,9)
 *  /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |       |       |       |       |       |       |       |       |       |
 *   |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 * (0,0)---(0,1)---(0,2)---(0,3)---(0,4)---(0,5)---(0,6)---(0,7)---(0,8)---(0,9)
 *  /       /       /       /       /       /       /       /       /       /
 *
 * Fibers are stored in a single-dimensional array. 
 * They are inserted into the array by row. 
 * Each node in the row has their OUT and RIGHT fiber added to the array.
 * Once that is done for the whole row, then the same row's UP fibers are added.
 * i.e.,
 *     /       /       /       /       /       /       /       /       /       /
 *    O---88--O---90--O---92--O---94--O---96--O---98--O--100--O--102--O--103--O
 *   /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 * 87 |    89 |    91 |    93 |    95 |    97 |    99 |   101 |   103 |   105 |
 *    77      78      79      80      81      82      83      84      85      86
 *    |       |       |       |       |       |       |       |       |       |
 *    |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *    O---59--O---61--O---63--O---65--O---67--O---69--O---71--O---73--O---75--O
 *   /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 * 58 |    60 |    62 |    64 |    66 |    68 |    70 |    72 |    74 |    76 |
 *    48      49      50      51      52      53      54      55      56      57
 *    |       |       |       |       |       |       |       |       |       |
 *    |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *    O---30--O---32--O---34--O---36--O---38--O---40--O---42--O---44--O---46--O
 *   /|      /|      /|      /|      /|      /|      /|      /|      /|      /|
 * 29 |    31 |    33 |    35 |    37 |    39 |    41 |    43 |    45 |    47 |
 *    19      20      21      22      23      24      25      26      27      28
 *    |       |       |       |       |       |       |       |       |       |
 *    |/      |/      |/      |/      |/      |/      |/      |/      |/      |/
 *    O---1---O---3---O---5---O---7---O---9---O---11--O---13--O---15--O---17--O
 *   /       /       /       /       /       /       /       /       /       /
 *  0       2       4       6       8      10      12      14      16      18
 *
 * The methods below convert between the node-and-direction and list 
 * representation of a fiber's location. That way all code can use the 
 * human-readable '(2,5,UP)' while computer memory uses the more efficient '82'. 
 ******************************************************************************/


/*
 * Returns the (single-dimensional) index of the fiber 
 * leaving node (nodeX, nodeY) in direction 'fiber'.
 * If an impossible request is made, the method returns '-1'.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberIndex(unsigned short nodeX, unsigned short nodeY, char direction) {
    if ((nodeX < 0) || (nodeY < 0)) 												// Checks that the node locations are non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (nodeX >= nodesInRow) 													// Checks that the node x-coordinate is within bounds
        throw invalid_argument("Index out of bounds.");
    else if (nodeY >= nodesInColumn) 												// Checks that the node y-coordinate is within bounds
        throw invalid_argument("Index out of bounds.");
    else
    	// This section only calculates the fiber index for UP, RIGHT, and OUT. 
    	// If one of the other directions are requested, it moves to a neighboring node
    	// and converts the request.
    	// i.e., if DOWN is requested, we recursively call the method again 
    	// asking for the UP fiber from the node below.
        switch (direction) { 														// Determine the direction of the fiber requested
            case DOWN:
                return (nodeY == 0) ? -1 : fiberIndex(nodeX, nodeY-1, UP); 			// Ensure that we are not at the bottom edge of the grid, then recurse
            case LEFT:                                                     
                return (nodeX == 0) ? -1 : fiberIndex(nodeX-1, nodeY, RIGHT);		// Ensure that we are not on the left edge of the grid, then recurse
            case IN:
                return fiberIndex(nodeX, nodeY, OUT); 								// Every node has an IN fiber, so just recurse
            case UP:
                return (nodeY == nodesInColumn - 1) ? -1 : 							// Ensure that we are not on the top edge of the grid
                    nodeY * fullRow + 													// The index is all fibers from all rows below, 
                     xzRow + 															// plus the x- and z-fibers from the current row,
                     nodeX;																// plus the y-fibers for all nodes to the left in this row.
            case RIGHT:
                return (nodeX == nodesInRow - 1) ? -1 : 							// Ensure that we are not on the right edge of the grid
                    nodeY * fullRow + 													// The index is all fibers from all rows below,
                     nodeX * 2 + 														// plus all the x- and z-fibers for all nodes to the left in this row,
                     1;																	// plus the z-fiber for the current node
            case OUT:
                return nodeY * fullRow + 											// The index is all fibers from all rows below,
                 nodeX * 2;															// plus all the x- and z-fibers for all nodes to the left in this row.
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
    if (fiberIndex < 0)  															// Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) 											// Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;   							// Find the fiber's position in its row 
    if (rowPosition >= xzRow)  														// If this fiber is in the second (y-fiber) part of its row
    	return rowPosition - xzRow;														// Just subtract the number of x- and z-fibers to get the node index
    else 																			// Else if this fiber is in the first (x-, z-fiber) part of its row
        return (unsigned short)(rowPosition / 2);										// Its node index is half its fiber index (rounded/cast) down
}

/*
 * Returns the y-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeY(unsigned short fiberIndex) {
    if (fiberIndex < 0)  															// Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) 											// Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    return (unsigned short)(fiberIndex / fullRow); 									// Determine the number of full rows preceding this fiber
}

/*
 * Returns the direction of a fiber given the 
 * (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberDirection(unsigned short fiberIndex) {
    if (fiberIndex < 0) 															// Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) 											// Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;								// Find the fiber's position in its row 
    if (rowPosition >= xzRow)   													// If this fiber is in the second (y-fiber) part of its row
    	return UP;																		// It must be an UP fiber
    else if (rowPosition % 2 == 0)													// Else, the even fibers are OUT and the odd ones are RIGHT
    	return OUT;
    else
    	return RIGHT; 
}

/*
 * Given the index of a fiber, this method returns the i'th neighbor of that fiber.
 * Each fiber has eight neighbors. The location of these neighbors depends which 
 * plane the fiber is in:
 * For an x-fiber the arrangement is as follows. 
 * Note that, if the fiber is at the top or bottom of the grid, the opposite fiber 
 * is substituted for the missing one.
 * i.e., if the fiber is at the top of the grid, neighbors 3 and 7 do not exist, so
 * fibers 0 and 4 occur as neighbors twice, once in their regular positions and once 
 * each in place of the non-existent fibers.
 *
 *        3  2    7  6
 *        | /     | /
 *        |/      |/
 *     ---O-------O---        
 *       /|      /|
 *      / |     / |
 *     1  0    5  4
 *
 * For a y-fiber the arrangement is as follows. 
 * Note that, if the fiber is at the left or right edge of the grid, the opposite fiber 
 * is substituted for the missing one. See the example for the x-fiber for more info.
 *
 *             6
 *          | /
 *          |/  
 *     4----O----7
 *         /|
 *        / |  2
 *       5  | /
 *          |/
 *     0----O----3
 *         /|
 *        / |
 *       1
 *
 * For a z-fiber the arrangement is as follows.
 * Since the current model only uses one layer of nodes, each neighbor occurs twice.
 * If the node is on an edge of the grid, the fibers are duplicated as in the x- and 
 * y-fiber cases above. This means that an individual fiber would occur 4 times as a neighbor.
 * i.e., If the node was at the top edge of the grid, the same DOWN fiber would be 
 * neighbors 4,5,6, & 7
 *
 *       6,7  
 *         | /
 *         |/   
 *   0,1---O---2,3
 *        /|    
 *       / |
 *         4,5
 *
 */
int getNeighbor(unsigned short fiber, unsigned short i){
	unsigned short x = nodeX(fiber); 												// Get the x-coordinate of the fiber's canonical node
    unsigned short y = nodeY(fiber); 												// Get the y-coordinate of the fiber's canonical node
    unsigned short dir = fiberDirection(fiber); 									// Get the direction of the fiber from its canonical node
	switch (dir) {
        case UP:
			switch(i){
				case 0:
				return (x != 0) ? fiberIndex(x,y,LEFT) : fiberIndex(x,y,RIGHT);						// Account for the left edge case
				case 1: 																			
				return fiberIndex(x, y, OUT); 
				case 2:
				return fiberIndex(x, y, IN);
				case 3: 
				return  (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT);		// Account for the right edge case
				case 4: 
				return (x != 0) ? fiberIndex(x,y+1,LEFT) : fiberIndex(x,y+1,RIGHT);
				case 5: 
				return  fiberIndex(x, y+1, OUT);
				case 6:
				return  fiberIndex(x, y+1, IN);
				case 7:
				return (x != nodesInRow-1) ? fiberIndex(x, y+1, RIGHT) : fiberIndex(x, y+1, LEFT); 
			}
		    case RIGHT:
			switch(i){
				case 0: 
				return (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); 					// Account for the bottom edge case
				case 1: 
				return fiberIndex(x, y, OUT);
				case 2:
				return fiberIndex(x, y, IN);
				case 3: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); 		// Account for the top edge case
				case 4: 
				return (y != 0) ? fiberIndex(x+1, y, DOWN) : fiberIndex(x+1, y, UP);
				case 5: 
				return fiberIndex(x, y, OUT);
				case 6:  
				return fiberIndex(x, y, IN);
				case 7: 
				return (y != nodesInColumn-1) ? fiberIndex(x+1, y, UP) : fiberIndex(x+1, y, DOWN);
			}
			 case OUT:
			switch(i){
				case 0: 
				return (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT); 				// Account for the left edge case
				case 1: 
				return (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT);
				case 2:
				return (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT); 		// Account for the right edge case
				case 3: 
				return (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT);
				case 4: 
				return  (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); 					// Account for the bottom edge case
				case 5: 
				return  (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP);
				case 6: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); 		// Account for the top edge case
				case 7: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN);
			}
			default: 																				// If something goes wrong with the method
            throw invalid_argument("Something went wrong while finding neighbors.");
		}
}

/*
 * Get the entire neighborhood of a fiber. This is returned in the 'neighbors' array
 */
void getNeighbors(unsigned short fiber, unsigned short neighbors[]) {
  for (int i = 0; i < 8; i++){														// Loop through the eight neighbors
		neighbors[i] = getNeighbor(fiber,i);										// Calculate each neighbor and store in the array
  }
}


/******************************************************************************
 ** Output methods 
 **
 ******************************************************************************/

const int black = 1;
const int red = 2;
const int green = 3;
const int yellow = 4;
const int blue = 5;
const int white = 6;

void initializeCurses() {
	initscr();		
	noecho();
 	curs_set(FALSE);
	start_color();	
	init_pair(black, COLOR_BLACK, COLOR_BLACK);
	init_pair(red, COLOR_RED, COLOR_BLACK);
	init_pair(green, COLOR_GREEN, COLOR_BLACK);
	init_pair(yellow, COLOR_YELLOW, COLOR_BLACK);
	init_pair(blue, COLOR_BLUE, COLOR_BLACK);
	init_pair(white, COLOR_WHITE, COLOR_BLACK);
}




string color(string in, string color) {
	// 
	// attron(COLOR_PAIR(1));
	// printw(in);
	// attroff(COLOR_PAIR(1));
	return in;
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
 * Prints the location of an individual molecule, its status, and the status of the fiber it is on
 */
void printMoleculeStatus(int testMolecule, unsigned int t) {
    cout << "Time: " << t * timeStep << " sec" << endl;
	cout << "Location of molecule #" << testMolecule << " is fiber #" << location[testMolecule] << endl;
	cout << "Molecule #" << testMolecule << " is " << (bound[testMolecule] ? "bound" : "unbound") << endl;
	cout << "Fiber #" << location[testMolecule] << " is " << (degraded[location[testMolecule]] ? "degraded" : "not degraded") << endl;
}

/*
 *
 */
void printNode(int x, int y) {
	int maxY = 0, maxX = 0;
	getmaxyx(stdscr, maxY, maxX);
	int baseY = maxY - y*(gridSize*2 + 2) - 5;
	int baseX = x*(gridSize*2 + 2) + 1;
	attron(COLOR_PAIR(white));
	mvprintw(baseY, baseX, "O");
	if (y >= firstFiberRow) {
		if (x != 0) {
			if (degraded[fiberIndex(x,y,LEFT)])
				attron(COLOR_PAIR(red));
			else
				attron(COLOR_PAIR(blue));
			for (int i = 0; i < gridSize+1; i++)
				mvprintw(baseY, baseX-1 - i, "-");
		}
		if (x != nodesInRow-1) {
			if (degraded[fiberIndex(x,y,RIGHT)])
				attron(COLOR_PAIR(red));
			else
				attron(COLOR_PAIR(blue));
			for (int i = 0; i < gridSize; i++)
				mvprintw(baseY, baseX+1 + i, "-");
		}
		if (y != firstFiberRow) {
			if (degraded[fiberIndex(x,y,DOWN)])
				attron(COLOR_PAIR(red));
			else
				attron(COLOR_PAIR(blue));
			for (int i = 0; i < gridSize+1; i++)
				mvprintw(baseY+1 + i, baseX, "|");
		}
		if (y != nodesInColumn-1) {
			if (degraded[fiberIndex(x,y,UP)])
				attron(COLOR_PAIR(red));
			else
				attron(COLOR_PAIR(blue));
			for (int i = 0; i < gridSize; i++)
				mvprintw(baseY-1 - i, baseX, "|");
		}
		if (degraded[fiberIndex(x,y,OUT)])
			attron(COLOR_PAIR(red));
		else
			attron(COLOR_PAIR(blue));
		mvprintw(baseY+1, baseX-1, "/");
		if (degraded[fiberIndex(x,y,IN)])
			attron(COLOR_PAIR(red));
		else
			attron(COLOR_PAIR(blue));
		mvprintw(baseY-1, baseX+1, "/");
	}
}

void printMolecule(int j) {
	int maxY = 0, maxX = 0;
	getmaxyx(stdscr, maxY, maxX);
	int x = nodeX(location[j]);
	int y = nodeY(location[j]);
	int baseY = maxY - y*(gridSize*2 + 2) - 5;
	int baseX = x*(gridSize*2 + 2) + 1;
	if (bound[j])
		attron(COLOR_PAIR(yellow));
	else
		attron(COLOR_PAIR(green));
	int direction = fiberDirection(location[j]);
	switch (direction) {
		case UP:
			mvprintw(baseY-2, baseX, "*");
			break;
		case RIGHT:
			mvprintw(baseY, baseX+2, "*");
			break;
		case OUT:
			mvprintw(baseY+1, baseX-1, "*");
			break;
		default:
			mvprintw(baseY, baseX, "X");
	}
}

/*
 *
 */
void printGrid(int t) {
	//clear();
	int maxY = 0, maxX = 0;
	getmaxyx(stdscr, maxY, maxX);
	maxX = nodesInRow*(gridSize*2 + 2) + 10;
	attron(COLOR_PAIR(red));
	mvprintw(10, maxX, "|");
	attron(COLOR_PAIR(blue));
	mvprintw(12, maxX, "|");
	attron(COLOR_PAIR(green));
	mvprintw(14, maxX, "*");
	attron(COLOR_PAIR(yellow));
	mvprintw(16, maxX, "*");
	attron(COLOR_PAIR(white));
	mvprintw(10, maxX+2, "Degraded Fiber");
	mvprintw(12, maxX+2, "Non-degraded Fiber");
	mvprintw(14, maxX+2, "Unbound Molecule");
	mvprintw(16, maxX+2, "Bound Molecule");
	mvprintw(0, 0, ("Time: " + to_string(t * timeStep) + " sec").c_str());
	for (int x = 0; x < nodesInRow; x++)
		for (int y = 0; y < nodesInColumn; y++)
			printNode(x,y);
	for (int j = 0; j < totalMolecules; j++)
		printMolecule(j);
	refresh();
//	getch();
}

/*
 *
 */
void printFiberStatus(int i, unsigned int t) {
    cout << "Time: " << t * timeStep << " sec" << endl;
    cout << "Fiber #" << i << " contains molecule(s) ";
    int currentMolecule = firstMolecule[i];
    while (currentMolecule != -1) {
    	if (bound[currentMolecule])
    		cout << to_string(currentMolecule) << ", ";
    	else
    		cout << to_string(currentMolecule) << ", ";
    	currentMolecule = nextMolecule[currentMolecule];
    }
    cout << endl;
}

/*
 * This method will save the data every 10 seconds by pushing the data into a vector until the loop
 * completes itself, then it will stop, and the vector can be viewed to see the gathered data.
 */
/*void saveData(int t) {
	// Include code from lines 866 - 878
	int tenSecond = ceil(10/timeStep);
		for (unsigned int t = 0; t < totalTimeSteps; t++){
	      if ((t % tenSeconds == 0) && verbose) {
		   dataVec.push_back(counterSet());
		   dataVec[t].locationCheck = location[totalMolecules];
		   dataVec[t].boundCheck = bound[totalMolecules];
		   dataVec[t].timeCheck = totalTime;
		   dataVec[t].degradedCheck = degraded[totalFibers];
		   //dataVec[t].totalBindingsCheck =
		  // dataVec[t].totalAttemptsCheck =
		  }
		}
}*/


/*
 * Collects the first undegraded edges in a column and puts them into a vector.
 */
/*void processData(unsigned short j, unsigned short p, ) {
	// ind is a vector containing the vertical planar edge numbers above node j
	
	for(int i = 0; i < xznodes; i++){
	  degreeVector[i] = getNeighbor(i,UP);
	//place(k) is the degradation state of each edge above node j
     degradeCheckVector[i].degCheck = degraded[degreeVector[i]];
	}
	// find the first undegraded vertical edge above node j
	for(int i = 0; i < xznodes; i++){
      if (degradeCheckVector][i].degcheck == false){
	   // Have some kind of array or vector keep track of undegraded edges	
	}
	
	
	//Find first inequal value 
	for(int i = 0;i < placeHolderLimit; i ++){
		    if (degreeVector[i].inequal){
	//If the first degree is zero than it is set to 1
	             if(degreeVector[i] == 0)
					 degreeVector[i] = 1;
			}
	placeHolderArray = degreeVector[i];
	//Otherwise it equals whatever number it can possibly be
	
	//Also store the data in an array or vector
	
	//Also find the first undegraded vertical edge above node j
	}
	
	//By now the successive y and x positions should be saved, and a bunch of data types need to be
	//saved for later use in Matlab.
	 
	 //Records all the successive points
	 //Lets you decide how many runs you want to save to make a movie
	 //Unsure how this will be implemented. Maybe something could be showen to the user before move processing begins.
	
	
	
	// INCLUDE FORTRAN CODE FROM LINES 883-981
}*/

/*
 * Processes data for a movie that can be produced with the data in this function after the fact.
 */
/*void movieProcessing() {
	
	Uses a grid of nodes, starting at the bottom left and moving right.
	Endpoints have the node numbers corresponding to the endpoints of the fiber (edge)
	
	Coordinates are set to zero on the grid
	If the next degree is equal to zero, then a counter known as countintact goes up then
	the system find the undegraded edge numbers and stores them somewhere
	
	Goes through the grid going through the endpoints of the fibers
	Starts with vertical edges
	Then moves to horizontal edges
	then it has a section for vertical edges
	
	
	if theres a horizontal edge it finds the y value at which the horizontal edge occurs
	if there is a vertical edge (or planar) it will find the x value at twhich the vertical edge occurs
	and also to find the bottom endpoint of the vertical edge
	
	Goes through and dots the location and boundedness of the grid, black dots if unbound, and green dots if bound.
	
	Method appears to go through the whole grid, checking on what kind of edges it is dealing with
	as it moves through the grid
	
	// INCLUDE FORTRAN CODE FROM LINES 982 - 1205 (1209 - 1382)
}*/

/*
 * Outputs all the processed data, then closes all the vectors???
 */
void outputData() {
	//Closes a bunch of units??? Nothing else apparently
	//Will print all the data out???
	
	// INCLUDE FORTRAN CODE FROM LINES 1403 - 1422)
	
}

/******************************************************************************
 ** Model running methods 
 **
 ******************************************************************************


/*
 * Sets the model state variables to their initial values, that is,
 * All fibers are non-degraded and have no degrade time set
 * All molecules are unbound, have no timer set, and their locations
 * are uniformly distributed on the ghost fibers
 */
void initializeVariables() {
    for (int i = 0; i < totalFibers; i++) { 										// Loop over all fibers 
        degraded[i] = false; 														// Set the fiber to non-degraded
        degradeTime[i] = totalTimeSteps + 1;										// Set the fiber to never degrade
        firstMolecule[i] = -1;														// Set the fiber to empty of molecules
    }
    for (int j = 0; j < totalMolecules; j++) { 										// Loop over all molecules
        bound[j] = false;															// Set the molecule to unbound
        location[j] = randomInt(lastGhostFiber + 1);								// Set the molecule's location to a random position outside the clot
        nextMolecule[j] = firstMolecule[location[j]];								// Insert this molecule at the beginning of the linked list for the fiber
        firstMolecule[location[j]] = j;
        bindUnbindTime[j] = totalTimeSteps + 1;										// Set the molecule to never bind
        //printFiberStatus(location[j], 0);
    }
}

/*
 * Checks all fibers and degrades those where the current time 't' is 
 * past their degradeTime. Also unbinds all molecules on those fibers
 */
void degradeFibers(unsigned int t) {
    for (int i = lastGhostFiber+1; i < totalFibers; i++) {							// Loop over all true fibers
        if (!degraded[i] && (degradeTime[i] <= t)) {								// If the fiber is not already degraded and its degradeTime has passed
            degraded[i] = true;														// Degrade the fiber
            int currentMolecule = firstMolecule[i];
            while (currentMolecule != -1) {											// Step through the linked list of molecules located on this fiber
                if (bound[currentMolecule]) {										// If the molecule is bound
                    bound[currentMolecule] = false;										// Unbind the molecule
                    bindUnbindTime[currentMolecule] = t;								// Set the molecule's unbind time to the current time
                }
                currentMolecule = nextMolecule[currentMolecule];
            }
            //printFiberStatus(i, t);
        }
    }
}

/*
 * Part of Bind. Finds the time to unbind, and returns colr2, for use in degradTime.
 */
int findUnbindTime(unsigned short j, unsigned int t, double r3) {
    int colr2 = ceil(r3*100);
	double slope = UnbindingTimeDistribution[colr2] - UnbindingTimeDistribution[colr2 - 1];
    double step = (r3*100) - colr2;
	double timeToUnbind = UnbindingTimeDistribution[colr2 -1] + (step*slope);
	bindUnbindTime[j] = t + timeToUnbind/timeStep;
	return colr2;
}

/**
 *
 */
int findBindingTime(unsigned short j, unsigned int t, double r) {
	bindUnbindTime[j] = t - log(r) / (bindingRate * bindingSites) - timeStep/2; //!! Ask Dr. Bannish about this line in the context of line 583 (586 in Fortran)
}

/*
 * Part of Bind. Finds the degradation time.
 * STILL NEEDS INPUT FROM DR BANNISH
 */
void setDegradeTime(unsigned short i, unsigned int t, int colr3, double r4) {
	int r400 = ceil(r4*100)+1;
	if (r400 <= lysesPerBlock[colr3-1]) {
		double slope = lysisTime[r400][colr3-1] - lysisTime[r400-1][colr3-1];
		double step = (r4*100) - (r400-1);
		double timeToDegrade = lysisTime[r400][colr3-1] + step*slope;
		degradeTime[i] = min(degradeTime[i], (int)(t + floor(timeToDegrade/timeStep)));
	}
}
    


/*
 * This method unbinds a molecule, while also finding the binding time.
 */
void unBind(unsigned short j, unsigned int t, double r) {
    bound[j] = false; //Sets bound to false, causing the point j to unbind.
    findBindingTime(j,t,r); //This will find the binding time of the molecule.
}

/*
 * This method will bind two molecules using the findUnbindTime and setDegradeTime methods.
 */
void bind(unsigned short j, unsigned short i, unsigned int t, double r1, double r2) {
    bound[j] = true; //Sets bound to true, which causes the point to bind.
    int colr2 = findUnbindTime(j,t,r1);
	setDegradeTime(i,t,colr2,r2);
}

/*
 * This method moves the molecule depending on moving probability, then decides which neighbor to move to.
 */
void moveMolecule(unsigned short j, unsigned int t, double r) {
    unsigned short neighbor = floor(8*((1-r) / movingProbability));
    //printMoleculeStatus(j, t);
    //printFiberStatus(location[j], t);
    if (firstMolecule[location[j]] == j) {
    	firstMolecule[location[j]] = nextMolecule[j];
    } else {
    	int currentMolecule = firstMolecule[location[j]];
    	while (nextMolecule[currentMolecule] != j) {
    		currentMolecule = nextMolecule[currentMolecule];
    		if (currentMolecule == -1)
    			throw logic_error("Broken Linked List.");
    	}
    	nextMolecule[currentMolecule] = nextMolecule[j];
    }
	location[j] = getNeighbor(location[j], neighbor);
	nextMolecule[j] = firstMolecule[location[j]];
	firstMolecule[location[j]] = j;
	findBindingTime(j, t, urcw1_());
}



/*  
 *
 */
void runModel() {
	int outputInterval = ceil(0.1/timeStep);
    for (unsigned int t = 0; t < totalTimeSteps; t++) {
    	//printGrid(t);
        if ((t % outputInterval == 0) && verbose) {
			//printMoleculeStatus(43, t);
			printGrid(t);
			
		}
        degradeFibers(t);
        for (unsigned int j = 0; j < totalMolecules; j++) {
            if (bound[j] && (bindUnbindTime[j] <= t)) // If the molecule is bound, should it unbind?
                unBind(j, t, urcw1_());
            if (!bound[j]) { // If the molecule is unbound, should it move?
				double r = urcw1_();
				// if (verbose) {
					// cout << "Move roll:" << r << endl;
				// }
                if (r <= 1-movingProbability) { // If the molecule does not move, can it bind?
                    if ((bindUnbindTime[j] <= t) && !degraded[location[j]]) { // It can bind.
                        bind(j, location[j], t, urcw1_(), urcw1_()); // If it can bind, bind it. If not, leave it as is.
                    }
                } else {
					double rn = urcw1_();
					if ((bindUnbindTime[j] <= t) && !degraded[location[j]]) { // If it can move, can it still bind?
						if(rn > (t - bindUnbindTime[j])/totalTimeSteps){
							moveMolecule(j, t, r);
						} else {
							bind(location[j], j, t, urcw1_(), urcw1_()); //If it is smaller, bind it.
						} 
					} else {// It could not bind, so move
						moveMolecule(j, t, r);
					}
                }
            }
        }
//		if ((t % tenSeconds) == 0)
//			saveData(t);
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
	if (verbose) cout << "Running Model..................." << endl;
	if (verbose) initializeCurses();
    runModel();
    if (verbose) getch();
    if (verbose) endwin();
    if (verbose) cout << "Run Complete." << endl;
    return 0;
}
