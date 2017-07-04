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
 ** Model Time Parameters
 ** i.e., those parameters selected to define 
 **       the running time and granularity of the model
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
unsigned int unbindingTime[totalMolecules]; // bind/t_leave

// Degradation state of each edge.
bool degraded[totalFibers]; // degrade

// Timestep when the fiber will degrade
int degradeTime[totalFibers]; // t_degrade


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
    else if (rowPosition % 2 == 0)													// Else, 
    	return OUT;
    else
    	return RIGHT; // Finds the direction of the fiber, and returns it
}






/*
 * Gets the neighbor surrounding the selected fiber, given the
 * index of a fiber.
 */
int getNeighbor(unsigned short fiber, unsigned short i){
	unsigned short x = nodeX(fiber); // Sets the x-coordinate to the fiber's x
    unsigned short y = nodeY(fiber); // Sets the y-coordinate to the fiber's x
    unsigned short dir = fiberDirection(fiber); // Sets the direction to the fiber's x
	switch (dir) {
        case UP:
			switch(i){
				case 0:
				return (x != 0) ? fiberIndex(x,y,LEFT) : fiberIndex(x,y,RIGHT);
				case 1: 
				return fiberIndex(x, y, OUT); 
				case 2:
				return fiberIndex(x, y, IN);
				case 3: 
				return  (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT);
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
				return (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // 
				case 1: 
				return fiberIndex(x, y, OUT); // 
				case 2:
				return fiberIndex(x, y, IN); // 
				case 3: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // 
				case 4: 
				return (y != 0) ? fiberIndex(x+1, y, DOWN) : fiberIndex(x+1, y, UP); // 
				case 5: 
				return fiberIndex(x, y, OUT); //  
				case 6:  
				return fiberIndex(x, y, IN); // 
				case 7: 
				return (y != nodesInColumn-1) ? fiberIndex(x+1, y, UP) : fiberIndex(x+1, y, DOWN);
			}
			 case OUT:
			switch(i){
				case 0: 
				return (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT); // z
				case 1: 
				return (x != 0) ? fiberIndex(x, y, LEFT) : fiberIndex(x, y, RIGHT); // z+1
				case 2:
				return (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT); // z
				case 3: 
				return (x != nodesInRow-1) ? fiberIndex(x, y, RIGHT) : fiberIndex(x, y, LEFT); // z+1
				case 4: 
				return  (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // z
				case 5: 
				return  (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // z+1
				case 6: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // z
				case 7: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // z+1
			}
			default: // If something goes wrongs with the method
            throw invalid_argument("Something went wrong while finding neighbors.");
		}
}

/*
 * Gets the neighbors surrounding the selected fiber, given the
 * index of a fiber.
 */
void getNeighbors(unsigned short fiber, unsigned short neighbors[]) {
  for (int i = 0; i < 8; i++){
		neighbors[i] = getNeighbor(fiber,i);
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
 * All fibers are non-degraded and have no degrade time set
 * All molecules are unbound, have no timer set, and their locations
 * are uniformly distributed on the ghost fibers
 */
void initializeVariables() {
    for (int i = 0; i < totalFibers; i++) { //Goes through all possible fibers 
        degraded[i] = false; 
        degradeTime[i] = totalTimeSteps + 1;
    }
    for (int i = 0; i < totalMolecules; i++) { //Goes through all possible fibers and unbinds them, while counting the unbind time.
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
 * Part of Bind. Finds the time to unbind, and returns colr2, for use in degradTime.
 */
int findUnbindTime(unsigned short j, unsigned int t, double r3) {
    int colr2 = ceil(r3*100);
	double slope = UnbindingTimeDistribution[colr2] - UnbindingTimeDistribution[colr2 - 1];
    double step = (r3*100) - colr2;
	double timeToUnbind = UnbindingTimeDistribution[colr2 -1] + (step*slope);
	unbindingTime[j] = t + timeToUnbind/timeStep;
	return colr2;
}

/**
 *
 */
int findBindingTime(unsigned short j, unsigned int t, double r) {
	unbindingTime[j] = t - log(r) / (bindingRate * bindingSites) - timeStep/2; //!! Ask Dr. Bannish about this line in the context of line 583 (586 in Fortran)
}

/*
 * Part of Bind. Finds the degradation time.
 * STILL NEEDS INPUT FROM DR BANNISH
 */
void setDegradeTime(unsigned short i, unsigned int t, int colr3, double r4) {
	int r400 = ceil(r4*100)+1;
	if(r400 >= lysesPerBlock[colr3-1]) {
		double slope = lysisTime[r400][colr3-1] - lysisTime[r400-1][colr3-1];
		double step = (r4*100)-(r400-1);
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
//	if (verbose) {
//		cout << "Finding neighbor " << neighbor << " for location " << location[j] << endl;
//	}
	location[j] = getNeighbor(location[j], neighbor);
	findBindingTime(j, t, urcw1_());
}










/*
 * This method will save the data every 10 seconds by pushing the data into a vector until the loop
 * completes itself, then it will stop, and the vector can be viewed to see the gathered data.
 */
//void saveData(int t) {
//	// Include code from lines 866 - 878
//	int tenSecond = ceil(10/timeStep);
//		for (unsigned int t = 0; t < totalTimeSteps; t++){
//	      if ((t % tenSeconds == 0) && verbose) {
//		   dataVec.push_back(counterSet());
//		   dataVec[t].locationCheck = location[totalMolecules];
//		   dataVec[t].boundCheck = bound[totalMolecules];
//		   dataVec[t].timeCheck = totalTime;
//		   dataVec[t].degradedCheck = degraded[totalFibers];
//		   //dataVec[t].totalBindingsCheck =
//		  // dataVec[t].totalAttemptsCheck =
//		  }
//		}
//}

/*  
 *
 */
void runModel() {
	int tenSeconds = ceil(10/timeStep);
    for (unsigned int t = 0; t < totalTimeSteps; t++) {
        if ((t % tenSeconds == 0) && verbose) {
			int testMolecule = 4364;
            cout << "Time: " << t * timeStep << " sec" << endl;
			cout << "Location of molecule #" << testMolecule << " is fiber #" << location[testMolecule] << endl;
			cout << "Molecule #" << testMolecule << " is " << (bound[testMolecule] ? "bound" : "unbound") << endl;
			cout << "Fiber #" << location[testMolecule] << " is " << (degraded[location[testMolecule]] ? "degraded" : "not degraded") << endl;
		}
        degradeFibers(t);
        for (unsigned int j = 0; j < totalMolecules; j++) {
            if (bound[j] && (unbindingTime[j] == t)) // If the molecule is bound, should it unbind?
                unBind(j, t, urcw1_());
            if (!bound[j]) { // If the molecule is unbound, should it move?
				double r = urcw1_();
				// if (verbose) {
					// cout << "Move roll:" << r << endl;
				// }
                if (r <= 1-movingProbability) { // If the molecule does not move, can it bind?
                    if ((unbindingTime[j] <= t) && !degraded[location[j]]) { // It can bind.
                        bind(location[j], j, t, urcw1_(), urcw1_()); // If it can bind, bind it. If not, leave it as is.
                    }
                } else {
					double rn = urcw1_();
					if ((unbindingTime[j] <= t) && !degraded[location[j]]) { // If it can move, can it still bind?
						if(rn > (t - unbindingTime[j])/totalTimeSteps){
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

/*
 * Collects the first undegraded edges in a column and puts them into a vector.
 */
//void processData(unsigned short j, unsigned short p, ) {
//	// ind is a vector containing the vertical planar edge numbers above node j
//	
//	for(int i = 0; i < xznodes; i++){
//	  degreeVector[i] = getNeighbor(i,UP);
//	//place(k) is the degradation state of each edge above node j
//      degradeCheckVector[i].degCheck = degraded[degreeVector[i]];
//	}
//	// find the first undegraded vertical edge above node j
//	for(int i = 0; i < xznodes; i++){
//       if (degradeCheckVector][i].degcheck == false){
//	   // Have some kind of array or vector keep track of undegraded edges	
//	}
//	
//	
//	//Find first inequal value 
//	for(int i = 0;i < placeHolderLimit; i ++){
//		    if (degreeVector[i].inequal){
//	//If the first degree is zero than it is set to 1
//	             if(degreeVector[i] == 0)
//					 degreeVector[i] = 1;
//			}
//	placeHolderArray = degreeVector[i];
//	//Otherwise it equals whatever number it can possibly be
//	
//	//Also store the data in an array or vector
//	
//	//Also find the first undegraded vertical edge above node j
//	}
//	
//	//By now the successive y and x positions should be saved, and a bunch of data types need to be
//	//saved for later use in Matlab.
//	 
//	 //Records all the successive points
//	 //Lets you decide how many runs you want to save to make a movie
//	 //Unsure how this will be implemented. Maybe something could be showen to the user before move processing begins.
//	
//	
//	
//	// INCLUDE FORTRAN CODE FROM LINES 883-981
//}

/*
 * Processes data for a movie that can be produced with the data in this function after the fact.
 */
void movieProcessing() {
	
	//Uses a grid of nodes, starting at the bottom left and moving right.
	//Endpoints have the node numbers corresponding to the endpoints of the fiber (edge)
	
	//Coordinates are set to zero on the grid
	// If the next degree is equal to zero, then a counter known as countintact goes up then
	// the system find the undegraded edge numbers and stores them somewhere
	
	//Goes through the grid going through the endpoints of the fibers
	//Starts with vertical edges
	//Then moves to horizontal edges
	//then it has a section for vertical edges
	
	
	//if theres a horizontal edge it finds the y value at which the horizontal edge occurs
	//if there is a vertical edge (or planar) it will find the x value at twhich the vertical edge occurs
	//and also to find the bottom endpoint of the vertical edge
	
	//Goes through and dots the location and boundedness of the grid, black dots if unbound, and green dots if bound.
	
	//Method appears to go through the whole grid, checking on what kind of edges it is dealing with
	//as it moves through the grid
	
	// INCLUDE FORTRAN CODE FROM LINES 982 - 1205 (1209 - 1382)
}

/*
 * Outputs all the processed data, then closes all the vectors???
 */
void outputData() {
	//Closes a bunch of units??? Nothing else apparently
	//Will print all the data out???
	
	// INCLUDE FORTRAN CODE FROM LINES 1403 - 1422)
	
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
	 if (verbose) cout << "Running Model...................";
    runModel();
    if (verbose) cout << "Done." << endl;
    return 0;
}
