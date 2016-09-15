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
 ** Output Data Storage Variables
 **
 ******************************************************************************/
 
//
int successfulBinds;

//
int unsuccessfulBinds;

//For the set of vectors for saveData.
struct counterSet{
	short locationCheck;
	bool boundCheck;
	int timeCheck;
	bool degradedCheck;
	int totalBindingsCheck;
	int totalaAttemptsCheck;
}
// A vector being used for processData.
struct dvSet{
	int x;
	int y;
	bool deg;
}

struct checkSet{
	bool degCheck;
}
//
vector<counterSet> dataVec;
//
vector<dvSet> degreeVector;
//
vector<checkSet> degradeCheckVector;

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
    success = readLysisTimeFromFile(); // Will read the file to check if it is a success
    if (!success)  // If it isn't a success
        return false; // Then readData will not return true
    success = readUnbindingTimeFromFile();
    if (!success)
        return false; // readData will not return true
    success = readLysesPerBlockFromFile();
    if (!success)
        return false; // readData will not return true
    return true; // If all are successful, then readData returns true
}

/*
 * Set up the random number generator 'kiss32' with the values from 'state'
 */
void initializeRandomGenerator() {
    set_kiss32_(state); // Each command initialize the random generators
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
    if ((nodeX < 0) || (nodeY < 0)) // Checks to see if the nodes are non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (nodeX >= nodesInRow) // Checks to see if the index is in bounds, row bound
        throw invalid_argument("Index out of bounds.");
    else if (nodeY >= nodesInColumn) // Checks to see if the index is in bounds, column bound
        throw invalid_argument("Index out of bounds.");
    else
        switch (direction) { // If the index is in bound, will get the index of the fiber depending on the direction given
            case DOWN:
                return (nodeY == 0) ? -1 : fiberIndex(nodeX, nodeY-1, UP); // Will check if node Y is zero, if it is, then -1 will be returned. Else, the function will be called-
            case LEFT:                                                     // -again but with parameters.
                return (nodeX == 0) ? -1 : fiberIndex(nodeX-1, nodeY, RIGHT);// Will check if node X is zero, if it is, then -1 will be returned. Else, the function will be called-
            case IN:
                return fiberIndex(nodeX, nodeY, OUT); // If the direction is IN, it will call the function again but with the direction being OUT
            case UP:
                return (nodeY == nodesInColumn - 1) ? -1 : //If node Y is equal to the nodesInColumn -1, it will return -1, else it will return the equation after the question mark
                    nodeY * fullRow + xzRow + nodeX;
            case RIGHT:
                return (nodeX == nodesInRow - 1) ? -1 : //If node X is equal to the nodesInRow -1, it will return -1, else it will return the equation after the question mark
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
    if (fiberIndex < 0)  // Checks if the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) // Checks if the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;    
    return (rowPosition >= xzRow) ? rowPosition - xzRow : // Finds the x-coordinate, and returns it
        (unsigned short)(rowPosition / 2);
}

/*
 * Returns the y-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeY(unsigned short fiberIndex) {
    if (fiberIndex < 0)  // Checks if the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) // Checks if the index is in bounds
        throw invalid_argument("Index out of bounds.");
    return (unsigned short)(fiberIndex / fullRow); // Finds the y-coordinate, and returns it
}

/*
 * Returns the direction of a fiber given the 
 * (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberDirection(unsigned short fiberIndex) {
    if (fiberIndex < 0) // Checks if the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers) // Checks if the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;
    return (rowPosition >= xzRow) ? UP : (rowPosition % 2 == 0 ? OUT : RIGHT); // Finds the direction of the fiber, and returns it
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
				return (y != 0) ? fiberIndex(x, y, DOWN) : fiberIndex(x, y, UP); // Goes through each touching element of the fiber selected
				case 1: 
				return fiberIndex(x, y, OUT); // Calls the fiberIndex function to get this neighbor
				case 2:
				return fiberIndex(x, y, IN); // Calls the fiberIndex function to get this neighbor
				case 3: 
				return (y != nodesInColumn-1) ? fiberIndex(x, y, UP) : fiberIndex(x, y, DOWN); // If y doesn't equal nodesInColumn-1, then it will call fiberIndex to UP. Else, DOWN
				case 4: 
				return (y != 0) ? fiberIndex(x+1, y, DOWN) : fiberIndex(x+1, y, UP); // If y doesn't equal nodesInColumn-1, then it will call fiberIndex to DOWN. Else, UP
				case 5: 
				return fiberIndex(x, y, OUT); // Calls the fiberIndex function to get this neighbor 
				case 6:  
				return fiberIndex(x, y, IN); // Calls the fiberIndex function to get this neighbor
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
 * All fibers are undegraded and have no degrade time set
 * All molecules are unbound, have no timer set, and their locations
 * are uniformly distributed on the ghost fibers
 */
void initializeVariables() {
    for (int i = 0; i < totalFibers; i++) { //Goes through all possible fibers and sets them all to false. Also adds to degreade time.
        degraded[i] = false; 
        degradeTime[i] = totalTimeSteps + 1;
    }
    for (int i = 0; i < totalMolecules; i++) {//Goes through all possible fibers and unbinds them, while counting the unbind time.
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
void saveData(int t) {
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
}

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
		if ((t % tenSeconds) == 0)
			saveData(t);
    }
}

/*
 * Collects the first undegraded edges in a column and puts them into a vector.
 */
void processData(unsigned short j, unsigned short p, ) {
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
	
	
	 
	
	
	
	// INCLUDE FORTRAN CODE FROM LINES 883-981
}

/*
 * Processes data for a movie that can be produced with the data in this function after the fact.
 */
void movieProcessing() {
	
	//Uses a grid of nodes, starting at the bottom left and moving right.
	//Endpoints have the node numbers corresponding to the endpoints of the fiber (edge)
	
	//Method appears to go through the whole grid, checking on what kind of edges it is dealing with
	//as it moves through the grid
	
	// INCLUDE FORTRAN CODE FROM LINES 981 - 1205 (1209 - 1382)
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
