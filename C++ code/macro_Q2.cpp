/* 
 * File:   macro_Q2.cpp
 * Author: Dr. Brittany Bannish
 * Converted to C++ by Dr. Bradley Paynter
 *
 * Created on December 18, 2015, 4:28 PM
 * 
 * Runs the macroscale model in a clot 
 * with 72.7 nm diameter fibers and 
 * pore size. 1.0135 uM. 
 * FB conc. = 8.8 uM
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <array>
#include <sstream>




//#include "kiss.o"

using namespace std;

const bool verbose = true;

// The number of lattice nodes in each (horizontal) row
const int nodesInRow = 93;

// The number of lattice nodes in each (vertical) column
const int nodesInColumn = 121;

/*
 * 1st node in vertical direction containing fibers. 
 * So if firstFiberRow=10, then rows 1-9 have no fibers, 
 * there's one more row of fiber-free planar vertical edges, 
 * and then the row starting with the firstFiberRow-th (e.g. 10th) vertical node 
 * is a full row of fibers
 */
const int firstFiberRow = 29;

// The number of independent trials to be run
const int numberOfTrials = 10;

// The total number of fibers in the model
const int numberOfFibers = (2*nodesInRow - 1)*nodesInColumn 
                                       + nodesInRow*(nodesInColumn - 1);

// The total number of tPA molecules: 
//      43074 is Colin's [tPA]=0.6 nM 
//      86148 is Colin's [tPA]=1.2 nM
const int totalMolecules = 43074;

// Total running time for model in seconds
const int totalTime = 10*60;

// The last edge number without fibrin
const int lastGhostFiber = (3*nodesInRow - 1)*(firstFiberRow - 1);

/*******************************************************************************
 * The following variables will hold the microscale data 
 * which will be read in from the named data files
 ******************************************************************************/
const string CDFtPAFile = "tPAleavePLG135_Q2.dat";
array<double, 101> CDFtPA;

const string tsec1File = "tsectPAPLG135_Q2.dat";
array<double, 101> tsec1;

/*
 * lysismat_PLG135_Q2.dat is a matrix with column corresponding to 
 * bin number (1-100) and with entries equal to the lysis times obtained 
 * in that bin. 
 * An entry of 6000 means lysis didn't happen.
 * i.e. the first column, lysismat(:,1), gives the lysis times 
 * for the first 100 tPA leaving times.
 */
const string lysisTimeFile = "lysismat_PLG135_Q2.dat";
array<array<double, 100>, 100> lysisTime;

/*
 * lenlysisvect_PLG135_Q2.dat saves the first row entry in each column of 
 * lysismat_PLG135_Q2.dat where lysis did not occur, 
 * i.e. the first entry there's a 6000
 */
const string lengthOfLysisFile = "lenlysisvect_PLG135_Q2.dat";
array<int,100> lengthOfLysis;




template <size_t x, size_t y>
bool read2DArrayFromFile(array<array<double, y>, x>& arr, std::string filename) {
    // Open data file for reading
    ifstream reader(filename);
    // Check if the file can be opened
    if (! reader) {
        // If not, return an error
        cerr << "Failed to open file " << filename << "." << endl;
        return false;
    } else {
        // If the file was successfully opened
        int i = 0;
        string line;
        while (!getline(reader, line)) {
            if (i >= arr.size())
                cerr << "Extra rows in " << filename << " discarded." << endl;
            else {
                istringstream stringReader(line);
                int j = 0;
                double a;
                // Read in the data one double at a time
                while (stringReader >> a) {
                    if (j >= arr.size())
                        cerr << "Extra entries in row " << i 
                                << " of " << filename << " discarded." << endl;
                    else
                        arr[i][j++] = a;
                }
                if (j < arr.size()) {
                    cerr << "Insufficient entries in row " << i 
                            << " of " << filename << "." << endl;
                    return false;
                }
            }
            i++;
        }
        // Close the file
        reader.close();
        if (i < arr.size()) { 
            cerr << "Insufficient rows in " << filename << "." << endl;
            return false;
        }
    }
    return true;
}

template <size_t size>
bool readArrayFromFile(array<double,size>& arr, string filename) {
    int i = 0;
    // Open data file for reading
    ifstream reader(filename);
    // Check if the file can be opened
    if (! reader) {
        // If not, return an error
        cerr << "Failed to open file " << filename << "." << endl;
        return false;
    } else {
        // If the file was successfully opened
        double a;
        // Read in the data one double at a time
        while (reader >> a) {
            if (i >= arr.size()) 
                cerr << "Extra entries in " << filename 
                        << " discarded." << endl;
            else
                arr[i++] = a;
        }
        // Close the file
        reader.close();
        if (i < arr.size()) {
            cerr << "Insufficient entries in " << filename << "." << endl;
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
    success = readArrayFromFile(CDFtPA, CDFtPAFile);
    if (!success)
        return false;
    success = readArrayFromFile(tsec1, tsec1File);
    if (!success)
        return false;
    return true;
}

/*
 * 
 */
int main(int argc, char** argv) {
    cout << readData() << endl;
    cout << CDFtPA[0] << "," << tsec1[0] << endl;
    cout << CDFtPA[28] << "," << tsec1[74] << endl;
    cout << CDFtPA[100] << "," << tsec1[100] << endl;
    cout << "I Worked!!" << endl;
    cout << "The last ghost fiber is in position " << lastGhostFiber << endl;
    return 0;
}

