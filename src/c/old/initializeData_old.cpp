#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "all.h"

/*
* The following three methods read in the microscale data
* from the appropriate files.
* Each method returns "false" if the data file cannot be read,
* or if it does not contain the correct amount of data.
* Each method returns "true" if the data is read in correctly.
*/

bool loadLysisTimeFromFile(const std::string &inLysisTimeFile, double** inLysisTime, const int inMaxLysesPerBlock,  const int inLysisBlocks) {

	using namespace std;

	ifstream reader(inLysisTimeFile);													// Open data file for reading
	if (!reader) {																	// Check if the file can be opened
		cerr << "Failed to open file " << inLysisTimeFile << "." << endl;
		return false;
	}
	else {																		// If the file was successfully opened
		int i = 0; 																		// Make a var to store the current line index
		string line; 																	// Make a var to hold the current line
		while (getline(reader, line)) { 												// Get the next line from the file
			if (i >= inMaxLysesPerBlock) 														// Check that we haven't exceeded the expected number of lines
				cerr << "Extra rows in " << inLysisTimeFile << " discarded." << endl;
			else {
				istringstream stringReader(line); 											// Set up the line for reading
				int j = 0; 																	// Make a var to store the current entry's index
				double a; 																	// Make a var to store the current entry
				while (stringReader >> a) { 												// Get the next entry from the current line
					if (j >= inLysisBlocks) 														// Check that we haven't exceeded the expected number of entries in this line
						cerr << "Extra entries in row " << i
						<< " of " << inLysisTimeFile << " discarded." << endl;
					else
						inLysisTime[i][j++] = a; 													// Store the data in the appropriate position in the array
				}
				if (j < inLysisBlocks) {														// Check if we've reached the end of the line with insufficient entries
					cerr << "Insufficient entries in row " << i
						<< " of " << inLysisTimeFile << "." << endl;
					return false;
				}
			}
			i++;																			// Increment the line counter
		}
		reader.close();																	// Close the file
		if (i < inMaxLysesPerBlock) {														// Check if we've reached the end of the file with insufficient rows
			cerr << "Insufficient rows in " << inLysisTimeFile << ":("
				<< i << ", " << inMaxLysesPerBlock << ")" << endl;
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
bool loadUnbindingTimeFromFile(const std::string &UnbindingTimeFile, double* UnbindingTimeDistribution, const int inUnbindingTimeDist) {

	using namespace std;

	ifstream reader(UnbindingTimeFile);												// Open data file for reading
	if (!reader) {																	// Check if the file can be opened
		cerr << "Failed to open file " << UnbindingTimeFile << "." << endl;				// If not, return an error
		return false;
	}
	else {																		// If the file was successfully opened
		int i = 0;																		// Make a var to store the current entry index
		double a;																		// Make a var to store the current entry
		while (reader >> a) {															// Read in the data one at a time
			if (i >= inUnbindingTimeDist)																	// Check if we've exceeded the expected number of entries
				cerr << "Extra entries in " << UnbindingTimeFile
				<< " discarded." << endl;
			else
				UnbindingTimeDistribution[i++] = a;											// Store the data in the appropriate position in the array
		}
		reader.close();																	// Close the file
		if (i < inUnbindingTimeDist) {
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
bool loadLysesPerBlockFromFile(const std::string totalLysesFile, int* inLysesPerBlock, const int &inLysisBlocks) {

	using namespace std;

	ifstream reader(totalLysesFile);												// Open data file for reading
	if (!reader) {																	// Check if the file can be opened
		cerr << "Failed to open file " << totalLysesFile << "." << endl;				// If not, return an error
		return false;
	}
	else {																		// If the file was successfully opened
		int i = 0;																		// Make a var to store the current entry index
		int a;																			// Make a var to store the current entry
		while (reader >> a) {															// Read in the data one entry at a time
			if (i >= inLysisBlocks)															// Check if we've exceeded the expected number of entries
				cerr << "Extra entries in " << totalLysesFile << " discarded." << endl;
			else
				inLysesPerBlock[i++] = a;														// Store the current entry in the array
		}
		reader.close();																	// Close the file
		if (i < inLysisBlocks) {															// Check if we've reached the end of the file with insufficient entries
			cerr << "Insufficient entries in " << totalLysesFile << "." << endl;
			return false;
		}
	}
	return true;																	// Everything was read in correctly
}

