#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <array>

using namespace std;

template <size_t x, size_t y>
bool read2DArrayFromFile(array<array<double, x>, y> arr, std::string filename) {
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
bool readArrayFromFile(array<double,size> arr, string filename) {
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
