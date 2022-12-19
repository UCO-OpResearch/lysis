/* 
 * File:   fileio.h
 * Author: bpaynter
 *
 * Created on December 19, 2015, 11:09 AM
 */

#ifndef FILEIO_H
#define	FILEIO_H
#include <string>
#include <array>

template <typename T, size_t size>
bool readArrayFromFile(std::array<T, size> arr, std::string filename);

template <typename T, size_t rowSize, size_t colSize>
bool read2DArrayFromFile(std::array<std::array<T, rowSize>, colSize> arr, std::string filename);


#endif	/* FILEIO_H */

