#pragma once

#include <iostream>
#include <algorithm>

#ifndef PADDING_H
#define PADDING_H

typedef enum {  // Padding mode def.
    constant = 1,
    edge,
    reflect,
    symmetric,
    maximum,
    minimum,
    mean,
    median
} PaddingMode;

class Padding {
public:
    /////////////////// Padding ///////////////////
    static int** zeroPadding(int** matrix, int numRows, int numCols, int padWidth, int padd_const);
    static int** edgePadding(int** matrix, int numRows, int numCols, int padWidth);
    static int** reflectPadding(int** image, int H, int W, int padWidth);
    static int** symmetricPadding(int** image, int H, int W, int padWidth);
    static int** maximumPadding(int** matrix, int H, int W, int padWidth);
    static int** minimumPadding(int** matrix, int H, int W, int padWidth);
    static int** meanPadding(int** matrix, int H, int W, int padWidth);
    static int** medianPadding(int** matrix, int H, int W, int padWidth);

    /////////////////// Dynamic Allocation ///////////////////
    static int** allocateMatrix(int rows, int cols);
    static void freeMatrix(int** matrix, int rows);
};

int** Pad(int matrix[MaxBMPSizeX][MaxBMPSizeY], int H, int W, int padWidth, PaddingMode mode, int padd_const);

#endif
