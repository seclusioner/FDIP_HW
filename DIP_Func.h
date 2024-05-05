#pragma once

#ifndef DIP_H
#define DIP_H
#endif

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>

#include "bmp.h"


#ifndef PI
#define PI 3.14159265358979323846
#endif

#define BIT 8 // Grayscale bits

using namespace std;


/****************************
 Functions Definitions
*****************************/
int** allocate2DArray(int rows, int cols);
float** allocate2DArray_f(int rows, int cols);

void free2DArray(int** arr, int rows);
void free2DArray_f(float** arr, int rows);

void selectionSort(int arr[], int n);
int Median(int* arr, int len);
void swap(int* xp, int* yp);

// Basic Functions
void Convolve2D(int** paddedMatrix, float** mask, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize);

int** ZeroPadding(int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernel_sz);
int** EdgePadding(int image[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernel_sz);

// Algorithm Functions
void Negative(int R[MaxBMPSizeX][MaxBMPSizeY], int G[MaxBMPSizeX][MaxBMPSizeY], int B[MaxBMPSizeX][MaxBMPSizeY],
              int r[MaxBMPSizeX][MaxBMPSizeY], int g[MaxBMPSizeX][MaxBMPSizeY], int b[MaxBMPSizeX][MaxBMPSizeY],
              int W, int H, bool GrayScale);

void Grayscale(int R[MaxBMPSizeX][MaxBMPSizeY], int G[MaxBMPSizeX][MaxBMPSizeY], int B[MaxBMPSizeX][MaxBMPSizeY],
              int r[MaxBMPSizeX][MaxBMPSizeY], int W, int H);

void Histo_Equal(int R[MaxBMPSizeX][MaxBMPSizeY], int r[MaxBMPSizeX][MaxBMPSizeY], int W, int H);

void Binarize(int R[MaxBMPSizeX][MaxBMPSizeY], int r[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int threshold);

void EdgeSharpen(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float alpha);

void ErrorDiffusion(int matrix[MaxBMPSizeX][MaxBMPSizeY], int output[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols);

void Gradient(int image[MaxBMPSizeX][MaxBMPSizeY], float gradient[MaxBMPSizeX][MaxBMPSizeY], float direction[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols);
void Canny(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int low_threshhold, int high_threshhold, int rows, int cols, float percentage);

int BilinearInterpolation(int image[MaxBMPSizeX][MaxBMPSizeY], float rotated_x, float rotated_y);
void Rotate_forward(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int width, int height, int &new_w, int &new_h, float deg);
void Rotate_backward(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int width, int height, int& new_w, int& new_h, float deg);


// Filter
void MedianFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize);
void SmoothFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize);
void  Bino_Filter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize, bool HPF);
void Sobel(int** paddedMatrix, int grad_x[MaxBMPSizeX][MaxBMPSizeY], int grad_y[MaxBMPSizeX][MaxBMPSizeY],
                  int sobel[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float percentage);
void Laplacian(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float percentage);
void GaussianFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize, float sigma);
