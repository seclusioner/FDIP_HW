/****************************************************
* Title: Digital Image Processing Basic Functions
* Author: D.S. (Jing Lu)
* Create: 2024/04/08
* Update: 2024/04/21
*****************************************************/

#include "DIP_Func.h"

int L = pow(2, BIT) - 1;

float gradient_magnitude[MaxBMPSizeX][MaxBMPSizeY] = { 0 };
float gradient_direction[MaxBMPSizeX][MaxBMPSizeY] = { 0 };

/**************************************************
Basic Functions
**************************************************/

/// <summary>
/// Allocates a 2D array for int type and initializes all elements to zero.
/// </summary>
/// <param name="rows">Number of rows in the array</param>
/// <param name="cols">Number of columns in the array</param>
/// <returns>Pointer to the allocated 2D array (int) </returns>
int** allocate2DArray(int rows, int cols) {
    int** array = new int* [rows];
    if(array == NULL){
        cout << "Allocation failed!\n";
        exit(1);
    }
    for (int i = 0; i < rows; ++i) {
        array[i] = new int[cols];
    }

    // Return to zero
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = 0;
        }
    }
    return array;
}

/// <summary>
/// Frees the memory allocated for a 2D array for int type.
/// </summary>
/// <param name="arr">Pointer to the 2D array</param>
/// <param name="rows">Number of rows in the array</param>
void free2DArray(int** arr, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] arr[i];
    }
    delete[] arr;
}

/* -------- Median Filter -------- */

/// <summary>
/// Swaps the values of two integers.
/// </summary>
/// <param name="xp">Pointer to the first integer</param>
/// <param name="yp">Pointer to the second integer</param>
void swap(int* xp, int* yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/// <summary>
/// Sorts an array in ascending order using the selection sort algorithm.
/// </summary>
/// <param name="arr">The array to be sorted</param>
/// <param name="n">The size of the array</param>
void selectionSort(int arr[], int n) {
    int i, j, min_idx;

    // One by one move boundary of unsorted subarray
    for (i = 0; i < n - 1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i + 1; j < n; j++)
            if (arr[j] < arr[min_idx])
                min_idx = j;

        // Swap the found minimum element with the first element
        if (min_idx != i)
            swap(&arr[min_idx], &arr[i]);
    }
}

/// <summary>
/// To calculate the median number for 1D array
/// </summary>
/// <param name="arr">Input 1D array to calculate median number</param>
/// <param name="len">The length of input array</param>
/// <returns></returns>
int Median(int* arr, int len) {
    sort(arr, arr + len);

    if (len % 2 != 0)
        return arr[len / 2];
    else
        return (arr[len / 2] + arr[len / 2 - 1]) / 2;
}

////////////// Binomial Calculation //////////////
float** allocate2DArray_f(int rows, int cols) {
    float** matrix = new float* [rows];
    if (matrix == NULL) {
        cout << "Allocation failed!\n";
        exit(1);
    }
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new float[cols];
    }

    // return to zero
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0;
        }
    }

    return matrix;
}

void free2DArray_f(float** arr, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] arr[i];
    }
    delete[] arr;
}

float Max(int** matrix, int rows, int cols) {
    int max = -999;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] > max)
                max = matrix[i][j];
        }
    }

    return max;
}

float Min(int** matrix, int rows, int cols) {
    int min = 999;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] < min)
                min = matrix[i][j];
        }
    }

    return min;
}

float Sum1dArray(float* arr, int sz) {
    float sum = 0.0;
    for (int i = 0; i < sz; i++) {
        sum += arr[i];
    }

    return sum;
}

int fac(int n) {
    if (n < 0) {
        cout << "n must >= 0." << "\n";
        exit(1);
    }

    if (n == 0 || n == 1)
        return 1;

    return n * fac(n - 1);
}

/// <summary>
/// Computes the coefficients of the binomial expansion for a given size.
/// </summary>
/// <param name="Size">Size of the binomial expansion</param>
/// <returns>Pointer to the array containing the coefficients</returns>
float* Bino_Coef(int Size) {
    float* coefficients = (float*)calloc(Size, sizeof(float));
    if(coefficients == NULL){
        cout << "Memory allocation failed.\n";
        exit(1);
    }
    for (int k = 0; k < Size; k++) {
        int coefficient = fac(Size - 1) / (fac(k) * fac(Size - 1 - k));
        coefficients[k] = coefficient;
    }

    // Normalize
    float sum = Sum1dArray(coefficients, Size);
    for (int i = 0; i < Size; i++) {
        coefficients[i] /= sum;
    }

    return coefficients;
}

/// <summary>
/// Create the binomial mask for binomial filter
/// </summary>
/// <param name="arr">input array (float)</param>
/// <param name="n">Size of the mask</param>
/// <returns>n x n binomial mask</returns>
float** Bino_mask(float* arr, int n) {
    float** mask = allocate2DArray_f(n, n);
    if(mask == NULL){
        cout << "Memory allocation failed.\n";
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mask[i][j] = arr[i] * arr[j];
        }
    }

    return mask;
}

/**************************************************
Padding
**************************************************/

/// <summary>
/// Padding with zeros (int type, single channel)
/// </summary>
/// <param name="matrix">Input matr</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernel_sz">Size of a mask</param>
/// <returns>Padded image</returns>
int** ZeroPadding(int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernel_sz) {
    int padding_sz = kernel_sz / 2;
    int** paddedMatrix = allocate2DArray(rows + padding_sz * 2, cols + padding_sz * 2);

    for (int i = padding_sz; i < rows + padding_sz; i++) {
        for (int j = padding_sz; j < cols + padding_sz; j++) {
            int x = i - (padding_sz);
            int y = j - (padding_sz);
            paddedMatrix[i][j] = matrix[x][y];
        }
    }

    return paddedMatrix;
}

/// <summary>
/// Padding with copy image's edge (int type, single channel)
/// </summary>
/// <param name="image">Input image</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernel_sz">Size of a mask</param>
/// <returns>Padded image</returns>
int** EdgePadding(int image[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernel_sz){
    int padding_size = kernel_sz / 2;
    int** padded_image = allocate2DArray(rows + padding_size * 2, cols + padding_size * 2);

    ////////////////////////// Fill //////////////////////////
    // Fill the center region with the original image
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            padded_image[i + padding_size][j + padding_size] = image[i][j];
        }
    }

    // Fill the top rows with the edge values
    for (int i = 0; i < padding_size; ++i) {
        for (int j = 0; j < cols + 2 * padding_size; ++j) {
            padded_image[i][j] = padded_image[padding_size][j];
        }
    }

    // Fill the bottom rows with the edge values
    for (int i = rows + padding_size; i < rows + 2 * padding_size; ++i) {
        for (int j = 0; j < cols + 2 * padding_size; ++j) {
            padded_image[i][j] = padded_image[rows + padding_size - 1][j];
        }
    }

    // Fill the left columns with the edge values
    for (int i = 0; i < rows + 2 * padding_size; ++i) {
        for (int j = 0; j < padding_size; ++j) {
            padded_image[i][j] = padded_image[i][padding_size];
        }
    }

    // Fill the right columns with the edge values
    for (int i = 0; i < rows + 2 * padding_size; ++i) {
        for (int j = cols + padding_size; j < cols + 2 * padding_size; ++j) {
            padded_image[i][j] = padded_image[i][cols + padding_size - 1];
        }
    }

    return padded_image;
}

/// <summary>
/// Padding with zeros (float type, single channel)
/// </summary>
/// <param name="matrix">Input matr</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernel_sz">Size of a mask</param>
/// <returns>Padded image</returns>
float** ZeroPadding_f(int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernel_sz) {
    int padding_sz = kernel_sz / 2;
    float** paddedMatrix = allocate2DArray_f(rows + padding_sz * 2, cols + padding_sz * 2);
    if(paddedMatrix == NULL){
        cout << "Allocation Failed!\n";
        exit(1);
    }
    for (int i = padding_sz; i < rows + padding_sz; i++) {
        for (int j = padding_sz; j < cols + padding_sz; j++) {
            int x = i - (padding_sz);
            int y = j - (padding_sz);
            paddedMatrix[i][j] = float(matrix[x][y]);
        }
    }

    return paddedMatrix;
}

/**************************************************
DIP Algorithms
**************************************************/

/// <summary>
/// Basic operation to apply a mask on the image.
/// </summary>
/// <param name="paddedMatrix">Padded image (input)</param>
/// <param name="mask">The kernel matrix(filter)</param>
/// <param name="matrix">The matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernelSize">Mask(Kernel)'s size</param>
void Convolve2D(int** paddedMatrix, float** mask, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize) {
    int paddingSize = kernelSize / 2;

    ///////////// 180 deg flip /////////////
    float** tmp = allocate2DArray_f(kernelSize, kernelSize);

    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            tmp[i][j] = mask[i][j];
        }
    }

    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            mask[i][j] = tmp[kernelSize - 1 - i][j];
        }
    }

    float sum = 0.0;
    for (int i = paddingSize; i < rows + paddingSize; i++) {
        for (int j = paddingSize; j < cols + paddingSize; j++) {
            sum = 0.0;
            for (int di = -(kernelSize / 2); di <= (kernelSize / 2); di++) {
                for (int dj = -(kernelSize / 2); dj <= (kernelSize / 2); dj++) {
                    int r = i + di;
                    int c = j + dj;

                    // Ensure boundary conditions
                    if (r >= 0 && r < rows + 2 * paddingSize && c >= 0 && c < cols + 2 * paddingSize) {
                        sum += float(paddedMatrix[r][c]) * mask[di + paddingSize][dj + paddingSize];
                    }
                }
            }
            matrix[i - paddingSize][j - paddingSize] = int(round(sum));
        }
    }

    free2DArray_f(tmp, kernelSize);
}

/// <summary>
/// Negative algorithm for the image. (3 channels)
/// if the input image is single channel, then you just have to input the same matrices
/// </summary>
/// <param name="R">Input  image's R channel</param>
/// <param name="G">Input  image's G channel</param>
/// <param name="B">Input  image's B channel</param>
/// <param name="r">Output  image's R channel</param>
/// <param name="g">Output  image's G channel</param>
/// <param name="b">Output  image's B channel</param>
/// <param name="W">The width of image</param>
/// <param name="H">The height of image</param>
/// <param name="GrayScale"></param>
void Negative(int R[MaxBMPSizeX][MaxBMPSizeY], int G[MaxBMPSizeX][MaxBMPSizeY], int B[MaxBMPSizeX][MaxBMPSizeY],
              int r[MaxBMPSizeX][MaxBMPSizeY], int g[MaxBMPSizeX][MaxBMPSizeY], int b[MaxBMPSizeX][MaxBMPSizeY],
              int W, int H, bool GrayScale){

    int x, y;

    // Example code (negative)
	for (y = H-1; y > -1; y--){ // typical image scanning order
	   for (x = 0; x < W; x++){
         if (GrayScale){
            r[x][y] = L - R[x][y];
         }else{
            r[x][y] = L - R[x][y];
            g[x][y] = L - G[x][y];
            b[x][y] = L - B[x][y];
         }
       } // r[0][0] is located at the bottom left corner.
    } // x: horizontal coordinate, y: vertical coordinate
}

/// <summary>
/// Change the RGB image to Grayscale, like BGR2GRAY
/// </summary>
/// <param name="R">Input  image's R channel</param>
/// <param name="G">Input  image's G channel</param>
/// <param name="B">Input  image's B channel</param>
/// <param name="r">Save the output grayscale result</param>
/// <param name="W">The width of input image</param>
/// <param name="H">The height of input image</param>
void Grayscale(int R[MaxBMPSizeX][MaxBMPSizeY], int G[MaxBMPSizeX][MaxBMPSizeY], int B[MaxBMPSizeX][MaxBMPSizeY],
              int r[MaxBMPSizeX][MaxBMPSizeY], int W, int H){
    int x, y;
	for (y = H-1; y > -1; y--){
	   for (x = 0; x < W; x++){
            //Luminace = 0.3086 * Red + 0.6094 * Green + 0.0820 * Blue
            //Luminace = 0.299 * Red + 0.587 * Green + 0.114 * Blue
            r[x][y] = int(round(0.299 * R[x][y] + 0.587 * G[x][y] + 0.114 * B[x][y]));
       }
	}
}

/// <summary>
/// Histogram Equalization for enhance the image's constract (single channel),
/// if you input the rgb image, then you have to do the grayscale first.
/// </summary>
/// <param name="R">Input image (grayscale)</param>
/// <param name="r">Save output result</param>
/// <param name="W">The width of input image</param>
/// <param name="H">The height of input image</param>
void Histo_Equal(int R[MaxBMPSizeX][MaxBMPSizeY], int r[MaxBMPSizeX][MaxBMPSizeY], int W, int H) {
    if(W <= 0 || H <= 0){
        cout << "Size is invalid!\n";
        exit(1);
    }

    // HE
    int L = pow(2, BIT);

    float* PMF = (float*)calloc(L+1, sizeof(float));
    float* CDF = (float*)calloc(L+1, sizeof(float));

    if(PMF==NULL || CDF == NULL){
        cout << "Allocation failed!\n";
        exit(1);
    }

    float sum = 0.0;
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            PMF[R[i][j]] += (1.0 / (W * H));
        }
    }

    for (int i = 0; i < L+1; i++) {
        sum += PMF[i];
        CDF[i] = sum;
    }

    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            if(R[i][j] > L)
                R[i][j] = L;
            else if(R[i][j] < 0)
                R[i][j] = 0;

            r[i][j] = int(L * (CDF[R[i][j]] - CDF[0]) / (1.0 - CDF[0]));
        }
    }

    free(PMF);
    free(CDF);
}

/// <summary>
/// Binarize the image according to the input threshold. (single channel)
/// </summary>
/// <param name="R">Input image</param>
/// <param name="r">Output image</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="threshold">Input threshold</param>
void Binarize(int R[MaxBMPSizeX][MaxBMPSizeY], int r[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int threshold){
    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            if(R[i][j] > threshold)
                r[i][j] = L;
            else
                r[i][j] = 0;
        }
    }

    return;
}

/// <summary>
/// To enhance the intensity of obejcets's edge (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded matrix</param>
/// <param name="matrix">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
void EdgeSharpen(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float alpha) {
    // Use Ex1 to be a formula
    int kernelSize = 3;
    int paddingSize = 1;

    float mask[3][3] = {
        {-alpha,  alpha-1, -alpha},
        {alpha-1, alpha+5, alpha-1},
        {-alpha,  alpha-1, -alpha}
    };

    for(int i=0;i<kernelSize;i++){
        for(int j=0;j<kernelSize;j++){
            mask[i][j] /= (alpha + 1);
        }
    }


    float sum = 0.0;
    for (int i = paddingSize; i < rows + paddingSize; i++) {
        for (int j = paddingSize; j < cols + paddingSize; j++) {
            sum = 0.0;
            for (int di = -(kernelSize / 2); di <= (kernelSize / 2); di++) {
                for (int dj = -(kernelSize / 2); dj <= (kernelSize / 2); dj++) {
                    int r = i + di;
                    int c = j + dj;

                    // Ensure boundary conditions
                    if (r >= 0 && r < rows + 2 * paddingSize && c >= 0 && c < cols + 2 * paddingSize) {
                        sum += float(paddedMatrix[r][c]) * mask[di + paddingSize][dj + paddingSize];
                    }
                }
            }

            int x = int(round(sum));
            if (x > L)
                x = L;
            else if (x < 0)
                x = 0;
            matrix[i - paddingSize][j - paddingSize] = x;
        }
    }

}

/// <summary>
/// Apply Floyd–Steinberg dithering algorithm on the image. (single channel)
/// </summary>
/// <param name="matrix">Input image</param>
/// <param name="output">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
void ErrorDiffusion(int matrix[MaxBMPSizeX][MaxBMPSizeY], int output[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols) {
    float threshold = (L / 2) + 1.0;

    float** img_dither = allocate2DArray_f(rows+1, cols+1);

    for (int i = 0; i < rows+1; i++) {
        for (int j = 0; j < cols+1; j++) {
            if (i < rows && j < cols)
                img_dither[i][j] = matrix[i][j];
            else
                img_dither[i][j] = 0;
        }
    }

    ////////////// Error Diffusion //////////////
    float old_pix, new_pix, quant_err;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            old_pix = img_dither[i][j];

            // Quantization(Binarize)
            if (img_dither[i][j] > threshold)
                new_pix = L;
            else
                new_pix = 0.0;

            img_dither[i][j] = new_pix;
            quant_err = old_pix - new_pix; // E

            if (j > 0)
                img_dither[i + 1][j - 1] = img_dither[i + 1][j - 1] + (quant_err * (3.0 / 16));
            img_dither[i + 1][j] = img_dither[i + 1][j] + (quant_err * (5.0 / 16));
            img_dither[i][j + 1] = img_dither[i][j + 1] + (quant_err * (7.0 / 16));
            img_dither[i + 1][j + 1] = img_dither[i + 1][j + 1] + (quant_err * (1.0 / 16));
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = int(round(img_dither[i][j]));
            output[i][j] = matrix[i][j];
        }
    }

    free2DArray_f(img_dither, rows+1);
}


///////////////////////////// Canny ////////////////////////////////
/// <summary>
/// Calculate the gradient (magnitude & phase) of image
/// </summary>
/// <param name="image">Input image</param>
/// <param name="gradient">Matrix used to save the gradient's magnitude</param>
/// <param name="direction">Matrix used to save the gradient's phase</param>
/// <param name="rows">The width fo image</param>
/// <param name="cols">The height of image</param>
void Gradient(int image[MaxBMPSizeX][MaxBMPSizeY], float gradient[MaxBMPSizeX][MaxBMPSizeY], float direction[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols) {
    int kenerl_size = 3;
    int paddingsize = kenerl_size / 2;

    // Define the vertical Sobel kernel
    float Vert[3][3] = {
        { -1, 0, 1 },
        { -2, 0, 2 },
        { -1, 0, 1 }
    };

    //Define the horizontal Sobel kernel
    float Hori[3][3] = {
        {  1,  2,  1 },
        {  0,  0,  0 },
        { -1, -2, -1 }
    };

    float** grad_x = allocate2DArray_f(rows, cols);
    float** grad_y = allocate2DArray_f(rows, cols);

    // Padding
    float** padded_x = ZeroPadding_f(image, rows, cols, kenerl_size);
    float** padded_y = ZeroPadding_f(image, rows, cols, kenerl_size);

    float sum_x, sum_y;
    for (int i = paddingsize; i < rows + paddingsize; i++) {
        for (int j = paddingsize; j < cols + paddingsize; j++) {
            sum_x = 0.0;
            sum_y = 0.0;
            for (int di = -paddingsize; di <= paddingsize; di++) {
                for (int dj = -paddingsize; dj <= paddingsize; dj++) {
                    sum_x += padded_x[i + di][j + dj] * Hori[di + paddingsize][dj + paddingsize];
                    sum_y += padded_y[i + di][j + dj] * Vert[di + paddingsize][dj + paddingsize];
                }
            }
            grad_x[i - paddingsize][j - paddingsize] = sum_x;
            grad_y[i - paddingsize][j - paddingsize] = sum_y;
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gradient[i][j] = sqrt((grad_x[i][j] * grad_x[i][j]) + (grad_y[i][j] * grad_y[i][j]));
            //float radian = grad_y[i][j] / grad_x[i][j];
            direction[i][j] = atan2(grad_y[i][j], grad_x[i][j]); //atan(radian);
        }
    }


    free2DArray_f(grad_x, rows);
    free2DArray_f(grad_y, rows);

    free2DArray_f(padded_x, rows + paddingsize * 2);
    free2DArray_f(padded_y, rows + paddingsize * 2);
}

float** NMS(float gradient_magnitude[MaxBMPSizeX][MaxBMPSizeY], float gradient_direction[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols) {
    float** nms_image = allocate2DArray_f(rows, cols);
    float neighbor1, neighbor2, mag, angle;
    mag = angle = 0;
    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            neighbor1 = neighbor2 = L;
            mag = gradient_magnitude[i][j];
            angle = gradient_direction[i][j] * 180.0 / PI;
            if (angle < 0)
                angle += 180.0;
            //////////////////////////////// Region judgement ////////////////////////////////
            // 0 deg
            if ((angle >= 0 && angle < 22.5) ||
                (angle >= 157.5 && angle <= 180)) {
                neighbor1 = gradient_magnitude[i][j + 1];
                neighbor2 = gradient_magnitude[i][j - 1];
            }

            // 45 deg
            else if( 22.5 <= angle && angle < 67.5){
                neighbor1 = gradient_magnitude[i + 1][j - 1];
                neighbor2 = gradient_magnitude[i - 1][j + 1];
            }
            // 90 deg
            else if( 67.5 <= angle && angle < 112.5){
                neighbor1 = gradient_magnitude[i + 1][j];
                neighbor2 = gradient_magnitude[i - 1][j];
            }
            // 135 deg
            else if(112.5 <= angle && angle < 157.5){
                neighbor1 = gradient_magnitude[i + 1][j + 1];
                neighbor2 = gradient_magnitude[i - 1][j - 1];
            }

            if (mag >= neighbor1 && mag >= neighbor2){
                nms_image[i][j] = mag;
            }
            else{
                nms_image[i][j] = 0;
            }
        }
    }

    return nms_image;
}

void Hysteresis_threshold(float** nms_image, float gradient_magnitude[MaxBMPSizeX][MaxBMPSizeY], int low_threshold, int high_threshhold, int rows, int cols) {
    float strong = L;
    float weak = L; //low_threshold * (1 - percentage) + percentage * high_threshhold;
    float background = 0;

    /////////////////////// Build Table ///////////////////////
    float** table = allocate2DArray_f(rows, cols); // to record the coordinate of boundary
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // three case: pixel > High_thr、Low_thr < pixel < High_thr、pixel < Low_thr
            float pixel = gradient_magnitude[i][j];
            // strong boundary
            if (pixel >= high_threshhold) {
                table[i][j] = strong;
            }
            // weak boundary
            else if (low_threshold < pixel && pixel < high_threshhold) {
                table[i][j] = weak;
            }
            // background
            else if (pixel <= low_threshold) {
                table[i][j] = background;
            }
        }
    }

    /////////////////////// Apply Table (for weak case) ///////////////////////
    /*
        |---------|---------|---------|
        | x-1,y-1 |  x-1,y  | x-1,y+1 |
        |---------|---------|---------|
        |  x,y-1  |  (x,y)  |  x,y+1  |
        |---------|---------|---------|
        | x+1,y-1 |  x+1,y  | x+1,y+1 |
        |---------|---------|---------|
    */

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            if (table[i][j] == weak) {
                if (table[i-1][j-1] == strong || table[i-1][j]   == strong ||
                    table[i-1][j+1] == strong || table[i][j-1]   == strong ||
                    table[i][j+1]   == strong || table[i+1][j-1] == strong ||
                    table[i+1][j]   == strong || table[i+1][j+1] == strong) {
                    nms_image[i][j] = weak;
                }
                else {
                    nms_image[i][j] = 0;
                }
            }
        }
    }

    free2DArray_f(table, rows);
}


/// <summary>
/// Edge detection algorithm, which includes 4 parts:
/// 1. GaussianFilter(): Eliminate the noise (keep the correct edge)
/// 2. Gradient(): Calculate the gradient to find the gradient magnitude & direction of image
/// 3. NMS(): Thining the edge and remove the unwanted edge
/// 4. Hysteresis_threshold(): Set the high & low threshold to divide strong edges and weak edges,
/// and connect the correct edge
/// 5. Canny()(): Ouput the result (single channel)
/// </summary>
/// <param name="image">Input image</param>
/// <param name="result">Matrix to save the result</param>
/// <param name="low_threshold">The low threshold for hysteresis threshold</param>
/// <param name="high_threshhold">The high threshold for hysteresis threshold</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="percentage">Number of pixels to keep (percentage), 0 < percentage < 1</param>
void Canny(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int low_threshold, int high_threshhold, int rows, int cols, float percentage){
    if(percentage > 1 || percentage < 0){
        cout << "% is < 1 and > 0!\n";
        return;
    }

    float bias;
    if(percentage == 1)
        bias = -0.001;
    else if(percentage == 0)
         bias = 0.001;
    else
        bias = 0;

    percentage += bias;

    low_threshold = low_threshold * percentage + (1 - percentage) * high_threshhold;

    //////////////////////////////////////////////////////////////
    int kernelsize = 5; // for gaussian filter
    int paddingsize = kernelsize / 2;

    int** paddedMatrix = EdgePadding(image, rows, cols, 5);
    GaussianFilter(paddedMatrix, image, rows, cols, 5, 1);

    free2DArray(paddedMatrix, rows + paddingsize * 2);

    Gradient(image, gradient_magnitude, gradient_direction, rows, cols);

    float** nms_img = NMS(gradient_magnitude, gradient_direction, rows, cols);

    Hysteresis_threshold(nms_img, gradient_magnitude, low_threshold, high_threshhold, rows, cols);


    //------------------------------------------------------------
    for (int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            result[i][j] = int(round(nms_img[i][j]));
        }
    }

    Binarize(result, result, rows, cols, low_threshold);

    free2DArray_f(nms_img, rows);
}

/////////////////////// Rotation ///////////////////////
/// <summary>
/// Use Bilinear Interpolation to interpolate the pixel value on the image and return
/// the interpolated pixel value
/// </summary>
/// <param name="image">The input image</param>
/// <param name="rotated_x">Input x-axis coordinate</param>
/// <param name="rotated_y">Input y-axis coordinate</param>
/// <returns>int type pixel</returns>
int BilinearInterpolation(int image[MaxBMPSizeX][MaxBMPSizeY], float rotated_x, float rotated_y) {
    /*
    (x_floor, y_ceil)       image(x_ceil, y_ceil)
                                               dy
                         X

                                              1-dy
    (x_floor, y_floor)       image(x_ceil, y_floor)
                     1-dx  dx
    */
    int x_floor = int(floor(rotated_x));
    int x_ceil  = x_floor + 1;
    int y_floor = int(floor(rotated_y));
    int y_ceil  = y_floor + 1;

    float dx = float(x_ceil) - rotated_x;
    float dy = float(y_ceil) - rotated_y;

    float tmp = image[x_floor][y_ceil] * dx * (1 - dy) + image[x_ceil][y_ceil] * (1 - dx) * (1 - dy) + image[x_floor][y_floor] * dx * dy + image[x_ceil][y_floor] * (1 - dx) * dy;

    int interpolated_pixel = int(round(tmp));

    return interpolated_pixel;
}

/// <summary>
/// Rotate the image with any angles (single channel, forward mapping)
/// </summary>
/// <param name="image">Input image (rotate)</param>
/// <param name="result">Matrix to save the rotated image</param>
/// <param name="width">The width of input image</param>
/// <param name="height">The height of input image</param>
/// <param name="new_w">The width after rotate (same as width)</param>
/// <param name="new_h">The height after rotate (same as height)</param>
/// <param name="deg">Dergree to rotate</param>
void Rotate_forward(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int width, int height, int &new_w, int &new_h, float deg){
    // input is degree, must change to radian for function.
    float theta = deg / 180.0 * PI;

    int new_width  = width;
    int new_height = height;

    //int new_width  = (int)ceil(fabs(height * cos(theta)) + fabs(width * sin(theta)));
    //int new_height = (int)ceil(fabs(width * cos(theta)) + fabs(height * sin(theta)));

    // Original center point
    float center_x = (float)width  / 2.0;
    float center_y = (float)height / 2.0;

    // New center point
    //float new_center_x = (float)new_width  / 2.0;
    //float new_center_y = (float)new_height / 2.0;

    int x, y;
    for (y = 0; y < new_height; y++) {  // Iterate over each pixel in the rotated image
        for (x = 0; x < new_width; x++) {
            // Translate the coord.to the center of the original image
            float x_translated = x - center_x;
            float y_translated = y - center_y;

            // Rotate the coordinates (R * [x, y].T) and shift to image's center
            float x_rotated = round(x_translated * cos(theta) - y_translated * sin(theta) + center_x);
            float y_rotated = round(x_translated * sin(theta) + y_translated * cos(theta) + center_y);

            // check if the coordinate is in the boundary
            if ((0 <= x_rotated && x_rotated < width) && (0 <= y_rotated && y_rotated < height)){
                /*
                // Perform bilinear interpolation
                int x_floor = int(floor(x_rotated));
                int x_ceil  = int(ceil(x_rotated));
                int y_floor = int(floor(y_rotated));
                int y_ceil  = int(ceil(y_rotated));

                float dx = x_rotated - x_floor;
                float dy = y_rotated - y_floor;

                // Interpolate the pixel values
                float interpolated_value = (1 - dx) * (1 - dy) * image[y_floor][x_floor] +
                                         dx * (1 - dy) * image[y_floor][x_ceil] +
                                         (1 - dx) * dy * image[y_ceil][x_floor] +
                                         dx * dy * image[y_ceil][x_ceil];

                // Assign the interpolated value to the corresponding pixel in the rotated image
                result[y][x] = int(round(interpolated_value));
                */
                // Fixed
                result[y][x] = BilinearInterpolation(image, y_rotated, x_rotated);
            }
        }
    }

    new_h = new_height;
    new_w = new_width;
}

//////////////////////////// Backward mapping ////////////////////////////
/// <summary>
/// Rotate the image with any angles (single channel, inverse/backward mapping)
/// </summary>
/// <param name="image">Input image (rotate)</param>
/// <param name="result">Matrix to save the rotated image</param>
/// <param name="width">The width of input image</param>
/// <param name="height">The height of input image</param>
/// <param name="new_w">The width after rotate (same as width)</param>
/// <param name="new_h">The height after rotate (same as height)</param>
/// <param name="deg">Dergree to rotate</param>
void Rotate_backward(int image[MaxBMPSizeX][MaxBMPSizeY], int result[MaxBMPSizeX][MaxBMPSizeY], int width, int height, int& new_w, int& new_h, float deg) {
    // input is degree, must change to radian.
    float theta = deg / 180.0 * PI;

    int new_width  = width;     //int(ceil(abs(width * cos(theta)) + abs(height * sin(theta))));
    int new_height = height;    //int(ceil(abs(height * cos(theta)) + abs(width * sin(theta))));

    new_h = new_height;
    new_w = new_width;
    // cout << "(" << new_width << ", " << new_height << ")" << endl;

    // Original center point
    float center_x = width / 2.0;
    float center_y = height / 2.0;

    // New center point
    //int new_center_x = new_width / 2;
    //int new_center_y = new_height / 2;

    int x, y;
    for (y = 0; y < new_height; y++) {  // Iterate over each pixel in the rotated image
        for (x = 0; x < new_width; x++) {
            // Translate the coordinates to the center of the original image
            float x_translated = x - center_x;
            float y_translated = y - center_y;

            // Rotate the coordinates
            int x_rotated = int(round(x_translated * cos(theta) + y_translated * sin(theta) + center_x));
            int y_rotated = int(round(-x_translated * sin(theta) + y_translated * cos(theta) + center_y));

            if ((0 <= x_rotated && x_rotated < width - 1) && (0 <= y_rotated && y_rotated < height - 1)) {
                result[x][y] = BilinearInterpolation(image, x_rotated, y_rotated);
            }
        }
    }
}


/////////////////////// Filter ///////////////////////
/// <summary>
/// Median filtering, you can input the kernel size to adjust the effect. (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="matrix">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernelSize">Size of a mask</param>
void MedianFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize) {
    // Single channel
    int paddingSize = kernelSize / 2;
    int* tmp = (int*)calloc(kernelSize * kernelSize, sizeof(int));
    int index = 0;

    if (kernelSize == 3) {
        for (int i = 1; i < rows + 1; i++) {
            for (int j = 1; j < cols + 1; j++) {
                // 3x3
                tmp[0] = paddedMatrix[i - 1][j - 1];
                tmp[1] = paddedMatrix[i - 1][j];
                tmp[2] = paddedMatrix[i - 1][j + 1];

                tmp[3] = paddedMatrix[i][j - 1];
                tmp[4] = paddedMatrix[i][j];
                tmp[5] = paddedMatrix[i][j + 1];

                tmp[6] = paddedMatrix[i + 1][j - 1];
                tmp[7] = paddedMatrix[i + 1][j];
                tmp[8] = paddedMatrix[i + 1][j + 1];

                sort(tmp, tmp + kernelSize * kernelSize);
                matrix[i - 1][j - 1] = tmp[kernelSize * kernelSize / 2]; //Median(tmp, kernelSize * kernelSize);
            }
        }
        return;
    }

    else{
        for (int i = paddingSize; i < rows + paddingSize; i++) {
            for (int j = paddingSize; j < cols + paddingSize; j++) {
                index = 0;
                for (int di = -(kernelSize / 2); di <= (kernelSize / 2); di++) {
                    for (int dj = -(kernelSize / 2); dj <= (kernelSize / 2); dj++) {
                        tmp[index++] = paddedMatrix[i + di][j + dj];
                    }
                }

                sort(tmp, tmp + kernelSize * kernelSize);
                //selectionSort(tmp, kernelSize * kernelSize);

                matrix[i - paddingSize][j - paddingSize] = Median(tmp, kernelSize * kernelSize);
            }
        }
    }

    free(tmp);

    return;
}

/// <summary>
/// Create the moving average mask to smooth the image (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="matrix">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernelSize">Size of a mask</param>
void SmoothFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize) {
    // Single channel
    int paddingSize = (kernelSize / 2);

    // Allocation
    float** smooth_mask = allocate2DArray_f(kernelSize, kernelSize);

    // Smooth mask
    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            smooth_mask[i][j]  = (1.0 / float(kernelSize * kernelSize));
        }
    }

    float sum = 0.0;
    for (int i = paddingSize; i < rows + paddingSize; i++) {
        for (int j = paddingSize; j < cols + paddingSize; j++) {
            sum = 0.0;
            for (int di = -(kernelSize / 2); di <= (kernelSize / 2); di++) {
                for (int dj = -(kernelSize / 2); dj <= (kernelSize / 2); dj++) {
                    int r = i + di;
                    int c = j + dj;

                    // Ensure boundary conditions
                    if (r >= 0 && r < rows + 2 * paddingSize && c >= 0 && c < cols + 2 * paddingSize) {
                        sum += float(paddedMatrix[r][c]) * smooth_mask[di + paddingSize][dj + paddingSize];
                    }
                }
            }
            matrix[i - paddingSize][j - paddingSize] = int(round(sum));
        }
    }

    // free mask
    free2DArray_f(smooth_mask, kernelSize);
}

/// <summary>
/// Create the Binomial coeiffient (1D) and transfer into mask (2D),
/// including high pass & low pass. (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="matrix">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernelSize">Size of a mask</param>
/// <param name="HPF">flag to decide if its LPF(false) or HPF (true)</param>
void Bino_Filter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize, bool HPF) {
    float* coefficients = Bino_Coef(kernelSize);
    float** mask = Bino_mask(coefficients, kernelSize);

    if (HPF) {
        for (int i = 0; i < kernelSize; i++) {
            for (int j = 0; j < kernelSize; j++) {
                if (i== kernelSize / 2 && j== kernelSize / 2) {
                    mask[i][j] = 1.0 - mask[i][j]; // 1 - center points
                }
                else {
                    mask[i][j] = - mask[i][j]; // 0 - other points
                }
            }
        }
    }

    Convolve2D(paddedMatrix, mask, matrix, rows, cols, kernelSize);


    if (HPF) {
        //Linear mapping
        /*
        //float max = Max(paddedMatrix, rows, cols);
        //float min = Min(paddedMatrix, rows, cols);
        for (int i = 0; i < rows + paddingSize * 2;i++) {
            for (int j = 0; j < cols + paddingSize * 2;j++) {
                int x = paddedMatrix[i][j];
                if( x > L )
                    paddedMatrix[i][j] = L;
                else if( x < 0)
                    paddedMatrix[i][j] = 0;
                //paddedMatrix[i][j] = int(round(L * (x - min) / (max - min)));
            }
        }
        */

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                int x = matrix[i][j];
                // Cilp
                if( x > L )
                    matrix[i][j] = L;
                else if( x < 0)
                    matrix[i][j] = 0;
            }
        }
    }

    // ----------------------------------------
    free(coefficients);
    free2DArray_f(mask, kernelSize);

    return;
}

/// <summary>
/// Edge detection using Sobel (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="grad_x">To save the gradiant magnitude of x-axis</param>
/// <param name="grad_y">To save the gradiant magnitude of y-axis</param>
/// <param name="sobel">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="percentage">Number of pixels to keep (percentage), 0 < percentage < 1</param>
void Sobel(int** paddedMatrix, int grad_x[MaxBMPSizeX][MaxBMPSizeY], int grad_y[MaxBMPSizeX][MaxBMPSizeY],
           int sobel[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float percentage) {

    /****************************
    Recommand：
    1. Use G channel as input, and output as single channel
    2. Threhold = 0.5
    ****************************/

    // sobel --> output result
    int kernel_size = 3;
    // Define the vertical Sobel kernel
    float Vert[3][3] = {
        { -1, 0, 1 },
        { -2, 0, 2 },
        { -1, 0, 1 }
    };

    //Define the horizontal Sobel kernel
    float Hori[3][3] = {
        {  1,  2,  1 },
        {  0,  0,  0 },
        { -1, -2, -1 }
    };

    float** tmp_mask = allocate2DArray_f(kernel_size, kernel_size);

    for (int i = 0; i < kernel_size; i++) {
        for (int j = 0; j < kernel_size; j++) {
            tmp_mask[i][j] = Vert[i][j];
        }
    }

    Convolve2D(paddedMatrix, tmp_mask, grad_y, rows, cols, kernel_size);

    for (int i = 0; i < kernel_size; i++) {
        for (int j = 0; j < kernel_size; j++) {
            tmp_mask[i][j] = Hori[i][j];
        }
    }

    Convolve2D(paddedMatrix, tmp_mask, grad_x, rows, cols, kernel_size);

    free2DArray_f(tmp_mask, kernel_size);

    //////////////// Sobel calculation ////////////////
    int** tmp_result = allocate2DArray(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            tmp_result[i][j] = int(round(sqrt(float(grad_x[i][j] * grad_x[i][j]) + float(grad_y[i][j] * grad_y[i][j]))));
        }
    }

    // Filtering (Threshold)
    float max = Max(tmp_result, rows, cols);
    float min = Min(tmp_result, rows, cols);

    if(percentage > 1.0 || percentage == 0){
        percentage = percentage / 100.0;
        float bias;
        if(percentage == 0.0)
            bias = 1e-5;
        else
            bias = 0;

        percentage += bias;
    }

    cout << "Processed Threshold: " << percentage * 100 << "%\n";
    float threshold = (L/2) -  (L/2) * percentage; // adjust the threshold

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int x = tmp_result[i][j];
            int pixel = L * (x - min) / (max - min); //linear mapping
            if (pixel >= threshold)
                sobel[i][j] = L;
            else
                sobel[i][j] = 0;
        }
    }

    free2DArray(tmp_result, rows);
}

/// <summary>
/// Edge detection using Laplacian (single channel)
/// The mask here used is:
///  { -1, -1, -1 }
///  { -1,  8, -1 }
///  { -1, -1, -1 }
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="matrix">Matrix to save the result</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="percentage">Number of pixels to keep (percentage), 0 < percentage < 1</param>
void Laplacian(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, float percentage) {
    int kernel_size = 3;

    // Define the vertical Sobel kernel
    float laplacian_kernel[3][3] = {
        { -1, -1, -1 },
        { -1,  8, -1 },
        { -1, -1, -1 }
    };

    float** tmp_mask = allocate2DArray_f(kernel_size, kernel_size);

    for (int i = 0; i < kernel_size; i++) {
        for (int j = 0; j < kernel_size; j++) {
            tmp_mask[i][j] = laplacian_kernel[i][j];
        }
    }

    Convolve2D(paddedMatrix, tmp_mask, matrix, rows, cols, kernel_size);

    free2DArray_f(tmp_mask, kernel_size);

    //////////////// Laplacian calculation ////////////////
    int** tmp_result = allocate2DArray(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            tmp_result[i][j] = matrix[i][j];
        }
    }

    // Filtering (Threshold)
    //float max = Max(tmp_result, rows, cols);
    //float min = Min(tmp_result, rows, cols);

    if (percentage > 1.0) {
        percentage = percentage / 100.0;
    }

    float threshold = (L/2) * (1 - percentage);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int x = tmp_result[i][j];
            int pixel = x;//L * (x - min) / (max - min); //linear mapping
            if (pixel >= threshold)
                matrix[i][j] = L;
            else
                matrix[i][j] = 0;
        }
    }

    free2DArray(tmp_result, rows);
}

/// <summary>
/// Create any Gaussian mask and filter the input image (single channel)
/// </summary>
/// <param name="paddedMatrix">Input padded image</param>
/// <param name="matrix">Save the image after gaussian filtering</param>
/// <param name="rows">The width of image</param>
/// <param name="cols">The height of image</param>
/// <param name="kernelSize">The size of gaussian kernel</param>
/// <param name="sigma">For gaussian formula</param>
void GaussianFilter(int** paddedMatrix, int matrix[MaxBMPSizeX][MaxBMPSizeY], int rows, int cols, int kernelSize, float sigma) {
    /////////// Mask generation ///////////
    if (kernelSize % 2 == 0)
        kernelSize++;

    float** x = allocate2DArray_f(kernelSize, kernelSize);
    float** y = allocate2DArray_f(kernelSize, kernelSize);
    float** gaussian_kernel = allocate2DArray_f(kernelSize, kernelSize);

    float* tmp = (float*)calloc(kernelSize, sizeof(float));
    float val = -kernelSize / 2;
    for (int i = 0; i < kernelSize; i++) {
        tmp[i] = val++;
    }

    int index = 0;
    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            x[i][j] = tmp[index];
            y[j][i] = tmp[index];
            index = (index + 1) % kernelSize;
        }
    }

    ////////////////// Gaussian kernel calculation //////////////////
    float sum = 0;
    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            gaussian_kernel[i][j] = exp(-(x[i][j] * x[i][j] + y[i][j] * y[i][j]) / (2 * sigma * sigma));
            sum += gaussian_kernel[i][j];
        }
    }

    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            gaussian_kernel[i][j] /= sum;
        }
    }

    ////////////////// Convolve in 2D //////////////////
    Convolve2D(paddedMatrix, gaussian_kernel, matrix, rows, cols, kernelSize);

    free2DArray_f(x, kernelSize);
    free2DArray_f(y, kernelSize);
    free2DArray_f(gaussian_kernel, kernelSize);
    free(tmp);
}
