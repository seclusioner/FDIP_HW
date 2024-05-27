#include "DIP_Func.h"
#include "Padding.h"

using namespace std;

///////////////////////// Padding Algorithms /////////////////////////
int** Padding::zeroPadding(int** matrix, int H, int W, int padWidth, int padd_const) {
    int newRows = H + 2 * padWidth;
    int newCols = W + 2 * padWidth;

    int** paddedMatrix = allocateMatrix(newRows, newCols);
    if (padd_const != 0) {
        for (int i = 0; i < newRows; ++i) {
            for (int j = 0; j < newCols; ++j) {
                paddedMatrix[i][j] = padd_const;
            }
        }
    }
    for (int i = padWidth; i < newRows - padWidth; ++i) {
        for (int j = padWidth; j < newCols - padWidth; ++j) {
            paddedMatrix[i][j] = matrix[i - padWidth][j - padWidth];
        }
    }

    return paddedMatrix;
}

int** Padding::edgePadding(int** matrix, int H, int W, int padWidth) {
    int newRows = H + 2 * padWidth;
    int newCols = W + 2 * padWidth;

    int** paddedMatrix = allocateMatrix(newRows, newCols);

    // Copy the original matrix
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            paddedMatrix[i + padWidth][j + padWidth] = matrix[i][j];
        }
    }

    // Fill the edge regions
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < padWidth; ++j) {
            paddedMatrix[i + padWidth][j] = matrix[i][0]; // left edge
            paddedMatrix[i + padWidth][newCols - 1 - j] = matrix[i][W - 1]; // right edge
        }
    }

    for (int j = 0; j < newCols; ++j) {
        for (int i = 0; i < padWidth; ++i) {
            paddedMatrix[i][j] = paddedMatrix[padWidth][j]; // top edge
            paddedMatrix[newRows - 1 - i][j] = paddedMatrix[newRows - 1 - padWidth][j]; // bottom edge
        }
    }

    return paddedMatrix;
}

int** Padding::reflectPadding(int** image, int H, int W, int padWidth) {
    int padding_size = padWidth;
    int new_H = H + 2 * padding_size;
    int new_W = W + 2 * padding_size;
    int** padded_image = allocateMatrix(new_H, new_W);

    ////////////////////////// Fill //////////////////////////
    // Fill the center region with the original image
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_image[i + padding_size][j + padding_size] = image[i][j];
        }
    }

    ////////////////////////// Padding //////////////////////////
    // Up & Down
    for (int i = 0; i < padding_size; i++) {
        for (int j = 0; j < W; j++) {
            padded_image[padding_size - 1 - i][padding_size + j] = image[i + 1][j]; // top
            padded_image[H + padding_size + i][padding_size + j] = image[H - 2 - i][j]; // bottom
        }
    }

    // Left & Right
    for (int i = 0; i < padding_size; i++) {
        for (int j = 0; j < H; j++) {
            padded_image[padding_size + j][padding_size - 1 - i] = image[j][i + 1]; // left
            padded_image[padding_size + j][W + padding_size + i] = image[j][W - 2 - i]; // right
        }
    }

    // Fill the four corners
    for (int i = 0; i < padding_size; i++) {
        for (int j = 0; j < padding_size; j++) {
            padded_image[padding_size - 1 - i][padding_size - 1 - j] = image[i + 1][j + 1];         // top-left
            padded_image[padding_size - 1 - i][W + padding_size + j] = image[i + 1][W - 2 - j];     // top-right
            padded_image[H + padding_size + i][padding_size - 1 - j] = image[H - 2 - i][j + 1];     // bottom-left
            padded_image[H + padding_size + i][W + padding_size + j] = image[H - 2 - i][W - 2 - j]; // bottom-right
        }
    }

    return padded_image;
}

int** Padding::symmetricPadding(int** image, int H, int W, int padWidth) {
    int new_H = H + 2 * padWidth;
    int new_W = W + 2 * padWidth;
    int** padded_image = allocateMatrix(new_H, new_W);

    ////////////////////////// Fill //////////////////////////
    // Fill the center region with the original image
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_image[i + padWidth][j + padWidth] = image[i][j];
        }
    }

    ////////////////////////// Padding //////////////////////////
    // Pad the top and bottom rows
    for (int i = 0; i < padWidth; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_image[i][j + padWidth] = image[padWidth - i - 1][j];
            padded_image[new_H - i - 1][j + padWidth] = image[H - padWidth + i][j];
        }
    }

    // Pad the left and right columns
    for (int i = 0; i < new_H; i++) {
        for (int j = 0; j < padWidth; j++) {
            padded_image[i][j] = padded_image[i][2 * padWidth - j - 1];
            padded_image[i][new_W - j - 1] = padded_image[i][new_W - 2 * padWidth + j];
        }
    }

    return padded_image;
}

int** Padding::maximumPadding(int** matrix, int H, int W, int pad_width) {
    int new_rows = H + 2 * pad_width;
    int new_cols = W + 2 * pad_width;

    int** padded_matrix = allocateMatrix(new_rows, new_cols);

    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_matrix[i + pad_width][j + pad_width] = matrix[i][j];
        }
    }

    int* col_max = new int[W];
    for (int j = 0; j < W; ++j) {
        col_max[j] = matrix[0][j];
        for (int i = 1; i < H; ++i) {
            col_max[j] = max(col_max[j], matrix[i][j]);
        }
    }

    /////////////////// Padding ///////////////////
    for (int i = 0; i < pad_width; ++i) {
        for (int j = pad_width; j < pad_width + W; ++j) {
            padded_matrix[i][j] = col_max[j - pad_width];
            padded_matrix[new_rows - i - 1][j] = col_max[j - pad_width];
        }
    }

    int* row_max = new int[H];
    for (int i = 0; i < H; ++i) {
        row_max[i] = matrix[i][0];
        for (int j = 1; j < W; ++j) {
            row_max[i] = max(row_max[i], matrix[i][j]);
        }
    }

    for (int i = pad_width; i < pad_width + H; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = row_max[i - pad_width];
            padded_matrix[i][new_cols - j - 1] = row_max[i - pad_width];
        }
    }

    int max_val = *max_element(col_max, col_max + W);
    for (int i = 0; i < pad_width; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = max_val;
            padded_matrix[i][new_cols - j - 1] = max_val;
            padded_matrix[new_rows - i - 1][j] = max_val;
            padded_matrix[new_rows - i - 1][new_cols - j - 1] = max_val;
        }
    }

    delete[] col_max;
    delete[] row_max;

    return padded_matrix;
}

int** Padding::minimumPadding(int** matrix, int H, int W, int pad_width) {
    int new_rows = H + 2 * pad_width;
    int new_cols = W + 2 * pad_width;

    int** padded_matrix = allocateMatrix(new_rows, new_cols);

    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_matrix[i + pad_width][j + pad_width] = matrix[i][j];
        }
    }

    int* col_min = new int[W];
    for (int j = 0; j < W; ++j) {
        col_min[j] = matrix[0][j];
        for (int i = 1; i < H; ++i) {
            col_min[j] = min(col_min[j], matrix[i][j]);
        }
    }

    for (int i = 0; i < pad_width; ++i) {
        for (int j = pad_width; j < pad_width + W; ++j) {
            padded_matrix[i][j] = col_min[j - pad_width];
            padded_matrix[new_rows - i - 1][j] = col_min[j - pad_width];
        }
    }

    int* row_min = new int[H];
    for (int i = 0; i < H; ++i) {
        row_min[i] = matrix[i][0];
        for (int j = 1; j < W; ++j) {
            row_min[i] = min(row_min[i], matrix[i][j]);
        }
    }

    for (int i = pad_width; i < pad_width + H; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = row_min[i - pad_width];
            padded_matrix[i][new_cols - j - 1] = row_min[i - pad_width];
        }
    }

    int min_val = *min_element(col_min, col_min + W);
    for (int i = 0; i < pad_width; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = min_val;
            padded_matrix[i][new_cols - j - 1] = min_val;
            padded_matrix[new_rows - i - 1][j] = min_val;
            padded_matrix[new_rows - i - 1][new_cols - j - 1] = min_val;
        }
    }

    delete[] col_min;
    delete[] row_min;

    return padded_matrix;
}

int** Padding::meanPadding(int** matrix, int H, int W, int pad_width) {
    int new_rows = H + 2 * pad_width;
    int new_cols = W + 2 * pad_width;

    int** padded_matrix = allocate2DArray(new_rows, new_cols);

    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_matrix[i + pad_width][j + pad_width] = matrix[i][j];
        }
    }

    double* col_mean = new double[W];
    for (int j = 0; j < W; ++j) {
        col_mean[j] = 0;
        for (int i = 0; i < H; ++i) {
            col_mean[j] += matrix[i][j];
        }
        col_mean[j] /= H;
    }

    for (int i = 0; i < pad_width; ++i) {
        for (int j = pad_width; j < pad_width + W; ++j) {
            padded_matrix[i][j] = round(col_mean[j - pad_width] * 1000) / 1000.0;
            padded_matrix[new_rows - i - 1][j] = round(col_mean[j - pad_width] * 1000) / 1000.0;
        }
    }

    double* row_mean = new double[H];
    for (int i = 0; i < H; ++i) {
        row_mean[i] = 0;
        for (int j = 0; j < W; ++j) {
            row_mean[i] += matrix[i][j];
        }
        row_mean[i] /= W;
    }

    for (int i = pad_width; i < pad_width + H; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = round(row_mean[i - pad_width] * 1000) / 1000.0;
            padded_matrix[i][new_cols - j - 1] = round(row_mean[i - pad_width] * 1000) / 1000.0;
        }
    }

    double total_mean = 0;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            total_mean += matrix[i][j];
        }
    }
    total_mean /= (H * W);

    for (int i = 0; i < pad_width; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = round(total_mean * 1000) / 1000.0;
            padded_matrix[i][new_cols - j - 1] = round(total_mean * 1000) / 1000.0;
            padded_matrix[new_rows - i - 1][j] = round(total_mean * 1000) / 1000.0;
            padded_matrix[new_rows - i - 1][new_cols - j - 1] = round(total_mean * 1000) / 1000.0;
        }
    }

    delete[] col_mean;
    delete[] row_mean;

    return padded_matrix;
}

int** Padding::medianPadding(int** matrix, int H, int W, int pad_width) {
    int new_rows = H + 2 * pad_width;
    int new_cols = W + 2 * pad_width;

    int** padded_matrix = allocate2DArray(new_rows, new_cols);

    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            padded_matrix[i + pad_width][j + pad_width] = matrix[i][j];
        }
    }

    int* column_values = new int[H];
    for (int j = 0; j < W; ++j) {
        for (int i = 0; i < H; ++i) {
            column_values[i] = matrix[i][j];
        }
        int median = Median(column_values, H);
        for (int i = 0; i < pad_width; ++i) {
            padded_matrix[i][j + pad_width] = median;
            padded_matrix[new_rows - i - 1][j + pad_width] = median;
        }
    }
    delete[] column_values;

    int* row_values = new int[W];
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            row_values[j] = matrix[i][j];
        }
        int median = Median(row_values, W);
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i + pad_width][j] = median;
            padded_matrix[i + pad_width][new_cols - j - 1] = median;
        }
    }
    delete[] row_values;

    int* all_values = new int[H * W];
    int idx = 0;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            all_values[idx++] = matrix[i][j];
        }
    }
    int total_median = Median(all_values, H * W);
    delete[] all_values;

    for (int i = 0; i < pad_width; ++i) {
        for (int j = 0; j < pad_width; ++j) {
            padded_matrix[i][j] = total_median;
            padded_matrix[i][new_cols - j - 1] = total_median;
            padded_matrix[new_rows - i - 1][j] = total_median;
            padded_matrix[new_rows - i - 1][new_cols - j - 1] = total_median;
        }
    }

    return padded_matrix;
}

//////////////////// Dynamic allocate ////////////////////
int** Padding::allocateMatrix(int rows, int cols) {
    int** matrix = (int**)std::calloc(rows, sizeof(int*));
    if (matrix == NULL) {
        cout << "Allocation failed!\n";
        exit(1);
    }
    for (int i = 0; i < rows; ++i) {
        matrix[i] = (int*)std::calloc(cols, sizeof(int));
    }
    return matrix;
}

void Padding::freeMatrix(int** matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        std::free(matrix[i]);
    }
    std::free(matrix);
}

int** Pad(int matrix[MaxBMPSizeX][MaxBMPSizeY], int H, int W, int padWidth, PaddingMode mode, int padd_const) {
    int** tmp = Padding::allocateMatrix(H, W);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            tmp[i][j] = matrix[i][j];
        }
    }

    ///////////////////// Padding /////////////////////
    int** padded_img;
    switch (mode)
    {
        case constant:
            padded_img = Padding::zeroPadding(tmp, H, W, padWidth, padd_const);
            return padded_img;
        case edge:
            padded_img = Padding::edgePadding(tmp, H, W, padWidth);
            return padded_img;
        case reflect:
            padded_img = Padding::reflectPadding(tmp, H, W, padWidth);
            return padded_img;
        case symmetric:
            padded_img = Padding::symmetricPadding(tmp, H, W, padWidth);
            return padded_img;
        case maximum:
            padded_img = Padding::maximumPadding(tmp, H, W, padWidth);
            return padded_img;
        case minimum:
            padded_img = Padding::minimumPadding(tmp, H, W, padWidth);
            return padded_img;
        case mean:
            padded_img = Padding::meanPadding(tmp, H, W, padWidth);
            return padded_img;
        case median:
            padded_img = Padding::medianPadding(tmp, H, W, padWidth);
            return padded_img;
    }

    Padding::freeMatrix(tmp, H);
}
