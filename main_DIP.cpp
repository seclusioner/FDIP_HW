/*===========================================================================
  Source¡G
  This demonstrative example is provided by the teaching assistant,
  Mr. Shih-Hung Liao (¹ù¥@§»), and modified by the instructor, Prof. Lan.

  (1) It can be compiled and executed correctly under the DEV-C++, and Visual C++
      environments.
  (2) In order to run this program, you should also have the "bmp.h" and
      "bmp.cpp" files installed in your current directory or whichever directory
      the C++ compiler is directed to search for.
  (3) The DEV-C++ is a free C++ development environment that is recommended for
      this course. Refer to http://www.bloodshed.net/dev/devcpp.html.

                             Apr. 3, 2006
  --------------------------------------------------------------------------------
  DIP Homework for 2024.03. ~ 2024.04.
  Author¡G D.S.

============================================================================*/

#include "bmp.h"
#include "DIP_Func.h"
#include "Padding.h"

#define MAXSTRING 100

int KERNEL;
int PADDING_SIZE;

using namespace std;

extern int L;

// input image(.bmp) RGB channels
int R[MaxBMPSizeX][MaxBMPSizeY]; // MaxBMPSizeX and MaxBMPSizeY are defined in "bmp.h"
int G[MaxBMPSizeX][MaxBMPSizeY];
int B[MaxBMPSizeX][MaxBMPSizeY];

// output image(.bmp) RGB channels
int r[MaxBMPSizeX][MaxBMPSizeY];
int g[MaxBMPSizeX][MaxBMPSizeY];
int b[MaxBMPSizeX][MaxBMPSizeY];

// tmp
//int grad_x[MaxBMPSizeX][MaxBMPSizeY];
//int grad_y[MaxBMPSizeX][MaxBMPSizeY];

// Register
int tmp_R[MaxBMPSizeX][MaxBMPSizeY];
//int tmp_G[MaxBMPSizeX][MaxBMPSizeY];
//int tmp_B[MaxBMPSizeX][MaxBMPSizeY];

//////////////////////// Functions ////////////////////////
void User();
int kernel_input();


int main(int argc, char *argv[])
{
    cout << "Start ...... \n";
    cout << "------------------------------------------------------------------------\n";

    //**********************************************************************
    cout << "Existing test file: \n" <<
    "test_images\\boat.bmp\n" <<
    "test_images\\framed_lena_color_256.bmp\n" <<
    "test_images\\house.bmp\n" <<
    "test_images\\lena.bmp\n" <<
    "test_images\\lena_pepper_and_salt_noise10%.bmp\n" <<
    "test_images\\lena_std.bmp\n" <<
    "test_images\\lighthouse.bmp\n" <<
    "test_images\\penguin.bmp\n" <<
    "test_images\\yuntech.bmp\n\n";
    cout << "------------------------------------------------------------------------\n";

    //*************************************************************************/
    User(); // user interface


    /* ------------------------------------------------------------------------ */

    cout << "\bFinished!\n";

	system("PAUSE"); /* so that the command window holds a while */
	return 0;
}

/////////////////////////// User Interface ///////////////////////////
int kernel_input() {
    int kernel;
    cout << "Input a kernel size: ";
    cin >> kernel;

    if (kernel == 0) {
        cout << "kernel size must > 0!\n";
        exit(1);
    }

    // temporary
    if (kernel % 2 == 0) kernel++;

    return kernel;
}

void User() {
    ///////////////// Open File /////////////////
    char filename[MAXSTRING];
    cout << "Enter the input file name: ";
    cin >> filename;

    char* input_image = (char*)filename;// "test_images\\lighthouse.bmp"
    cout << "Open \"" << input_image << "\" .\n";

    // Input if the image is gray scale
    int Gray_flag = 0;
    char GrayScale;
    cout << "Input image is Grayscale? [Y/N] : ";
    cin >> GrayScale;

    // Open .bmp file function ...
    int width, height;
    if (GrayScale == 'Y') {
        Gray_flag = 1;
        cout << "Input image is single channel.\n";
        open_bmp(input_image, R, R, R, width, height);
    }
    else if(GrayScale == 'N') {
        Gray_flag = 0;
        cout << "Input image is rgb channel.\n";
        open_bmp(input_image, R, G, B, width, height);
    }
    else {
        cout << "Input is invalid\n";
        exit(1);
    }


    // Choose according to the Problem (Algorithm)
    int problem;

    //////////////////////// Choose Problem ////////////////////////
    cout << "----------------------------------------------------------\n";
    cout << "| 01                    Median Filter                    |\n";
    cout << "| 02  Binomial Lowpass Filtering and Highpass Filtering  |\n";
    cout << "| 03                Histogram Equalization               |\n";
    cout << "| 04                    Edge Sharpening                  |\n";
    cout << "| 05                    Error Diffusion                  |\n";
    cout << "| 06                       Rotation                      |\n";
    cout << "| 07                  Basic Edge Detection               |\n";
    cout << "----------------------------------------------------------\n";

    cout << "Choose a problem: ";
    cin >> problem;

    cout << "Problem 0" << problem << ": \n";

    /*********************************************************************
    Main program (Digital Image Processing)
    **********************************************************************/
    switch (problem)
    {
        case 1:
            KERNEL = kernel_input();
            PADDING_SIZE = KERNEL / 2;

            cout << "Kernel size: " << KERNEL << "\n";
            cout << "Apply " << KERNEL << " x " << KERNEL << " median filter ..... \n";
            if (Gray_flag) { // Single channel
                int** paddedMatrix = Pad(R, height, width, PADDING_SIZE, edge, 0);
                MedianFilter(paddedMatrix, r, width, height, KERNEL);
                free2DArray(paddedMatrix, width + PADDING_SIZE * 2);
            }
            else { // RGB channel
                int** paddedMatrix_R = Pad(R, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_G = Pad(G, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_B = Pad(B, height, width, PADDING_SIZE, edge, 0);

                MedianFilter(paddedMatrix_R, r, width, height, KERNEL);
                MedianFilter(paddedMatrix_G, g, width, height, KERNEL);
                MedianFilter(paddedMatrix_B, b, width, height, KERNEL);

                free2DArray(paddedMatrix_R, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_G, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_B, width + PADDING_SIZE * 2);
            }
            break;
        case 2:
            // Decide LPF / HPF
            bool HPF;
            char filter;
            cout << "Choose LPF / HPF [L/H]: ";
            cin >> filter;

            // Input a kernel size (3x3) & padding size (kernel size / 2)
            KERNEL = kernel_input();
            PADDING_SIZE = KERNEL / 2;

            if (filter == 'L') {
                HPF = false;
                cout << "Apply Binomial " << KERNEL << " x " << KERNEL<< " Low Pass Filter ..... \n";
            }
            else if(filter == 'H'){
                HPF = true;
                cout << "Apply Binomial " << KERNEL << " x " << KERNEL << " HIGH Pass Filter ..... \n";
            }
            else {
                cout << "Input is invalid!\n";
                exit(1);
            }

            if (Gray_flag) {
                int** paddedMatrix = Pad(R, height, width, PADDING_SIZE, edge, 0);
                Bino_Filter(paddedMatrix, r, width, height, KERNEL, HPF);
                free2DArray(paddedMatrix, width + PADDING_SIZE * 2);
            }
            else{
                int** paddedMatrix_R = Pad(R, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_G = Pad(G, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_B = Pad(B, height, width, PADDING_SIZE, edge, 0);

                Bino_Filter(paddedMatrix_R, r, width, height, KERNEL, HPF);
                Bino_Filter(paddedMatrix_G, g, width, height, KERNEL, HPF);
                Bino_Filter(paddedMatrix_B, b, width, height, KERNEL, HPF);

                free2DArray(paddedMatrix_R, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_G, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_B, width + PADDING_SIZE * 2);
            }
            break;

        case 3:
            cout << "Apply Histogram Equalization ..... \n";
            if(Gray_flag){
                Histo_Equal(R, r, width, height);
            }
            else {
                Grayscale(R, G, B, tmp_R, width, height);
                Gray_flag = 1;
                Histo_Equal(R, r, width, height);
            }

            break;
        case 4: // Edge Sharpening
            cout << "Apply Edge Sharpening ..... \n";

            KERNEL = 3;
            PADDING_SIZE = KERNEL / 2;

            if (Gray_flag) {
                int** paddedMatrix = Pad(R, height, width, PADDING_SIZE, constant, 0);
                EdgeSharpen(paddedMatrix, r, width, height, 0.2);
                free2DArray(paddedMatrix, width + PADDING_SIZE * 2);
            }
            else {
                int** paddedMatrix_R = Pad(R, height, width, PADDING_SIZE, constant, 0);
                int** paddedMatrix_G = Pad(G, height, width, PADDING_SIZE, constant, 0);
                int** paddedMatrix_B = Pad(B, height, width, PADDING_SIZE, constant, 0);

                EdgeSharpen(paddedMatrix_R, r, width, height, 0.2);
                EdgeSharpen(paddedMatrix_G, g, width, height, 0.2);
                EdgeSharpen(paddedMatrix_B, b, width, height, 0.2);

                free2DArray(paddedMatrix_R, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_G, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_B, width + PADDING_SIZE * 2);
            }
            break;
        case 5: // Error Diffusion (OK)
            cout << "Apply Error Diffusion ..... \n";
            if(Gray_flag){
                ErrorDiffusion(R, r, width, height);
            }
            else{
                ErrorDiffusion(R, r, width, height);
                ErrorDiffusion(G, g, width, height);
                ErrorDiffusion(B, b, width, height);
            }
            break;
        case 6: // Rotation
            int new_w, new_h;

            float degree;
            cout << "Enter a degree to rotate: ";
            cin >> degree;
            cout << "Rotate \"" << input_image << "\" " << degree << " degs. \n";

            if(Gray_flag){
                Rotate_forward(R, r, width, height, new_w, new_h, degree);
            }
            else{
                Rotate_forward(R, r, width, height, new_w, new_h, degree);
                Rotate_forward(G, g, width, height, new_w, new_h, degree);
                Rotate_forward(B, b, width, height, new_w, new_h, degree);
            }
            break;
        case 7: // Edge Detection
            float percentage;
            cout << "Enter the percentage of retain pixels (0 < % < 1): ";
            cin >> percentage;

            cout << "Apply Edge Detection (Canny) ..... \n";
            if (Gray_flag) {
                Canny(R, r, 120, 200, width, height, percentage);
            }
            else {
                Grayscale(R, G, B, tmp_R, width, height);
                Canny(tmp_R, r, 120, 200, width, height, percentage);
                Gray_flag = 1;
            }
            break;
        default:  // This region is for test / record
            cout << "Others ... \n"; // "No this problem! must be 1 ~ 7\n";

            // Settings param .
            KERNEL = 13;
            float sigma = 2;
            PADDING_SIZE = KERNEL / 2;

            if(Gray_flag){ // single channel
                //int** paddedMatrix = EdgePadding(R, width, height, KERNEL);
                int** paddedMatrix = Pad(R, height, width, PADDING_SIZE, constant, 0);
                GaussianFilter(paddedMatrix, r, width, height, KERNEL, sigma);
                free2DArray(paddedMatrix, width + PADDING_SIZE * 2);

            }else{ // rgb channel
                // Padding
                int** paddedMatrix_R = Pad(R, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_G = Pad(G, height, width, PADDING_SIZE, edge, 0);
                int** paddedMatrix_B = Pad(B, height, width, PADDING_SIZE, edge, 0);

                // Processed
                GaussianFilter(paddedMatrix_R, r, width, height, KERNEL, sigma);
                GaussianFilter(paddedMatrix_G, g, width, height, KERNEL, sigma);
                GaussianFilter(paddedMatrix_B, b, width, height, KERNEL, sigma);

                // Free
                free2DArray(paddedMatrix_R, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_G, width + PADDING_SIZE * 2);
                free2DArray(paddedMatrix_B, width + PADDING_SIZE * 2);
            }
            break;
    }

    cout << "OK!\n";
    cout << "-----------------------------------\n";

    // Output the result as .bmp file
    cout << "Enter the output file name: ";
    cin >> filename;

    char* output_image = (char*)filename;

    cout << "Save \"" << output_image << "\" ";
    if (Gray_flag) {
        cout << "as single channel.\n";
        save_bmp(output_image, r, r, r); // for gray images (one channel)
    }

    else {  // for true color images
        cout << "as RGB channel.\n";
        save_bmp(output_image, r, g, b);
    }

    // Close .bmp file
    close_bmp();
}
