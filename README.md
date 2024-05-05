# ✨FDIP_HW✨

## Introduction

This repository contains homework assignments from the FDIP class in 2024. I've implemented various basic algorithms for Digital Image Processing (DIP) and will continue to optimize, expand, and update them.

``` bash
git clone https://github.com/seclusioner/FDIP_HW.git
```

After download, you have to create the project to add all `.h`/`.cpp` files under the same directory, Otherwise you must copy the program files to the project's directory.

If you use another IDE, then you may have to **change the settings of compiler**.

## Implementation

| Item     | Description                   |
|----------|-------------------------------|
| IDE      | CodeBlock ver 20.03           |
| Language | C / C++                       |

## What I've Implemented

### Image File Processing (.bmp)

- Open / Save
- Split 3 channels (RGB)
- Note: Image height (H) and width (W) must be the same so far.

### Padding

- Zeropadding
- EdgePadding

### Filters

- Smoothing
- Gaussian
- Binomial (LPF/HPF)
- Median
- Sobel
- Laplacian

### Algorithms

- Rotation
- Grayscale
- Negative
- Histogram Equalization
- Binarization
- Edge Sharpening
- Error Diffusion
- Canny Edge Detection

### Others

- Bilinear Interpolation
- 2D Convolution
> **Note:** For the even kernel size, there exists another method, but here is just limited to odd kernel sizes so far.
