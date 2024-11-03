#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define NX 4  // Size in the x-direction
#define NY 4  // Size in the y-direction

// Function to generate some sample data (e.g., a sinusoidal wave in 2D)
void generate_real_space_data(double **data, int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            data[i][j] = sin(2 * M_PI * i / nx) * cos(2 * M_PI * j / ny); // Example: 2D sine-cosine wave
        }
    }
}

// Function to allocate a 2D array using malloc for fftw_complex
fftw_complex** allocate_2d_complex_array(int nx, int ny) {
    fftw_complex**  array = malloc(nx * sizeof(fftw_complex*));
    for (int i = 0; i < nx; i++) {
        array[i] = malloc(ny * sizeof(fftw_complex));
    }
    return array;
}

// Function to allocate a 2D array using malloc for double
double** allocate_2d_array(int nx, int ny) {
    double **array = (double**) malloc(nx * sizeof(double*));
    for (int i = 0; i < nx; i++) {
        array[i] = (double*) malloc(ny * sizeof(double));
    }
    return array;
}

// Function to free a 2D array of doubles
void free_2d_array(double **array, int nx) {
    for (int i = 0; i < nx; i++) {
        free(array[i]);
    }
    free(array);
}

// Function to free a 2D array of fftw_complex
void free_2d_complex_array(fftw_complex **array, int nx) {
    for (int i = 0; i < nx; i++) {
        free(array[i]);
    }
    free(array);
}

int main() {
    // Allocate 2D array for real space data
    double **real_data = allocate_2d_array(NX, NY);
    
    // Allocate 2D array for k-space (Fourier space) data using malloc
    fftw_complex **k_space_data = allocate_2d_complex_array(NX, NY / 2 + 1);

    // Create plan for forward FFT (real to complex)
    fftw_plan plan_forward = fftw_plan_dft_r2c_2d(NX, NY, &real_data[0][0], &k_space_data[0][0], FFTW_MEASURE);

    // Generate 2D real space data
    generate_real_space_data(real_data, NX, NY);

    // Print the real space data
    printf("Real space data:\n");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            printf("%6.3f ", real_data[i][j]);
        }
        printf("\n");
    }

    // Execute the forward FFT
    fftw_execute(plan_forward);

    // Print the k-space (Fourier space) data
    printf("\nK-space data:\n");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY / 2 + 1; j++) {
            printf("k[%d][%d] = %6.3f + %6.3fi\n", i, j, k_space_data[i][j][0], k_space_data[i][j][1]);
        }
    }

    // Free allocated memory
    fftw_destroy_plan(plan_forward);
    free_2d_complex_array(k_space_data, NX);
    free_2d_array(real_data, NX);

    return 0;
}

