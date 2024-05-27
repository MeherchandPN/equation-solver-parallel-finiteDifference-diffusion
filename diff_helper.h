#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

// Function to create a 2D array of doubles
// Parameters:
//   r: number of rows
//   c: number of columns
// Returns:
//   A pointer to a 2D array of doubles
double** array2D(int r, int c) {
    // Allocate memory for both array of pointers and 2D array in a single block
    double** A = (double**)malloc(sizeof(double*) * r + sizeof(double) * r * c);
    if (A == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Set pointers to rows
    double* data = (double*)(A + r);
    for (int i = 0; i < r; i++)
        A[i] = data + i * c;

    return A;
}

// Function to create a directory if it does not exist
void make_dir() {
    struct stat st = {0};
    if (stat("diff_out", &st) == -1) {
        if (mkdir("diff_out", 0700) == -1) {
            perror("mkdir failed");
            exit(EXIT_FAILURE);
        }
    }
}

// Function to save a 2D array of doubles to a file
// Parameters:
//   a: pointer to the 2D array
//   r: number of rows
//   c: number of columns
//   n: file index
void save_array(double** a, int r, int c, int n) {
    char loc[50];
    // Construct file path with index
    snprintf(loc, sizeof(loc), "diff_out/data.%06d", n);
    FILE* file = fopen(loc, "w");
    if (file == NULL) {
        perror("fopen failed");
        exit(EXIT_FAILURE);
    }

    // Write array data to file
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            fprintf(file, "%f ", a[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

// Function to save a 2D array of doubles to a file with time label
// Parameters:
//   a: pointer to the 2D array
//   r: number of rows
//   c: number of columns
//   t: time label
void save_array_timelabel(double** a, int r, int c, double t) {
    char loc[50];
    // Construct file path with time label
    snprintf(loc, sizeof(loc), "diff_out/data.t=%.6f", t);
    FILE* file = fopen(loc, "w");
    if (file == NULL) {
        perror("fopen failed");
        exit(EXIT_FAILURE);
    }

    // Write array data to file
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            fprintf(file, "%f ", a[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
}
