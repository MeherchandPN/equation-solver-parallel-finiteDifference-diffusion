#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diff_helper.h"

#define MASTER	0

// Function to initialize the computational domain and boundary conditions
double** Initialize(int r, int c, double* x, double* y) {
    double** u = array2D(r, c); // Allocate memory for the 2D array

    // Set initial and boundary conditions
    double radius;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            radius = sqrt(x[j] * x[j] + y[i] * y[i]);
            u[i][j] = (radius >= 0.5 && radius <= 0.8) ? 1.0 : 0.0;
        }
    }

    make_dir(); // Create directory for output data

    return u;
}

int main(int argc, char* argv[]) {
    int quit = 0, r = atoi(argv[1]), c = atoi(argv[2]); // Dimensions of the grid
    double** u;
    int numtasks, taskid;
    MPI_Status status;

    /***** MPI Initialization *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if (r % numtasks != 0 || c % numtasks != 0) {
        if (taskid == MASTER)
            printf("Quitting...\nNumber of rows & columns tasks must be divisible by the number of processors deployed.\n");
        MPI_Abort(MPI_COMM_WORLD, quit);
        exit(0);
    }

    int rend = r / numtasks; // Rows per task for domain decomposition

    /*Setting up the Computational Domain*/
    double xmin = -1, xmax = 1, ymin = -1, ymax = 1;
    double* x = linspace(xmin, xmax, c);
    double* y = linspace(ymin, ymax, r);

    double** u_old;
    double** u_new;
    if (taskid == MASTER) {
        u = Initialize(r, c, x, y); // Initialize the domain in MASTER process
        save_array(u, r, c, 0); // Save initial state

        u_old = array2D(rend + 1, c); // Allocate memory for receiving buffer
        u_new = array2D(rend + 1, c); // Allocate memory for computation
        // Send initial domain decomposition to other processes
        for (int splitter = 1; splitter < numtasks; splitter++) {
            int seek = rend * splitter - 1;
            MPI_Send(&u[seek][0], (rend + 1) * c, MPI_DOUBLE, splitter, 75 + splitter, MPI_COMM_WORLD);
        }
    } else {
        if (taskid == numtasks - 1)
            u_old = array2D(rend + 1, c); // Allocate memory for receiving buffer
        else
            u_old = array2D(rend + 2, c); // Allocate memory for receiving buffer
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Proc %d: Computation Domain successfully set.\n", taskid);
    MPI_Barrier(MPI_COMM_WORLD);

    // Computation starts here
    double uxx, uyy, a = 1.0;
    double tstop = 1.0;
    double delx = x[1] - x[0], dely = y[1] - y[0];
    double dt = 0.5 * (1 / (2 * a)) * (pow(delx, 2) * pow(dely, 2)) / (pow(delx, 2) + pow(dely, 2));
    int nt = (int)(1 / dt), stoprow, startrow;
    int collectfreq = (int)(nt / 10);
    stoprow = rend + 1;
    startrow = 1;
    if (taskid == MASTER || taskid == numtasks - 1)
        stoprow = rend;
    int n;
    double t;
    char pcnt = '%';
    double disp;

    for (n = 1, t = 0.0; n <= nt; n++, t = t + dt) {
        for (int i = startrow; i < stoprow; i++) {
            for (int j = 1; j < c - 1; j++) {
                disp = u_old[i][j];
                uxx = (u_old[i + 1][j] - 2 * u_old[i][j] + u_old[i - 1][j]) / pow((x[1] - x[0]), 2);
                uyy = (u_old[i][j + 1] - 2 * u_old[i][j] + u_old[i][j - 1]) / pow((y[1] - y[0]), 2);
                disp = u_old[i][j + 1];
                u_new[i][j] = dt * a * (uxx + uyy) + u_old[i][j];
                disp = u_old[i + 1][j];
            }
        }

        /*Swap*/
        for (int i = startrow; i < stoprow; i++) {
            for (int j = 1; j < c - 1; j++) {
                disp = u_old[i][j];
                u_old[i][j] = u_old[i][j] + u_new[i][j];
                u_new[i][j] = u_old[i][j] - u_new[i][j];
                u_old[i][j] = u_old[i][j] - u_new[i][j];
                disp = u_old[i][j];
                if (disp > 1)
                    printf("Error! Check Convergence!");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /*Communicate Domain edge values*/
        if (taskid == MASTER && numtasks != 1) {
            MPI_Send(&u_old[rend - 1][0], c, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD);
            MPI_Recv(&u_old[rend][0], c, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD, &status);
        } else if (taskid == numtasks - 1 && numtasks != 1) {
            MPI_Send(&u_old[1][0], c, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD);
            MPI_Recv(&u_old[0][0], c, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD, &status);
        } else if (taskid != MASTER && taskid != numtasks - 1 && numtasks != 1) {
            MPI_Send(&u_old[1][0], c, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD);
            MPI_Send(&u_old[rend][0], c, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD);
            MPI_Recv(&u_old[0][0], c, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD, &status);
            MPI_Recv(&u_old[rend + 1][0], c, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /*Collect computed result to MASTER*/
        if (n % collectfreq == 0 || n == nt) {
            if (taskid == MASTER)
                printf("t=%f (%.2f%c done)\n", t, (t / tstop) * 100, pcnt);
            if (taskid == MASTER) {
                if (numtasks != 1) {
                    for (int receiver = 1; receiver < numtasks; receiver++)
                        MPI_Recv(&u[receiver * rend][0], rend * c, MPI_DOUBLE, receiver, 111 + receiver, MPI_COMM_WORLD, &status);
                }
                for (int i = 1; i < rend; i++) {
                    for (int j = 1; j < c - 1; j++)
                        u[i][j] = u_old[i][j];
                }
            } else if (numtasks != 1) {
                MPI_Send(&u_old[1][0], rend * c, MPI_DOUBLE, MASTER, 111 + taskid, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (taskid == MASTER)
                save_array(u, r, c, n);
        }
    }

    MPI_Finalize();
    return 0;
}
