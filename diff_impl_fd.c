#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diff_helper.h"
#define MASTER		0

double** Initialize(int rows, int cols, double* x, double* y)
{
    /* Sets up the Initial and Boundary Conditions */
    double** u = array2D(rows, cols);
    double radius;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            radius = pow(pow(x[j], 2) + pow(y[i], 2), 0.5);
            if (radius >= 0.5 && radius <= 0.8) 
                u[i][j] = 1.0;
            else 
                u[i][j] = 0.0;
        }
    }
    make_dir();

    return u;
}

int main(int argc, char* argv[])
{
    int quit = 0, rows = atoi(argv[1]), cols = atoi(argv[2]); // rows X cols discrete points

    double** u;
    int numtasks, taskid;
    MPI_Status status;

    /***** MPI Initialization *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if (rows % numtasks != 0 || cols % numtasks != 0)
    {
        if (taskid == MASTER)
            printf("Quitting...\nNumber of rows & columns tasks must be divisible by the number of processors deployed.\n");
        MPI_Abort(MPI_COMM_WORLD, quit);
        exit(0);
    }

    int rows_per_task = (int)rows / numtasks; // For Domain Decomposition
    if (numtasks == 1)
        rows_per_task--;
    // printf("rows_per_task=%d\n", rows_per_task);

    /* Setting up the Computational Domain */
    double xmin = -1, xmax = 1, ymin = -1, ymax = 1;
    double* x = linspace(xmin, xmax, cols);
    double* y = linspace(ymin, ymax, rows);
    double delx = x[1] - x[0], dely = y[1] - y[0];
    // double xyd = 0.;
    // if (taskid == MASTER)
    // for (int m = 0; m < rows; m++) { for (int n = 0; n < cols; n++) { printf("(%f,%f)", x[n], y[m]); } printf("\n"); }
    printf("Proc %d: Computation Domain successfully set.\n", taskid);
    MPI_Barrier(MPI_COMM_WORLD);

    double** u_old;
    double** u_new;
    double** resd;
    MPI_Barrier(MPI_COMM_WORLD);

    if (taskid == MASTER)
    {
        u = Initialize(rows, cols, x, y); // Initialize and set the Boundary for the entire computational domain (only stored in MASTER)
        save_array_timelabel(u, rows, cols, 0.);

        /* Vertical (row) domain decomposition scheme (send to worker) */
        int seek = 0;
        if (numtasks != 1)
            printf("Proc %d: Sending Domain Decomposition....\n", taskid);
        for (int splitter = 0; splitter < numtasks; splitter++)
        {
            seek = rows_per_task * splitter - 1;
            if (splitter == 0)
            {
                u_old = array2D(rows_per_task + 1, cols);
                u_new = array2D(rows_per_task + 1, cols);
                resd = array2D(rows_per_task + 1, cols);
                for (int i = 0; i < (rows_per_task + 1); i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        u_old[i][j] = u[i][j];
                        u_new[i][j] = u[i][j];
                        resd[i][j] = 0.;
                    }
                }
            }
            else if (splitter == numtasks - 1)
                MPI_Send(&u[seek][0], (rows_per_task + 1) * cols, MPI_DOUBLE, splitter, 75 + splitter, MPI_COMM_WORLD);
            else if (splitter != MASTER && splitter != numtasks - 1)
                MPI_Send(&u[seek][0], (rows_per_task + 2) * cols, MPI_DOUBLE, splitter, splitter + 75, MPI_COMM_WORLD);
        }
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    /* Initialize u for each process and Vertical (row) domain decomposition scheme (received from master) */
    else if (taskid == numtasks - 1)
    {
        u_old = array2D(rows_per_task + 1, cols);
        u_new = array2D(rows_per_task + 1, cols);
        resd = array2D(rows_per_task + 1, cols);
        MPI_Recv(&u_old[0][0], (rows_per_task + 1) * cols, MPI_DOUBLE, MASTER, 75 + taskid, MPI_COMM_WORLD, &status);
        for (int i = 0; i < (rows_per_task + 1); i++)
            for (int j = 0; j < cols; j++)
            {
                u_new[i][j] = u_old[i][j];
                resd[i][j] = 0.;
            }
        printf("Proc %d: Receiving Domain Decomposition....\n", taskid);
    }
    else if (taskid != MASTER && taskid != numtasks - 1)
    {
        u_old = array2D(rows_per_task + 2, cols);
        u_new = array2D(rows_per_task + 2, cols);
        resd = array2D(rows_per_task + 2, cols);
        MPI_Recv(&u_old[0][0], (rows_per_task + 2) * cols, MPI_DOUBLE, MASTER, 75 + taskid, MPI_COMM_WORLD, &status);
        for (int i = 0; i < (rows_per_task + 2); i++)
            for (int j = 0; j < cols; j++)
            {
                u_new[i][j] = u_old[i][j];
                resd[i][j] = 0.;
            }
        printf("Proc %d: Receiving Domain Decomposition....\n", taskid);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Computation starts from here */
    printf("Proc %d: Computation starts...\n", taskid);
    double uxx, uyy, bij, a = 1.0;
    double tstop = 1.0;
    double theta = 0.5, omega = 1.5;
    double dt = 1e-6 * tstop;
    double alpha = (theta * dt * a) / pow(delx, 2), beta = (theta * dt * a) / pow(dely, 2);
    int nt = 250, stoprow, startrow;
    // int collectfreq = (int) (nt/10);
    stoprow = rows_per_task + 1;
    startrow = 1;
    if (taskid == MASTER || taskid == numtasks - 1)
        stoprow = rows_per_task;
    int n;
    double t;
    char pcnt = '%';
    int count = 0, frq = 10;

    for (double t = dt; t <= tstop; t = t + dt, count++)
    {

        for (n = 1; n <= nt; n++)
        {
            for (int i = startrow; i < stoprow; i++)
            {
                for (int j = 1; j < cols - 1; j++)
                {
                    uxx = (u_old[i + 1][j] - 2 * u_old[i][j] + u_old[i - 1][j]) / pow((x[1] - x[0]), 2);
                    uyy = (u_old[i][j + 1] - 2 * u_old[i][j] + u_old[i][j - 1]) / pow((y[1] - y[0]), 2);
                    bij = u_old[i][j] + (1 - theta) * a * dt * (uxx + uyy);
                    resd[i][j] = bij - (1 + 2 * (alpha + beta)) * u_old[i][j] + alpha * (u_old[i + 1][j] + u_old[i - 1][j]) + beta * (u_old[i][j + 1] + u_old[i][j - 1]);
                    u_new[i][j] = u_old[i][j] + (omega / (1 + 2 * (alpha + beta))) * resd[i][j];
                }
            }
            /* Swap */
            for (int i = startrow; i < stoprow; i++)
            {
                for (int j = 1; j < cols - 1; j++)
                {
                    u_old[i][j] = u_old[i][j] + u_new[i][j];
                    u_new[i][j] = u_old[i][j] - u_new[i][j];
                    u_old[i][j] = u_old[i][j] - u_new[i][j];
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
            /* Communicate Domain edge values */
            if (taskid == MASTER && numtasks != 1)
            {
                MPI_Send(&u_old[rows_per_task - 1][0], cols, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD);
                MPI_Recv(&u_old[rows_per_task][0], cols, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD, &status);
            }
            else if (taskid == numtasks - 1 && numtasks != 1)
            {
                MPI_Send(&u_old[1][0], cols, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD);
                MPI_Recv(&u_old[0][0], cols, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD, &status);
            }
            else if (taskid != MASTER && taskid != numtasks - 1 && numtasks != 1)
            {
                MPI_Send(&u_old[1][0], cols, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD);
                MPI_Send(&u_old[rows_per_task][0], cols, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD);
                MPI_Recv(&u_old[0][0], cols, MPI_DOUBLE, taskid - 1, 30, MPI_COMM_WORLD, &status);
                MPI_Recv(&u_old[rows_per_task + 1][0], cols, MPI_DOUBLE, taskid + 1, 30, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        /* Collect computed result to MASTER */
        if (taskid == MASTER)
            printf("t=%f (%.2f%c done)\n", t, (t / tstop) * 100, pcnt);
        if (taskid == MASTER)
        {
            if (numtasks != 1)
            {
                for (int receiver = 1; receiver < numtasks; receiver++)
                    MPI_Recv(&u[receiver * rows_per_task][0], rows_per_task * cols, MPI_DOUBLE, receiver, 111 + receiver, MPI_COMM_WORLD, &status);
            }
            for (int i = 1; i < rows_per_task; i++)
                for (int j = 1; j < cols - 1; j++)
                    u[i][j] = u_old[i][j];
        }
        else if (numtasks != 1)
        {
            MPI_Send(&u_old[1][0], rows_per_task * cols, MPI_DOUBLE, MASTER, 111 + taskid, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (taskid == MASTER && count % frq == 0)
            save_array_timelabel(u, rows, cols, t);
    }

    MPI_Finalize();
    return 0;
}
