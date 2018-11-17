#pragma once
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include "math.h"

int MPI_Gather_m(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, int root, MPI_Comm comm);