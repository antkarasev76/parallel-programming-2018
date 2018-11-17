// PPlab-2.cpp: определяет точку входа для консольного приложения.
//
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include "gather_m.h"

int main(int argc, char* argv[])
{
	int ProcRank, ProcNum, count = 1, root_proc = 0, error;
	long double time1, time2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status stat;
	time1 = MPI_Wtime();

	int *sbuf = new int[count]; 
	int *rbuf = new int[count * ProcNum];

	for (int i = 0; i < count; ++i)
	{
		sbuf[i] = ProcRank * 100 + i;
		std::cout << sbuf[i] << " ";
	}
	std::cout << std::endl;

	error = MPI_Gather_m(sbuf, count, MPI_INT, rbuf, count, MPI_INT, root_proc, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == root_proc)
	{
		for (int i = 0; i < ProcNum * count; ++i)
		{
			std::cout << rbuf[i] << " ";
		}
		std::cout << std::endl;
		time2 = MPI_Wtime();
		std::cout << "The operating time of parallel algorithm: " << time2 - time1 << std::endl;
		//std::cout << error << std::endl;
	}
	MPI_Finalize();
	return 0;
}
