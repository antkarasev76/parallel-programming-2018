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
	int ProcRank, ProcNum, count = 1, root_proc = 1, error;
	long double time1, time2, time3, time4;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int *sbuf = new int[count];
	int *rbuf = new int[count * ProcNum], *crbuf = new int[count* ProcNum];

	for (int i = 0; i < count; ++i)
	{
		sbuf[i] = ProcRank * 100 + i;
		std::cout << sbuf[i] << " ";
	}
	std::cout << std::endl;

	MPI_Barrier(MPI_COMM_WORLD);

	time1 = MPI_Wtime();

	MPI_Gather(sbuf, count, MPI_INT, crbuf, count, MPI_INT, root_proc, MPI_COMM_WORLD);

	time2 = MPI_Wtime();

	MPI_Group mpi_group_world, reindex_gr;
	MPI_Comm rcomm;

	int *proc_ranks = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		proc_ranks[i] = (i + root_proc) % ProcNum;
	}

	MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
	MPI_Group_incl(mpi_group_world, ProcNum, proc_ranks, &reindex_gr);
	MPI_Comm_create(MPI_COMM_WORLD, reindex_gr, &rcomm);

	MPI_Barrier(MPI_COMM_WORLD);

	time3 = MPI_Wtime();

	error = MPI_Gather_m(sbuf, count, MPI_INT, rbuf, count, MPI_INT, root_proc, rcomm);

	time4 = MPI_Wtime();

	MPI_Comm_free(&rcomm);
	MPI_Group_free(&reindex_gr);
	MPI_Group_free(&mpi_group_world);

	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == root_proc)
	{
		for (int i = 0; i < ProcNum * count; i++)
		{
			std::cout << crbuf[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < ProcNum * count; ++i)
		{
			std::cout << rbuf[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "The operating time of library Gather: " << time2 - time1 << std::endl;
		std::cout << "The operating time of embodied Gather: " << time4 - time3 << std::endl;
		std::cout << "Asceleration: " << (time2 - time1)/(time4 - time3) << std::endl;
	}
	//std::cout << error << std::endl;
	MPI_Finalize();
	return 0;
}
