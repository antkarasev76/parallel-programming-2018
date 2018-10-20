#include <mpi.h>
#include <iostream>
#include "accessory.h"
#include <vector>
#include <string>
#include <time.h>

int main(int argc, char* argv[])
{
	int size_st = 0;
	char *work_str = CreateRandStr(&size_st);
	int ProcRank, ProcNum, coeff_prop;
	long double time1, time2, time3;
	int count_chars[256] = { 0 };
	int rem = 0, crem = 0, pos = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		time3 = MPI_Wtime();
		uniproc_alg(work_str, size_st, count_chars);
		time3 = MPI_Wtime() - time3;
		for (int i = 0; i < 256; i++)
		{
			count_chars[i] = 0;
		}
		time1 = MPI_Wtime();
		std::cout << "Randomly generated string: " << std::endl;
		for (int i = 0; i < size_st; i++)
		{
			std::cout << work_str[i];
		}
		std::cout << std::endl;
		std::cout << "Her length: " << size_st << std::endl;
	}

	coeff_prop = size_st / ProcNum;
	rem = size_st % ProcNum;
	crem = size_st % ProcNum;

	int *sendc = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		if (crem > 0)
		{
			sendc[i] = coeff_prop + 1;
			crem--;
		}
		else
		{
			sendc[i] = coeff_prop;
		}
	}

	int *displs = new int[ProcNum];
	displs[0] = 0;
	for (int i = 0; i < ProcNum - 1; i++)
	{
		pos += sendc[i];
		displs[i + 1] = pos;
	}

	char *recvb = new char[sendc[ProcRank]];
	for (int i = 0; i < sendc[ProcRank]; i++)
	{
		recvb[i] = 'a';
	}


	MPI_Scatterv(work_str, sendc, displs, MPI_CHAR, recvb, sendc[ProcRank], MPI_CHAR, 0, MPI_COMM_WORLD);

	int part_count_chars[256] = { 0 };

	for (int i = 0; i < sendc[ProcRank]; i++)
	{
		std::cout << "The process number " << ProcRank << "works" << std::endl;
		int asc2 = static_cast<int>((unsigned char)recvb[i]);
		part_count_chars[asc2]++;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(part_count_chars, count_chars, 256, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		time2 = MPI_Wtime();
		std::cout << "The operating time of uniprocessor algorithm:" << time3 << std::endl;
		std::cout << "The operating time of parallel algorithm: " << time2 - time1 << std::endl;
		std::cout << "Speedup T1/Tp: " << time3 / (time1 - time2) << std::endl;

		for (int i = 0; i < 256; i++)
		{
			std::cout << count_chars[i] << std::endl;
		}
	}

	MPI_Finalize();

	char t = getchar();
	return 0;
}
