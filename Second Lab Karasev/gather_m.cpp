#include "gather_m.h"

int MPI_Gather_m(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, int root, MPI_Comm comm)
{
	int ProcRank = 0, ProcNum = 0, reindex_Rank = 0, copy_Rank = 0, stock = 0, height_tree = 0, error = MPI_SUCCESS;
	MPI_Aint d_length, ex;
	MPI_Status stat;

	MPI_Comm_size(comm, &ProcNum);
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Type_get_extent(rtype, &ex, &d_length);

	stock = ProcNum;
	copy_Rank = ProcRank;

	if (copy_Rank % 2 == 1)
	{
		stock = 1;
	}
	else if (copy_Rank % 2 == 0 && copy_Rank != 0)
	{
		int i = 1, k = 1;
		while (ProcRank % (1 << k) == 0)
		{
			i *= 2;
			k++;
		}
		stock = i;
		if (ProcRank == ProcNum - 1)
		{
			stock = 1;
		}
		if (ProcRank == 0)
		{
			stock = ProcNum;
		}
	}
	else
	{
		stock = ProcNum;
	}


	void *job_mas = calloc(stock*scount, d_length);
	int destination_proc = ProcRank - (ProcRank % 2);

	height_tree = 1 + ceil(log2(ProcNum));

	std::cout << stock << std::endl;

	if (ProcRank != destination_proc)
	{
		error = MPI_Send(sbuf, scount, stype, destination_proc, 0, comm);
	}
	else
	{
		memcpy(job_mas, sbuf, (size_t)(scount * d_length));
	}

	if (height_tree == 1)
	{
		memcpy(rbuf, job_mas, (size_t)(scount * d_length));
		return error;
	}

	for (int i = 1; i < height_tree; ++i)
	{
		destination_proc = ProcRank - (ProcRank % (1 << (i + 1)));
		int count_provise = scount * (1 << i);

		if (ProcRank % (1 << i) == 0)
		{
			if (i != height_tree - 1)
			{
				if ((ProcRank + (1 << (i - 1))) < ProcNum)
				{
					error = MPI_Recv((int *)(job_mas)+count_provise / 2, count_provise, rtype, ProcRank + (1 << (i - 1)), 0, comm, &stat);
				}

				if (ProcRank != destination_proc)
				{
					error = MPI_Send(job_mas, count_provise, stype, destination_proc, 0, comm);
				}

			}
			else
			{
				memcpy(rbuf, job_mas, count_provise * d_length);
				error = MPI_Recv((int *)(rbuf)+count_provise / 2, count_provise, rtype, ProcRank + (1 << (i - 1)), 0, comm, &stat);			
			}
		}
	}

	MPI_Barrier(comm);

	if (ProcRank == 0)
	{
		if (root != 0)
		{
			error = MPI_Send(rbuf, scount * ProcNum, rtype, root, 0, comm);
		}
	}
	if (ProcRank == root)
	{
		if (root != 0)
		{
			error = MPI_Recv(rbuf, scount * ProcNum, rtype, 0, 0, comm, &stat);
		}
	}
	return error;
}