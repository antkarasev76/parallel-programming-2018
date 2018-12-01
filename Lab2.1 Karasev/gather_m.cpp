#include "gather_m.h"

int MPI_Gather_m(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, int root, MPI_Comm comm)
{

	int ProcNum = 0, ProcRank = 0, stock = 0, height_tree = 0, error = MPI_SUCCESS, cur_count = 0;
	MPI_Aint d_length, ex;
	MPI_Status stat;

	MPI_Comm_size(comm, &ProcNum);
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Type_get_extent(rtype, &ex, &d_length);

	stock = ProcNum;

	if (ProcRank % 2 == 1)
	{
		stock = 1;
	}
	else if (ProcRank % 2 == 0 && ProcRank != 0)
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
	}
	else
	{
		stock = ProcNum;
	}

	if ((ProcRank + stock - 1) > ProcNum - 1)
	{
		stock = ProcNum - ProcRank;
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
		cur_count += 1;
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
					if ((count_provise / 2) > ((stock - cur_count) * scount))
					{
						error = MPI_Recv((int *)(job_mas)+(count_provise / 2), (stock - cur_count) * scount, rtype, ProcRank + (1 << (i - 1)), 0, comm, &stat);
						cur_count += (stock - cur_count);
					}
					else
					{
						error = MPI_Recv((int *)(job_mas)+(count_provise / 2), (count_provise / 2), rtype, ProcRank + (1 << (i - 1)), 0, comm, &stat);
						cur_count += (count_provise / 2 * scount);
					}
				}

				if (ProcRank != destination_proc)
				{
					if (count_provise > scount * stock)
					{
						error = MPI_Send(job_mas, scount * stock, stype, destination_proc, 0, comm);
					}
					else
					{
						error = MPI_Send(job_mas, count_provise, stype, destination_proc, 0, comm);
					}
				}

			}
			else
			{
				memcpy(rbuf, job_mas, cur_count * scount * d_length);
				error = MPI_Recv((int *)(rbuf) + (count_provise / 2), (stock  - cur_count) * scount, rtype, ProcRank + (1 << (i - 1)), 0, comm, &stat);
			}
		}
	}

	MPI_Barrier(comm);
	return error;
}