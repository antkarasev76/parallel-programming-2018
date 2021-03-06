// PPlab-3.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include <mpi.h>
#include <cstdlib>

int compare_int(const void* a, const void* b)
{
	const int *a_c = (const int*)a;
	const int *b_c = (const int*)b;


	if (*a_c < *b_c)
	{
		return -1;
	}
	else if (*a_c > *b_c)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int size_part_tree(int ProcRank, int ProcNum, int *sendc)
{
	int stock = ProcNum; 
	int real_count = 0;

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

	for (int i = ProcRank; i < ProcRank + stock; ++i)
	{
		real_count += sendc[i];

	}
	return real_count;
}

int BinSearching(int *mas, int l, int r, int x) 
{
	if (l == r)
	{
		return l;
	}
	if (l + 1 == r)
	{
		if (x < mas[l])
		{
			return l;
		}
		else
		{
			return r;
		}
	}
	int m = (l + r) / 2; 
	if (x < mas[m])
	{
		r = m;
	} 
	else if (x > mas[m])
	{
		l = m;
	} 
	else
	{
		return m;
	} 
	return BinSearching(mas, l, r, x);
}

int* simple_merge(int* mas1, int size_m1, int* mas2, int size_m2)
{
	int *tmp = new int[size_m1 + size_m2];
	int ind1 = 0, ind2 = 0, ind_r = 0;

	while ((ind1 != size_m1) && (ind2 != size_m2))
	{ 
		if (mas1[ind1] <= mas2[ind2])
		{
			tmp[ind_r] = mas1[ind1];
			ind1++; 
		} 
		else 
		{ 
			tmp[ind_r] = mas2[ind2]; 
			ind2++; 
		} 
		ind_r++; 
	} 
	if (ind1 == size_m1) 
	{
		int j = ind2; 
		for (; j < size_m2; j++, ind_r++)
		{
			tmp[ind_r] = mas2[j];
		}
	}
	else
	{
		int j = ind1;
		for (; j < size_m1; j++, ind_r++)
		{
			tmp[ind_r] = mas1[j];
		}
	}

	if (size_m2 > 0)
	{
		delete[] mas2;
	}

	return tmp;
}



//int ind(int *mas, int l, int r) {
//	int l_old = l, r_old = r;
//	int r_el = mas[l];
//	while (l < r)
//	{
//		while ((mas[r] >= r_el) && (l < r))
//		{
//			r--;
//		}
//		if (l != r)
//		{
//			mas[l] = mas[r];
//			l++; 
//		}
//		while ((mas[l] <= r_el) && (l < r))
//		{
//			l++;
//		}
//		if (l != r)
//		{
//			mas[r] = mas[l];
//			r--;
//		}
//	}
//	mas[l] = r_el;
//	r_el = l;
//	l = l_old;
//	r = r_old;
//
//	return r_el;
//}
//
//void quick_sort(int* mas, int l, int r) 
//{
//	int index = ind(mas, l, r);
//	std::cout << index << std::endl;
//	if (l < index)
//		quick_sort(mas, l, index - 1);
//	if (index < r)
//		quick_sort(mas, index + 1, r);
//}

void transfer_d(MPI_Datatype type, int sch, int d_proc, int *job_mas, int count, MPI_Aint d_length, int real_count)
{
	MPI_Status stat;
	int ProcNum = 0, ProcRank = 0;
	int *sort_mas, *s_elem_size, *t_elem_size, ind = 0;
	//int *transfer_mas;
	int size_d_proc = 0, size_t = 0;

	s_elem_size = new int[2];
	t_elem_size = new int[2];

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	//std::cout << real_count << std::endl;

	if (sch == 0)
	{
		ind = count - count / 2 - 1;
		t_elem_size[0] = job_mas[ind];
		t_elem_size[1] = count - ind;

		MPI_Send(t_elem_size, 2, MPI_INT, d_proc, 0, MPI_COMM_WORLD);
		MPI_Recv(s_elem_size, 2, MPI_INT, d_proc, 0, MPI_COMM_WORLD, &stat);

		//transfer_mas = new int[s_elem_size[1]];
		sort_mas = new int[s_elem_size[1]];

		if (s_elem_size[1] > 0)
		{
			MPI_Recv(sort_mas, s_elem_size[1], MPI_INT, d_proc, 0, MPI_COMM_WORLD, &stat);
		}
		else
		{
			delete[] sort_mas;
			sort_mas = nullptr;
		}

		if (t_elem_size[1] > 0)
		{
			MPI_Send(job_mas + ind, t_elem_size[1], MPI_INT, d_proc, 0, MPI_COMM_WORLD);
		}

		/*for (int i = 0; i < s_elem_size[1]; i++)
		{
			std::cout << transfer_mas[i] << " ";
		}
		std::cout << std::endl;

		for (int i = 0; i < count - t_elem_size[1]; i++)
		{
			std::cout << job_mas[i] << " ";
		}
		std::cout << std::endl;*/

		/*sort_mas = simple_merge(job_mas, count - t_elem_size[1], sort_mas, s_elem_size[1]);

		for (int i = 0; i < count - t_elem_size[1] + s_elem_size[1]; i++)
		{
			std::cout << sort_mas[i] << " ";
		}
		std::cout << std::endl;*/

		/*for (int i = 0; i < count - t_elem_size[1] + s_elem_size[1]; i++)
		{
			job_mas[i] = sort_mas[i];
		}*/

		/*for (int i = 0; i < count - t_elem_size[1] + s_elem_size[1]; i++)
		{
			std::cout << job_mas[i] << " ";
		}
		std::cout << std::endl;*/

		MPI_Recv(&size_d_proc, 1, MPI_INT, d_proc, 0, MPI_COMM_WORLD, &stat);

		size_t = count + size_d_proc - (count - t_elem_size[1] + s_elem_size[1]);

		if (size_t > 0)
		{
			MPI_Recv(job_mas + count - t_elem_size[1] + s_elem_size[1], size_t, type, d_proc, 0, MPI_COMM_WORLD, &stat);
		}
		/*if (d_proc == 1)
		{
			for (int i = 0; i < size_t; i++)
			{
				std::cout << job_mas[i] << std::endl;
			}
		}*/
	}
	else
	{
		MPI_Recv(s_elem_size, 2, MPI_INT, d_proc, 0, MPI_COMM_WORLD, &stat);

		ind = BinSearching(job_mas, 0, count - 1, s_elem_size[0]);

		while (job_mas[ind] >= s_elem_size[0])
		{
			ind--;
		}

		t_elem_size[0] = job_mas[ind];
		t_elem_size[1] = ind + 1;

		MPI_Send(t_elem_size, 2, MPI_INT, d_proc, 0, MPI_COMM_WORLD);

		if (t_elem_size[1] > 0)
		{
			MPI_Send(job_mas, t_elem_size[1], MPI_INT, d_proc, 0, MPI_COMM_WORLD);
		}

		sort_mas = new int[s_elem_size[1]];

		if (s_elem_size[1] > 0)
		{
			MPI_Recv(sort_mas, s_elem_size[1], MPI_INT, d_proc, 0, MPI_COMM_WORLD, &stat);
		}
		else
		{
			delete[] sort_mas;
			sort_mas = nullptr;
		}

		/*for (int i = 0; i < s_elem_size[1]; i++)
		{
			std::cout << transfer_mas[i] << " ";
		}
		std::cout << std::endl;

		for (int i = t_elem_size[1]; i < count; i++)
		{
			std::cout << job_mas[i] << " ";
		}
		std::cout << std::endl;*/

		sort_mas = simple_merge((job_mas + t_elem_size[1]), count - t_elem_size[1], sort_mas, s_elem_size[1]);

		/*for (int i = 0; i < count - t_elem_size[1] + s_elem_size[1]; i++)
		{
			std::cout << sort_mas[i] << " ";
		}
		std::cout << std::endl;*/

		MPI_Send(&count, 1, MPI_INT, d_proc, 0, MPI_COMM_WORLD);

		if (count - t_elem_size[1] + s_elem_size[1] > 0)
		{
			MPI_Send(sort_mas, count - t_elem_size[1] + s_elem_size[1], MPI_INT, d_proc, 0, MPI_COMM_WORLD);
		}
	}

	//delete[] transfer_mas;
	delete[] sort_mas;
	delete[] s_elem_size;
	delete[] t_elem_size;
	sort_mas = s_elem_size = t_elem_size = nullptr;
}

void tree_g(int *sendc, int *recvb, MPI_Datatype type, int *rand_mas)
{
	int ProcNum = 0, ProcRank = 0, height_tree = 0, cur_count = 0, real_count = 0;
	MPI_Aint d_length, ex;
	MPI_Status stat;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Type_get_extent(type, &ex, &d_length);

	real_count = size_part_tree(ProcRank, ProcNum, sendc);

	int *job_mas = new int[real_count];

	int destination_proc = ProcRank - (ProcRank % 2);

	height_tree = 1 + ceil(log2(ProcNum));

	if (ProcRank != destination_proc)
	{
		cur_count += sendc[ProcRank];
		if (ProcNum != 1)
		{
			memcpy(job_mas, recvb, (size_t)(sendc[ProcRank] *d_length));
			transfer_d(MPI_INT, 1, destination_proc, job_mas, sendc[ProcRank], d_length, real_count);
		}
	}
	else
	{
		memcpy(job_mas, recvb, (size_t)(sendc[ProcRank] * d_length));
		cur_count += sendc[ProcRank];
	}

	if (height_tree == 1)
	{
		rand_mas = job_mas;
		return;
	}

	for (int i = 1; i < height_tree; ++i)
	{
		destination_proc = ProcRank - (ProcRank % (1 << (i + 1)));

		if (ProcRank % (1 << i) == 0)
		{
			if (i != height_tree - 1)
			{
				if ((ProcRank + (1 << (i - 1))) < ProcNum)
				{
					transfer_d(MPI_INT, 0, ProcRank + (1 << (i - 1)), job_mas, cur_count, d_length, real_count);
					int size_s = size_part_tree(ProcRank + (1 << (i - 1)), ProcNum, sendc);
					cur_count += size_s;
				}

				//MPI_Barrier(MPI_COMM_WORLD);

				if (ProcRank != destination_proc)
				{
					transfer_d(MPI_INT, 1, destination_proc, job_mas, cur_count, d_length, real_count);
				}

			}
			else
			{
				transfer_d(MPI_INT, 0, ProcRank + (1 << (i - 1)), job_mas, cur_count, d_length, real_count);
				int size_s = size_part_tree(ProcRank + (1 << (i - 1)), ProcNum, sendc);
				cur_count += size_s;
				std::cout << cur_count;
				/*for (int i = 0; i < cur_count; i++)
				{
					std::cout << job_mas[i] << " ";
				}
				std::cout << std::endl;*/
				if (ProcRank == 0)
				{
					rand_mas = job_mas;
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


int* paral_sort(int *sendc, int *recvb, int *rand_mas)
{
	int ProcRank = 0, ProcNum = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	tree_g(sendc, recvb, MPI_INT, rand_mas);

	int *res = rand_mas;
	return res;
}

int main(int argc, char* argv[])
{
	int *rand_mas = nullptr, ProcRank, ProcNum, coeff_prop, ressidue_c, count = 100, pos = 0;
	long double time1, time2, time3, time4;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	time1 = MPI_Wtime();

	if (ProcRank == 0)
	{
		rand_mas = new int[count];

		for (int i = 0; i < count; ++i)
		{
			rand_mas[i] = 1 + rand() % 10;
		}
		/*for (int i = 0; i < count; i++)
		{
			std::cout << rand_mas[i] << " ";
		}
		std::cout << std::endl;*/
	}

	MPI_Barrier(MPI_COMM_WORLD);

	coeff_prop = count / ProcNum;
	ressidue_c = count % ProcNum;

	int *sendc = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		if (ressidue_c > 0)
		{
			sendc[i] = coeff_prop + 1;
			--ressidue_c;
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

	int *recvb = new int[sendc[ProcRank]];
	for (int i = 0; i < sendc[ProcRank]; i++)
	{
		recvb[i] = 0;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Scatterv(rand_mas, sendc, displs, MPI_INT, recvb, sendc[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

	if (rand_mas != nullptr)
	{
		delete[] rand_mas;
		rand_mas = nullptr;
	}

	qsort(recvb, sendc[ProcRank], sizeof(int), compare_int);
	//quick_sort(recvb, 0, sendc[ProcRank] - 1);

	/*if (ProcRank == 1)
	{
		for (int i = 0; i < sendc[ProcRank]; i++)
		{
			std::cout << recvb[i] << " ";
		}
		std::cout << std::endl;
	}*/

	MPI_Barrier(MPI_COMM_WORLD);

	int* res = nullptr;

	res = paral_sort(sendc, recvb, res);

	if (ProcRank == 0)
	{
		if (res != nullptr)
		{
			for (int i = 0; i < count; i++)
			{
			std::cout << res[i] << " ";
			}
			std::cout << std::endl;
		}
		time2 = MPI_Wtime();
		std::cout << "The operating time of parallel algorithm: " << time2 - time1 << std::endl;
	}
	MPI_Finalize();
	return 0;
}
