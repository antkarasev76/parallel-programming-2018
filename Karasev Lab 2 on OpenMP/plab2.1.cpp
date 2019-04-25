#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <locale.h>
#include <Windows.h>
#include <omp.h>
#include <random>
#include <iostream>
#include <ctime>
#include <string>

#define Min 200;

int **offset = nullptr;
std::vector<int> *tmp;

inline int get_byte(int number, int nbyte)
{
	return (number >> (nbyte * 8)) & 0xFF;
}

void show_mas(int *rand_mas, int length, int i = 0)
{
	std::cout << i << " ";
	for (int j = 0; j < length; j++)
	{
		std::cout << rand_mas[j] << " ";
	}
	std::cout << std::endl;
}

void get_sort_parts(int *rand_mas, int &length, int count_of_threads)
{
	tmp = new std::vector<int>[count_of_threads];

	for (int i = 0; i < count_of_threads; i++)
	{
		for (int j = offset[i][0]; j < offset[i][0] + offset[i][1]; j++)
		{
			tmp[i].push_back(rand_mas[j]);
		}
	}
}

void valid_of_sort(int *rand_mas, int &length, int idth)
{
	int  valid = 0;
	int i = 0;

	while (tmp[idth].size() != 0)
	{	
		int e = rand_mas[i];
		auto it = std::find(tmp[idth].begin(), tmp[idth].end(), e);

		if (it == tmp[idth].end())
		{
			valid = 1;
			break;
		}
		else
		{
			tmp[i].erase(it);
		}
	}

	offset[idth][2] = valid;
}

void counting_sort(int *rand_mas, int *temp, int &length, int &nbyte)
{
	int count[257] = { 0 };

	for (int i = 0; i < length; i++)
	{
		count[get_byte(rand_mas[i], nbyte) + 1]++;
	}

	for (int i = 1; i < 256; i++)
	{
		count[i] += count[i - 1];
	}

	for (int i = 0; i < length; i++)
	{
		temp[count[get_byte(rand_mas[i], nbyte)]++] = rand_mas[i];
	}

	for (int i = 0; i < length; i++)
	{
		rand_mas[i] = temp[i];
	}
}

void radix_sort(int *rand_mas, int &length)
{
	unsigned int *job = (unsigned int *)rand_mas;
	int count = sizeof(int);
	int *temp = new int[length];

	for (int i = 0; i < length; i++)
	{
		job[i] ^= INT_MIN;
	}

	for (int i = 0; i < count; i++)
	{
		counting_sort(rand_mas, temp, length, i);
	}

	for (int i = 0; i < length; i++)
	{
		job[i] ^= INT_MIN;
	}
	delete temp;
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

int* simple_merge(int *mas1, int size_m1, int *mas2, int size_m2)
{
	int* res = new int[size_m1 + size_m2];

	int ind1 = 0, ind2 = 0, ind_r = 0;

	while ((ind1 != size_m1) && (ind2 != size_m2))
	{
		if (mas1[ind1] <= mas2[ind2])
		{
			res[ind_r] = mas1[ind1];
			ind1++;
		}
		else
		{
			res[ind_r] = mas2[ind2];
			ind2++;
		}
		ind_r++;
	}
	if (ind1 == size_m1)
	{
		int j = ind2;
		for (; j < size_m2; j++, ind_r++)
		{
			res[ind_r] = mas2[j];
		}
	}
	else
	{
		int j = ind1;
		for (; j < size_m1; j++, ind_r++)
		{
			res[ind_r] = mas1[j];
		}
	}

	return res;
}

void transfer_pair_threads2(int *rand_mas, int thread_l, int thread_r)
{
	int idth;
    #pragma omp parallel shared(rand_mas) private(idth) num_threads(2)
	{
		idth = omp_get_thread_num();
		int ind1, ind2, elem1, elem2;
		int *res = nullptr;

		if (idth == 0)
		{
			elem1 = rand_mas[(offset[thread_l][0] + offset[thread_l][0] - offset[thread_l][0] / 2 - 1)];
			ind1 = offset[thread_l][0] + offset[thread_l][0] - offset[thread_l][0] / 2 - 1;

			ind2 = BinSearching(rand_mas, offset[thread_r][0], offset[thread_r][0] + offset[thread_r][1] - 1, elem1);

			while (rand_mas[ind2] >= elem1)
			{
				ind2--;
			}
			if (ind2 < offset[thread_r][0])
			{
				ind2 = offset[thread_r][0] - 1;
			}

			elem2 = rand_mas[ind2];

			res = simple_merge(rand_mas + offset[thread_l][0], ind1 + 1 - offset[thread_l][0], rand_mas + offset[thread_r][0], ind2 + 1 - offset[thread_r][0]);

			int length = ind1 + 1 - offset[thread_l][0] + ind2 + 1 - offset[thread_r][0], k = 0;

			for (int i = offset[thread_l][0]; i < offset[thread_l][0] + length; i++)
			{
				rand_mas[i] = res[k];
				k++;
			}
		}
		else
		{
			ind2 = offset[thread_l][1] - offset[thread_l][1] / 2 - 1 + offset[thread_l][0];
			elem2 = rand_mas[ind2];

			ind1 = BinSearching(rand_mas, offset[thread_r][0], offset[thread_r][0] + offset[thread_r][1] - 1, elem2);

			while (rand_mas[ind1] >= elem2)
			{
				ind1--;
			}
			if (ind1 < offset[thread_r][0])
			{
				ind1 = offset[thread_r][0] - 1;
			}

			elem1 = rand_mas[ind1];

			res = simple_merge(rand_mas + ind1 + 1, offset[thread_r][1] - (ind1 - offset[thread_r][0] + 1), rand_mas + ind2 + 1, offset[thread_l][1] - (ind2 + 1 - offset[thread_l][0]));

			int length = offset[thread_r][1] - (ind1 - offset[thread_r][0] + 1) + offset[thread_l][1] - (ind2 + 1 - offset[thread_l][0]);
			int k = 0;
			for (int i = (ind2 + 1 - offset[thread_l][0] + ind1 + 1 - offset[thread_r][0]); i < ind2 + 1 - offset[thread_l][0] + ind1 + 1 - offset[thread_r][0] + length; i++)
			{
				rand_mas[i] = res[k];
				k++;
			}
		}

		delete res;
	}
}

void thread_tree(int *rand_mas, int &count_of_threads, int& length)
{
	int height_tree = 1 + ceil(log2(count_of_threads)), pos = 0;

	for (int i = 1; i < height_tree; i++)
	{
		/*for (int j = 0; j < count_of_threads; j += (1 << (i - 1)))
		{
			int count = 0;

			transfer_pair_threads2(rand_mas, j, j + (1 << (i - 1)));

			pos = 0;
			for (int l = 0; l < count_of_threads; l += (1 << (i - 1)))
			{
				pos += offset[l][1];
				offset[l + (1 << (i - 1))][0] = pos;
			}
		}*/

		int count = 0;

		for (int j = 0; j < count_of_threads; j += (1 << i))
		{
			count++;
		}

		/*int *mas_rangs = new int[count];

		for (int j = 0; j < count; j += (1 << i))
		{
			offset[j][2] = j;
			offset[j + (1 << (i - 1))][2] = j + (1 << (i - 1));
		}

        #pragma omp parallel shared(rand_mas) num_threads(count)
		{
			transfer_pair_threads2(rand_mas, j, j + (1 << (i - 1)));
		}*/

		if (count > 1)
		{
			omp_set_nested(1);

            #pragma omp parallel shared(rand_mas) num_threads(count)
			{
				int idth = omp_get_thread_num(), l = idth;
				if (idth != 0)
				{
					l *= (1 << i);
				}

				int r = l + (1 << (i - 1));

				transfer_pair_threads2(rand_mas, l, r);
			}
		}
		else
		{
			transfer_pair_threads2(rand_mas, 0, 1);
		}

		//#pragma omp barrier;

		pos = 0;
		for (int l = 0; l < count_of_threads; l += (1 << (i - 1)))
		{
			pos += offset[l][1];
			offset[l + (1 << (i - 1))][0] = pos;
		}
		
	}
}

void parallel_radix_sort(int *rand_mas, int &length, int &count_of_threads)
{
	int i, idth;
    #pragma omp parallel shared(rand_mas) private(i, idth) num_threads(count_of_threads)
	{
		int idth = omp_get_thread_num();
		radix_sort(rand_mas + offset[idth][0], offset[idth][1]);

		//valid_of_sort(rand_mas, count_of_elem, idth);

        #pragma omp barrier
	}

	thread_tree(rand_mas, count_of_threads, length);
}

int main(int argc, char* argv[])
{
	int count_of_elem = 100, i, count_of_threads = 2;
	double time1 = 0, time2 = 0, time3 = 0, time4 = 0;
	bool valid;

	std::cout << "Input count of threads: ";
	std::cin >> count_of_threads;

	if (count_of_threads > omp_get_max_threads())
	{
		count_of_threads = omp_get_max_threads();
	}

	std::cout << "Input count of element in massiv: ";
	std::cin >> count_of_elem;

	int *rand_mas = new int[count_of_elem];
	offset = new int*[count_of_threads];

	for (int i = 0; i < count_of_elem; i++)
	{
		offset[i] = new int[3];
	}


	int coeff_prop = count_of_elem / count_of_threads;
	int ressidue_c = count_of_elem % count_of_threads;
	int pos = 0;

	for (int i = 0; i < count_of_threads; i++)
	{
		if (ressidue_c > 0)
		{
			offset[i][1] = coeff_prop + 1;
			--ressidue_c;
		}
		else
		{
			offset[i][1] = coeff_prop;
		}
	}

	offset[0][0] = 0;

	for (int i = 0; i < count_of_threads - 1; i++)
	{
		pos += offset[i][1];
		offset[i + 1][0] = pos;
	}


	std::default_random_engine gen(std::time(nullptr));
	std::uniform_int_distribution<> d(-100, 100);

	if (count_of_elem > 2)
	{
             #pragma omp parallel shared(rand_mas) private(i) num_threads(count_of_threads)
		     {
                #pragma omp for schedule(static)

			    for (i = 0; i < count_of_elem; ++i)
			    {
				    rand_mas[i] = d(gen);
			    }

                #pragma omp barrier
		     } 
	}
	else
	{
		for (i = 0; i < count_of_elem; ++i)
		{
			rand_mas[i] = d(gen);
		}
	}

	time1 = omp_get_wtime();

	radix_sort(rand_mas, count_of_elem);

	time2 = omp_get_wtime();

	time3 = time2 - time1;

	//get_sort_parts(rand_mas, count_of_elem, 1);

	//valid_of_sort(rand_mas, count_of_elem, 0);

	/*if (offset[0][2] == 1)
	{
		std::cout << "Massiv wasn't sorted" << std::endl;
	}*/


	if (count_of_elem > 2)
	{
        #pragma omp parallel shared(rand_mas) private(i) num_threads(count_of_threads)
		{
            #pragma omp for schedule(static)

			for (i = 0; i < count_of_elem; ++i)
			{
				rand_mas[i] = d(gen);
			}

            #pragma omp barrier
		}
	}
	else
	{
		for (i = 0; i < count_of_elem; ++i)
		{
			rand_mas[i] = d(gen);
		}
	}

	//get_sort_parts(rand_mas, count_of_elem, 1);

	time1 = omp_get_wtime();

	parallel_radix_sort(rand_mas, count_of_elem, count_of_threads);

	time2 = omp_get_wtime();

	time4 = time2 - time1;

	/*bool valid = true;

	for (int i = 0; i < count_of_threads; i++)
	{
		if (offset[i][2] == 1)
		{
			valid = false;
		}
	}*/

	/*if (valid == false)
	{
	std::cout << "Massiv wasn't sorted" << std::endl;
	}*/

	std::cout << "The operating time of serial algorithm: " << time3 << std::endl;
	std::cout << "The operating time of parallel algorithm: " << time4 << std::endl;

	system("pause");

	return 0;
}