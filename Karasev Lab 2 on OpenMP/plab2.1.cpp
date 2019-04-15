#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <locale.h>
#include <Windows.h>
#include <omp.h>
#include <random>
#include <iostream>
#include <ctime>

#define Min 200;

int *tmp = nullptr, **offset = nullptr, *vp = nullptr;

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

bool valid_of_sort(int *rand_mas, int &length)
{
	bool valid = true;

	for (int i = 0; i < length - 1; i++)
	{
		if (rand_mas[i] > rand_mas[i + 1])
		{
			valid = false;
		}
	}

	return valid;
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

int size_part_tree(int thread_num, int thread_count, int& length)
{
	int stock = thread_count, real_count = 0;

	if (thread_num % 2 == 1)
	{
		stock = 1;
	}
	else if (thread_num % 2 == 0 && thread_num != 0)
	{
		int i = 1, k = 1;
		while (thread_num % (1 << k) == 0)
		{
			i *= 2;
			k++;
		}
		stock = i;
		if (thread_num == thread_count - 1)
		{
			stock = 1;
		}
	}
	else
	{
		stock = thread_count;
	}

	if ((thread_num + stock) > thread_count)
	{
		stock = thread_count - thread_num;
	}

	for (int i = thread_num; i < thread_num + stock; ++i)
	{
		real_count += offset[i][2];
	}
	return real_count;
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

void transfer_pair_threads(int own_size, int dest_thread, int size_s, int sch, int count_of_threads, int size_elem, int *rand_mas)
{
	int idth = omp_get_thread_num(), ind1, ind2, elem1, elem2, ind_n = 0;
	int count = own_size, *res = nullptr;

	/*offset[idth][0] = count;

	int k = 0;

    #pragma omp barrier

    #pragma omp master
	{
		offset[0][1] = offse
		offset[0][0] = 0;
		
		for (int i = 1; i < count_of_threads; i++)
		{
			offset[i][1] = offset[i][0];
			offset[i][0] = offset[i - 1][0] + offset[i - 1][1];
		}
	}*/

	//count = offset[idth][0];

	if (sch == 0)
	{
		elem1 = rand_mas[(offset[idth][0] + count - count / 2 - 1)];
		ind1 = offset[idth][0] + count - count / 2 - 1;

		//size_d = (dest_thread == count_of_threads - 1) ? (size_elem - offset[dest_thread][0]) : (offset[dest_thread][0] - offset[dest_thread - 1][0] - 1);

		ind2 = BinSearching(rand_mas, offset[dest_thread][0], offset[dest_thread][0] + size_s - 1, elem1);

		while (rand_mas[ind2] >= elem1)
		{
			ind2--;
		}
		/*if (st_ind > ind2)
		{
			ind2 = st_ind;
		}*/

		elem2 = rand_mas[ind2];

		ind_n = offset[idth][0];

		res = simple_merge(rand_mas + offset[idth][0], ind1 + 1 - offset[idth][0], rand_mas + offset[dest_thread][0], ind2 + 2 - offset[dest_thread][0]);
	}
	else
	{
		//ind2 = (offset[dest_thread + 1][0] - offset[dest_thread][0]) - ((offset[dest_thread + 1][0] - offset[dest_thread][0]) / 2) - 1;
		ind2 = offset[dest_thread][2] - offset[dest_thread][2] / 2 - 1 + offset[dest_thread][2];
		elem2 = rand_mas[ind2];

		//size_d = (dest_thread == count_of_threads - 1) ? (size_elem - offset[dest_thread][0]) : (offset[dest_thread + 1][0] - offset[dest_thread][0] - 1);

		ind1 = BinSearching(rand_mas, offset[dest_thread][0], offset[dest_thread][0] + offset[dest_thread][2] - 1, elem2);

		while (rand_mas[ind1] >= elem2)
		{
			ind1--;
		}
		/*if (st_ind > ind1)
		{
			ind1 = st_ind;
		}*/

		elem1 = rand_mas[ind1];

		ind_n = offset[dest_thread][0] + own_size - (ind1 - offset[idth][0] + 1) + offset[idth][0] + offset[idth][2] - ind2 - 1;

		res = simple_merge(rand_mas + offset[idth][0] + ind1 + 1, own_size - (ind1 - offset[idth][0] + 1), rand_mas + offset[dest_thread][0] + ind2 + 1, offset[idth][0] + offset[idth][2] - ind2 - 1);
	}

    #pragma omp barrier

	if (sch == 0)
	{
		int length = ind1 + 1 - offset[idth][0] + ind2 + 2 - offset[dest_thread][0];
		int k = 0;
		for (int i = ind_n; i < ind2 + length; i++)
		{
			rand_mas[i] = res[k];
			k++;
		}
	}
	else
	{
		int length = own_size - (ind1 - offset[idth][0] + 1) + offset[idth][0] + offset[idth][2] - ind2 - 1;
		int k = 0;
		for (int i = ind_n; i < ind2 + length; i++)
		{
			rand_mas[i] = res[k];
			k++;
		}
	}


	delete res;
    #pragma omp barrier
}

void thread_tree(int *rand_mas, int &count_of_threads, int& length)
{
	int idth = omp_get_thread_num();

	int height_tree = 0, cur_count = 0, real_count = 0, size_s = 0, size_d;

	real_count = size_part_tree(idth, count_of_threads, length);

	int dest_thread = idth - (idth % 2);
	
	cur_count = offset[idth][2];

	size_s = size_part_tree(dest_thread, count_of_threads, length);

	if (idth != dest_thread)
	{
		size_d = size_part_tree(dest_thread, count_of_threads, length);
		transfer_pair_threads(cur_count, dest_thread, size_s, 1, count_of_threads, length, rand_mas);
	}

	height_tree = 1 + ceil(log2(count_of_threads));

	for (int i = 1; i < height_tree; ++i)
	{
		dest_thread = idth - (idth % (1 << (i + 1)));

		if (idth % (1 << i) == 0)
		{
			if (i != height_tree - 1)
			{
				if ((idth + (1 << (i - 1))) < count_of_threads)
				{
					size_s = size_part_tree(idth + (1 << (i - 1)), count_of_threads, length);
				    transfer_pair_threads(cur_count, idth + (1 << (i - 1)), size_s, 0, count_of_threads, length, rand_mas);
					cur_count += size_s;
				}

                #pragma omp barrier

				if (idth != dest_thread)
				{
					offset[i][2] = cur_count;
				}

                #pragma omp master
				{
					int pos = 0;
					for (int j = 0; j < count_of_threads; j += (1 << (i - 1)))
					{
						pos += offset[j][2];
						offset[j + (1 << (i - 1))][0] = pos;
					}
				}

				if (idth != dest_thread)
				{
					size_s = size_part_tree(idth + (1 << (i - 1)), count_of_threads, length);
					transfer_pair_threads(cur_count, dest_thread, size_s, 1, count_of_threads, length, rand_mas);
					cur_count += size_s;
				}

			}
			else
			{
				size_s = size_part_tree(idth + (1 << (i - 1)), count_of_threads, length);
				transfer_pair_threads(cur_count, idth + (1 << (i - 1)), size_s, 0, count_of_threads, length, rand_mas);
				cur_count += size_s;
			}
		}
	}

    #pragma omp barrier

	/*vp = rand_mas;
	rand_mas = tmp;
	tmp = vp;*/
}

void parallel_radix_sort(int *rand_mas, int &length, int &count_of_threads)
{
	int i, idth;
    #pragma omp parallel shared(rand_mas) private(i, idth) num_threads(count_of_threads)
	{
		int idth = omp_get_thread_num();
		radix_sort(rand_mas + offset[idth][2], offset[idth][0]);

        #pragma omp barrier

		/*if (idth == 0)
		{
			show_mas(rand_mas, length, 0);
		}*/

		thread_tree(rand_mas, count_of_threads, length);

        #pragma omp barrier

		/*if (idth == 0)
		{
			show_mas(rand_mas, length, 0);
		}*/
	}
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
	tmp = new int[count_of_elem];
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
			offset[i][2] = coeff_prop + 1;
			--ressidue_c;
		}
		else
		{
			offset[i][2] = coeff_prop;
		}
	}

	offset[0][0] = 0;

	for (int i = 0; i < count_of_threads - 1; i++)
	{
		pos += offset[i][2];
		offset[i + 1][0] = pos;
	}


	std::default_random_engine gen(std::time(nullptr));
	std::uniform_int_distribution<> d(-100, 100);

	//std::cout << "Random massiv: ";

	if (count_of_elem > 2)
	{
             #pragma omp parallel shared(rand_mas) private(i) num_threads(count_of_threads)
		     {
                #pragma omp for schedule(static)

			    for (i = 0; i < count_of_elem; ++i)
			    {
				    rand_mas[i] = d(gen);
					/*if (omp_get_thread_num() == 0)
					{
						std::cout << rand_mas[i] << " ";
					}*/
				    
			    }

                #pragma omp barrier
		     } 
	}
	else
	{
		for (i = 0; i < count_of_elem; ++i)
		{
			rand_mas[i] = d(gen);
			//std::cout << rand_mas[i] << " ";
		}
	}

	//show_mas(rand_mas, count_of_elem);

	//std::cout << std::endl;

	time1 = omp_get_wtime();

	radix_sort(rand_mas, count_of_elem);

	time2 = omp_get_wtime();

	//show_mas(rand_mas, count_of_elem);

	time3 = time2 - time1;

	valid = valid_of_sort(rand_mas, count_of_elem);

	if (valid == false)
	{
		std::cout << "Massiv wasn't sorted" << std::endl;
	}

	//std::cout << "The operating time of serial algorithm: " << time3 << std::endl;

	//std::cout << "Sorting massiv: ";

	/*for (int i = 0; i < count_of_elem; ++i)
	{
		std::cout << rand_mas[i] << " ";
	}
	std::cout << std::endl;*/

	if (count_of_elem > 2)
	{
        #pragma omp parallel shared(rand_mas) private(i) num_threads(count_of_threads)
		{
            #pragma omp for schedule(static)

			for (i = 0; i < count_of_elem; ++i)
			{
				rand_mas[i] = d(gen);
				//std::cout << rand_mas[i] << " ";
			}

            #pragma omp barrier
		}
	}
	else
	{
		for (i = 0; i < count_of_elem; ++i)
		{
			rand_mas[i] = d(gen);
			//std::cout << rand_mas[i] << " ";
		}
	}

	//show_mas(rand_mas, count_of_elem);

	time1 = omp_get_wtime();

	parallel_radix_sort(rand_mas, count_of_elem, count_of_threads);

	time2 = omp_get_wtime();

	time4 = time2 - time1;

	//show_mas(rand_mas, count_of_elem);


	valid = valid_of_sort(rand_mas, count_of_elem);

	if (valid == false)
	{
		std::cout << "Massiv wasn't sorted" << std::endl;
	}

	std::cout << "The operating time of serial algorithm: " << time3 << std::endl;
	std::cout << "The operating time of parallel algorithm: " << time4 << std::endl;

	//free(tmp);

	return 0;
}