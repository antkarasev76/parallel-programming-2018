#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <random>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"
#include "tbb/parallel_sort.h"

using namespace tbb;
using namespace std;

inline int GetByte(int number_, int nbyte_)
{
	return (number_ >> (nbyte_ * 8)) & 0xFF;
}

void SequentialCountingSort(int *rand_mas_, int *temp_, int &length_, int &nbyte_)
{
	int count[257] = { 0 };

	for (int i = 0; i < length_; i++)
	{
		count[GetByte(rand_mas_[i], nbyte_) + 1]++;
	}

	for (int i = 1; i < 256; i++)
	{
		count[i] += count[i - 1];
	}

	for (int i = 0; i < length_; i++)
	{
		temp_[count[GetByte(rand_mas_[i], nbyte_)]++] = rand_mas_[i];
	}

	for (int i = 0; i < length_; i++)
	{
		rand_mas_[i] = temp_[i];
	}
}

void SequentialRadixSort(int *rand_mas_, int &length_)
{
	unsigned int *job = (unsigned int *)rand_mas_;
	int count = sizeof(int);
	int *temp = new int[length_];

	for (int i = 0; i < length_; i++)
	{
		job[i] ^= INT_MIN;
	}

	for (int i = 0; i < count; i++)
	{
		SequentialCountingSort(rand_mas_, temp, length_, i);
	}

	for (int i = 0; i < length_; i++)
	{
		job[i] ^= INT_MIN;
	}
	delete[] temp;
}

int* GetCopyNotSortMas(int* rand_mas_, int size_)
{
	int *copy_ = new int[size_];

	for (int i = 0; i < size_; i++)
	{
		copy_[i] = rand_mas_[i];
	}

	return copy_;
}

void GetCopySortMas(int* rand_mas_, int size_, int count_of_threads_)
{
	task_scheduler_init init(count_of_threads_);

	tbb::parallel_sort(rand_mas_, rand_mas_ + size_);

	init.terminate();
}

bool IsMasCorSorted(int *rand_mas_, int *copy_rand_mas_, int size_)
{
	bool valid = true;

	for (int i = 0; i < size_; i++)
	{
		if (rand_mas_[i] != copy_rand_mas_[i])
		{
			valid = false;
		}
	}

	return valid;
}

void ShowMas(int *rand_mas_, int length_)
{
	for (int i = 0; i < length_; i++)
	{
		std::cout << rand_mas_[i] << " ";
	}
	std::cout << std::endl;
}

class Merger:public task
{
private:
	int *masL_, *masR_, *res_;
	int sizeL_, sizeR_;
public:
	Merger(int *masL, int *masR, int *res, int sizeL, int sizeR) :
		masL_(masL), masR_(masR), res_(res), sizeL_(sizeL), sizeR_(sizeR) {}

	task* execute()
	{
		int ind1 = 0, ind2 = 0, ind_r = 0;

		while ((ind1 != sizeL_) && (ind2 != sizeR_))
		{
			if (masL_[ind1] <= masR_[ind2])
			{
				res_[ind_r] = masL_[ind1];
				ind1++;
			}
			else
			{
				res_[ind_r] = masR_[ind2];
				ind2++;
			}
			ind_r++;
		}
		if (ind1 == sizeL_)
		{
			int j = ind2;
			for (; j < sizeR_; j++, ind_r++)
			{
				res_[ind_r] = masR_[j];
			}
		}
		else
		{
			int j = ind1;
			for (; j < sizeL_; j++, ind_r++)
			{
				res_[ind_r] = masL_[j];
			}
		}

		return NULL;
	}
};

class ParallelRadixSort:public task
{
private:
	int *mas_, *res_;
	int size_, part_;
	int count_of_threads_;

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

public:
	ParallelRadixSort(int *mas, int *res, int size, int part, int count_of_threads) :
		mas_(mas), res_(res), size_(size), part_(part), count_of_threads_(count_of_threads) {}

	task* execute()
	{
		if (size_ <= part_)
		{
			SequentialRadixSort(mas_, size_);
		}
		else
		{
			ParallelRadixSort &s1 = *new (allocate_child()) 
				ParallelRadixSort(mas_, res_, size_ / 2, part_, count_of_threads_ / 2);

			ParallelRadixSort &s2 = *new (allocate_child()) 
				ParallelRadixSort(mas_ + size_ / 2, res_ + size_ / 2, size_ - size_ / 2, part_, count_of_threads_ / 2); 
			
			set_ref_count(3); 
			
			spawn(s1); 
			
			spawn_and_wait_for_all(s2); 
			
			Merger **m = new Merger*[count_of_threads_ - 1]; 
			
			int tmp_size_ = size_ / 2; 
			tmp_size_ = tmp_size_ / count_of_threads_;

			int l1 = 0, r1 = tmp_size_, l2 = 0, r2 = 0; 
			
			for (int i = 0; i < count_of_threads_ - 1; i++)
			{
				int x = mas_[r1]; 
				r2 = BinSearching(mas_ + size_ / 2, 0, size_ - size_ / 2, x);

				m[i] = new (allocate_child())
					Merger(mas_ + l1, mas_ + size_ / 2 + l2, res_ + l1 + l2, r1 - l1, r2 - l2);

				l1 += tmp_size_; 
				r1 += tmp_size_; 
				l2 = r2;
			}

			Merger &spl = *new (allocate_child()) 
				Merger(mas_ + l1, mas_ + size_ / 2 + l2, res_ + l1 + l2, size_ / 2 - l1, size_ - size_ / 2 - l2); 
			
			set_ref_count(count_of_threads_ + 1);

			for (int i = 0; i < count_of_threads_ - 1; i++)
			{
				spawn(*(m[i]));
			}

			spawn_and_wait_for_all(spl);
			
			for (int i = 0; i < size_; i++)
			{
				mas_[i] = res_[i];
			} 
			
			delete[] m;
		} 
		
		return NULL;
	}
};

void RunParallelRadixSort(int *source_, int size_, int count_of_threads_) 
{
	int *res_ = new int[size_]; 
	int part_ = size_ / count_of_threads_; 
	
	if (size_ % count_of_threads_ != 0)
	{
		part_++;
	}	

	ParallelRadixSort& s = *new (task::allocate_root()) 
		ParallelRadixSort(source_, res_, size_, part_, count_of_threads_); 
	
	task::spawn_root_and_wait(s); 
	
	delete[] res_;
}

int main(int argc, char **argv)
{
	int count_of_elem = 100, i, count_of_threads = 2;
	tick_count time1, time2;
	double time3 = 0, time4 = 0;
	bool valid;

	std::cout << "Input count of threads: ";
	std::cin >> count_of_threads;

	std::cout << "Input count of element in massiv: ";
	std::cin >> count_of_elem;

	int *rand_mas = new int[count_of_elem], *copy_rand_mas;

	std::default_random_engine gen(std::time(nullptr));
	std::uniform_int_distribution<> d(-100, 100);

	for (i = 0; i < count_of_elem; ++i)
	{
		rand_mas[i] = d(gen);
	}

	//ShowMas(rand_mas, count_of_elem);

	copy_rand_mas = GetCopyNotSortMas(rand_mas, count_of_elem);

	GetCopySortMas(copy_rand_mas, count_of_elem, count_of_threads);

	time1 = tick_count::now();

	SequentialRadixSort(rand_mas, count_of_elem);

	time2 = tick_count::now();

	time3 = (time2 - time1).seconds();

	//ShowMas(rand_mas, count_of_elem);

	valid = IsMasCorSorted(rand_mas, copy_rand_mas, count_of_elem);

	if (valid != true)
	{
		std::cout << "Massiv  not sorted correctly";
	}
	else
	{
		std::cout << "Massiv  sorted correctly";
	}
	std::cout << std::endl;

	for (i = 0; i < count_of_elem; ++i)
	{
		rand_mas[i] = d(gen);
	}

	delete[] copy_rand_mas;

	copy_rand_mas = GetCopyNotSortMas(rand_mas, count_of_elem);

	GetCopySortMas(copy_rand_mas, count_of_elem, count_of_threads);

	//ShowMas(rand_mas, count_of_elem);

	time1 = tick_count::now();

	RunParallelRadixSort(rand_mas, count_of_elem, count_of_threads);

	time2 = tick_count::now();

	time4 = (time2 - time1).seconds();

	//ShowMas(rand_mas, count_of_elem);

	valid = IsMasCorSorted(rand_mas, copy_rand_mas, count_of_elem);

	if (valid != true)
	{
		std::cout << "Massiv not sorted correctly";
	}
	else
	{
		std::cout << "Massiv  sorted correctly";
	}
	std::cout << std::endl;

	//std::cout << "The operating time of serial algorithm: " << time3 << std::endl;
	std::cout << "The operating time of parallel TBB algorithm: " << time4 << std::endl;

	system("pause");

	return 0;
}
