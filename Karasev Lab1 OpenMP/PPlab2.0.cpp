#include "pch.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <random>
#include <omp.h>

int *tmp = nullptr;

inline int get_byte(int number, int nbyte)
{
	return (number >> (nbyte * 8)) & 0xFF;
}

void counting_sort(int *rand_mas, int &length, int &nbyte)
{
	int count[257] = { 0 };

	for (int i = 0; i < length; i++)
	{
		count[get_byte(rand_mas[i], nbyte) + 1]++;
	}

	for (int i = 1; i < 256; i++)
	{
		count[i] += count[i-1];
	}

	for (int i = 0; i < length; i++)
	{
		tmp[count[get_byte(rand_mas[i], nbyte)]++] = rand_mas[i];
	}

	for (int i = 0; i < length; i++)
	{
		rand_mas[i] = tmp[i];
	}
}

void radix_sort(int *rand_mas, int &length)
{
	unsigned int *job = (unsigned int *)rand_mas;
	int count = sizeof(int);

	for (int i = 0; i < length; i++)
	{
		job[i] ^= INT_MIN;
	}

	for (int i = 0; i < count; i++)
	{
		counting_sort(rand_mas, length, i);
	}

	for (int i = 0; i < length; i++)
	{
		job[i] ^= INT_MIN;
	}
}

int main(int argc, char* argv[])
{
	int count = 100;
	long long time = 0;
	double time1, time2;

	std::cout << "Input count of element in massiv: ";
	std::cin >> count;

	int *rand_mas = new int[count];
	tmp = new int[count];

	std::default_random_engine gen(std::time(nullptr));
	std::uniform_int_distribution<> d(-100, 100);

	std::cout << "Random massiv: ";

	for (int i = 0; i < count; ++i)
	{
		rand_mas[i] = d(gen);
		std::cout << rand_mas[i] << " ";
	}
	std::cout << std::endl;

	time1 = omp_get_wtime();

	radix_sort(rand_mas, count);

	time2 = omp_get_wtime() - time1;

	std::cout << "Sorting massiv: ";

	for (int i = 0; i < count; ++i)
	{
		std::cout << rand_mas[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "The operating time of serial algorithm: " << time2 << std::endl;

	free(tmp);

	return 0;
}