#include "accessory.h"
#include <time.h>
#include <algorithm>
#include <vector>
#include <time.h>

char* CreateRandStr(int *size)
{
	srand(time(NULL));
	int len = rand() % 2 + 3;
	char* rand_str = new char[len];
	for (int i = 0; i < len; i++)
	{
		rand_str[i] = (char)(rand());
	}
	*size = len;
	return rand_str;
}

void uniproc_alg(char *work_str, int size, int *mas_char_count)
{
	for (int i = 0; i < size; i++)
	{
		int asc2 = static_cast<int>((unsigned char)work_str[i]);
		mas_char_count[asc2]++;
	}
}