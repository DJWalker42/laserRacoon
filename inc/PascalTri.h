#ifndef PASCALTRI_HPP
#define PASCALTRI_HPP

#include <vector>

namespace{
	typedef struct tPTRI{
		int row;
		int* pos;
	} PTRI;
}

namespace Pascal{

	const size_t n_max = 21;

	int num0[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	int num1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	int num2[] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190};
	int num3[] = {1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560, 680, 816, 969, 1140};
	int num4[] = {1, 5, 15, 35, 70, 126, 210, 330, 495, 715, 1001, 1365, 1820, 2380, 3060, 3876, 4845};
	int num5[] = {1, 6, 21, 56, 126, 252, 462, 792, 1287, 2002, 3003, 4368, 6188, 8568, 11628, 15504};
	int num6[] = {1, 7, 28, 84, 210, 462, 924, 1716, 3003, 5005, 8008, 12376, 18564, 27132, 38760};
	int num7[] = {1, 8, 36, 120, 330, 792, 1716, 3432, 6435, 11440, 19448, 31824, 50388, 77520};
	int num8[] = {1, 9, 45, 165, 495, 1287, 3003, 6435, 12870, 24310, 43758, 75582, 125970};
	int num9[] = {1, 10, 55, 220, 715, 2002, 5005, 11440, 24310, 48620, 92378, 167960};
	int num10[] = {1, 11, 66, 286, 1001, 3003, 8008, 19448, 43758, 92378, 184756};
	int num11[] = {1, 12, 78, 364, 1365, 4368, 12376, 31824, 75582, 167960};
	int num12[] = {1, 13, 91, 455, 1820, 6188, 18564, 50388, 125970};
	int num13[] = {1, 14, 105, 560, 2380, 8568, 27132, 77520};
	int num14[] = {1, 15, 120, 680, 3060, 11628, 38760};
	int num15[] = {1, 16, 136, 816, 3876, 15504};
	int num16[] = {1, 17, 153, 969, 4845};
	int num17[] = {1, 18, 171, 1140};
	int num18[] = {1, 19, 190};
	int num19[] = {1, 20};
	int num20[] = {1};

	PTRI triangle [] = {
								{0, num0},
								{1, num1},
								{2, num2},
								{3, num3},
								{4, num4},
								{5, num5},
								{6, num6},
								{7, num7},
								{8, num8},
								{9, num9},
								{10, num10},
								{11, num11},
								{12, num12},
								{13, num13},
								{14, num14},
								{15, num15},
								{16, num16},
								{17, num17},
								{18, num18},
								{19, num19},
								{20, num20}
	};

	const int triangle_size = sizeof(triangle)/sizeof(triangle[0]);

	std::vector<int> compute_triangle(size_t n)
	{		
		std::vector<int> row;
		std::vector<std::vector<int>> tri(n, row);

		//load pre-calculated values in triangle
		for(size_t i = 0; i < n_max; ++i)
			for(size_t j = 0; j < n_max-i; ++j)
				tri[i].push_back(triangle[i].pos[j]);

		//compute new values -- if n <= nmax this will be skipped
		for(size_t i = n_max; i < n; i++)
			tri[i][0] = tri[0][i] = 1;
			
		for(size_t i = 1; i < n-1; i++)
			for(size_t j = n_max; j < n-i; j++)
				tri[i][j] = tri[i-1][j] + tri[i][j-1];

		std::vector<int> retval(n);

		//extract the binomial coefficients
		for(size_t i = 0; i < n; ++i)
			retval[i] = tri[n-i-1][i]; 

		return retval;
	}

}

#endif