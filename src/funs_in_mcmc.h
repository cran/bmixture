// ------------------------------------------------------------------------------------------------|
//     Copyright (C) 2017 - 2018 Reza Mohammadi                                                    |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     "bmixture" is free software: you can redistribute it and/or modify it under                 |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// ------------------------------------------------------------------------------------------------|

#ifndef matrix_H
#define matrix_H

#include <R.h>
#include <Rmath.h>
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <math.h>        // pow 
#include <algorithm>     // for sort function

using namespace std;

extern "C" {

	class Sorter {
	public:
		bool operator()( pair<int, double> const &a, pair<int, double> const &b )
			{
				return a.second < b.second;
			}
	};
   
	void order_vec( double vec[], int order[], int *size_vec );
	
	void sample_c( double prob[], int *index_selected, int *k_c );

}

#endif
