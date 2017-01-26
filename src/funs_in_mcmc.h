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
