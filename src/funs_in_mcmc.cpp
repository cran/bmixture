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

#include "funs_in_mcmc.h"

// function to find the order of a vector
// vec : input/output as a vector of double values
// order: output which is a order of vector "vec" as a vector of int
void order_vec( double vec[], int order[], int *size_vec )
{
    int n = *size_vec;

    vector<pair<int,double> > pair_vec( n );
    for( int i = 0; i < n; i++ )
        pair_vec[i] = pair<int,double>( i, vec[i] );

    Sorter sorter;
    sort( pair_vec.begin(), pair_vec.end(), sorter );

    for( int i = 0; i < n; i++ )
    {
        order[ pair_vec[i].first ] = i;
        vec[i] = pair_vec[i].second;
    }
}
   
// sample function 
// prob:              input as a vector with size k
// selected : input/output 
void sample_c( double prob[], int *index_selected, int *k_c )
{
	//GetRNGstate();
	int k_star = *k_c;

	for( int i = 1; i < k_star; i++ )
		prob[i] += prob[ i - 1 ];

	double random_value = prob[k_star - 1] * unif_rand();

	for( int i = 0; i < k_star; i++ )
	{
		if( prob[i] > random_value )
		{
			*index_selected = i;
			break;
		}
	}
	//PutRNGstate();
} 
    
