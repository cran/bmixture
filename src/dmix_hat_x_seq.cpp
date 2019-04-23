// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2017 - 2019 Reza Mohammadi                                                    |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     "bmixture" is free software: you can redistribute it and/or modify it under                 |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#include "util_bmixture.h"
#include "funs_in_mcmc.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

extern "C" {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixnorm_hat_x_seq_unknow_k( double x_seq[], double f_hat_x_seq[], 
					 double pi_sample[], double mu_sample[], double sig_sample[],
					 int k[], int *sweep_r, int *size_x_seq_r )
{
	int sweep      = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			int k_j = k[ j ];
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < k_j; i++ )
			{
				int ji = i * sweep + j;
				//f_mix_t = f_mix_t + pi_r[ j, i ] * dnorm( x_seq[t], mu[ j, i ], sqrt( sig[ j, i ] ) )
				f_hat_x_t += pi_sample[ ji ] * dnorm( x_seq[ t ], mu_sample[ ji ], sqrt( sig_sample[ ji ] ), 0 );
			}
			
			//f_hat_x_seq[t] = f_hat_x_seq[t] + f_mix_t
			f_hat_x_seq[ t ] += f_hat_x_t;
			//f_hat_x_seq[t] += ( f_hat_x_t * weights[j] );
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixnorm_hat_x_seq_fixed_k( double x_seq[], double f_hat_x_seq[], 
					 double pi_sample[], double mu_sample[], double sig_sample[],
					 int *size_mix, int *sweep_r, int *size_x_seq_r )
{
	int sweep      = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	int size_mix_c = *size_mix;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < size_mix_c; i++ )
			{
				int ji = i * sweep + j;
				//f_mix_t = f_mix_t + pi_r[ j, i ] * dnorm( x_seq[t], mu[ j, i ], sqrt( sig[ j, i ] ) )
				f_hat_x_t += pi_sample[ ji ] * dnorm( x_seq[ t ], mu_sample[ ji ], sqrt( sig_sample[ ji ] ), 0 );
			}
			
			//f_hat_x_seq[t] = f_hat_x_seq[t] + f_mix_t
			f_hat_x_seq[ t ] += f_hat_x_t;
		}
	}
}

// ----- For mixture of t-distribution ----------------------------------------|

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixt_hat_x_seq_fixed_k( double x_seq[], double f_hat_x_seq[], int *df_t, 
					 double pi_sample[], double mu_sample[], double sig_sample[],
					 int *size_mix, int *sweep_r, int *size_x_seq_r )
{
	int sweep      = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	int size_mix_c = *size_mix;
	int df         = *df_t;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < size_mix_c; i++ )
			{
				int ji = i * sweep + j;
				f_hat_x_t += pi_sample[ ji ] * Rf_dt( ( x_seq[ t ] - mu_sample[ ji ] ) / sqrt( sig_sample[ ji ] ), df, 0 );
			}
			
			f_hat_x_seq[ t ] += f_hat_x_t;
		}
	}
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixt_hat_x_seq_unknow_k( double x_seq[], double f_hat_x_seq[], int *df_t, 
					 double pi_sample[], double mu_sample[], double sig_sample[],
					 int k[], int *sweep_r, int *size_x_seq_r )
{
	int sweep = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	int df         = *df_t;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			int k_j = k[ j ];
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < k_j; i++ )
			{
				int ji = i * sweep + j;
				f_hat_x_t += pi_sample[ ji ] * Rf_dt( ( x_seq[ t ] - mu_sample[ ji ] ) / sqrt( sig_sample[ ji ] ), df, 0 );
			}
			
			f_hat_x_seq[ t ] += f_hat_x_t;
		}
	}
}
     
// ----- For mixture of gamma distribution ------------------------------------|

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixgamma_hat_x_seq_fixed_k( double x_seq[], double f_hat_x_seq[], 
					 double pi_sample[], double alpha_sample[], double beta_sample[],
					 int *size_mix, int *sweep_r, int *size_x_seq_r )
{
	int sweep      = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	int size_mix_c = *size_mix;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < size_mix_c; i++ )
			{
				int ji = i * sweep + j;
				//f_mix_t = f_mix_t + pi_r[ j, i ] * dnorm( x_seq[t], mu[ j, i ], sqrt( sig[ j, i ] ) )
				f_hat_x_t += pi_sample[ ji ] * dgamma( x_seq[ t ], alpha_sample[ ji ], 1.0 / beta_sample[ ji ], 0 );
			}
			
			f_hat_x_seq[ t ] += f_hat_x_t;
		}
	}
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void dmixgamma_hat_x_seq_unknow_k( double x_seq[], double f_hat_x_seq[], 
					 double pi_sample[], double alpha_sample[], double beta_sample[],
					 int k[], int *sweep_r, int *size_x_seq_r )
{
	int sweep      = *sweep_r;
	int size_x_seq = *size_x_seq_r;
	
	for( int t = 0; t < size_x_seq; t++ )
	{
		for( int j = 0; j < sweep; j++ )
		{
			int k_j = k[ j ];
			double f_hat_x_t = 0.0;
			
			for( int i = 0; i < k_j; i++ )
			{
				int ji = i * sweep + j;
				f_hat_x_t += pi_sample[ ji ] * dgamma( x_seq[ t ], alpha_sample[ ji ], 1.0 / beta_sample[ ji ], 0 );
			}
			
			f_hat_x_seq[ t ] += f_hat_x_t;
		}
	}
}
     
} // End of exturn "C"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

