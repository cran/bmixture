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
// update z and calculate n_i based on z
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_z_gamma( int z[], int n_i[], 
			double data_c[], int *n_c, int *k_c, 
			double alpha_c[], double beta_c[], double pi_c[] 
			)
{
	//GetRNGstate();
	vector<double> prob_z( *k_c ); 
	// Sample group memberships from multinominal distribution
	for( int j = 0; j < *n_c; j++ )
	{
		// prob_z[i] = pi_r[i] * dgamma( data[j], alpha[i], beta[i] )
		for( int i = 0; i < *k_c; i++ ) 
			prob_z[ i ] = pi_c[ i ] * Rf_dgamma( data_c[ j ], alpha_c[ i ], 1.0 / beta_c[ i ], 0 ); 
		
		int selected_i;
		sample_c( &prob_z[0], &selected_i, k_c );
		
		for( int i = 0; i < *k_c; i++ )
			z[ j * *k_c + i ] = 0;
			
		z[ j * *k_c + selected_i ] = 1;
	}
	
	//-- Computing n_i :   n_i = apply( z, 1, sum ) -------------------------------|
	int sum_n;
	for( int i = 0; i < *k_c; i++ )
	{	
		sum_n = 0;
		
		for( int j = 0; j < *n_c; j++ )
			sum_n += z[ j * *k_c + i ];
		
		n_i[ i ] = sum_n;
	}		
	//PutRNGstate();
}	
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// update pi : Sample mixtures proportions (pi) from a dirichlet distribution
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_pi_gamma( double pi_c[], int n_i[], int *n_c, int *k_c )
{
	//GetRNGstate();
	// for( i in 1:k ) pi[i] = rgamma( n = 1, shape = 1 + n_i[i], scale = 1 ) 
	for( int i = 0; i < *k_c; i++ ) 
		pi_c[ i ] = rgamma( 1.0 + n_i[ i ], 1.0 ); 
	
	// Normalise to dirichelet
	// pi = pi / sum( pi )
	double sum_pi = 0.0;
	for( int i = 0; i < *k_c; i++ ) sum_pi    += pi_c[ i ];
	for( int i = 0; i < *k_c; i++ ) pi_c[ i ] /= sum_pi;
	//PutRNGstate();
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// update alpha and beta
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_alpha_beta( double alpha_c[], double beta_c[],
			double data_c[], int z[], int n_i[], int *n_c, int *k_c, 
			double *mu_c, double *nu_c, double *kesi_c, double *tau_c
			)
{
	//GetRNGstate();
	for( int i = 0; i < *k_c; i++ ) 
	{
		// step 3.1 : updating beta
		// beta[i] <- rgamma( 1, kesi + n_i[i] * alpha[i], tau + sum( data * z[i,] ) ) 
		double sum_data_z = 0.0;
		for( int j = 0; j < *n_c; j++ )
			sum_data_z += z[ j * *k_c + i ] * data_c[ j ];
		
		beta_c[ i ] = rgamma( *kesi_c + n_i[ i ] * alpha_c[ i ], 1.0 / ( *tau_c + sum_data_z ) ); 
		
		// step 3.2 : updating alpha using a Metropolis-Hastings
		double proposed_alpha = rgamma( *mu_c, 1.0 / *nu_c ); // proposed_alpha <- rgamma( 1, mu, nu )

		double sum_log_data_z = 0.0;
		for( int j = 0; j < *n_c; j++ ) 
		{
			int ij = j * *k_c + i;
			if( z[ ij ] == 1 ) sum_log_data_z += log( beta_c[ i ] * data_c[ j ] * z[ ij ] );
		}
				
		//     log_accept_alpha = n_i[i] * ( lgamma( alpha[i] ) - lgamma( proposed_alpha ) ) + ( proposed_alpha - alpha[i] ) * sum_log_data_z 
		double log_accept_alpha = n_i[ i ] * ( lgammafn( alpha_c[ i ] ) - lgammafn( proposed_alpha ) ) + ( proposed_alpha - alpha_c[ i ] ) * sum_log_data_z; 
				
		//if( log_accept_alpha > log( runif( 1 ) ) ) alpha[i] = proposed_alpha 
		if( log_accept_alpha > log( unif_rand() ) ) alpha_c[ i ] = proposed_alpha; 
	}
	//PutRNGstate();
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// update all the parameters
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_mcmc_bmixgamma(
			double data_c[], int *n_c, int *k_c, 
			double *mu_c, double *nu_c, double *kesi_c, double *tau_c,
			double alpha_c[], double beta_c[], double pi_c[]
			 )
{
	//-- STEP 1: sample latend variables (matrix z) -------------------------------|
  
	vector<int> z( *k_c * *n_c );         // z  = matrix( 0, nrow = k, ncol = n ) 
	vector<int> n_i( *k_c );
	update_z_gamma( &z[0], &n_i[0], data_c, n_c, k_c, alpha_c, beta_c, pi_c );

	//-- STEP 2: updating pi ------------------------------------------------------|

	update_pi_gamma( pi_c, &n_i[0], n_c, k_c );
		
	//-- STEP 3: updating alpha and beta ------------------------------------------|
   	
	update_alpha_beta( alpha_c, beta_c, data_c, &z[0], &n_i[0], n_c, k_c, mu_c, nu_c, kesi_c, tau_c );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sort all parameters based on pi
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sort_sample_gamma( double alpha_c[], double beta_c[], double pi_c[], int *k_c )
{
	if( *k_c > 1 )   
	{ 		
		vector<int> order_pi( *k_c );
		order_vec( pi_c, &order_pi[0], k_c );

		vector<double> alpha_copy( *k_c );
		vector<double> beta_copy(  *k_c );
		memcpy( &alpha_copy[0], &alpha_c[0], sizeof( double ) * *k_c );
		memcpy( &beta_copy[0] , &beta_c[0] , sizeof( double ) * *k_c );
	
		for( int i = 0; i < *k_c; i++ )
		{
			alpha_c[ i ] = alpha_copy[ order_pi[ i ] ];
			beta_c[  i ] = beta_copy[  order_pi[ i ] ];
		}
	}	
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// MCMC for finite mixture of Normal distribution  
// for case the number of component (k) is known
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void bmix_gamma_unknown_k( double data_r[], int *n, int *k, int *k_max_r, int *iter, int *burnin, double *lambda_r,
						double pi_sample[], double alpha_sample[], double beta_sample[],
						int all_k[], double all_weights[],
			            double *mu, double *nu, double *kesi, double *tau,
						double alpha[], double beta[], double pi_r[] )
{
	GetRNGstate();
	int n_c = *n, k_c = *k, iteration = *iter, burn_in = *burnin, sweep = iteration - burn_in;
	int i, j, ij;
	
	double mu_c = *mu, nu_c = *nu, kesi_c = *kesi, tau_c = *tau;
	double lambda_c = *lambda_r, k_max_c = *k_max_r;

	vector<double> pi_c(    k_c ); 
	vector<double> alpha_c( k_c ); 
	vector<double> beta_c(  k_c ); 
	
	memcpy( &pi_c[0],    pi_r,  sizeof( double ) * k_c );
	memcpy( &alpha_c[0], alpha, sizeof( double ) * k_c );
	memcpy( &beta_c[0],  beta,  sizeof( double ) * k_c );
	
	vector<double> data_c( n_c ); 
	memcpy( &data_c[0], data_r, sizeof( double ) * n_c );
	
	int counter = 0;
	double max_numeric_limits_ld = numeric_limits<double>::max() / 10000;
	double min_numeric_limits_ld = numeric_limits<double>::min() * 10000;
	
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 100 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
	
// ------------  Birth-death step ---------------------------------------------|
		
		// first obtain death_rates (vector of death rates for each component)	
		vector<double> death_rates( k_c );  // zero death rate as a default
		if( k_c > 1 )   
		{ 		
			for( j = 0; j < k_c; j++ )
			{
				double log_death_rates = 0.0;
				for( int t = 0; t < n_c; t++ )
				{
					// likelhood = sum( pi_r * dgamma( data[t], alpha, beta ) )
					double likelhood = 0.0;
					for( i = 0; i < k_c; i++ ) 
						likelhood += pi_c[ i ] * dgamma( data_c[ t ], alpha_c[ i ], 1.0 / beta_c[ i ], 0 ); 
					
					if( likelhood == 0 ) likelhood = min_numeric_limits_ld; 
					
					//     d_data_t_j = dgamma( data[t], alpha[j], beta[j] )
					double d_data_t_j = dgamma( data_c[ t ], alpha_c[ j ], 1.0 / beta_c[ j ], 0 );
					
					//log_death_rates = log_death_rates + log( 1 - ( pi_r[j] * d_data_t_j ) / likelhood ) - log( 1 - pi_r[j] )
					//log_death_rates += log( 1.0 - ( pi_c[j] * d_data_t_j ) / likelhood ) - log( 1.0 - pi_c[j] );
					log_death_rates += log1p( - pi_c[ j ] * d_data_t_j / likelhood ) - log1p( - pi_c[ j ] );
				}
				
				if( log_death_rates == R_NegInf ) 
				{ 
					death_rates[ j ] = min_numeric_limits_ld; 
				}else{
					death_rates[ j ] = exp( log_death_rates );
				}
				
				if( death_rates[ j ] == R_PosInf ) 
					death_rates[ j ] = max_numeric_limits_ld;			
			}	
		}
				
		// simulate the time s to the next jump, from an exponential distribution
		//~ weight = 1 / ( lambda + sum( death_rates ) )
		double sum_death_rates = 0.0;
		for( i = 0; i < k_c; i++ ) 
			sum_death_rates += death_rates[ i ];
		
		double weight = 1.0 / ( lambda_c + sum_death_rates );
		
		// first Birth and death run and modify the parameters
		double birthp = lambda_c * weight;
		
		if( ( unif_rand() < birthp ) && ( k_c < k_max_c ) )
		{
			// BIRTH EVENT --------------
			
			double pi_new = rbeta( 1.0, k_c );
						
			//~ pi_r   <- c( pi_r * ( 1 - pi_new ), pi_new )
			for( i = 0; i < k_c; i++ )
				pi_c[ i ] *= ( 1 - pi_new );
				
			pi_c.push_back( pi_new );
			
			// alpha_new = rgamma( 1, mu  , nu  )
			double alpha_new = rgamma( mu_c, 1.0 / nu_c );
			alpha_c.push_back( alpha_new );
		
			// beta_new = rgamma( 1, kesi, tau )
			double beta_new = rgamma( kesi_c, 1.0 / tau_c );
			beta_c.push_back( beta_new );
			
			++k_c;	
		}else{
			// DEATH EVENT --------------
			
			//~ j <- sample( 1:k, 1, prob = death_rates )
			int selected_j;
			sample_c( &death_rates[0], &selected_j, &k_c );
			
			//~ pi_r   <- pi_r[-j] / ( 1 - pi_r[j] )
			double pi_c_j = 1 - pi_c[ selected_j ];
			pi_c.erase( pi_c.begin() + selected_j );     // pi  <- pi[-j]
			for( i = 0; i < k_c; i++ )
				pi_c[ i ] /= pi_c_j;
			
			alpha_c.erase( alpha_c.begin() + selected_j );   // alpha  <- alpha[-j]
			beta_c.erase(  beta_c.begin()  + selected_j );   // sig <- sig[-j]
			
			--k_c;
		}

// ------------  Update all parameters ----------------------------------------|
  
        update_mcmc_bmixgamma(
			&data_c[0], &n_c, &k_c, 
			&mu_c, &nu_c, &kesi_c, &tau_c,
			&alpha_c[0], &beta_c[0], &pi_c[0] );

// ----------- Sort parameters based on pi ------------------------------------|

		sort_sample_gamma( &alpha_c[0], &beta_c[0], &pi_c[0], &k_c );	

// ---------- saving result ---------------------------------------------------|	


		all_k[       i_mcmc ] = k_c;    // saving all values of k for chicking convergency
		all_weights[ i_mcmc ] = weight;

		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < k_c; i++ ) 
			{
				ij = i * sweep + counter;
				
				pi_sample[    ij ] = pi_c[    i ];
				alpha_sample[ ij ] = alpha_c[ i ];
				beta_sample[  ij ] = beta_c[  i ];
			}
				
			++counter;
		}

    } // End main MCMC for loop
	PutRNGstate();
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// MCMC for finite mixture of Normal distribution  
// for case the number of component (k) is unknown
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void bmix_gamma_fixed_k( 
			double data_r[], int *n, int *k, int *iter, int *burnin, 
			double pi_sample[], double alpha_sample[], double beta_sample[],
			double *mu, double *nu, double *kesi, double *tau,
			double alpha[], double beta[], double pi_r[] )
{
	GetRNGstate();
	int n_c = *n, k_c = *k, iteration = *iter, burn_in = *burnin, sweep = iteration - burn_in;
	int i, ij;
	
	double mu_c = *mu, nu_c = *nu, kesi_c = *kesi, tau_c = *tau;

	vector<double> pi_c(    k_c ); 
	vector<double> alpha_c( k_c ); 
	vector<double> beta_c(  k_c ); 
	
	memcpy( &pi_c[0],    pi_r,  sizeof( double ) * k_c );
	memcpy( &alpha_c[0], alpha, sizeof( double ) * k_c );
	memcpy( &beta_c[0],  beta,  sizeof( double ) * k_c );
	
	vector<double> data_c( n_c ); 
	memcpy( &data_c[0], data_r, sizeof( double ) * n_c );
		
	int counter = 0;

	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 100 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

// ------------  Update all parameters ----------------------------------------|

        update_mcmc_bmixgamma(
			&data_c[0], &n_c, &k_c, 
			&mu_c, &nu_c, &kesi_c, &tau_c,
			&alpha_c[0], &beta_c[0], &pi_c[0] );

// ----------- Sort parameters based on pi ------------------------------------|

		sort_sample_gamma( &alpha_c[0], &beta_c[0], &pi_c[0], &k_c );	

// ---------- saving result ---------------------------------------------------|	

		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < k_c; i++ ) 
			{
				ij = i * sweep + counter;
				
				pi_sample[    ij ] = pi_c[    i ];
				alpha_sample[ ij ] = alpha_c[ i ];
				beta_sample[  ij ] = beta_c[  i ];
			}
				
			++counter;
		}

    } // End main MCMC for loop
	PutRNGstate();
}
    

} // End of exturn "C"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

