#include <R.h>
#include <Rmath.h>
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <math.h>        // pow 
#include <algorithm>     // for sort function
#include <limits>        // for numeric_limits<long double>::max()
#include <R_ext/Arith.h> // for R_PosInf, R_NegInf, R_NaReal (more commonly used as NA_REAL)

#include "funs_in_mcmc.h"

using namespace std;

extern "C" {
	
// update z and calculate n_i based on z
void update_z( int z[], int n_i[], 
			double data_c[], int *n_c, int *k_c, 
			double mu_c[], double sig_c[], double pi_c[] 
			)
{
	//GetRNGstate();
	vector<double> prob_z( *k_c ); 
	// Sample group memberships from multinominal distribution
	// for( j in 1:n ) z[, j] = rmultinom( n = 1, size = 1, prob = pi * dnorm( data[j], mean = mu, sd = sqrt( sig ) ) )	
	for( int j = 0; j < *n_c; j++ )
	{
		for( int i = 0; i < *k_c; i++ ) 
			prob_z[i] = pi_c[i] * dnorm( data_c[j], mu_c[i], sqrt( sig_c[i] ), 0 );
		
		int selected_i;
		sample_c( &prob_z[0], &selected_i, k_c );
		
		for( int i = 0; i < *k_c; i++ )
			z[j * *k_c + i] = 0;
			
		z[j * *k_c + selected_i] = 1;
	}
	
//-- Computing n_i :   n_i = apply( z, 1, sum ) -------------------------------|
	int sum_n;
	for( int i = 0; i < *k_c; i++ )
	{	
		sum_n = 0;
		for( int j = 0; j < *n_c; j++ )
			sum_n += z[ j * *k_c + i ];
		n_i[i] = sum_n;
	}		
	//PutRNGstate();
}	

// update beta
void update_beta( 
			double *beta_new, int *n_c, int *k_c, 
			double *alpha_c, double *g_c, double *h_c,
			double sig_c[] 
			)
{
	//GetRNGstate();
	// beta_r = rgamma( 1, g + k * alpha, h + sum( 1 / sig ) )
	double sum_inv_sig = 0.0;
	for( int i = 0; i < *k_c; i++ ) 
		sum_inv_sig += 1.0 / sig_c[i];
	
	*beta_new = rgamma( *g_c + *k_c * *alpha_c, 1.0 / ( *h_c + sum_inv_sig ) );	
	//PutRNGstate();
}
    
// update pi : Sample mixtures proportions (pi) from a dirichlet distribution
void update_pi( double pi_c[], int n_i[], int *n_c, int *k_c )
{
	//GetRNGstate();
	// for( i in 1:k ) pi[i] = rgamma( n = 1, shape = 1 + n_i[i], scale = 1 ) 
	for( int i = 0; i < *k_c; i++ ) 
		pi_c[i] = rgamma( 1 + n_i[i], 1 ); 
	
	// Normalise to dirichelet
	// pi = pi / sum( pi )
	double sum_pi = 0.0;
	for( int i = 0; i < *k_c; i++ ) sum_pi  += pi_c[i];
	for( int i = 0; i < *k_c; i++ ) pi_c[i] /= sum_pi;
	//PutRNGstate();
}
    
// update mu
void update_mu( 
			double data_c[], int z[], int n_i[], int *n_c, int *k_c, 
			double *epsilon_c, double *kappa_c,
			double mu_c[], double sig_c[] 
			)
{
	//GetRNGstate();
	for( int i = 0; i < *k_c; i++ ) 
	{
		// Set mixture sample mean
		// sum_data = sum( z[i,] * data )
		double sum_data = 0.0;
		for( int j = 0; j < *n_c; j++ )
			sum_data += z[ j * *k_c + i ] * data_c[j];
		
		// Set posterior mean and standard deviation
		// mu_sig  = 1 / ( n_i[i] / sig[i] + kappa_r ) 
		double mu_sig  = 1.0 / ( n_i[i] / sig_c[i] + *kappa_c ); 
		//     mu_mean = ( sum_data / sig[i] + kappa * epsilon ) * mu_sig
		double mu_mean = ( sum_data / sig_c[i] + *kappa_c * *epsilon_c ) * mu_sig;

		// Sample new mean from normal distribution
		mu_c[i] = rnorm( mu_mean, sqrt( mu_sig ) );
	}
	//PutRNGstate();
}
    
// update sig
void update_sig( double *beta_new, 
			double data_c[], int z[], int n_i[], int *n_c, int *k_c, 
			double *alpha_c,
			double mu_c[], double sig_c[] 
			)
{
	//GetRNGstate();
	for( int i = 0; i < *k_c; i++ ) 
	{
		// Set inverse gamma posterior parameters
		double alpha_sig = *alpha_c + n_i[i] / 2;
		
		// beta_sig  = beta_r  + ( sum( z[i,] * ( data - mu[i] ) ^ 2 ) ) / 2 
		double sum_z_datamu = 0.0;
		for( int j = 0; j < *n_c; j++ ) 
		{
			double data_mu  = data_c[j] - mu_c[i];
			double z_datamu = z[ j * *k_c + i ] * pow( data_mu, 2.0 );
			
			sum_z_datamu += z_datamu ;
		}
		
		double beta_sig  = *beta_new + sum_z_datamu / 2.0; 
		
		// Sample variance from gamma random variable and invert it
		// sig[i] = 1 / rgamma( n = 1, alpha_sig, beta_sig )
		sig_c[i] = 1.0 / rgamma( alpha_sig, 1.0 / beta_sig );
	}	
	//PutRNGstate();
}
    
// update all the parameters
void update_parameters_bmixnorm(
			double data_c[], int *n_c, int *k_c, 
			double *epsilon_c, double *kappa_c, double *alpha_c, double *g_c, double *h_c,
			double mu_c[], double sig_c[], double pi_c[] 
			)
{
//-- STEP 1: sample latend variables (matrix z) -------------------------------|
  
	vector<int> z( *k_c * *n_c );         // z  = matrix( 0, nrow = k, ncol = n ) 
	vector<int> n_i( *k_c );
	update_z( &z[0], &n_i[0], data_c, n_c, k_c, mu_c, sig_c, pi_c );

//-- STEP 2: updating beta ----------------------------------------------------|
   
	double beta_new;
	update_beta( &beta_new, n_c, k_c, alpha_c, g_c, h_c, sig_c );

//-- STEP 3: updating pi_c ----------------------------------------------------|
   			
	update_pi( pi_c, &n_i[0], n_c, k_c );
		
//-- STEP 4: updating mu_c ----------------------------------------------------|
   
	update_mu( data_c, &z[0], &n_i[0], n_c, k_c, epsilon_c, kappa_c, mu_c, sig_c	);
	
//-- STEP 5: updating sig_c ---------------------------------------------------|
   
	update_sig( &beta_new, data_c, &z[0], &n_i[0], n_c, k_c, alpha_c, mu_c, sig_c );
}

// sort all parameters based on pi
void sort_sample( double mu_c[], double sig_c[], double pi_c[], int *k_c )
{
	if( *k_c > 1 )   
	{ 		
		vector<int> order_pi( *k_c );
		order_vec( pi_c, &order_pi[0], k_c );

		vector<double> mu_copy( *k_c );
		vector<double> sig_copy( *k_c );
		memcpy( &mu_copy[0],  &mu_c[0],  sizeof( double ) * *k_c );
		memcpy( &sig_copy[0], &sig_c[0], sizeof( double ) * *k_c );
	
		for( int i = 0; i < *k_c; i++ )
		{
			mu_c[i]  = mu_copy[ order_pi[i] ];
			sig_c[i] = sig_copy[ order_pi[i] ];
		}
	}	
}

/*
 * MCMC for finite mixture of Normal distribution  
 * for case the number of component (k) is known
*/
void bmix_norm_k_unknown( double data_r[], int *n, int *k, int *k_max_r, int *iter, int *burnin, double *lambda_r,
						double pi_sample[], double mu_sample[], double sig_sample[],
						int all_k[], double all_weights[],
						double *epsilon, double *kappa_r, double *alpha, double *beta_r, double *g, double *h,
						double mu[], double sig[], double pi_r[] )
{
	GetRNGstate();
	int n_c = *n, k_c = *k, iteration = *iter, burn_in = *burnin, sweep = iteration - burn_in;
	int i, j, ij;
	
	double epsilon_c = *epsilon, kappa_c = *kappa_r, alpha_c = *alpha, g_c = *g, h_c = *h;
	double lambda_c = *lambda_r, k_max_c = *k_max_r, beta_c = *beta_r;

	vector<double> pi_c(  k_c ); 
	vector<double> mu_c(  k_c ); 
	vector<double> sig_c( k_c ); 
	
	memcpy( &pi_c[0],  pi_r,  sizeof( double ) * k_c );
	memcpy( &mu_c[0],  mu,    sizeof( double ) * k_c );
	memcpy( &sig_c[0], sig,   sizeof( double ) * k_c );
	
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
					// likelhood = sum( pi_r * dnorm( data[t], mu, sqrt( sig ) ) )
					double likelhood = 0.0;
					for( i = 0; i < k_c; i++ ) 
						likelhood += pi_c[i] * dnorm( data_c[t], mu_c[i], sqrt( sig_c[i] ), 0 ); 
					
					if( likelhood == 0 ) likelhood = 1e-300; 
					
					//d_data_t_j = dnorm( data[t], mu[j], sqrt( sig[j] ) )
					double d_data_t_j = dnorm( data_c[t], mu_c[j], sqrt( sig_c[j] ), 0 );
					
					//log_death_rates = log_death_rates + log( 1 - ( pi_r[j] * d_data_t_j ) / likelhood ) - log( 1 - pi_r[j] )
					//log_death_rates += log( 1.0 - ( pi_c[j] * d_data_t_j ) / likelhood ) - log( 1.0 - pi_c[j] );
					log_death_rates += log1p( - pi_c[j] * d_data_t_j / likelhood ) - log1p( - pi_c[j] );
				}
				
				if( log_death_rates == R_NegInf ) 
				{ 
					death_rates[j] = min_numeric_limits_ld; 
				}else{
					death_rates[j] = exp( log_death_rates );
				}
				
				if( death_rates[j] == R_PosInf ) 
					death_rates[j] = max_numeric_limits_ld;			
			}	
		}
				
		// simulate the time s to the next jump, from an exponential distribution
		//~ weight = 1 / ( lambda + sum( death_rates ) )
		double sum_death_rates = 0.0;
		for( i = 0; i < k_c; i++ ) 
			sum_death_rates += death_rates[i];
		
		double weight = 1.0 / ( lambda_c + sum_death_rates );
		
		// first Birth and death run and modify the parameters
		double birthp = lambda_c * weight;
		
		if( ( unif_rand() < birthp ) && ( k_c < k_max_c ) )
		{
			// BIRTH EVENT --------------
			
			double pi_new = rbeta( 1.0, k_c );
						
			//~ pi_r   <- c( pi_r * ( 1 - pi_new ), pi_new )
			for( i = 0; i < k_c; i++ )
				pi_c[i] *= ( 1. - pi_new );
				
			pi_c.push_back( pi_new );
			
			//~ mu <- c( mu,  rnorm( 1, epsilon, sqrt( 1 / kappa_r ) ) )
			double mu_new = rnorm( epsilon_c, sqrt( 1.0 / kappa_c ) );
			mu_c.push_back( mu_new );
		
			//~ sig <- c( sig, 1 / rgamma( 1, alpha, beta_r ) )
			double sig_new = 1. / rgamma( alpha_c, 1 / beta_c );
			sig_c.push_back( sig_new );
			
			++k_c;	
		}else{
			// DEATH EVENT --------------
			
			//~ j <- sample( 1:k, 1, prob = death_rates )
			int selected_j;
			sample_c( &death_rates[0], &selected_j, &k_c );
			
			//~ pi_r   <- pi_r[-j] / ( 1 - pi_r[j] )
			double pi_c_j = 1. - pi_c[ selected_j ];
			pi_c.erase( pi_c.begin() + selected_j );     // pi  <- pi[-j]
			for( i = 0; i < k_c; i++ )
				pi_c[i] /= pi_c_j;
			
			mu_c.erase(  mu_c.begin()  + selected_j );   // mu  <- mu[-j]
			sig_c.erase( sig_c.begin() + selected_j );   // sig <- sig[-j]
			
			--k_c;
		}

// ------------  Update all parameters ----------------------------------------|
  
		// updating mu and sig 
        update_parameters_bmixnorm(
			&data_c[0], &n_c, &k_c, 
			&epsilon_c, &kappa_c, &alpha_c, &g_c, &h_c,
			&mu_c[0], &sig_c[0], &pi_c[0] 
			);

// ----------- Sort parameters based on pi ------------------------------------|

		sort_sample( &mu_c[0], &sig_c[0], &pi_c[0], &k_c );	

// ---------- saving result ---------------------------------------------------|	

		all_k[i_mcmc]       = k_c;    // saving all values of k for chicking convergency
		all_weights[i_mcmc] = weight;

		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < k_c; i++ ) 
			{
				ij = i * sweep + counter;
				
				pi_sample[  ij ] = pi_c[i];
				mu_sample[  ij ] = mu_c[i];
				sig_sample[ ij ] = sig_c[i];
			}
				
			++counter;
		}
    } // End main MCMC for loop
	PutRNGstate();
}
    
/*
 * MCMC for finite mixture of Normal distribution  
 * for case the number of component (k) is unknown
*/
void bmix_norm_k_fixed( 
			double data_r[], int *n, int *k, int *iter, int *burnin, 
			double pi_sample[], double mu_sample[], double sig_sample[],
			double *epsilon, double *kappa_r, double *alpha, double *g, double *h,
			double mu[], double sig[], double pi_r[] 
			)
{
	GetRNGstate();
	int n_c = *n, k_c = *k, iteration = *iter, burn_in = *burnin, sweep = iteration - burn_in;
	int i, ij;
	
	double epsilon_c = *epsilon, kappa_c = *kappa_r, alpha_c = *alpha, g_c = *g, h_c = *h;

	vector<double> pi_c(  k_c ); 
	vector<double> mu_c(  k_c ); 
	vector<double> sig_c( k_c ); 
	
	memcpy( &pi_c[0],  pi_r,  sizeof( double ) * k_c );
	memcpy( &mu_c[0],  mu,    sizeof( double ) * k_c );
	memcpy( &sig_c[0], sig,   sizeof( double ) * k_c );
	
	vector<double> data_c( n_c ); 
	memcpy( &data_c[0], data_r, sizeof( double ) * n_c );
		
	int counter = 0;

	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 100 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

// ------------  Update all parameters ----------------------------------------|

        update_parameters_bmixnorm(
			&data_c[0], &n_c, &k_c, 
			&epsilon_c, &kappa_c, &alpha_c, &g_c, &h_c,
			&mu_c[0], &sig_c[0], &pi_c[0] 
			);

// ----------- Sort parameters based on pi ------------------------------------|

		sort_sample( &mu_c[0], &sig_c[0], &pi_c[0], &k_c );	

// ---------- saving result ---------------------------------------------------|	

		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < k_c; i++ ) 
			{
				ij = i * sweep + counter;
				
				pi_sample[  ij ] = pi_c[i];
				mu_sample[  ij ] = mu_c[i];
				sig_sample[ ij ] = sig_c[i];
			}
				
			++counter;
		}

    } // End main MCMC for loop
	PutRNGstate();
}
    

} // End of exturn "C"
