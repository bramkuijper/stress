// the evolution of the stress response in a fluctuating envt
// Bram Kuijper
//  

// important libraries
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>

// random number generation, via the GNU scientific library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// various functions, such as unique filename creation
// which helps when running multiple instances of simulation
// in the same folder
#include "auxiliary.h"

//#define NDEBUG
//
// the debug sign should only be turned on when one wants
// to assess the complete distribution of phenotypes
// 

using namespace std;


// gnu scientific library random number generator initialization
// http://www.gnu.org/software/gsl/ 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_r; // gnu scientific rng 


// number of generations
const long int NumGen = 50000;

// population size
const int Npop = 5000; 

// number of generations to skip when outputting data
const long int skip = 10;

// track population sizes in each of the environments
int numP = 0;
int numNP = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
long int generation = 0;
long int t_change = 0;

// mutation rate parameters
double mu_feedback = 0.0;
double mu_stress_influx = 0.0;
double mu_influx = 0.0;

// environmental switching parameters
// indices here indicate world 1
// and world 2, to include development
// if development is absent, values for
// switch rates are identical
//
// here alpha means switch from P to NP
// and beta from NP to P
double s_PN2P[2] = {0.0,0.0};
double s_NP2P[2] = {0.0,0.0};

// per generation 
// switching probability from world 1 to 2
double gamma12 = 0.0;
double gamma21  = 0.0;

// baseline survival rate
double surv_baseline = 0;

// store the random seed
// to 'replay' the simulation afterwards
unsigned seed = 0;


// the individual struct
struct Individual
{
    // diploid loci specifying the evolving traits
    //
    // self-dependent increase/decrease in hormone
    double feedback[2];

    // influx of new hormone when encountering stress
    double stress_influx[2];

    // stress independent hormone influx
    double influx[2];

    // components of the individual's state
    double hormone; // current hormone level
    double damage; // current damage level

    // whether individual is now living in P or not
    // this is merely for error checking
    bool envt_is_P;
};

// allocate a population and a population of survivors
typedef Individual Population[Npop];

// specify populations living in predator (P)
// or non-predator (NP) patches
Population P, NP;

// generate a unique filename for the output file
string filename("sim_stress");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

// initialize simulations from command line arguments
void init_arguments(int argc, char *argv[])
{
}

// mutation according to a continuum of alleles model
void mutate(double &G, double mu, double sdmu)
{
    G += gsl_rng_uniform(rng_r) < mu ? gsl_ran_gaussian(rng_r, sdmu) : 0;
}

// same but with a lower bound of 0 in the resulting trait
void mutate0(double &G, double mu, double sdmu)
{
    G += gsl_rng_uniform(rng_r) < mu ? gsl_ran_gaussian(rng_r, sdmu) : 0;

    // set the bound
    if (G < 0)
    {
        G = 0;
    }
}

// write the parameters 
// (typically at the end of the output file)
void write_parameters()
{
	DataFile << endl
		<< endl
		<< "seed;" << seed << ";"<< endl;
}



// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void init_population()
{
    // get the timestamp (using nanosecs,
    // so that two concurrent simulations still have
    // different seeds)
    // to initialize the random seed
	seed = get_nanoseconds();
    
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng_r = gsl_rng_alloc(T);
    gsl_rng_set(rng_r, seed);

	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        // randomly allocate over P and NP
        if (gsl_rng_uniform(rng_r) < 0.5)
        {
            P[

        Pop[i].phen = intercept;

        for (int j = 0; j < n_alleles_g; ++j)
        {
            Pop[i].g[j] = intercept == 0 ? 0 : intercept/n_alleles_g;
        }

        for (int j = 0; j < n_alleles_b; ++j)
        {
            Pop[i].b[j] = 0;
        }

        for (int j = 0; j < n_alleles_m; ++j)
        {
            Pop[i].m_m[j] = 0.0;
            Pop[i].m_e[j] = 0.0;
            Pop[i].m_g[j] = 0.0;
            Pop[i].m_m_b[j] = 0.0;
            Pop[i].m_e_b[j] = 0.0;
            Pop[i].m_g_b[j] = 0.0;
        }

        if (ampl * sin(rate * i) > max_sin)
        {
            max_sin = ampl * sin(rate * i);
        } 
	}
}

// create an offspring
void create_offspring(int mother, int father, Individual &kid)
{
    double sum_g = 0; // sum over all the breeding values of the offspring coding for the actual phenotype
    double sum_b = 0; // sum over all the breeding values of the offspring coding for the norm of reaction
    double sum_m_m = 0; // sum over all the breeding values of the offspring coding for the maternal effect
    double sum_m_g = 0; // sum over all the breeding values of the offspring coding for the maternal effect
    double sum_m_e = 0; // sum over all the breeding values of the offspring coding for the maternal effect

    double sum_m_m_b = 0; 
    double sum_m_g_b = 0; 
    double sum_m_e_b = 0; 
    
    // we assume all loci are unlinked 
    for (int i = 0; i < n_alleles_g;++i)
    {
        kid.g[i] = i % 2 == 0 ? 
            Survivors[mother].g[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].g[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.g[i], mu_g, sdmu);
        sum_g += kid.g[i];
    }

    for (int i = 0; i < n_alleles_b; ++i)
    {
        kid.b[i] = i % 2 == 0 ? 
            Survivors[mother].b[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].b[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.b[i], mu_b, sdmu);

        sum_b += kid.b[i];
    }

    for (int i = 0; i < n_alleles_m; ++i)
    {
        kid.m_m[i] = i % 2 == 0 ? 
            Survivors[mother].m_m[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_m[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.m_m[i], mu_m_m, sdmu);
        sum_m_m += kid.m_m[i];
        
        kid.m_e[i] = i % 2 == 0 ? 
            Survivors[mother].m_e[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_e[i - 1 + gsl_rng_uniform_int(rng_r, 2)];


        mutate(kid.m_e[i], mu_m_e, sdmu);
        sum_m_e += kid.m_e[i];
        
        kid.m_g[i] = i % 2 == 0 ? 
            Survivors[mother].m_g[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_g[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.m_g[i], mu_m_g, sdmu);
        sum_m_g += kid.m_g[i];

        // interaction coefficients
        kid.m_m_b[i] = i % 2 == 0 ? 
            Survivors[mother].m_m_b[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_m_b[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.m_m_b[i], mu_m_m_b, sdmu);
        sum_m_m_b += kid.m_m_b[i];


        
        kid.m_e_b[i] = i % 2 == 0 ? 
            Survivors[mother].m_e_b[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_e_b[i - 1 + gsl_rng_uniform_int(rng_r, 2)];


        mutate(kid.m_e_b[i], mu_m_e_b, sdmu);
        sum_m_e_b += kid.m_e_b[i];
       


        kid.m_g_b[i] = i % 2 == 0 ? 
            Survivors[mother].m_g_b[i + gsl_rng_uniform_int(rng_r, 2)] 
            : 
            Survivors[father].m_g_b[i - 1 + gsl_rng_uniform_int(rng_r, 2)];

        mutate(kid.m_g_b[i], mu_m_g_b, sdmu);
        sum_m_g_b += kid.m_g_b[i];
    }

    kid.phen_m_m = sum_m_m;
    kid.phen_m_g = sum_m_g;
    kid.phen_m_e = sum_m_e;

    kid.phen_m_m_b = sum_m_m_b;
    kid.phen_m_g_b = sum_m_g_b;
    kid.phen_m_e_b = sum_m_e_b;

    kid.phen_g = sum_g;
    kid.phen_b = sum_b;

    // complete maternal control
    kid.phen = 
        kid.phen_g // elevation
        + gsl_ran_gaussian(rng_r,sigma_e)  // developmental noise
        + kid.phen_b * epsilon_sens  // plasticity
        + kid.phen_m_m * Survivors[mother].phen // maternal phenotypic effect
        + kid.phen_m_e * Survivors[mother].envt // maternal phenotypic effect
        + kid.phen_m_g * Survivors[mother].phen_g // maternal genetic effect
        + kid.phen_m_m_b * Survivors[mother].phen * epsilon_sens 
        + kid.phen_m_e_b * Survivors[mother].envt * epsilon_sens 
        + kid.phen_m_g_b * Survivors[mother].phen_g * epsilon_sens; 

    kid.envt = epsilon_sens;

    assert(isnan(kid.phen) == 0);
}


// Survival of juveniles to reproductive adults
void survive()
{
    double W;

    double theta = A + B * epsilon;

    if (!envt_is_changed)
    {
        if (generation > t_change 
                && 
                ampl * sin(rate * generation) > max_sin - 0.1
           )
        {
            rate = rateptb;
            intercept = intptb;
            ampl = amplptb;

            envt_is_changed = true;
        }
    }


    NSurv = 0;

    for (int i = 0; i < Npop; ++i)
    {
        W = wmin + (1.0 -  wmin) * exp(-.5 * (
                                            pow((Pop[i].phen - theta),2.0)/omega2 
                                            + pow(Pop[i].phen_b,2.0)/omega_b_2
                                            + pow(Pop[i].phen_m_m,2.0)/omega_m_m_2
                                            + pow(Pop[i].phen_m_g,2.0)/omega_m_g_2
                                            + pow(Pop[i].phen_m_e,2.0)/omega_m_e_2
                                            + pow(Pop[i].phen_m_m_b,2.0)/omega_m_m_2
                                            + pow(Pop[i].phen_m_g_b,2.0)/omega_m_g_2
                                            + pow(Pop[i].phen_m_e_b,2.0)/omega_m_e_2
                                        )
                                    );


        assert(isnan(W) == 0);

        // let individual survive or not
        if (gsl_rng_uniform(rng_r) < W)
        {
            Survivors[NSurv++] = Pop[i];
        }
    }

    if (NSurv == 0)
    {
        write_parameters();
        exit(1);
    }

    // update the environment for the next generation
    // as an autocorrelated gaussian random variable
    ksi = rho_t*ksi + gsl_ran_gaussian(rng_r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);

    epsilon = intercept + ampl * sin(rate * (generation+1)) + ksi;

    // in the likely case there is a developmental timelag, tau,
    // update the environment for a number of 'sub' timesteps
    // to achieve a 'sensed' (rather than real) value of epsilon
    if (tau > 0)
    {
        int timesteps = rint(1.0 / tau);

        for (int time_i = 0; time_i < timesteps; ++time_i)
        {
            ksi = rho_t*ksi + gsl_ran_gaussian(rng_r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);
        }
    }

    // update the value of the sensed environment
    epsilon_sens = intercept + ampl * sin(rate * (generation-tau+1)) + ksi;

    // finally, let the survivors produce N offspring 
    for (int i = 0; i < Npop; ++i)
    {
        Individual Kid;

        create_offspring(
                gsl_rng_uniform_int(rng_r,NSurv), 
                gsl_rng_uniform_int(rng_r,NSurv), 
                Kid);

        Pop[i] = Kid;
    }
}


// write down summary statistics
void write_data()
{
    double meanphen = 0;
    double meang = 0;
    double ssg = 0;
    double meanm_m = 0;
    double meanm_g = 0;
    double meanm_e = 0;
    double meanm_m_b = 0;
    double meanm_g_b = 0;
    double meanm_e_b = 0;
    double ssm_m = 0;
    double ssm_g = 0;
    double ssm_e = 0;
    double ssm_m_b = 0;
    double ssm_g_b = 0;
    double ssm_e_b = 0;
    double meanb = 0;
    double ssb = 0;

    // get stats from the population
    for (int i =  0; i < Npop; ++i)
    {
        // stats for m
        meang += Pop[i].phen_g;
        ssg += Pop[i].phen_g * Pop[i].phen_g;

        // stats for m
        meanm_m += Pop[i].phen_m_m;
        ssm_m += Pop[i].phen_m_m * Pop[i].phen_m_m;
        
        meanm_g += Pop[i].phen_m_g;
        ssm_g += Pop[i].phen_m_g * Pop[i].phen_m_g;
        
        meanm_e += Pop[i].phen_m_e;
        ssm_e += Pop[i].phen_m_e * Pop[i].phen_m_e;
        
        // stats for m
        meanm_m_b += Pop[i].phen_m_m_b;
        ssm_m_b += Pop[i].phen_m_m_b * Pop[i].phen_m_m_b;
        
        meanm_g_b += Pop[i].phen_m_g_b;
        ssm_g_b += Pop[i].phen_m_g_b * Pop[i].phen_m_g_b;
        
        meanm_e_b += Pop[i].phen_m_e_b;
        ssm_e_b += Pop[i].phen_m_e_b * Pop[i].phen_m_e_b;

        meanb += Pop[i].phen_b;
        ssb += Pop[i].phen_b * Pop[i].phen_b;

        meanphen += Pop[i].phen;
    }

    DataFile << generation << ";" << epsilon << ";" << NSurv << ";" << ksi << ";";

    DataFile 
            << (meanphen/Npop) << ";"
            << (meang/(Npop)) << ";"
            << (ssg/(Npop) - pow(meang/Npop,2.0)) << ";"
            << (meanm_m/(Npop)) << ";"
            << (ssm_m/Npop - pow(meanm_m/Npop,2.0)) << ";" 
            << (meanm_g/(Npop)) << ";"
            << (ssm_g/Npop - pow(meanm_g/Npop,2.0)) << ";" 
            << (meanm_e/(Npop)) << ";"
            << (ssm_e/Npop - pow(meanm_e/Npop,2.0)) << ";" 
            
            << (meanm_m_b/(Npop)) << ";"
            << (ssm_m_b/Npop - pow(meanm_m_b/Npop,2.0)) << ";" 
            << (meanm_g_b/(Npop)) << ";"
            << (ssm_g_b/Npop - pow(meanm_g_b/Npop,2.0)) << ";" 
            << (meanm_e_b/(Npop)) << ";"
            << (ssm_e_b/Npop - pow(meanm_e_b/Npop,2.0)) << ";" 
           
            << (meanb/Npop) << ";"
            << (ssb/Npop - pow(meanb/Npop,2.0)) << ";"  << endl;
}

// write the headers of a datafile
void write_data_headers()
{
    DataFile << "generation;" 
        << "epsilon;"
        << "nsurv;"
        << "ksi;"
        << "meanz;"
        << "meang;"
        << "varg;"
        << "meanm_m;"
        << "varm_m;"
        << "meanm_g;"
        << "varm_g;"
        << "meanm_e;"
        << "varm_e;"
        << "meanm_m_b;"
        << "varm_m_b;"
        << "meanm_g_b;"
        << "varm_g_b;"
        << "meanm_e_b;"
        << "varm_e_b;"
        << "meanb;"
        << "varb;" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	write_data_headers();
	init_population();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		survive();

        if (do_stats)
		{
			write_data();
		}
	}

	write_parameters();
}
