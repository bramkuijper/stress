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
const long int NumGen = 10000000;

// population size
const int Npop = 5000; 

// number of generations to skip when outputting data
const long int skip = 1000;

// track population sizes in each of the environments
int numP = 0;
int numNP = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
long int generation = 0;

// mutation rate parameters
double mu_feedback = 0.0;
double mu_stress_influx = 0.0;
double mu_influx = 0.0;
double sdmu = 0.0;

// environmental switching parameters
// indices here indicate world 1
// and world 2, to include development
// if development is absent, values for
// switch rates are identical
//
// switch rate from P to NP
double s_P_2_NP[2] = {0.0,0.0};
// switch rate from NP to P
double s_NP_2_P[2] = {0.0,0.0};

// per generation 
// switching probability from world 1 to 2
// and from world 2 to 1
double s_12[2] = {0.0, 0.0};

// initial values for the evolving traits
double init_feedback = 0.0;
double init_stress_influx = 0.0;
double init_influx = 0.0;

// cue probabilities
double cue_P = 0.0;
double cue_NP = 0.0;

// baseline survival rate
double s0 = 0;

// store the random seed
// to 'replay' the simulation afterwards
unsigned seed = 0;

double ad = 0;
double aP = 0;
double dmax = 0;
double zmax = 0;
double r = 0;
double u = 0;

// which of the two world currently applies
bool current_world = false;


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
// also make two stacks for the new populations after 
// environmental change
Population P, NP, Pnew, NPnew;

// generate a unique filename for the output file
string filename("sim_stress");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

// initialize simulations from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu_feedback  = atof(argv[1]);
    mu_stress_influx  = atof(argv[2]);
    mu_influx  = atof(argv[3]);
    sdmu = atof(argv[4]);
    s_P_2_NP[0]  = atof(argv[5]);
    s_P_2_NP[1]  = atof(argv[6]);
    s_NP_2_P[0]  = atof(argv[7]);
    s_NP_2_P[1]  = atof(argv[8]);
    s_12[0]  = atof(argv[9]);
    s_12[1]  = atof(argv[10]);
    init_feedback  = atof(argv[11]);
    init_stress_influx  = atof(argv[12]);
    init_influx  = atof(argv[13]);
    cue_P  = atof(argv[14]);
    cue_NP  = atof(argv[15]);
    s0  = atof(argv[16]);
    ad  = atof(argv[17]);
    aP  = atof(argv[18]);
    dmax  = atof(argv[19]);
    zmax  = atof(argv[20]);
    r = atof(argv[21]);
    u = atof(argv[22]);
}

// mutation according to a continuum of alleles model
void mutate(
        double &G, 
        double const mu, 
        double const sdmu)
{
    G += gsl_rng_uniform(rng_r) < mu ? gsl_ran_gaussian(rng_r, sdmu) : 0;
}

// mutation but with a lower bound
void mutate_bound_lower(
        double &G, 
        double const mu, 
        double const sdmu, 
        double const bound_lower)
{
    G += gsl_rng_uniform(rng_r) < mu ? gsl_ran_gaussian(rng_r, sdmu) : 0;

    // set the bound
    if (G < bound_lower)
    {
        G = bound_lower;
    }
}

// mutation but within lower and upper bounds
void mutate_bound(
        double &G, 
        double const mu, 
        double const sdmu, 
        double const bound_lower,
        double const bound_upper)
{
    G += gsl_rng_uniform(rng_r) < mu ? gsl_ran_gaussian(rng_r, sdmu) : 0;

    // set the bound
    if (G < bound_lower)
    {
        G = bound_lower;
    }
    else if (G > bound_upper)
    {
        G = bound_upper;
    }
}

// write the parameters 
// (typically at the end of the output file)
void write_parameters()
{
	DataFile << endl
		<< endl
		<< "seed;" << seed << ";"<< endl
		<< "Npop;" << Npop << ";"<< endl
		<< "mu_feedback;" << mu_feedback << ";"<< endl
		<< "mu_stress_influx;" << mu_stress_influx << ";"<< endl
		<< "mu_influx;" << mu_influx << ";"<< endl
		<< "sP2NP_1;" << s_P_2_NP[0] << ";"<< endl
		<< "sP2NP_2;" << s_P_2_NP[1] << ";"<< endl
		<< "sNP2P_1;" << s_NP_2_P[0] << ";"<< endl
		<< "sNP2P_2;" << s_NP_2_P[1] << ";"<< endl
		<< "gamma1;" << s_12[0] << ";"<< endl
		<< "gamma2;" << s_12[1] << ";"<< endl
		<< "init_feedback;" << init_feedback << ";"<< endl
		<< "init_stress_influx;" << init_stress_influx << ";"<< endl
		<< "init_influx;" << init_influx << ";"<< endl
		<< "cue_P;" << cue_P << ";"<< endl
		<< "cue_NP;" << cue_NP << ";"<< endl
		<< "s0;" << s0 << ";"<< endl
		<< "ad;" << ad << ";"<< endl
		<< "aP;" << aP << ";"<< endl
		<< "dmax;" << dmax << ";"<< endl
		<< "zmax;" << zmax << ";"<< endl
		<< "r;" << r << ";"<< endl
		<< "u;" << u << ";"<< endl;
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

    // set counters to 0
    numP = 0;
    numNP = 0;

	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        Individual newInd;

        // set alleles
        for (int allele_i = 0; allele_i < 2; ++allele_i)
        {
            newInd.feedback[allele_i] = init_feedback;
            newInd.stress_influx[allele_i] = init_stress_influx;
            newInd.influx[allele_i] = init_influx;
        }

        // initialize current hormone level and damage as follows:
        newInd.hormone = 0.0;
        newInd.damage = 0.0;

        // put individuals in environment P or NP accordingt to their
        // overall frequency in the first world. This does not really matter
        // too much, as I expect frequencies of individuals in either
        // environment to rapidly attain their ecological equilibria
        newInd.envt_is_P = gsl_rng_uniform(rng_r) < s_P_2_NP[0] / 
            (s_NP_2_P[0] + s_P_2_NP[0]);

        if (newInd.envt_is_P)
        {
            P[numP++] = newInd;
        }
        else
        {
            NP[numNP++] = newInd;
        }
    } 

    current_world = false;
} // end init_population

// create an offspring
void create_offspring(
        Individual &mother, 
        Individual &father, 
        Individual &kid)
{
    // inherit the different gene loci
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        kid.feedback[allele_i] = 
            allele_i < 1 ? 
            mother.feedback[gsl_rng_uniform_int(rng_r,2)] // mother inherits 
            :
            father.feedback[gsl_rng_uniform_int(rng_r,2)]; // father inherits

        mutate_bound(
                kid.feedback[allele_i],
                mu_feedback,
                sdmu,
                0.0,
                1.0);

        kid.stress_influx[allele_i] = 
            allele_i < 1 ? 
            mother.stress_influx[gsl_rng_uniform_int(rng_r,2)] 
            :
            father.stress_influx[gsl_rng_uniform_int(rng_r,2)];

        mutate_bound_lower(
                kid.stress_influx[allele_i],
                mu_stress_influx, 
                sdmu,
                0.0
                );

        kid.influx[allele_i] = 
            allele_i < 1 ? 
            mother.influx[gsl_rng_uniform_int(rng_r,2)] 
            :
            father.influx[gsl_rng_uniform_int(rng_r,2)];

        mutate_bound_lower(
                kid.influx[allele_i],
                mu_influx, 
                sdmu,
                0.0);
    } // inheritance done

    // set hormone level and damage to their baseline values
    kid.hormone = 0.0;
    kid.damage = 0.0;
}

// calculate survival probability
// dependent on 
// - whether or not individual is in environment P
// - its level of damage
// - its current hormone level
double survival(
        bool const is_in_envt_P, 
        double const damage,
        double const hormone_level)
{
    assert(damage >= 0);
    assert(damage <= dmax);
    double p_survive = s0 + (1.0 - s0) * (1.0 - pow(damage/dmax, ad));

    // survive dependent on predator
    if (is_in_envt_P)
    {
        p_survive *= pow(hormone_level/zmax, aP);
    }

    return(p_survive);
}

// survival of individuals
void survive()
{
    // counters for individuals that move to a different
    // patch due to environmental change, with switch rates
    // s_P_2_NP and s_NP_2_P
    int numNPnew = 0;

    // survival in the P population
    for (int ind_i = 0; ind_i < numP; ++ind_i)
    {
        // standard influx and outflux in hormone level
        // z(t+1) = influx + (1 - feedback) * z(t)
        P[ind_i].hormone = 
            0.5 * (P[ind_i].influx[0] + P[ind_i].influx[1])
             + (1.0 - 0.5 * (P[ind_i].feedback[0] + P[ind_i].feedback[1]))
             * P[ind_i].hormone;

        // individual gets predator cue 
        if (gsl_rng_uniform(rng_r) < cue_P)
        {
            // gets cue, spike the hormone level
            P[ind_i].hormone += 
                0.5 * (P[ind_i].stress_influx[0] + P[ind_i].stress_influx[1]);
        }

        // set boundaries to hormone level
        if (P[ind_i].hormone > zmax)
        {
            P[ind_i].hormone = 
                P[ind_i].hormone > zmax ? 
                zmax 
                : 
                P[ind_i].hormone < 0.0 ? 
                    0.0 
                    : 
                    P[ind_i].hormone;
        }
       
        // update damage levels
        P[ind_i].damage = (1.0 - r) * P[ind_i].damage + u * P[ind_i].hormone;

        // set boundaries to damage level
        assert(P[ind_i].damage >= 0);

        if (P[ind_i].damage > dmax)
        {
            P[ind_i].damage = dmax;
        }


        // OK, did not survive dependent on the level of stress & damage
        if (gsl_rng_uniform(rng_r) > 
                survival(true, 
                    P[ind_i].damage,
                    P[ind_i].hormone))
        {
            // delete individual from stack
            P[ind_i] = P[--numP];
            --ind_i;

            assert(numP >= 0);
            assert(numP <= Npop);
            assert(numNP >= 0);
            assert(numNP <= Npop);
        }
        else // Individual survived. Check for potential for envt'al change
        {
            // individual switches to a different environment
            if (gsl_rng_uniform(rng_r)  < s_P_2_NP[current_world])
            {
                // copy individual to NPnew stack which is a 
                // placeholder for all individuals who go to NP
                //
                // we cannot directly transfer individuals to the NP
                // stack itself, as the NP stack still needs to undergo
                // survival selection.
                NPnew[numNPnew++] = P[ind_i];
                P[ind_i] = P[--numP];
                --ind_i;
                assert(numP >= 0);
                assert(numP <= Npop);
                assert(numNP >= 0);
                assert(numNP <= Npop);
            }
        }
    }

    // survival in the NP population
    for (int ind_i = 0; ind_i < numNP; ++ind_i)
    {
        // standard influx and outflux in hormone level
        // z(t+1) = influx + (1 - feedback) * z(t)
        NP[ind_i].hormone = 
            0.5 * (NP[ind_i].influx[0] + NP[ind_i].influx[1])
             + (1.0 - 0.5 * (NP[ind_i].feedback[0] + NP[ind_i].feedback[1]))
             * NP[ind_i].hormone;

        // individual gets predator cue 
        if (gsl_rng_uniform(rng_r) < cue_NP)
        {
            // gets cue, spike the hormone level
            NP[ind_i].hormone += 
                0.5 * (NP[ind_i].stress_influx[0] + NP[ind_i].stress_influx[1]);
        }

        // set boundaries to hormone level
        if (NP[ind_i].hormone > zmax)
        {
            NP[ind_i].hormone = 
                NP[ind_i].hormone > zmax ? 
                zmax 
                : 
                NP[ind_i].hormone < 0.0 ? 
                    0.0 
                    : 
                    NP[ind_i].hormone;
        }

        // update damage levels
        NP[ind_i].damage = (1.0 - r) * NP[ind_i].damage + u * NP[ind_i].hormone;

        // set boundaries to damage level
        assert(NP[ind_i].damage >= 0);

        if (NP[ind_i].damage > dmax)
        {
            NP[ind_i].damage = dmax;
        }

        // OK, did not survive dependent on the level of stress
        if (gsl_rng_uniform(rng_r) > 
                survival(false, 
                    NP[ind_i].damage,
                    NP[ind_i].hormone))
        {
            // delete individual from stack
            NP[ind_i] = NP[--numNP];
            --ind_i;

            assert(numP >= 0);
            assert(numP <= Npop);
            assert(numNP >= 0);
            assert(numNP <= Npop);
        }
        else // Individual survived. Check for potential for envt'al change
        {
            // individual switches to a different environment
            if (gsl_rng_uniform(rng_r)  < s_P_2_NP[current_world])
            {
                // as all survival is done we can directly copy individual to
                // P stack and delete it here
                P[numP++] = NP[ind_i];
                NP[ind_i] = NP[--numNP];
                --ind_i;
                
                assert(numP >= 0);
                assert(numP <= Npop);
                assert(numNP >= 0);
                assert(numNP <= Npop);
            }
        }
    }

    if (numP + numNP == 0)
    {
        write_parameters();
        exit(1);
    }

    // see if we need to switch worlds
    if (gsl_rng_uniform(rng_r) < s_12[current_world])
    {
        current_world = !current_world;
    }
}

// now reproduce
void reproduce()
{
    // number of offspring to be made
    int Noffspring = Npop - numP - numNP;

    assert(Noffspring >= 0);
    assert(Noffspring <= Npop);

    if (Noffspring == 0)
    {
        return;
    }

    Individual kids[Noffspring];

    // calculate fraction of adults in P
    double fraction_in_P = (double) numP / (numP + numNP);

    // go through all the offspring that need to
    // made
    for (int kid_i = 0; kid_i < Noffspring; ++kid_i)
    {
        // sample father
        Individual father;

        if (gsl_rng_uniform(rng_r) < fraction_in_P)
        {
            father = P[gsl_rng_uniform_int(rng_r, numP)];
        }
        else
        {
            father = NP[gsl_rng_uniform_int(rng_r, numNP)];
        }

        // sample mother
        Individual mother;

        if (gsl_rng_uniform(rng_r) < fraction_in_P)
        {
            mother = P[gsl_rng_uniform_int(rng_r, numP)];
        }
        else
        {
            mother = NP[gsl_rng_uniform_int(rng_r, numNP)];
        }

        // create the offspring
        create_offspring(mother, father, kids[kid_i]);

    }// end for kid_i

    // now distribute kids over the environments
    
    // calculate ratio of the P envt relative to the other envts
    double ratio_P_NP_envt = s_NP_2_P[current_world] / 
        (s_P_2_NP[current_world] + s_NP_2_P[current_world]);

    assert(ratio_P_NP_envt >= 0);
    assert(ratio_P_NP_envt <= 1);

    // now start redistributing kids
    for (int kid_i = 0; kid_i < Noffspring; ++kid_i)
    {
        // kid to envt P
        if (gsl_rng_uniform(rng_r) < ratio_P_NP_envt)
        {
            P[numP++] = kids[kid_i];
        }
        else
        {
            NP[numNP++] = kids[kid_i];
        }
    }

    assert(numP + numNP == Npop);
}


// write down summary statistics
void write_data()
{
    double mean_feedback = 0;
    double ss_feedback = 0;
    double mean_stress_influx = 0;
    double ss_stress_influx = 0;
    double mean_influx = 0;
    double ss_influx = 0;
    double mean_hormone = 0;
    double ss_hormone = 0;
    double mean_damage = 0;
    double ss_damage = 0;

    double freq_P = (double) numP / (numP + numNP);

    // get stats from the population
    for (int ind_i =  0; ind_i < numP; ++ind_i)
    {
        // feedback
        double feedback = 0.5 * (P[ind_i].feedback[0] + P[ind_i].feedback[1]);
        mean_feedback += feedback;
        ss_feedback += feedback * feedback;
        
        // stress_influx
        double stress_influx = 
            0.5 * (P[ind_i].stress_influx[0] + P[ind_i].stress_influx[1]);

        mean_stress_influx += stress_influx;
        ss_stress_influx += stress_influx * stress_influx;

        // influx
        double influx = 0.5 * (P[ind_i].influx[0] + P[ind_i].influx[1]);
        mean_influx += influx;
        ss_influx += influx * influx;
        
        double hormone = P[ind_i].hormone;
        mean_hormone += hormone;
        ss_hormone += hormone * hormone;
        
        double damage = P[ind_i].damage;
        mean_damage += damage;
        ss_damage += damage * damage;
    }
    
    // get stats from the population
    for (int ind_i =  0; ind_i < numNP; ++ind_i)
    {
        // feedback
        double feedback = 0.5 * (NP[ind_i].feedback[0] + NP[ind_i].feedback[1]);
        mean_feedback += feedback;
        ss_feedback += feedback * feedback;
        
        // stress_influx
        double stress_influx = 
            0.5 * (NP[ind_i].stress_influx[0] + NP[ind_i].stress_influx[1]);

        mean_stress_influx += stress_influx;
        ss_stress_influx += stress_influx * stress_influx;

        // influx
        double influx = 0.5 * (NP[ind_i].influx[0] + NP[ind_i].influx[1]);
        mean_influx += influx;
        ss_influx += influx * influx;
        
        double hormone = NP[ind_i].hormone;
        mean_hormone += hormone;
        ss_hormone += hormone * hormone;
        
        double damage = NP[ind_i].damage;
        mean_damage += damage;
        ss_damage += damage * damage;
    }

    mean_feedback /= numP + numNP;
    mean_stress_influx /= numP + numNP;
    mean_influx /= numP + numNP;
    mean_hormone /= numP + numNP;
    mean_damage /= numP + numNP;

    double var_feedback = ss_feedback / (numP + numNP) - pow(mean_feedback,2.0);
    double var_stress_influx = ss_stress_influx / (numP + numNP) - pow(mean_stress_influx,2.0);
    double var_influx = ss_influx / (numP + numNP) - pow(mean_influx,2.0);
    double var_hormone = ss_hormone / (numP + numNP) - pow(mean_hormone,2.0);
    double var_damage = ss_damage / (numP + numNP) - pow(mean_damage,2.0);

    DataFile << generation << ";" 
        << freq_P << ";"
        << mean_feedback << ";"
        << mean_stress_influx << ";"
        << mean_influx << ";"
        << mean_hormone << ";"
        << mean_damage << ";"
        << var_feedback << ";"
        << var_stress_influx << ";"
        << var_influx << ";"
        << var_hormone << ";"
        << var_damage << ";" << endl;
}

// write the headers of a datafile
void write_data_headers()
{
    DataFile << "generation;" 
        << "freq_P" << ";"
        << "mean_feedback" << ";"
        << "mean_stress_influx" << ";"
        << "mean_influx" << ";"
        << "mean_hormone" << ";"
        << "mean_damage" << ";"
        << "var_feedback" << ";"
        << "var_stress_influx" << ";"
        << "var_influx" << ";"
        << "var_hormone" << ";"
        << "var_damage" << ";" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	init_population();

    // write headers to the datafile
	write_data_headers();

    // the main part of the code
	for (generation = 0; generation <= NumGen; ++generation)
	{

		survive();

		reproduce();
        
        do_stats = generation % skip == 0;

        // print statistics every nth generation
        if (do_stats)
		{
			write_data();
		}
	}

    // finally write some params
	write_parameters();
}
