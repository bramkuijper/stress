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
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <vector>
#include <random>


// various functions, such as unique filename creation
// which helps when running multiple instances of simulation
// in the same folder
#include "auxiliary.hpp"

//#define NDEBUG
//
// the debug sign should only be turned on when one wants
// to assess the complete distribution of phenotypes
// 

using namespace std;


// set up the random number generator
// set random seed etc
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
bernoulli_distribution random_allele(0.5);


// number of generations
const long int NumGen = 100000;

// population size
const int Npop = 5000; 

// number of generations to skip when outputting data
const long int skip = 100;

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

// background mortality
double mort_background = 0.0;

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

double ad = 0;
double aP = 0;
double dmax = 0;
double zmax = 0;
double r = 0;
double u = 0;
int death_t = 0;

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

// vector with doubles on cumulative damage levels
double fecundity_cumul[Npop];
double sum_fecundity;

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
    mort_background = atof(argv[23]);
}

// apply boundaries to a certain value
void clamp(double &val, double const min, double const max)
{
    if (val > max)
    {
        val = max;
    } 
    else if (val < min)
    {
        val = min;
    }
}

// mutation according to a continuum of alleles model
void mutate(
        double &G, 
        double const mu, 
        double const sdmu)
{

    if (uniform(rng_r) < mu)
    {
        normal_distribution<> allele_dist(0.0, sdmu);
        G += allele_dist(rng_r);
    }
}

// write the parameters 
// (typically at the end of the output file)
void write_parameters(ofstream &DataFile)
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
		<< "mort_background;" << mort_background << ";"<< endl
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
        newInd.envt_is_P = uniform(rng_r) < s_P_2_NP[0] / 
            (s_NP_2_P[0] + s_P_2_NP[0]);

        if (newInd.envt_is_P)
        {
            P[numP++] = newInd;
        }
        else
        {
            NP[numNP++] = newInd;
        }

        assert(numP + numNP <= Npop);
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
            mother.feedback[random_allele(rng_r)] // mother inherits 
            :
            father.feedback[random_allele(rng_r)]; // father inherits

        assert(kid.feedback[allele_i] >= 0.0);
        assert(kid.feedback[allele_i] <= 1.0);

        mutate(
                kid.feedback[allele_i],
                mu_feedback,
                sdmu);

        clamp(kid.feedback[allele_i], 0.0, 1.0);

        kid.stress_influx[allele_i] = 
            allele_i < 1 ? 
            mother.stress_influx[random_allele(rng_r)] 
            :
            father.stress_influx[random_allele(rng_r)];

        mutate(
                kid.stress_influx[allele_i],
                mu_stress_influx, 
                sdmu
                );

        clamp(kid.stress_influx[allele_i], 0.0, DBL_MAX);

        kid.influx[allele_i] = 
            allele_i < 1 ? 
            mother.influx[random_allele(rng_r)] 
            :
            father.influx[random_allele(rng_r)];

        mutate(
                kid.influx[allele_i],
                mu_influx, 
                sdmu);
        
        clamp(kid.influx[allele_i], 0.0, DBL_MAX);
    } // inheritance done

    // set hormone level and damage to their baseline values
    kid.hormone = 0.0;
    kid.damage = 0.0;
}

// calculate the reduction in fecundity due to damage
double fecundity_damage(double const damage)
{
    double val = 1.0 - pow(damage/dmax,ad);

    if (val <= 0.0)
    {
        val = 0.0;
    }

    return(val);
}

// calculate mortality probability
// dependent on 
// - whether or not individual is in environment P
// - its current hormone level
double pkill(double const hormone_level, bool envt_is_P)
{
    double kill_prob = mort_background;
   
    // in predator environment, mortality is relative to hormone level
    if (envt_is_P)
    {
        kill_prob += (1.0 - mort_background) * (1.0 - pow(hormone_level/zmax, aP));
    }

    assert(kill_prob >= 0);
    assert(kill_prob <= 1.0);

    return(kill_prob);
}

// move between environments
void environmental_switching()
{
    // counters for individuals that move to a different
    // patch due to environmental change, with switch rates
    // s_P_2_NP and s_NP_2_P
    int numNPnew = 0;

    for (int ind_i = 0; ind_i < numP; ++ind_i)
    {
        if (uniform(rng_r)  < s_P_2_NP[current_world])
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
            assert(numNPnew >= 0);
            assert(numNP <= Npop);
        }
    }

    for (int ind_i = 0; ind_i < numNP; ++ind_i)
    {
        if (uniform(rng_r)  < s_NP_2_P[current_world])
        {
            // copy individual to NPnew stack which is a 
            // placeholder for all individuals who go to NP
            //
            // we cannot directly transfer individuals to the NP
            // stack itself, as the NP stack still needs to undergo
            // survival selection.
            P[numP++] = NP[ind_i];
            NP[ind_i] = NP[--numNP];
            --ind_i;
            assert(numP >= 0);
            assert(numP <= Npop);
            assert(numNP >= 0);
            assert(numNP <= Npop);
        }
    }

    // now copy NPnew to NP
    for (int ind_i = 0; ind_i < numNPnew; ++ind_i)
    {
        NP[numNP++] = NPnew[ind_i];

        assert(numNP + numP >= 0);
        assert(numNP + numP <= Npop);
    }
} //end void environmental_switching()


// survival of individuals
void survive(ofstream &datafile)
{
    sum_fecundity = 0.0;

    death_t = 0;

    assert(numP >= 0);
    assert(numNP >= 0);
    assert(numP + numNP <= Npop);

    // survival in the P population
    for (int ind_i = 0; ind_i < numP; ++ind_i)
    {
        // standard influx and outflux in hormone level
        // z(t+1) = influx + (1 - feedback) * z(t)
        P[ind_i].hormone = 
            0.5 * (P[ind_i].influx[0] + P[ind_i].influx[1])
             + (1.0 - 0.5 * (P[ind_i].feedback[0] + P[ind_i].feedback[1]))
             * P[ind_i].hormone;

        assert(P[ind_i].influx[0] >= 0.0);        
        assert(P[ind_i].influx[1] >= 0.0);        
        
        assert(P[ind_i].stress_influx[0] >= 0.0);        
        assert(P[ind_i].stress_influx[1] >= 0.0);        

        // individual gets predator cue 
        if (uniform(rng_r) < cue_P)
        {
            // gets cue, spike the hormone level
            P[ind_i].hormone += 
                0.5 * (P[ind_i].stress_influx[0] + P[ind_i].stress_influx[1]);
        }

        clamp(P[ind_i].hormone, 0.0, zmax);
        
        // OK, did not survive dependent on the level of stress & damage
        if (uniform(rng_r) < 
                pkill(P[ind_i].hormone, true))
        {
            ++death_t;
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
            // update damage levels
            P[ind_i].damage = (1.0 - r) * P[ind_i].damage + u * P[ind_i].hormone;

            // set boundaries to damage level
            assert(P[ind_i].damage >= 0);

            clamp(P[ind_i].damage,0.0,dmax);

            // add to cumulative distribution of damage
            fecundity_cumul[ind_i] = sum_fecundity + fecundity_damage(P[ind_i].damage);
            sum_fecundity = fecundity_cumul[ind_i];
        }
    } // end for (int ind_i = 0; ind_i < numP; ++ind_i)

    // survival in the NP population
    for (int ind_i = 0; ind_i < numNP; ++ind_i)
    {
        // standard influx and outflux in hormone level
        // z(t+1) = influx + (1 - feedback) * z(t)
        NP[ind_i].hormone = 
            0.5 * (NP[ind_i].influx[0] + NP[ind_i].influx[1])
             + (1.0 - 0.5 * (NP[ind_i].feedback[0] + NP[ind_i].feedback[1]))
             * NP[ind_i].hormone;

        assert(NP[ind_i].influx[0] >= 0.0);        
        assert(NP[ind_i].influx[1] >= 0.0);        
        
        assert(NP[ind_i].stress_influx[0] >= 0.0);        
        assert(NP[ind_i].stress_influx[1] >= 0.0);        

        // individual gets predator cue 
        if (uniform(rng_r) < cue_NP)
        {
            // gets cue, spike the hormone level
            NP[ind_i].hormone += 
                0.5 * (NP[ind_i].stress_influx[0] + NP[ind_i].stress_influx[1]);
        }

        clamp(P[ind_i].hormone, 0.0, zmax);

        // individual did not survive 
        if (uniform(rng_r) < 
                pkill(NP[ind_i].hormone, false))
        {
            ++death_t;
            // delete individual from stack
            NP[ind_i] = NP[--numNP];
            --ind_i;

            assert(numNP >= 0);
            assert(numNP <= Npop);
        }
        else // Individual survived. Check for potential for envt'al change
        {
            // update damage levels
            NP[ind_i].damage = (1.0 - r) * NP[ind_i].damage + u * NP[ind_i].hormone;

            // set boundaries to damage level
            assert(NP[ind_i].damage >= 0);

            clamp(NP[ind_i].damage, 0.0, dmax);

            if (numP > 0)
            {
                assert(fecundity_cumul[numP + ind_i - 1] == sum_fecundity);
            }

            // add to cumulative distribution of damage
            fecundity_cumul[numP + ind_i] = sum_fecundity + fecundity_damage(NP[ind_i].damage);
            sum_fecundity = fecundity_cumul[numP + ind_i];
        }
    } // end for for (int ind_i = 0; ind_i < numNP; ++ind_i)

    assert(numP >= 0);
    assert(numNP >= 0);
    assert(numP + numNP <= Npop);

    if (numP + numNP == 0)
    {
        cout << "extinct" << endl;
        write_parameters(datafile);
        exit(1);
    }

    // see if we need to switch worlds
    if (uniform(rng_r) < s_12[current_world])
    {
        current_world = !current_world;
    }
}


// check reproduction
void reproduce_check(ofstream &datafile)
{
    assert(numP >= 0);
    assert(numNP >= 0);
    assert(numP + numNP > 0);
    assert(numP + numNP <= Npop);

    // number of offspring to be made
    int Noffspring = Npop - numP - numNP;

    assert(Noffspring >= 0);
    assert(Noffspring <= Npop);

    if (Noffspring == 0)
    {
        return;
    }

    Individual kids[Noffspring];

    double cumul_deviate = 0;

    int father, mother;

    if (numP + numNP < 1)
    {
        write_parameters(datafile);
        exit(1);
    }   

    uniform_int_distribution<int> random_parent(0, numP + numNP - 1);

    for (int kid_i = 0; kid_i < Noffspring; ++kid_i)
    {
        father = random_parent(rng_r);
        mother = random_parent(rng_r);

        assert(father < Npop);
        assert(father >= 0);
        assert(mother < Npop);
        assert(mother >= 0);

        // first obtain father from cumul dist
        cumul_deviate = uniform(rng_r) * sum_fecundity;

        for (int ind_i = 0; ind_i < numP + numNP; ++ind_i)
        {
            assert(ind_i >= 0);
            assert(ind_i < numP + numNP);
            assert(ind_i < Npop);

            if (cumul_deviate < fecundity_cumul[ind_i])
            {
                father = ind_i;
                break;
            }
        }
        
        // then obtain father from cumul dist
        cumul_deviate = uniform(rng_r) * sum_fecundity;

        assert(cumul_deviate >= 0);
        assert(cumul_deviate <= sum_fecundity);

        for (int ind_i = numP; ind_i < numP + numNP; ++ind_i)
        {
            assert(ind_i >= 0);
            assert(ind_i < numP + numNP);
            assert(ind_i < Npop);

            if (cumul_deviate < fecundity_cumul[ind_i])
            {
                mother = ind_i;
                break;
            }
        }

        assert(father >= 0);
        assert(father < numP + numNP);
        
        assert(mother >= 0);
        assert(mother < numP + numNP);
        
        Individual father_ind;

        if (father < numP)
        {
            father_ind = P[father];
        }
        else
        {
            assert(father - numP >= 0);
            assert(father - numP < numNP);
            father_ind = NP[father - numP];
        }
        
        
        Individual mother_ind;

        if (mother < numP)
        {
            mother_ind = P[mother];
        }
        else
        {
            assert(mother - numP >= 0);
            assert(mother - numP < numNP);
            mother_ind = NP[mother - numP];
        }
        
        // create the offspring
        create_offspring(mother_ind, father_ind, kids[kid_i]);
    } // end for (int kid_i = 0; kid_i < Noffspring; ++kid_i)
    
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
        if (uniform(rng_r) < ratio_P_NP_envt)
        {
            P[numP++] = kids[kid_i];
        }
        else
        {
            NP[numNP++] = kids[kid_i];
        }
    } // end for (int kid_i = 0; kid_i < Noffspring; ++kid_i)

    assert(numP + numNP == Npop);
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

    int Nparents = 2 * Noffspring;

    // allocate array with random samples from cumulative distribution
    vector <double> cumul_dist_samples(Nparents, 0.0);

    // store the parents that result from the sample (because random shuffle)
    vector <int> parents(Nparents, 0);

    // go through all the offspring that need to
    // made
    for (int parent_i = 0; parent_i < Nparents; ++parent_i)
    {
        // now add random numbers from the cumulative distribution to the list
        // of samples that need to be found in the cumulative distribution
        cumul_dist_samples[parent_i] = uniform(rng_r) * sum_fecundity;
    }

    cout << "american airpowert made the difference" << endl;

    // sort cumulative samples from low to high
    sort(cumul_dist_samples.begin(), cumul_dist_samples.end());

    int cumul_counter = 0;

    // associate random numbers with population
    for (int ind_i = 0; ind_i < numP + numNP; ++ind_i)
    {
        for (; cumul_counter < Nparents; ++cumul_counter)
        {
            // ok random deviate lower than current percentile
            if (cumul_dist_samples[cumul_counter] <= fecundity_cumul[ind_i])
            {
                parents[cumul_counter] = ind_i;
            }
            else
            {
                // break cumul counter for loop
                // and iterate to the next deviate of the cumulative distribution
                break;
            }
        }

        if (cumul_counter >= Nparents)
        {
            break;
        }
    } // end for (int ind_i = 0

    assert(cumul_counter >= Nparents);

    // now shuffle the female dispersers
    shuffle(
            parents.begin(),
            parents.end(),
            rng_r);

    for (int kid_i = 0; kid_i < Noffspring; ++kid_i)
    {
        // get father
        assert(parents[2 * kid_i] >= 0);
        assert(parents[2 * kid_i] <= Npop);

        Individual father;

        if (parents[2*kid_i] < numP)
        {
            father = P[parents[2*kid_i]];
        }
        else
        {
            father = NP[parents[2*kid_i] - numP];
        }
        
        
        // get father
        assert(parents[2 * kid_i + 1] >= 0);
        assert(parents[2 * kid_i + 1] <= Npop);


        // sample mother
        Individual mother;
        
        if (parents[2*kid_i + 1] < numP)
        {
            mother = P[parents[2*kid_i + 1]];
        }
        else
        {
            mother = NP[parents[2*kid_i + 1] - numP];
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
        if (uniform(rng_r) < ratio_P_NP_envt)
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
void write_data(ofstream &DataFile)
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
        << var_damage << ";" 
        << (double)death_t/Npop << ";" 
        << endl;
}

// write the headers of a datafile
void write_data_headers(ofstream &DataFile)
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
        << "var_damage" << ";" 
        << "prop_dead" << ";"
        << endl;
}

// iterate individuals for tmax timesteps 
// to plot the stress response curve for
// different individuals
void write_simple_iter(ofstream &IterFile)
{
    // number of individuals
    int nrep = 100;
    int tmax = 500;
    int tstress = tmax - 100;

    double stress, stress_tplus1;

    Individual ind;

    // calculate freq of individuals in 
    // predator vs nonpredator patch
    double freq_p = (double) numP / (numP + numNP);

    IterFile << "time;individual;stress;" << endl;

    uniform_int_distribution<int> random_from_P(0, numP - 1);
    uniform_int_distribution<int> random_from_NP(0, numNP - 1);

    // sample individuals
    for (int ind_i = 0; ind_i < nrep; ++ind_i)
    {
        if (uniform(rng_r) < freq_p)
        {
            ind = P[random_from_P(rng_r)];
        }
        else
        {
            ind = NP[random_from_NP(rng_r)];
        }

        assert(ind.feedback[0] >= 0.0);
        assert(ind.feedback[1] <= 1.0);

        stress = 0.0;
        stress_tplus1 = 0.0;

        // iterate the stress response for this individual
        for (int timestep = 0; timestep < tmax; ++timestep)
        {
            stress_tplus1 = 0.5 * (ind.influx[0] + ind.influx[1]) 
                                + (1.0 - 0.5 * (ind.feedback[0] + ind.feedback[1]))
                                * stress;
            
            if (timestep == tmax - tstress)
            {
                stress_tplus1 += 0.5 * (ind.stress_influx[0] + ind.stress_influx[1]);
            }

            stress = stress_tplus1;

            IterFile << timestep << ";" << ind_i << ";" << stress << ";" << endl;
        }
    }
}


// the guts of the code
int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	init_population();

    // generate a unique filename for the output file
    string filename("sim_stress");
    create_filename(filename);
    ofstream DataFile(filename.c_str());  

    string filename_iters(filename + "iters");
    ofstream IterFile(filename_iters.c_str());  

    
    // write some params
	write_parameters(DataFile);

    // write headers to the datafile
	write_data_headers(DataFile);

    // the main part of the code
	for (generation = 0; generation <= NumGen; ++generation)
	{
        environmental_switching();

		survive(DataFile);

		reproduce_check(DataFile);
        
        do_stats = generation % skip == 0;

        // print statistics every nth generation
        if (do_stats)
		{
			write_data(DataFile);
		}
	}

    // iterate the stress response curves for a subset of individuals
    write_simple_iter(IterFile);
//    // finally write some params
//	write_parameters(DataFile);
}
