// the evolution of the stress response in a fluctuating envt
// Bram Kuijper
//

// The NDEBUG macro should be defined (before including <cassert>) if we DON'T
// want to do debugging (using the assert macro)
#define NDEBUG

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


// Note from Olle: using namespace std is often regarded as bad practice, but I
// have not changed it here
using namespace std;


// set up the random number generator using a good way of getting random seed
random_device rd;
unsigned seed = rd();
mt19937 rng_r(seed);
uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
bernoulli_distribution random_allele(0.5);


// the NumGen parameter is actually the number of time steps with overlapping
// generations; it is read from the command line
int NumGen = 1000;

// population size
const int Npop = 5000;

// number of generations to skip when outputting data
const int skip = 500;

// track population sizes in each of the environments; these variables are
// mainly for bookkeeping
int numP = 0;
int numNP = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
int generation = 0;

// mutation rate parameters
double mu_clearance = 0.0;
double mu_stress_influx = 0.0;
double mu_h1_S = 0.0;
double mu_h1_I = 0.0;
double mu_influx = 0.0;
double mu_hstart = 0.0;
double sdmu_clearance = 0.0;
double sdmu_stress_influx = 0.0;
double sdmu_h1_S = 0.0;
double sdmu_h1_I = 0.0;
double sdmu_influx = 0.0;
double sdmu_hstart = 0.0;

// background mortality
double mort_background = 0.0;

double opt_baseline = 0.3;

// switch rate from P to NP
double s_P_2_NP = 0.0;

// switch rate from NP to P
double s_NP_2_P = 0.0;

// probability of attack if predator is nearby
double p_att = 1.0;

// for convenience (and speed), introduce  (equilibrium) probabilities of the
// environment being P
double pr_envt_is_P = 0.0;

// the allelic values of clearance, stress_influx, influx and hstart are
// restricted to the unit interval; for this reason, zmax and dmax should be
// equal to 1.0

// initial values for the evolving traits
double init_clearance = 0.0;
double init_h1_S = 0.0;
double init_h1_I = 0.0;
double init_stress_influx = 0.0;
double init_influx = 0.0;
double init_hstart = 0.0;

// cue probabilities; this is the cue that is given before a predation
// potentially happens if the cue is equal in environments P or NP then it is
// uninformative
double cue_P = 0.0;
double cue_NP = 0.0;

// baseline survival rate
double s0 = 0;

double ad = 1.0;
double aP = 1.0;
double dmax = 1.0;
double zmax = 1.0;
double stress_influx_max = 0.25;
double min_clearance = 1.0;
double damage_clearance = 1.0; // strength of damage feedback
double u = 1.0; // how much damage each new load of hormone causes
int tmax_stress_influx = 5; // the maximum duration post stressor that a stress influx takes


// parameters for input/output of pop to text file; they are read from the
// command line
int ioind = 0; // 0 = none; 1 = output to file; 2 = both in/output
string base_name;

// mortalities across the two environments in the order NP, P
int Nmort_stats[2] = {0,0};

// keep track of number of living individuals in pop (for bookkeeping)
int num_alive = Npop;

// vector with doubles on cumulative damage levels
double fecundity_stats[2];


// the individual struct
struct Individual
{
    // diploid loci specifying the evolving traits
    // self-dependent increase/decrease in hormone
    double clearance[2];

    double stress_influx[2];

    // strength of neg feedback on stress influx
    double h1_S[2];
    
    // strength of neg feedback on normal influx
    double h1_I[2];

    // stress independent hormone influx
    double influx[2];

    // starting hormone level
    double hstart[2];

    // control
    double hx;

    // components of the individual's state
    double hormone; // current hormone level
    double damage;  // current damage level

    // whether individual is now living in P or not; this indicates the "true
    // situation", which need not be known by the individual
    bool envt_is_P;

    // indicator if the individual is alive
    bool alive;

    int time_since_last_stressor;
};


// storage of individuals is a single vector; each individual has indicators of
// their environment and if they are alive
std::vector<Individual> pop(Npop);
// fecundity values are stored in this variable
std::vector<double> fecundity(Npop);

// initialize simulations from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu_clearance = atof(argv[1]);
    mu_stress_influx  = atof(argv[2]);
    mu_h1_S  = atof(argv[3]);
    mu_h1_I  = atof(argv[4]);
    mu_influx  = atof(argv[5]);
    mu_hstart  = atof(argv[6]);

    sdmu_clearance = atof(argv[7]);
    sdmu_stress_influx = atof(argv[8]);
    sdmu_h1_S = atof(argv[9]);
    sdmu_h1_I = atof(argv[10]);
    sdmu_influx = atof(argv[11]);
    sdmu_hstart = atof(argv[12]);

    s_P_2_NP  = atof(argv[13]);
    s_NP_2_P  = atof(argv[14]);

    init_clearance = atof(argv[15]);
    init_stress_influx  = atof(argv[16]);
    init_h1_S  = atof(argv[17]);
    init_h1_I  = atof(argv[18]);
    init_influx  = atof(argv[19]);
    init_hstart  = atof(argv[20]);

    cue_P  = atof(argv[21]);
    cue_NP  = atof(argv[22]);

    s0  = atof(argv[23]);
    ad  = atof(argv[24]);
    aP  = atof(argv[25]);
    opt_baseline = atof(argv[26]);
    dmax  = atof(argv[27]);
    zmax  = atof(argv[28]);
    stress_influx_max = atof(argv[29]);
    min_clearance = atof(argv[30]);

    damage_clearance = atof(argv[31]);
    u = atof(argv[32]);

    mort_background = atof(argv[33]);
    p_att = atof(argv[34]);
    NumGen = atoi(argv[35]);
    tmax_stress_influx = atoi(argv[36]);
    ioind = atoi(argv[37]);
    base_name = argv[38];

    // equilibrium probabilities
    pr_envt_is_P = s_NP_2_P / (s_NP_2_P + s_P_2_NP);

} // end init_arguments


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

double g(double const h
        ,double const h1)
{
    if (h > h1)
    {
        return(0);
    }

    return(1.0 - h/h1);
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
		<< "mu_clearance;" << mu_clearance << ";"<< endl
		<< "mu_stress_influx;" << mu_stress_influx << ";"<< endl
		<< "mu_h1_S;" << mu_h1_S << ";"<< endl
		<< "mu_h1_I;" << mu_h1_I << ";"<< endl
        << "mu_influx;" << mu_influx << ";"<< endl
		<< "mu_hstart;" << mu_hstart << ";"<< endl
		<< "sP2NP_1;" << s_P_2_NP << ";"<< endl
		<< "sNP2P_1;" << s_NP_2_P << ";"<< endl
		<< "init_clearance;" << init_clearance << ";"<< endl
		<< "init_stress_influx;" << init_stress_influx << ";"<< endl
		<< "init_h1_S;" << init_h1_S << ";"<< endl
		<< "init_h1_I;" << init_h1_I << ";"<< endl
        << "init_influx;" << init_influx << ";"<< endl
		<< "init_hstart;" << init_hstart << ";"<< endl
        << "tmax_stress_influx;" << tmax_stress_influx << ";" << endl
		<< "cue_P;" << cue_P << ";"<< endl
		<< "cue_NP;" << cue_NP << ";"<< endl
        << "autocorr;" << 1.0 - (s_P_2_NP + s_NP_2_P) << ";" << endl
        << "risk;" << s_NP_2_P/(s_P_2_NP + s_NP_2_P) <<  ";" << endl
		<< "s0;" << s0 << ";"<< endl
		<< "ad;" << ad << ";"<< endl
		<< "aP;" << aP << ";"<< endl
		<< "opt_baseline;" << opt_baseline << ";"<< endl
		<< "mort_background;" << mort_background << ";"<< endl
		<< "dmax;" << dmax << ";"<< endl
		<< "zmax;" << zmax << ";"<< endl
		<< "stress_influx_max;" << stress_influx_max << ";"<< endl
		<< "min_clearance;" << min_clearance << ";"<< endl
		<< "damage_clearance;" << damage_clearance << ";"<< endl
        << "u;" << u << ";"<< endl
		<< "p_att;" << p_att << ";"<< endl
        << "NumGen;" << NumGen << ";"<< endl
        << "base_name;" << base_name << ";"<< endl;
}


// initialize the simulation by giving all the individuals genotypic values
void init_population()
{
    // set counters to 0
    numP = 0;
    numNP = 0;

	// initialize the whole population
	for (int i = 0; i < Npop; ++i)
	{
        Individual newInd;

        // set alleles
        for (int allele_i = 0; allele_i < 2; ++allele_i)
        {
            newInd.clearance[allele_i] = init_clearance;
            newInd.stress_influx[allele_i] = init_stress_influx;
            newInd.h1_S[allele_i] = init_h1_S;
            newInd.h1_I[allele_i] = init_h1_I;
            newInd.influx[allele_i] = init_influx;
            newInd.hstart[allele_i] = init_hstart;
        }

        // initialize hormone level and damage
        newInd.hormone = 0.5*(newInd.hstart[0] + newInd.hstart[1]);
        clamp(newInd.hormone, 0.0, zmax);
        newInd.damage = (damage_clearance > 0) ? u * newInd.hormone/damage_clearance : 1.0;
        clamp(newInd.damage, 0.0, dmax);

        newInd.hx = 0;

        // use equilibrium probability for P and NP
        newInd.envt_is_P = uniform(rng_r) < pr_envt_is_P;

        // make individuals alive and put in pop
        newInd.alive = true;
        pop[i] = newInd;

        newInd.time_since_last_stressor = tmax_stress_influx;

        if (newInd.envt_is_P)
        {
            ++numP;
        }
        else
        {
            ++numNP;
        }

        assert(numP + numNP <= Npop);
    }
} // end init_population


// create an offspring
void create_offspring(
        Individual &mother,
        Individual &father,
        Individual &kid)
{
    // assume allele 0 is from maternal gamete and allele 1 from paternal
    // gamete; the gametes are formed using free recombination of parental loci
    double mat_clearance = mother.clearance[random_allele(rng_r)];
    double pat_clearance = father.clearance[random_allele(rng_r)];
    double mat_stress_influx = mother.stress_influx[random_allele(rng_r)];
    double pat_stress_influx = father.stress_influx[random_allele(rng_r)];

    double mat_h1_S = mother.h1_S[random_allele(rng_r)];
    double pat_h1_S = father.h1_S[random_allele(rng_r)];

    double mat_h1_I = mother.h1_I[random_allele(rng_r)];
    double pat_h1_I = father.h1_I[random_allele(rng_r)];

    double mat_influx = mother.influx[random_allele(rng_r)];
    double pat_influx = father.influx[random_allele(rng_r)];
    double mat_hstart = mother.hstart[random_allele(rng_r)];
    double pat_hstart = father.hstart[random_allele(rng_r)];

    kid.clearance[0] = mat_clearance;
    kid.clearance[1] = pat_clearance;
    kid.stress_influx[0] = mat_stress_influx;
    kid.stress_influx[1] = pat_stress_influx;
    kid.h1_S[0] = mat_h1_S;
    kid.h1_S[1] = pat_h1_S;
    kid.h1_I[0] = mat_h1_I;
    kid.h1_I[1] = pat_h1_I;
    kid.influx[0] = mat_influx;
    kid.influx[1] = pat_influx;
    kid.hstart[0] = mat_hstart;
    kid.hstart[1] = pat_hstart;

    kid.time_since_last_stressor = tmax_stress_influx;

    // take into account mutation and note that allelic values are on logit
    // scale (and need not be clamped)
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        mutate(
                kid.clearance[allele_i],
                mu_clearance,
                sdmu_clearance);
        clamp(kid.clearance[allele_i], min_clearance, 1.0);

        mutate(
                kid.stress_influx[allele_i],
                mu_stress_influx,
                sdmu_stress_influx
                );
        clamp(kid.stress_influx[allele_i], 0.0, 1.0);
        
        mutate(
                kid.h1_S[allele_i],
                mu_h1_S,
                sdmu_h1_S
                );
        clamp(kid.h1_S[allele_i], 0.0, 1.0);
        
        mutate(
                kid.h1_I[allele_i],
                mu_h1_I,
                sdmu_h1_I
                );
        clamp(kid.h1_I[allele_i], 0.0, 1.0);

        mutate(
                kid.influx[allele_i],
                mu_influx,
                sdmu_influx);
        clamp(kid.influx[allele_i], 0.0, 1.0);

        mutate(
                kid.hstart[allele_i],
                mu_hstart,
                sdmu_hstart);
        clamp(kid.hstart[allele_i], 0.0, 1.0);
    }

    // set hormone level and damage
    kid.hormone = 0.5*(kid.hstart[0] + kid.hstart[1]);
    clamp(kid.hormone, 0.0, zmax);
    kid.damage = u * kid.hormone;
    clamp(kid.damage, 0.0, dmax);

    kid.hx = 0;

    // let kid start its life in "random" environment
    kid.envt_is_P = uniform(rng_r) < pr_envt_is_P;

    // set kid alive
    kid.alive = true;

    // update numP and numNP
    if (kid.envt_is_P) {
        ++numP;
    } else {
        ++numNP;
    }
}


// switch between environments
void environmental_switching()
{
    // just update envt_is_P field for each alive individual
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (pop[ind_i].alive) { // all individuals ought to be alive here
            if (pop[ind_i].envt_is_P) {
                if (uniform(rng_r) < s_P_2_NP) { // switch
                    pop[ind_i].envt_is_P = false;
                    --numP;
                    ++numNP;
                } // else no change
            } else {
                if (uniform(rng_r) < s_NP_2_P) { // switch
                    pop[ind_i].envt_is_P = true;
                    ++numP;
                    --numNP;
                } // else no change
            }
        }
    }
} //end void environmental_switching()


// survival of individuals
void survive(ofstream &datafile)
{
    // reset fecundity averages
    // that are tracked for stats
    fecundity_stats[0] = 0.0;
    fecundity_stats[1] = 0.0;

    // reset numbers of dead individuals
    // that are tracked for stats
    Nmort_stats[0] = 0;
    Nmort_stats[1] = 0;

    assert(numP >= 0);
    assert(numNP >= 0);
    assert(numP + numNP <= Npop);


    // run through alive individuals in pop and implement per time step
    // predation and background mortality
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (pop[ind_i].alive) { // all individuals ought to be alive here

            // increase time since last stressor counter
            ++pop[ind_i].time_since_last_stressor;

            // update hormone level from "background" influx and outflux
            double clearance = 0.5*(pop[ind_i].clearance[0] +
                              pop[ind_i].clearance[1]);

            double stress_influx = 0.5*(pop[ind_i].stress_influx[0] +
                                   pop[ind_i].stress_influx[1]);

            double h1_S = 0.5 * (pop[ind_i].h1_S[0] + pop[ind_i].h1_S[1]);
            
            double h1_I = 0.5 * (pop[ind_i].h1_I[0] + pop[ind_i].h1_I[1]);

            double influx = 0.5*(pop[ind_i].influx[0] +
                            pop[ind_i].influx[1]);


            pop[ind_i].hormone = (1.0 - clearance)*pop[ind_i].hormone + influx;

            double p_cue = pop[ind_i].envt_is_P ? cue_P : cue_NP;

            if (uniform(rng_r) < p_cue) {
                // individual gets predator cue
                //pop[ind_i].hormone += stress_influx;
                pop[ind_i].time_since_last_stressor = 0;
            }

            clamp(pop[ind_i].hormone, 0.0, zmax);

            // take into account possible predator attack
            if (pop[ind_i].envt_is_P) { // attack only possible if P
                if (uniform(rng_r) < p_att) { // predator attacks

                    pop[ind_i].time_since_last_stressor = 0;

                    double kill_prob = 1.0 - pow(pop[ind_i].hormone/zmax, aP);

                    if (uniform(rng_r) < kill_prob) {
                        // individual is killed by predator
                        pop[ind_i].alive = false;
                        --num_alive;
                        --numP;
                        ++Nmort_stats[1];
                        if (num_alive == 0)
                        {
                            cout << "extinct" << endl;
                            write_parameters(datafile);
                            exit(1);
                        }
                    } 
//                    else {
//                        // individual survives and gets hormone spike
//                        pop[ind_i].hormone += stress_influx;
//                        clamp(pop[ind_i].hormone, 0.0, zmax);
//                    }
                } // end if (uniform(rng_r) < p_att
            } // end if (pop[ind_i].envt_is_P

            if (pop[ind_i].time_since_last_stressor < tmax_stress_influx)
            {
                // if individuals are within the time interval [0,tmax_stress_influx)
                // then add stress influx and hormone-dependent stress influx
                pop[ind_i].hormone += stress_influx_max * g(pop[ind_i].hx, h1_S);
                clamp(pop[ind_i].hormone, 0.0, zmax);
                
                if (pop[ind_i].hx < pop[ind_i].hormone)
                {
                    pop[ind_i].hx = pop[ind_i].hormone;
                }
            }
        } // if (pop[ind_i].alive)


        if (pop[ind_i].alive) {
            // now take into account background mortality
            if (uniform(rng_r) < mort_background) {
                pop[ind_i].alive = false;
                --num_alive;
                if (pop[ind_i].envt_is_P) {
                    --numP;
                    ++Nmort_stats[1];
                } else {
                    --numNP;
                    ++Nmort_stats[0];
                }
                // dead individual, no fecundity
                fecundity[ind_i] = 0.0;
                if (num_alive == 0)
                {
                    cout << "extinct" << endl;
                    write_parameters(datafile);
                    exit(1);
                }
            } else {
                // update damage levels
                pop[ind_i].damage = (1.0 - damage_clearance) * pop[ind_i].damage +
                    u * pop[ind_i].hormone;

                clamp(pop[ind_i].damage,0.0,dmax);
                // damage-dependent fecundity
                fecundity[ind_i] = 1.0 - pow(pop[ind_i].damage/dmax,ad);
                fecundity_stats[pop[ind_i].envt_is_P] += fecundity[ind_i];
            }
        } else { // dead individual, no fecundity
            fecundity[ind_i] = 0.0;
        }
    }
    // make fecundity stats be per capita
    fecundity_stats[0] /= numNP;
    fecundity_stats[1] /= numP;
}


// fail-safe function to do reproduction
void reproduce_check(ofstream &datafile)
{
    // run through population and replace all dead individuals with offspring;
    // the parents for each offspring are randomly selected from the entire
    // population (this is a bit unrealistic); to select random parents, we use
    // a std::discrete_distribution<int> with fecundity vector as weights (note:
    // fecundity is zero for dead individuals)
    discrete_distribution<int> par_distr(fecundity.begin(), fecundity.end());
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (!pop[ind_i].alive) { // replace this dead individual
            int father = par_distr(rng_r);
            int mother = par_distr(rng_r);
            // create the offspring in the right position of pop
            create_offspring(pop[mother], pop[father], pop[ind_i]);
        }
    }
    num_alive = Npop; // all should now be alive
    assert(numP + numNP == Npop);
}


// write summary statistics
void write_data(ofstream &DataFile)
{
    double mean_clearance = 0;
    double ss_clearance = 0;
    double mean_stress_influx = 0;
    double ss_stress_influx = 0;
    double mean_h1_S = 0;
    double ss_h1_S = 0;
    double mean_h1_I = 0;
    double ss_h1_I = 0;
    double mean_influx = 0;
    double ss_influx = 0;
    double mean_hormone = 0;
    double ss_hormone = 0;
    double mean_damage = 0;
    double ss_damage = 0;

    double freq_P = (double) numP / (numP + numNP);

    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (pop[ind_i].alive) { // all should be alive at this point
            // Note from Olle: phenotypes are assumed to be mean of allelic
            // values (often sum of allelic values are assumed)
            // clearance
            double clearance = 0.5*(pop[ind_i].clearance[0] +
                              pop[ind_i].clearance[1]);

            mean_clearance += clearance;
            ss_clearance += clearance * clearance;
            // stress_influx
            double stress_influx = 0.5*(pop[ind_i].stress_influx[0] +
                                   pop[ind_i].stress_influx[1]);
            mean_stress_influx += stress_influx;
            ss_stress_influx += stress_influx * stress_influx;

            double h1_S = 0.5*(pop[ind_i].h1_S[0] +
                                   pop[ind_i].h1_S[1]);

            mean_h1_S += h1_S;
            ss_h1_S += h1_S * h1_S;
            
            double h1_I = 0.5*(pop[ind_i].h1_I[0] +
                                   pop[ind_i].h1_I[1]);

            mean_h1_I += h1_I;
            ss_h1_I += h1_I * h1_I;

            // influx
            double influx = 0.5*(pop[ind_i].influx[0] +
                            pop[ind_i].influx[1]);
            mean_influx += influx;
            ss_influx += influx * influx;
            // hormone
            double hormone = pop[ind_i].hormone;
            mean_hormone += hormone;
            ss_hormone += hormone * hormone;
            // damage
            double damage = pop[ind_i].damage;
            mean_damage += damage;
            ss_damage += damage * damage;
        }
    }
    mean_clearance /= Npop;
    mean_stress_influx /= Npop;
    mean_h1_S /= Npop;
    mean_h1_I /= Npop;
    mean_influx /= Npop;
    mean_hormone /= Npop;
    mean_damage /= Npop;

    double sd_clearance = sqrt(ss_clearance / Npop - pow(mean_clearance,2.0));
    double sd_stress_influx = sqrt(ss_stress_influx / Npop -
                                   pow(mean_stress_influx,2.0));
    double sd_h1_S = sqrt(ss_h1_S / Npop -
                                   pow(mean_h1_S,2.0));
    double sd_h1_I = sqrt(ss_h1_I / Npop -
                                   pow(mean_h1_I,2.0));
    double sd_influx = sqrt(ss_influx / Npop - pow(mean_influx,2.0));
    double sd_hormone = sqrt(ss_hormone / Npop - pow(mean_hormone,2.0));
    double sd_damage = sqrt(ss_damage / Npop - pow(mean_damage,2.0));

    DataFile << generation << ";"
        << freq_P << ";"
        << mean_clearance << ";"
        << mean_stress_influx << ";"
        << mean_h1_S << ";"
        << mean_h1_I << ";"
        << mean_influx << ";"
        << mean_hormone << ";"
        << mean_damage << ";"
        << sd_clearance << ";"
        << sd_stress_influx << ";"
        << sd_h1_S << ";"
        << sd_h1_I << ";"
        << sd_influx << ";"
        << sd_hormone << ";"
        << sd_damage << ";"
        << (double)Nmort_stats[1]/numP << ";"
        << (double)Nmort_stats[0]/numNP << ";"
        << (double)fecundity_stats[1] << ";"
        << (double)fecundity_stats[0] << ";"
        << endl;
}

// write the headers of a datafile
void write_data_headers(ofstream &DataFile)
{
    DataFile << "generation;"
        << "freq_P" << ";"
        << "mean_clearance " << ";"
        << "mean_stress_influx" << ";"
        << "mean_h1_S" << ";"
        << "mean_h1_I" << ";"
        << "mean_influx" << ";"
        << "mean_hormone" << ";"
        << "mean_damage" << ";"
        << "sd_clearance" << ";"
        << "sd_stress_influx" << ";"
        << "sd_h1_S" << ";"
        << "sd_h1_I" << ";"
        << "sd_influx" << ";"
        << "sd_hormone" << ";"
        << "sd_damage" << ";"
        << "prop_dead_P" << ";"
        << "prop_dead_NP" << ";"
        << "mean_fecundity_P" << ";"
        << "mean_fecundity_NP" << ";"
        << endl;
}

// iterate individuals for tmax timesteps
// to plot the stress response curve for
// different individuals
void write_simple_iter(ofstream &IterFile)
{
    // number of individuals
    int nrep = 20;
    int tmax = 500;
    int tstress = 100;

    // use an individual's genotype to compute its
    // hormone response to a single stressor event at time step tstress

    double hormone, hormone_tplus1;

    Individual ind;

    IterFile << "time;individual;hormone;" << endl;

    uniform_int_distribution<int> rint(0, Npop - 1);

    // sample individuals
    for (int ind_i = 0; ind_i < nrep; ++ind_i)
    {
        ind = pop[rint(rng_r)];

        // get clearance, stress_influx, influx and hstart from genotype
        double clearance = 0.5*(ind.clearance[0] + ind.clearance[1]);
        double stress_influx = 0.5*(ind.stress_influx[0] +
                               ind.stress_influx[1]);
        double h1_S = 0.5*(ind.h1_S[0] +
                               ind.h1_S[1]);
        double h1_I = 0.5*(ind.h1_I[0] +
                               ind.h1_I[1]);
        double influx = 0.5*(ind.influx[0] + ind.influx[1]);
        double hstart = 0.5*(ind.hstart[0] + ind.hstart[1]);

        // hormone baseline level
        hormone = (clearance > 0) ? influx/clearance : 1.0;
        hormone = hstart;
        clamp(hormone, 0.0, zmax);
        hormone_tplus1 = 0.0;

        IterFile << 0 << ";" << ind_i << ";" << hormone << ";" << endl;

        // iterate the stress response for this individual
        for (int timestep = 0; timestep < tmax; ++timestep)
        {
            hormone_tplus1 = (1.0 - clearance)*hormone + influx;

            if (timestep >= tstress && timestep < tstress + tmax_stress_influx)
            {
                hormone_tplus1 += stress_influx_max * g(hormone_tplus1, h1_S);
                
                clamp(hormone_tplus1, 0.0, zmax);
                
                if (pop[ind_i].hx < pop[ind_i].hormone)
                {
                    pop[ind_i].hx = pop[ind_i].hormone;
                }
            }

            clamp(hormone_tplus1, 0.0, zmax);
            hormone = hormone_tplus1;

            IterFile << timestep << ";" << (ind_i + 1) << ";"
                     << hormone << ";" << endl;
        }
    }
}

// code to read population from tab-separated text file
void read_pop_from_file(string infilename)
{
    ifstream infile(infilename.c_str());
    if (!infile) {
        std::cerr << "Could not open file " << infilename << '\n';
    } else {
        // skip first line in file (contains headers)
        char c = '\0';
        while (c != '\n' && infile) infile.get(c);
        // count individuals
        int ind_i = 0;
        Individual indi;
        double clearance; // dummy variable
        double stress_influx; // dummy variable
        double h1_S; // dummy variable
        double h1_I; // dummy variable
        double influx; // dummy variable
        double hstart; // dummy variable
        int envtP = 0;
        int alv = 0;
        while (infile && ind_i < Npop) {
            infile >> indi.clearance[0];
            infile >> indi.clearance[1];
            infile >> indi.stress_influx[0];
            infile >> indi.stress_influx[1];
            infile >> indi.h1_S[0];
            infile >> indi.h1_S[1];
            infile >> indi.h1_I[0];
            infile >> indi.h1_I[1];
            infile >> indi.influx[0];
            infile >> indi.influx[1];
            infile >> indi.hstart[0];
            infile >> indi.hstart[1];
            infile >> clearance;
            infile >> stress_influx;
            infile >> h1_S;
            infile >> h1_I;
            infile >> influx;
            infile >> hstart;
            infile >> indi.hormone;
            infile >> indi.damage;
            infile >> envtP;
            indi.envt_is_P = envtP > 0;
            infile >> alv;
            indi.alive = alv > 0;
            pop[ind_i] = indi;
            ++ind_i;
        }
        if (ind_i < Npop) {
            std::cerr << "Failed to read pop from " << infilename << "!\n";
        } else { // compute numP and numNP
            numP = 0;
            numNP = 0;
            for (int i = 0; i < Npop; ++i) {
                if (pop[i].envt_is_P) {
                    ++numP;
                } else {
                    ++numNP;
                }
            }
        }
        infile.close();
    }
}

// code to write population to tab-separated text file
void write_pop_to_file(ofstream &PopFile)
{
    // first write headers
    PopFile << "clearance_1" << "\t" << "clearance_2" << "\t"
        << "stress_influx_1" << "\t" << "stress_influx_2" << "\t"
        << "h1_S_1" << "\t" << "h1_S_2" << "\t"
        << "h1_I_1" << "\t" << "h1_I_2" << "\t"
        << "influx_1" << "\t" << "influx_2" << "\t"
        << "hstart_1" << "\t" << "hstart_2" << "\t"
        << "clearance" << "\t"
        << "stress_influx" << "\t"
        << "h1_S" << "\t"
        << "h1_I" << "\t"
        << "influx" << "\t"
        << "hstart" << "\t"
        << "hormone" << "\t"
        << "damage" << "\t"
        << "envt_is_P" << "\t"
        << "alive" << "\n";
    // then write population
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        double clearance = 0.5*(pop[ind_i].clearance[0] +
                          pop[ind_i].clearance[1]);
        double stress_influx = 0.5*(pop[ind_i].stress_influx[0] +
                               pop[ind_i].stress_influx[1]);
        double h1_S = 0.5*(pop[ind_i].h1_S[0] +
                               pop[ind_i].h1_S[1]);
        double h1_I = 0.5*(pop[ind_i].h1_I[0] +
                               pop[ind_i].h1_I[1]);
        double influx = 0.5*(pop[ind_i].influx[0] +
                        pop[ind_i].influx[1]);
        double hstart = 0.5*(pop[ind_i].hstart[0] +
                        pop[ind_i].hstart[1]);
        PopFile << pop[ind_i].clearance[0] << "\t"
                << pop[ind_i].clearance[1] << "\t"
                << pop[ind_i].stress_influx[0] << "\t"
                << pop[ind_i].stress_influx[1] << "\t"
                << pop[ind_i].h1_S[0] << "\t"
                << pop[ind_i].h1_S[1] << "\t"
                << pop[ind_i].h1_I[0] << "\t"
                << pop[ind_i].h1_I[1] << "\t"
                << pop[ind_i].influx[0] << "\t"
                << pop[ind_i].influx[1] << "\t"
                << pop[ind_i].hstart[0] << "\t"
                << pop[ind_i].hstart[1] << "\t"
                << clearance << "\t"
                << stress_influx << "\t"
                << h1_S << "\t"
                << h1_I << "\t"
                << influx << "\t"
                << hstart << "\t"
                << pop[ind_i].hormone << "\t"
                << pop[ind_i].damage << "\t"
                << pop[ind_i].envt_is_P << "\t"
                << pop[ind_i].alive << "\n";
    }
}

// the guts of the code
int main(int argc, char ** argv)
{
	init_arguments(argc, argv);

    // this code makes it possible to read pop from file
    if (ioind > 1) {
        // read pop from file
        read_pop_from_file(base_name + "pop.txt");
    } else {
        // initialize as before
        init_population();
    }

    // base_name for files
    ofstream DataFile(base_name.c_str());
    ofstream IterFile((base_name + "iters.csv").c_str());

    // write some params
	write_parameters(DataFile);

    // write headers to the datafile
	write_data_headers(DataFile);

    // the main part of the code
	for (generation = 0; generation <= NumGen; ++generation)
	{
        // generation is really time step, with overlapping generations
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

    if (ioind > 0) {
        // write pop to file
        ofstream PopFile((base_name + "pop.txt").c_str());
        write_pop_to_file(PopFile);
    }

    // iterate the stress response curves for a subset of individuals
    write_simple_iter(IterFile);
//    // finally write some params
//	write_parameters(DataFile);
}
