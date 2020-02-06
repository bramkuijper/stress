// the evolution of the stress response in a fluctuating envt
// Bram Kuijper
//

// Note from Olle: The main idea in this version of the program is to implement
// the traits in a different way: one trait is the same as the previous feedback
// (can be interpreted as the strength of a process breaking down hormone), the
// next trait is the proportion of the total influx capacity that is produced
// for the baseline (pr_baseline), and the third trait is the total influx
// capacity (max_influx).

// Note from Olle:
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


// set up the random number generator
// Note from Olle: below is a reasonable way of getting random seed
random_device rd;
unsigned seed = rd();
mt19937 rng_r(seed);
uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
bernoulli_distribution random_allele(0.5);


// Note from Olle: the NumGen parameter is actually the number of time steps
// with overlapping generations; it is now read from the command line
int NumGen = 1000;

// population size
const int Npop = 5000;

// number of generations to skip when outputting data
const int skip = 500;

// track population sizes in each of the environments
// Note from Olle: with my implementation, these variables are mainly for
// bookkeeping
int numP = 0;
int numNP = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
int generation = 0;

// mutation rate parameters
double mu_feedback = 0.0;
double mu_pr_baseline = 0.0;
double mu_max_influx = 0.0;
double sdmu = 0.0;

// background mortality
double mort_background = 0.0;

// switch rate from P to NP
double s_P_2_NP = 0.0;

// switch rate from NP to P
double s_NP_2_P = 0.0;

// probability of attack if predator is nearby
double p_att = 1.0;

// Note from Olle: for convenience (and speed), introduce  (equilibrium)
// probabilities of the environment being P
double pr_envt_is_P = 0.0;

// Note from Olle: one possibly useful change I did is to let the allelic
// values of feedback, pr_baseline and max_influx be on a logit scale; this
// restricts the values to the (0, 1) interval and might "help evolution" find
// equilibria where the actual values are close to 0 or close to 1; for this
// reason, zmax and dmax should be equal to 1.0

// initial values for the evolving traits
double init_lfeedback = 0.0;
double init_lpr_baseline = 0.0;
double init_lmax_influx = 0.0;

// cue probabilities
// this is the cue that is given before
// a predation potentially happens
// if the cue is equal in environments P or NP
// then it is uninformative
double cue_P = 0.0;
double cue_NP = 0.0;

// baseline survival rate
double s0 = 0;

double ad = 1.0;
double aP = 1.0;
double dmax = 1.0;
double zmax = 1.0;
double r = 0;
double u = 0;

// Note from Olle: here are parameters for input/output of pop to text file;
// they are read from the command line
int ioind = 0; // 0 = none; 1 = output to file; 2 = both in/output
string base_name;

// mortalities across the two environments
// Note from Olle: I took them in the order NP, P
int Nmort_stats[2] = {0,0};

// Note from Olle: keep track of number of living individuals in pop
// (for bookkeeping)
int num_alive = Npop;

// vector with doubles on cumulative damage levels
double fecundity_stats[2];


// the individual struct
struct Individual
{
    // assume alleles are on logit scale

    // diploid loci specifying the evolving traits
    // self-dependent increase/decrease in hormone
    // double feedback[2];
    double lfeedback[2];
    // proportion of max total influx produced for baseline
    double lpr_baseline[2];
    // maximum total hormone influx
    double lmax_influx[2];

    // components of the individual's state
    double hormone; // current hormone level
    double damage; // current damage level

    // whether individual is now living in P or not; this indicates the "true
    // situation", which need not be known by the individual
    bool envt_is_P;

    // indicator if the individual is alive
    bool alive;
};

// storage of individuals is a single vector; each individual has indicators of
// their environment and if they are alive
std::vector<Individual> pop(Npop);
// fecundity values are stored in this variable
std::vector<double> fecundity(Npop);

// initialize simulations from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu_feedback  = atof(argv[1]);
    mu_pr_baseline  = atof(argv[2]);
    mu_max_influx  = atof(argv[3]);
    sdmu = atof(argv[4]);
    s_P_2_NP  = atof(argv[5]);
    s_NP_2_P  = atof(argv[6]);
    init_lfeedback  = atof(argv[7]);
    init_lpr_baseline  = atof(argv[8]);
    init_lmax_influx  = atof(argv[9]);
    cue_P  = atof(argv[10]);
    cue_NP  = atof(argv[11]);
    s0  = atof(argv[12]);
    ad  = atof(argv[13]);
    aP  = atof(argv[14]);
    dmax  = atof(argv[15]);
    zmax  = atof(argv[16]);
    r = atof(argv[17]);
    u = atof(argv[18]);
    mort_background = atof(argv[19]);
    // Note from Olle: added these parameters
    p_att = atof(argv[20]);
    NumGen = atoi(argv[21]);
    ioind = atoi(argv[22]);
    base_name = argv[23];
    // equilibrium probabilities
    pr_envt_is_P = s_NP_2_P / (s_NP_2_P + s_P_2_NP);
}

// logistic function
double logistic(double val)
{
    return 1.0/(1.0 + std::exp(-val));
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
		<< "mu_stress_influx;" << mu_pr_baseline << ";"<< endl
		<< "mu_influx;" << mu_max_influx << ";"<< endl
		<< "sP2NP_1;" << s_P_2_NP << ";"<< endl
		<< "sNP2P_1;" << s_NP_2_P << ";"<< endl
		<< "init_lfeedback;" << init_lfeedback << ";"<< endl
		<< "init_lpr_baseline;" << init_lpr_baseline << ";"<< endl
		<< "init_lmax_influx;" << init_lmax_influx << ";"<< endl
		<< "cue_P;" << cue_P << ";"<< endl
		<< "cue_NP;" << cue_NP << ";"<< endl
		<< "s0;" << s0 << ";"<< endl
		<< "ad;" << ad << ";"<< endl
		<< "aP;" << aP << ";"<< endl
		<< "mort_background;" << mort_background << ";"<< endl
		<< "dmax;" << dmax << ";"<< endl
		<< "zmax;" << zmax << ";"<< endl
		<< "r;" << r << ";"<< endl
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
            newInd.lfeedback[allele_i] = init_lfeedback;
            newInd.lpr_baseline[allele_i] = init_lpr_baseline;
            newInd.lmax_influx[allele_i] = init_lmax_influx;
        }

        // initialize hormone level and damage to baseline values:
        double feedback =
            logistic(0.5*(newInd.lfeedback[0] + newInd.lfeedback[1]));
        double pr_baseline =
            logistic(0.5*(newInd.lpr_baseline[0] + newInd.lpr_baseline[1]));
        double max_influx =
            logistic(0.5*(newInd.lmax_influx[0] + newInd.lmax_influx[1]));
        double influx = pr_baseline*max_influx;
        newInd.hormone = (feedback > 0) ? influx/feedback : 1.0;
        clamp(newInd.hormone, 0.0, zmax);
        newInd.damage = (r > 0) ? u * newInd.hormone/r : 1.0;
        clamp(newInd.damage, 0.0, dmax);

        // use equilibrium probability for P and NP
        newInd.envt_is_P = uniform(rng_r) < pr_envt_is_P;

        // make individuals alive and put in pop
        newInd.alive = true;
        pop[i] = newInd;

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
    double mat_lfeedback = mother.lfeedback[random_allele(rng_r)];
    double pat_lfeedback = father.lfeedback[random_allele(rng_r)];
    double mat_lpr_baseline = mother.lpr_baseline[random_allele(rng_r)];
    double pat_lpr_baseline = father.lpr_baseline[random_allele(rng_r)];
    double mat_lmax_influx = mother.lmax_influx[random_allele(rng_r)];
    double pat_lmax_influx = father.lmax_influx[random_allele(rng_r)];

    kid.lfeedback[0] = mat_lfeedback;
    kid.lfeedback[1] = pat_lfeedback;
    kid.lpr_baseline[0] = mat_lpr_baseline;
    kid.lpr_baseline[1] = pat_lpr_baseline;
    kid.lmax_influx[0] = mat_lmax_influx;
    kid.lmax_influx[1] = pat_lmax_influx;

    // take into account mutation and note that allelic values are on logit
    // scale (and need not be clamped)
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        mutate(
                kid.lfeedback[allele_i],
                mu_feedback,
                sdmu);

        mutate(
                kid.lpr_baseline[allele_i],
                mu_pr_baseline,
                sdmu
                );

        mutate(
                kid.lmax_influx[allele_i],
                mu_max_influx,
                sdmu);
    }

    // set hormone level and damage to their baseline values
    double feedback =
        logistic(0.5*(kid.lfeedback[0] + kid.lfeedback[1]));
    double pr_baseline =
        logistic(0.5*(kid.lpr_baseline[0] + kid.lpr_baseline[1]));
    double max_influx =
        logistic(0.5*(kid.lmax_influx[0] + kid.lmax_influx[1]));
    double influx = pr_baseline*max_influx;
    kid.hormone = (feedback > 0) ? influx/feedback : 1.0;
    clamp(kid.hormone, 0.0, zmax);
    kid.damage = u * kid.hormone;
    clamp(kid.damage, 0.0, dmax);

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
    // sum_fecundity = 0.0;

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
            // update hormone level from "background" influx and outflux
            double feedback =
                logistic(0.5*(pop[ind_i].lfeedback[0] +
                              pop[ind_i].lfeedback[1]));
            double pr_baseline =
                logistic(0.5*(pop[ind_i].lpr_baseline[0] +
                              pop[ind_i].lpr_baseline[1]));
            double max_influx =
                logistic(0.5*(pop[ind_i].lmax_influx[0] +
                              pop[ind_i].lmax_influx[1]));
            double influx = pr_baseline*max_influx;
            double stress_influx = (1.0 - pr_baseline)*max_influx;
            pop[ind_i].hormone = (1.0 - feedback)*pop[ind_i].hormone + influx;
            double p_cue = pop[ind_i].envt_is_P ? cue_P : cue_NP;
            if (uniform(rng_r) < p_cue) {
                // individual gets predator cue
                pop[ind_i].hormone += stress_influx;
            }
            clamp(pop[ind_i].hormone, 0.0, zmax);
            // take into account possible predator attack
            if (pop[ind_i].envt_is_P) { // attack only possible if P
                if (uniform(rng_r) < p_att) { // predator attacks
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
                    } else {
                        // individual survives and gets hormone spike
                        pop[ind_i].hormone += stress_influx;
                        clamp(pop[ind_i].hormone, 0.0, zmax);
                    }
                }
            }
        }
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
                pop[ind_i].damage = (1.0 - r) * pop[ind_i].damage +
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

    // Note from Olle: changed to run through pop
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (pop[ind_i].alive) { // all should be alive at this point
            // Note from Olle: phenotypes are assumed to be mean of allelic
            // values (often sum of allelic values are assumed)
            // feedback
            double feedback =
                logistic(0.5*(pop[ind_i].lfeedback[0] +
                              pop[ind_i].lfeedback[1]));
            double pr_baseline =
                logistic(0.5*(pop[ind_i].lpr_baseline[0] +
                              pop[ind_i].lpr_baseline[1]));
            double max_influx =
                logistic(0.5*(pop[ind_i].lmax_influx[0] +
                              pop[ind_i].lmax_influx[1]));
            mean_feedback += feedback;
            ss_feedback += feedback * feedback;
            // stress_influx
            double stress_influx = (1.0 - pr_baseline)*max_influx;
            mean_stress_influx += stress_influx;
            ss_stress_influx += stress_influx * stress_influx;
            // influx
            double influx = pr_baseline*max_influx;
            mean_influx += influx;
            ss_influx += influx * influx;
            // feedback
            double hormone = pop[ind_i].hormone;
            mean_hormone += hormone;
            ss_hormone += hormone * hormone;
            // damage
            double damage = pop[ind_i].damage;
            mean_damage += damage;
            ss_damage += damage * damage;
        }
    }
    mean_feedback /= Npop;
    mean_stress_influx /= Npop;
    mean_influx /= Npop;
    mean_hormone /= Npop;
    mean_damage /= Npop;

    double sd_feedback = sqrt(ss_feedback / Npop - pow(mean_feedback,2.0));
    double sd_stress_influx = sqrt(ss_stress_influx / Npop -
                                   pow(mean_stress_influx,2.0));
    double sd_influx = sqrt(ss_influx / Npop - pow(mean_influx,2.0));
    double sd_hormone = sqrt(ss_hormone / Npop - pow(mean_hormone,2.0));
    double sd_damage = sqrt(ss_damage / Npop - pow(mean_damage,2.0));

    DataFile << generation << ";"
        << freq_P << ";"
        << mean_feedback << ";"
        << mean_stress_influx << ";"
        << mean_influx << ";"
        << mean_hormone << ";"
        << mean_damage << ";"
        << sd_feedback << ";"
        << sd_stress_influx << ";"
        << sd_influx << ";"
        << sd_hormone << ";"
        << sd_damage << ";"
        << (double)Nmort_stats[1]/numP << ";"
        << (double)Nmort_stats[0]/numNP << ";"
        << (double)fecundity_stats[1] << ";"
        << (double)fecundity_stats[0] << ";"
        << endl;
    // Note from Olle: changed interpretation of fecundity stats (seems I got
    // it wrong somehow)
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
        << "sd_feedback" << ";"
        << "sd_stress_influx" << ";"
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
    int nrep = 100;
    int tmax = 500;
    int tstress = 100;

    // Note from Olle: I only use an individual's genotype to compute its
    // hormone response to a single stressor event at time step tstress; also,
    // I assume the hormone is at the baseline equilibrium level at the start

    double hormone, hormone_tplus1;

    Individual ind;

    // calculate freq of individuals in
    // predator vs non-predator patch
    // double freq_p = (double) numP / (numP + numNP);

    IterFile << "time;individual;hormone;" << endl;

    // uniform_int_distribution<int> random_from_P(0, numP - 1);
    // uniform_int_distribution<int> random_from_NP(0, numNP - 1);
    uniform_int_distribution<int> rint(0, Npop - 1);

    // sample individuals
    for (int ind_i = 0; ind_i < nrep; ++ind_i)
    {
        // if (uniform(rng_r) < freq_p)
        // {
        //     ind = P[random_from_P(rng_r)];
        // }
        // else
        // {
        //     ind = NP[random_from_NP(rng_r)];
        // }

        ind = pop[rint(rng_r)];

        // assert(ind.feedback[0] >= 0.0);
        // assert(ind.feedback[1] <= 1.0);

        // get feedback, pr_baseline and max_influx from genotype
        double feedback =
            logistic(0.5*(ind.lfeedback[0] +
                          ind.lfeedback[1]));
        double pr_baseline =
            logistic(0.5*(ind.lpr_baseline[0] +
                          ind.lpr_baseline[1]));
        double max_influx =
            logistic(0.5*(ind.lmax_influx[0] +
                          ind.lmax_influx[1]));

        double influx = pr_baseline*max_influx;
        double stress_influx = (1.0 - pr_baseline)*max_influx;

        // hormone baseline level
        hormone = (feedback > 0) ? influx/feedback : 1.0;
        clamp(hormone, 0.0, zmax);
        hormone_tplus1 = 0.0;

        IterFile << 0 << ";" << ind_i << ";" << hormone << ";" << endl;

        // iterate the stress response for this individual
        for (int timestep = 0; timestep < tmax; ++timestep)
        {
            hormone_tplus1 = (1.0 - feedback)*hormone + influx;
            if (timestep == tstress)
            {
                hormone_tplus1 += stress_influx;
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
        double feedback; // dummy variable
        double stress_influx; // dummy variable
        double influx; // dummy variable
        int envtP = 0;
        int alv = 0;
        while (infile && ind_i < Npop) {
            infile >> indi.lfeedback[0];
            infile >> indi.lfeedback[1];
            infile >> indi.lpr_baseline[0];
            infile >> indi.lpr_baseline[1];
            infile >> indi.lmax_influx[0];
            infile >> indi.lmax_influx[1];
            infile >> feedback;
            infile >> stress_influx;
            infile >> influx;
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
    PopFile << "lfeedback1" << "\t" << "lfeedback2" << "\t"
        << "lmax_influx1" << "\t" << "lmax_influx2" << "\t"
        << "linflux1" << "\t" << "linflux2" << "\t"
        << "feedback" << "\t"
        << "stress_influx" << "\t"
        << "influx" << "\t"
        << "hormone" << "\t"
        << "damage" << "\t"
        << "envt_is_P" << "\t"
        << "alive" << "\n";
    // then write population
    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        double feedback =
            logistic(0.5*(pop[ind_i].lfeedback[0] +
                          pop[ind_i].lfeedback[1]));
        double pr_baseline =
            logistic(0.5*(pop[ind_i].lpr_baseline[0] +
                          pop[ind_i].lpr_baseline[1]));
        double max_influx =
            logistic(0.5*(pop[ind_i].lmax_influx[0] +
                          pop[ind_i].lmax_influx[1]));
        double stress_influx = (1.0 - pr_baseline)*max_influx;
        double influx = pr_baseline*max_influx;
        PopFile << pop[ind_i].lfeedback[0] << "\t"
                << pop[ind_i].lfeedback[1] << "\t"
                << pop[ind_i].lpr_baseline[0] << "\t"
                << pop[ind_i].lpr_baseline[1] << "\t"
                << pop[ind_i].lmax_influx[0] << "\t"
                << pop[ind_i].lmax_influx[1] << "\t"
                << feedback << "\t"
                << stress_influx << "\t"
                << influx << "\t"
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
        // Note from Olle: read pop from file
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
