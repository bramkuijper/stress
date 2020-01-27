// the evolution of the stress response in a fluctuating envt
// Bram Kuijper
//

// Note from Olle: in this version I removed many commented away sections and
// introduced separate stress_influx (response to actual predator attack) and
// cue_influx (response to predator cue)

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

std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);
std::uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
std::bernoulli_distribution random_allele(0.5);


// maximum number of timesteps of the simulation
int maxtime = 30;

// population size
const int Npop = 5000;

// number of generations to skip when outputting data
const int skip = 500;

// track population sizes in each of the environments (bookkeeping)
int numP = 0;
int numNP = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

// keep track of the current generation number
int generation = 0;

// mutation rate parameters
double mu_feedback = 0.0;
double mu_cue_influx = 0.0;
double mu_stress_influx = 0.0;
double mu_influx = 0.0;
double sdmu = 0.0;

// background mortality
double mort_background = 0.0;

// environmental switching parameters
// switch rate from P to NP
double s_P_2_NP = 0.0;
// switch rate from NP to P
double s_NP_2_P = 0.0;

// probability of attack if predator is nearby
double p_att = 1.0;

// for convenience (and speed), introduce (equilibrium) probabilities of the
// environment being P
double pr_envt_is_P = 0.0;

// initial values for the evolving traits
double init_lfeedback = 0.0;
double init_lcue_influx = 0.0;
double init_lstress_influx = 0.0;
double init_linflux = 0.0;

// cue probabilities
// this is the cue that is given before
// a predation potentially happens
// if the cue is equal in environments P or NP
// then it is uninformative
double cue_P = 0.0;
double cue_NP = 0.0;

double ad = 1.0;
double aP = 1.0;
double r = 0;
double u = 0;

// parameters for input/output of pop to text file; they are read from the
// command line
int ioind = 0; // 0 = none; 1 = output to file; 2 = both in/output
std::string base_name;

// mortalities across the two environments in the order NP, P
int Nmort_stats[2] = {0,0};

// keep track of number of living individuals in pop (for bookkeeping)
int num_alive = Npop;

double fecundity_stats[2];


// the individual struct
struct Individual
{
    // assume alleles are on logit scale

    // diploid loci specifying the evolving traits
    //
    // self-dependent increase/decrease in hormone
    double lfeedback[2];

    // influx of new hormone when encountering predator cue
    double lcue_influx[2];

    // influx of new hormone when encountering predator attack
    double lstress_influx[2];

    // stress independent hormone influx
    double linflux[2];

    // components of the individual's state
    double hormone; // current hormone level
    double damage;  // current damage level

    // whether individual is now living in P or not; used to indicate the "true
    // situation", which need not be known by the individual
    bool envt_is_P;

    // indicator if the individual is alive
    bool alive;
};

// the population is stored in this variable
std::vector<Individual> pop(Npop);
// fecundity values are stored in this variable
std::vector<double> fecundity(Npop);


// initialize simulations from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu_feedback  = atof(argv[1]);
    mu_cue_influx  = atof(argv[2]);
    mu_stress_influx  = atof(argv[3]);
    mu_influx  = atof(argv[4]);
    sdmu = atof(argv[5]);
    s_P_2_NP  = atof(argv[6]);
    s_NP_2_P  = atof(argv[7]);
    init_lfeedback  = atof(argv[8]);
    init_lcue_influx  = atof(argv[9]);
    init_lstress_influx  = atof(argv[10]);
    init_linflux  = atof(argv[11]);
    cue_P  = atof(argv[12]);
    cue_NP  = atof(argv[13]);
    ad  = atof(argv[14]);
    aP  = atof(argv[15]);
    r = atof(argv[16]);
    u = atof(argv[17]);
    mort_background = atof(argv[18]);
    p_att = atof(argv[19]);
    maxtime = atoi(argv[20]);
    ioind = atoi(argv[21]);
    base_name = argv[22];

    // set equilibrium probabilities
    pr_envt_is_P = s_NP_2_P / (s_NP_2_P + s_P_2_NP);
}


// logistic function (to go from genotypic to phenotypic values)
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
        std::normal_distribution<> allele_dist(0.0, sdmu);
        G += allele_dist(rng_r);
    }
}


// write the parameters
void write_parameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< "seed;" << seed << ";"<< std::endl
		<< "Npop;" << Npop << ";"<< std::endl
		<< "mu_feedback;" << mu_feedback << ";"<< std::endl
        << "mu_cue_influx;" << mu_cue_influx << ";"<< std::endl
		<< "mu_stress_influx;" << mu_stress_influx << ";"<< std::endl
		<< "mu_influx;" << mu_influx << ";"<< std::endl
		<< "sP2NP;" << s_P_2_NP << ";"<< std::endl
		<< "sNP2P;" << s_NP_2_P << ";"<< std::endl
		<< "init_lfeedback;" << init_lfeedback << ";"<< std::endl
        << "init_lcue_influx;" << init_lcue_influx << ";"<< std::endl
		<< "init_lstress_influx;" << init_lstress_influx << ";"<< std::endl
		<< "init_linflux;" << init_linflux << ";"<< std::endl
		<< "cue_P;" << cue_P << ";"<< std::endl
		<< "cue_NP;" << cue_NP << ";"<< std::endl
		<< "ad;" << ad << ";"<< std::endl
		<< "aP;" << aP << ";"<< std::endl
		<< "r;" << r << ";"<< std::endl
        << "u;" << u << ";"<< std::endl
        << "mort_background;" << mort_background << ";"<< std::endl
		<< "p_att;" << p_att << ";"<< std::endl
        << "maxtime;" << maxtime << ";"<< std::endl
        << "base_name;" << base_name << ";"<< std::endl;
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
            newInd.lcue_influx[allele_i] = init_lcue_influx;
            newInd.lstress_influx[allele_i] = init_lstress_influx;
            newInd.linflux[allele_i] = init_linflux;
        }

        // initialize hormone level and damage to baseline values:
        // Note from Olle: tried to do this below
        double feedback =
            logistic(0.5*(newInd.lfeedback[0] + newInd.lfeedback[1]));
        double influx =
            logistic(0.5*(newInd.linflux[0] + newInd.linflux[1]));
        newInd.hormone = (feedback > 0) ? influx/feedback : 1.0;
        clamp(newInd.hormone, 0.0, 1.0);
        newInd.damage = (r > 0) ? u * newInd.hormone/r : 1.0;
        clamp(newInd.damage, 0.0, 1.0);

        // use of equilibrium probabilities for initial environments
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
    // assume allele 0 to be from maternal gamete and allele 1 from paternal
    // gamete; the gametes are formed using free recombination of parental loci
    double mat_lfeedback = mother.lfeedback[random_allele(rng_r)];
    double pat_lfeedback = father.lfeedback[random_allele(rng_r)];
    double mat_lcue_influx = mother.lcue_influx[random_allele(rng_r)];
    double pat_lcue_influx = father.lcue_influx[random_allele(rng_r)];
    double mat_lstress_influx = mother.lstress_influx[random_allele(rng_r)];
    double pat_lstress_influx = father.lstress_influx[random_allele(rng_r)];
    double mat_linflux = mother.linflux[random_allele(rng_r)];
    double pat_linflux = father.linflux[random_allele(rng_r)];

    kid.lfeedback[0] = mat_lfeedback;
    kid.lfeedback[1] = pat_lfeedback;
    kid.lcue_influx[0] = mat_lcue_influx;
    kid.lcue_influx[1] = pat_lcue_influx;
    kid.lstress_influx[0] = mat_lstress_influx;
    kid.lstress_influx[1] = pat_lstress_influx;
    kid.linflux[0] = mat_linflux;
    kid.linflux[1] = pat_linflux;

    // take into account mutation and note that allelic values are on logit
    // scale (and need not be clamped)
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        mutate(kid.lfeedback[allele_i], mu_feedback, sdmu);
        mutate(kid.lcue_influx[allele_i], mu_cue_influx, sdmu);
        mutate(kid.lstress_influx[allele_i], mu_stress_influx, sdmu);
        mutate(kid.linflux[allele_i], mu_influx, sdmu);
    }

    // set hormone level and damage to their baseline values
    double feedback =
        logistic(0.5*(kid.lfeedback[0] + kid.lfeedback[1]));
    double influx =
        logistic(0.5*(kid.linflux[0] + kid.linflux[1]));
    kid.hormone = (feedback > 0) ? influx/feedback : 1.0;
    clamp(kid.hormone, 0.0, 1.0);
    kid.damage = (r > 0) ? u * kid.hormone/r : 1.0;
    clamp(kid.damage, 0.0, 1.0);

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
void survive(std::ofstream &datafile)
{
    // reset fecundity averages that are tracked for stats
    fecundity_stats[0] = 0.0;
    fecundity_stats[1] = 0.0;

    // reset numbers of dead individuals that are tracked for stats
    Nmort_stats[0] = 0;
    Nmort_stats[1] = 0;

    assert(numP >= 0);
    assert(numNP >= 0);
    assert(numP + numNP <= Npop);

    for (int ind_i = 0; ind_i < Npop; ++ind_i)
    {
        if (pop[ind_i].alive) { // all individuals ought to be alive here
            // update hormone level from "background" influx and outflux
            double feedback =
                logistic(0.5*(pop[ind_i].lfeedback[0] +
                              pop[ind_i].lfeedback[1]));
            double cue_influx =
                logistic(0.5*(pop[ind_i].lcue_influx[0] +
                              pop[ind_i].lcue_influx[1]));
            double stress_influx =
                logistic(0.5*(pop[ind_i].lstress_influx[0] +
                              pop[ind_i].lstress_influx[1]));
            double influx =
                logistic(0.5*(pop[ind_i].linflux[0] +
                              pop[ind_i].linflux[1]));
            pop[ind_i].hormone = (1.0 - feedback)*pop[ind_i].hormone + influx;
            double p_cue = pop[ind_i].envt_is_P ? cue_P : cue_NP;
    
            if (uniform(rng_r) < p_cue) {
                // individual gets predator cue
                pop[ind_i].hormone += cue_influx;
            }

            clamp(pop[ind_i].hormone, 0.0, 1.0);
            // take into account possible predator attack
            if (pop[ind_i].envt_is_P) { // attack only possible if P
                if (uniform(rng_r) < p_att) { // predator attacks
                    double kill_prob = 1.0 - pow(pop[ind_i].hormone, aP);
                    if (uniform(rng_r) < kill_prob) {
                        // individual is killed by predator
                        pop[ind_i].alive = false;
                        --num_alive;
                        --numP;
                        ++Nmort_stats[1];
                        if (num_alive == 0)
                        {
                            std::cout << "extinct" << std::endl;
                            write_parameters(datafile);
                            exit(1);
                        }
                    } else {
                        // individual survives and gets hormone spike
                        pop[ind_i].hormone += stress_influx;
                        clamp(pop[ind_i].hormone, 0.0, 1.0);
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
                    std::cout << "extinct" << std::endl;
                    write_parameters(datafile);
                    exit(1);
                }
            } else {
                // update damage levels
                pop[ind_i].damage = (1.0 - r) * pop[ind_i].damage +
                    u * pop[ind_i].hormone;
                clamp(pop[ind_i].damage,0.0, 1.0);
                // damage-dependent fecundity
                fecundity[ind_i] = 1.0 - pow(pop[ind_i].damage,ad);
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
void reproduce_check(std::ofstream &datafile)
{
    // run through population and replace all dead individuals with offspring;
    // the parents for each offspring are randomly selected from the entire
    // population (this is a bit unrealistic); to select random parents, we use
    // a std::discrete_distribution<int> with fecundity vector as weights
    // (note: fecundity is zero for dead individuals)
    std::discrete_distribution<int> par_distr(fecundity.begin(),
                                              fecundity.end());
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
void write_data(std::ofstream &DataFile)
{
    double mean_feedback = 0;
    double ss_feedback = 0;
    double mean_cue_influx = 0;
    double ss_cue_influx = 0;
    double mean_stress_influx = 0;
    double ss_stress_influx = 0;
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
            // feedback
            double feedback =
                logistic(0.5*(pop[ind_i].lfeedback[0] +
                              pop[ind_i].lfeedback[1]));
            mean_feedback += feedback;
            ss_feedback += feedback * feedback;
            // cue_influx
            double cue_influx =
                logistic(0.5*(pop[ind_i].lcue_influx[0] +
                              pop[ind_i].lcue_influx[1]));
            mean_cue_influx += cue_influx;
            ss_cue_influx += cue_influx * cue_influx;
            // stress_influx
            double stress_influx =
                logistic(0.5*(pop[ind_i].lstress_influx[0] +
                              pop[ind_i].lstress_influx[1]));
            mean_stress_influx += stress_influx;
            ss_stress_influx += stress_influx * stress_influx;
            // influx
            double influx =
                logistic(0.5*(pop[ind_i].linflux[0] +
                              pop[ind_i].linflux[1]));
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
    mean_feedback /= Npop;
    mean_cue_influx /= Npop;
    mean_stress_influx /= Npop;
    mean_influx /= Npop;
    mean_hormone /= Npop;
    mean_damage /= Npop;

    double sd_feedback = sqrt(ss_feedback / Npop - pow(mean_feedback,2.0));
    double sd_cue_influx = sqrt(ss_cue_influx / Npop -
                                pow(mean_cue_influx,2.0));
    double sd_stress_influx = sqrt(ss_stress_influx / Npop -
                                   pow(mean_stress_influx,2.0));
    double sd_influx = sqrt(ss_influx / Npop - pow(mean_influx,2.0));
    double sd_hormone = sqrt(ss_hormone / Npop - pow(mean_hormone,2.0));
    double sd_damage = sqrt(ss_damage / Npop - pow(mean_damage,2.0));

    DataFile << generation << ";"
        << freq_P << ";"
        << mean_feedback << ";"
        << mean_cue_influx << ";"
        << mean_stress_influx << ";"
        << mean_influx << ";"
        << mean_hormone << ";"
        << mean_damage << ";"
        << sd_feedback << ";"
        << sd_cue_influx << ";"
        << sd_stress_influx << ";"
        << sd_influx << ";"
        << sd_hormone << ";"
        << sd_damage << ";"
        << (double)Nmort_stats[1]/numP << ";"
        << (double)Nmort_stats[0]/numNP << ";"
        << (double)fecundity_stats[1] << ";"
        << (double)fecundity_stats[0] << ";"
        << std::endl;
}

// write the headers of a datafile
void write_data_headers(std::ofstream &DataFile)
{
    DataFile << "generation;"
        << "freq_P" << ";"
        << "mean_feedback" << ";"
        << "mean_cue_influx" << ";"
        << "mean_stress_influx" << ";"
        << "mean_influx" << ";"
        << "mean_hormone" << ";"
        << "mean_damage" << ";"
        << "sd_feedback" << ";"
        << "sd_cue_influx" << ";"
        << "sd_stress_influx" << ";"
        << "sd_influx" << ";"
        << "sd_hormone" << ";"
        << "sd_damage" << ";"
        << "prop_dead_P" << ";"
        << "prop_dead_NP" << ";"
        << "mean_fecundity_P" << ";"
        << "mean_fecundity_NP" << ";"
        << std::endl;
}

// iterate individuals for tmax timesteps to plot the stress response curve for
// different individuals
void write_simple_iter(std::ofstream &IterFile)
{
    // number of individuals
    int nrep = 50;
    int tmax = 500;
    int tstress = 100;

    // use an individual's genotype to compute its hormone response to, first,
    // a single predator cue at time step tstress, and second, a predator
    // attack at time step tstress; assume the hormone is at the background
    // equilibrium level at the start

    Individual ind;

    IterFile << "time;individual;hormone;" << std::endl;

    std::uniform_int_distribution<int> rint(0, Npop - 1);

    // sample individuals
    for (int ind_i = 0; ind_i < nrep; ++ind_i)
    {
        ind = pop[rint(rng_r)];
        // get feedback, cue_influx, stress_influx and influx from genotype
        double feedback =
            logistic(0.5*(ind.lfeedback[0] +
                          ind.lfeedback[1]));
        double cue_influx =
            logistic(0.5*(ind.lcue_influx[0] +
                          ind.lcue_influx[1]));
        double stress_influx =
            logistic(0.5*(ind.lstress_influx[0] +
                          ind.lstress_influx[1]));
        double influx =
            logistic(0.5*(ind.linflux[0] +
                          ind.linflux[1]));

        // hormone background level
        double hormone0 = (feedback > 0) ? influx/feedback : 1.0;
        clamp(hormone0, 0.0, 1.0);
        double hormone_cue = hormone0;
        double hormone_att = hormone0;

        IterFile << "cue" << ";" << 0 << ";" << ind_i << ";"
                 << hormone_cue << ";" << std::endl;
        IterFile << "att" << ";" << 0 << ";" << ind_i << ";"
                 << hormone_att << ";" << std::endl;

        // iterate the stress responses for this individual
        for (int timestep = 0; timestep < tmax; ++timestep)
        {
            hormone_cue = (1.0 - feedback)*hormone_cue + influx;
            hormone_att = (1.0 - feedback)*hormone_att + influx;

            if (timestep == tstress)
            {
                hormone_cue += cue_influx;
                hormone_att += stress_influx;
            }

            clamp(hormone_cue, 0.0, 1.0);
            clamp(hormone_att, 0.0, 1.0);

            IterFile << "cue" << ";" << timestep << ";" << (ind_i + 1) << ";"
                     << hormone_cue << ";" << std::endl;
            IterFile << "att" << ";" << timestep << ";" << (ind_i + 1) << ";"
                     << hormone_att << ";" << std::endl;
        }
    }
}


// read population from tab-separated text file
void read_pop_from_file(std::string infilename)
{
    std::ifstream infile(infilename.c_str());
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
        double cue_influx; // dummy variable
        double stress_influx; // dummy variable
        double influx; // dummy variable
        int envtP = 0;
        int alv = 0;
        while (infile && ind_i < Npop) {
            infile >> indi.lfeedback[0];
            infile >> indi.lfeedback[1];
            infile >> indi.lcue_influx[0];
            infile >> indi.lcue_influx[1];
            infile >> indi.lstress_influx[0];
            infile >> indi.lstress_influx[1];
            infile >> indi.linflux[0];
            infile >> indi.linflux[1];
            infile >> feedback;
            infile >> cue_influx;
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

// write population to tab-separated text file
void write_pop_to_file(std::ofstream &PopFile)
{
    // first write headers
    PopFile << "lfeedback1" << "\t" << "lfeedback2" << "\t"
        << "lcue_influx1" << "\t" << "lcue_influx2" << "\t"
        << "lstress_influx1" << "\t" << "lstress_influx2" << "\t"
        << "linflux1" << "\t" << "linflux2" << "\t"
        << "feedback" << "\t"
        << "cue_influx" << "\t"
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
        double cue_influx =
            logistic(0.5*(pop[ind_i].lcue_influx[0] +
                          pop[ind_i].lcue_influx[1]));
        double stress_influx =
            logistic(0.5*(pop[ind_i].lstress_influx[0] +
                          pop[ind_i].lstress_influx[1]));
        double influx =
            logistic(0.5*(pop[ind_i].linflux[0] +
                          pop[ind_i].linflux[1]));
        PopFile << pop[ind_i].lfeedback[0] << "\t"
                << pop[ind_i].lfeedback[1] << "\t"
                << pop[ind_i].lcue_influx[0] << "\t"
                << pop[ind_i].lcue_influx[1] << "\t"
                << pop[ind_i].lstress_influx[0] << "\t"
                << pop[ind_i].lstress_influx[1] << "\t"
                << pop[ind_i].linflux[0] << "\t"
                << pop[ind_i].linflux[1] << "\t"
                << feedback << "\t"
                << cue_influx << "\t"
                << stress_influx << "\t"
                << influx << "\t"
                << pop[ind_i].hormone << "\t"
                << pop[ind_i].damage << "\t"
                << pop[ind_i].envt_is_P << "\t"
                << pop[ind_i].alive << "\n";
    }
}


int main(int argc, char ** argv)
{
	init_arguments(argc, argv);

    if (ioind > 1) {
        // read pop from file
        read_pop_from_file(base_name + "pop.txt");
    } else {
        // initialize from input parameters
        init_population();
    }

    std::ofstream DataFile(base_name.c_str());
    std::ofstream IterFile((base_name + "iters.csv").c_str());

    // write some params
	write_parameters(DataFile);

    // write headers to the datafile
	write_data_headers(DataFile);

    // generation is really time step, with overlapping generations
	for (generation = 0; generation <= maxtime; ++generation)
	{
        environmental_switching();

		survive(DataFile);

		reproduce_check(DataFile);

        do_stats = generation % skip == 0;

        // print statistics every skip generation
        if (do_stats)
		{
			write_data(DataFile);
		}
	}

    if (ioind > 0) {
        // write pop to file
        std::ofstream PopFile((base_name + "pop.txt").c_str());
        write_pop_to_file(PopFile);
    }

    // iterate the stress response curves for a subset of individuals
    write_simple_iter(IterFile);
}
