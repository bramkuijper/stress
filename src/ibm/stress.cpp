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



// set up the random number generator using a good way of getting random seed
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);
std::uniform_real_distribution<> uniform(0.0,1.0);

// for meiotic segregation
std::bernoulli_distribution random_allele(0.5);


// the NumGen parameter is actually the number of time steps with overlapping
// generations; it is read from the command line
int NumGen = 200000;

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
double mu_influx = 0.0;
double mu_hstart = 0.0;
double sdmu_clearance = 0.0;
double sdmu_stress_influx = 0.0;
double sdmu_influx = 0.0;
double sdmu_hstart = 0.0;

// background mortality
double mort_background = 0.0;

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
// restricted to the unit interval; for this reason, zmax and dopt should be
// equal to 1.0

// initial values for the evolving traits
double init_clearance = 0.0;
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
double dopt = 1.0;
double sel_width = 4.0;
double zmax = 1.0;
double min_clearance = 1.0;
double r = 1.0;
double u = 1.0;



// parameters for input/output of pop to text file; they are read from the
// command line
int ioind = 0; // 0 = none; 1 = output to file; 2 = both in/output
std::string base_name;

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
    // influx of new hormone when encountering stress
    double stress_influx[2];
    
    // stress independent hormone influx
    double influx[2];

    // starting hormone level
    double hstart[2];

    // components of the individual's state
    double hormone; // current hormone level
    double damage;  // current damage level

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
    mu_clearance = atof(argv[1]);
    mu_stress_influx  = atof(argv[2]);
    mu_influx  = atof(argv[3]);
    mu_hstart  = atof(argv[4]);

    sdmu_clearance = atof(argv[5]);
    sdmu_stress_influx = atof(argv[6]);
    sdmu_influx = atof(argv[7]);
    sdmu_hstart = atof(argv[8]);

    s_P_2_NP  = atof(argv[9]);
    s_NP_2_P  = atof(argv[10]);

    init_clearance = atof(argv[11]);
    init_stress_influx  = atof(argv[12]);
    init_influx  = atof(argv[13]);
    init_hstart  = atof(argv[14]);

    cue_P  = atof(argv[15]);
    cue_NP  = atof(argv[16]);

    s0  = atof(argv[17]);
    ad  = atof(argv[18]);
    aP  = atof(argv[19]);
    dopt  = atof(argv[20]);
    sel_width = atof(argv[21]);
    zmax  = atof(argv[22]);
    min_clearance = atof(argv[23]);

    r = atof(argv[24]);
    u = atof(argv[25]);

    mort_background = atof(argv[26]);
    p_att = atof(argv[27]);
    NumGen = atoi(argv[28]);
    ioind = atoi(argv[29]);
    base_name = argv[30];
    // equilibrium probabilities
    pr_envt_is_P = s_NP_2_P / (s_NP_2_P + s_P_2_NP);
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
// (typically at the end of the output file)
void write_parameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< std::endl
		<< "seed;" << seed << ";"<< std::endl
		<< "Npop;" << Npop << ";"<< std::endl
		<< "mu_clearance;" << mu_clearance << ";"<< std::endl
		<< "mu_stress_influx;" << mu_stress_influx << ";"<< std::endl
        << "mu_influx;" << mu_influx << ";"<< std::endl
		<< "mu_hstart;" << mu_hstart << ";"<< std::endl
		<< "sP2NP_1;" << s_P_2_NP << ";"<< std::endl
		<< "sNP2P_1;" << s_NP_2_P << ";"<< std::endl
		<< "init_clearance;" << init_clearance << ";"<< std::endl
		<< "init_stress_influx;" << init_stress_influx << ";"<< std::endl
        << "init_influx;" << init_influx << ";"<< std::endl
		<< "init_hstart;" << init_hstart << ";"<< std::endl
		<< "cue_P;" << cue_P << ";"<< std::endl
		<< "cue_NP;" << cue_NP << ";"<< std::endl
		<< "s0;" << s0 << ";"<< std::endl
		<< "ad;" << ad << ";"<< std::endl
		<< "aP;" << aP << ";"<< std::endl
        << "risk;" << s_NP_2_P / (s_NP_2_P + s_P_2_NP)  << ";" << std::endl
        << "autocorr;" << 1.0 - s_NP_2_P - s_P_2_NP  << ";" << std::endl
		<< "mort_background;" << mort_background << ";"<< std::endl
		<< "dopt;" << dopt << ";"<< std::endl
		<< "sel_width;" << sel_width << ";"<< std::endl
		<< "zmax;" << zmax << ";"<< std::endl
		<< "min_clearance;" << min_clearance << ";"<< std::endl
		<< "r;" << r << ";"<< std::endl
        << "u;" << u << ";"<< std::endl
		<< "p_att;" << p_att << ";"<< std::endl
        << "NumGen;" << NumGen << ";"<< std::endl
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
            newInd.clearance[allele_i] = init_clearance;
            newInd.stress_influx[allele_i] = init_stress_influx;
            newInd.influx[allele_i] = init_influx;
            newInd.hstart[allele_i] = init_hstart;
        }

        // initialize hormone level and damage
        newInd.hormone = 0.5*(newInd.hstart[0] + newInd.hstart[1]);
        newInd.hormone = std::clamp(newInd.hormone, 0.0, zmax);
        newInd.damage = (r > 0) ? u * newInd.hormone/r : 1.0;
        newInd.damage = std::clamp(newInd.damage, 0.0, 1.0);

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
    double mat_clearance = mother.clearance[random_allele(rng_r)];
    double pat_clearance = father.clearance[random_allele(rng_r)];
    double mat_stress_influx = mother.stress_influx[random_allele(rng_r)];
    double pat_stress_influx = father.stress_influx[random_allele(rng_r)];
    double mat_influx = mother.influx[random_allele(rng_r)];
    double pat_influx = father.influx[random_allele(rng_r)];
    double mat_hstart = mother.hstart[random_allele(rng_r)];
    double pat_hstart = father.hstart[random_allele(rng_r)];

    kid.clearance[0] = mat_clearance;
    kid.clearance[1] = pat_clearance;
    kid.stress_influx[0] = mat_stress_influx;
    kid.stress_influx[1] = pat_stress_influx;
    kid.influx[0] = mat_influx;
    kid.influx[1] = pat_influx;
    kid.hstart[0] = mat_hstart;
    kid.hstart[1] = pat_hstart;

    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        mutate(
                kid.clearance[allele_i],
                mu_clearance,
                sdmu_clearance);

        kid.clearance[allele_i] = std::clamp(
                kid.clearance[allele_i], min_clearance, 1.0);

        mutate(
                kid.stress_influx[allele_i],
                mu_stress_influx,
                sdmu_stress_influx
                );

        kid.stress_influx[allele_i] = std::clamp(kid.stress_influx[allele_i], 0.0, 1.0);

        mutate(
                kid.influx[allele_i],
                mu_influx,
                sdmu_influx);
        kid.influx[allele_i] = std::clamp(kid.influx[allele_i], 0.0, 1.0);

        mutate(
                kid.hstart[allele_i],
                mu_hstart,
                sdmu_hstart);
        kid.hstart[allele_i] = std::clamp(kid.hstart[allele_i], 0.0, 1.0);
    }

    // set hormone level and damage
    kid.hormone = 0.5*(kid.hstart[0] + kid.hstart[1]);
    kid.hormone = std::clamp(kid.hormone, 0.0, zmax);

    kid.damage = u * kid.hormone;
    kid.damage = std::clamp(kid.damage, 0.0, 1.0);

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
void survive(std::ofstream &datafile)
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
            // update hormone level from "background" influx and outflux
            double clearance = 0.5*(pop[ind_i].clearance[0] +
                              pop[ind_i].clearance[1]);
            double stress_influx = 0.5*(pop[ind_i].stress_influx[0] +
                                   pop[ind_i].stress_influx[1]);
            double influx = 0.5*(pop[ind_i].influx[0] +
                            pop[ind_i].influx[1]);
            pop[ind_i].hormone = (1.0 - clearance)*pop[ind_i].hormone + influx;
            double p_cue = pop[ind_i].envt_is_P ? cue_P : cue_NP;
            if (uniform(rng_r) < p_cue) {
                // individual gets predator cue
                pop[ind_i].hormone += stress_influx;
            }
            pop[ind_i].hormone = std::clamp(pop[ind_i].hormone, 0.0, zmax);
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
                            std::cout << "extinct" << std::endl;
                            write_parameters(datafile);
                            exit(1);
                        }
                    } else {
                        // individual survives and gets hormone spike
                        pop[ind_i].hormone += stress_influx;
                        pop[ind_i].hormone = std::clamp(pop[ind_i].hormone, 0.0, zmax);
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
                pop[ind_i].damage = std::clamp(pop[ind_i].damage,0.0,1.0);
                // damage-dependent fecundity
                fecundity[ind_i] = 1.0 - sel_width * pow(pop[ind_i].damage - dopt,ad);  //1.0 - pow(pop[ind_i].damage/dmax,ad);
                fecundity[ind_i] = std::clamp(fecundity[ind_i],0.0,1.0);
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
    // a std::discrete_distribution<int> with fecundity vector as weights (note:
    // fecundity is zero for dead individuals)
    std::discrete_distribution<int> par_distr(fecundity.begin(), fecundity.end());
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
    double mean_clearance = 0;
    double ss_clearance = 0;

    double mean_stress_influx = 0;
    double ss_stress_influx = 0;

    double mean_influx = 0;
    double ss_influx = 0;

    double mean_hormone = 0;
    double ss_hormone = 0;

    double mean_damage = 0;
    double ss_damage = 0;

    double freq_P = (double) numP / (numP + numNP);

    int n_alive = 0;

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

            ++n_alive;
        }
    }

    mean_clearance /= n_alive;
    mean_stress_influx /= n_alive;
    mean_influx /= n_alive;
    mean_hormone /= n_alive;
    mean_damage /= n_alive;

    double sd_clearance = sqrt(ss_clearance / n_alive - 
                                        pow(mean_clearance,2.0));

    double sd_stress_influx = sqrt(ss_stress_influx / n_alive -
                                   pow(mean_stress_influx,2.0));

    double sd_influx = sqrt(ss_influx / n_alive - pow(mean_influx,2.0));
    double sd_hormone = sqrt(ss_hormone / n_alive - pow(mean_hormone,2.0));
    double sd_damage = sqrt(ss_damage / n_alive - pow(mean_damage,2.0));

    DataFile << generation << ";"
        << freq_P << ";"
        << mean_clearance << ";"
        << mean_stress_influx << ";"
        << mean_influx << ";"
        << mean_hormone << ";"
        << mean_damage << ";"
        << sd_clearance << ";"
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
        << "mean_clearance" << ";"
        << "mean_stress_influx" << ";"
        << "mean_influx" << ";"
        << "mean_hormone" << ";"
        << "mean_damage" << ";"
        << "sd_clearance" << ";"
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

// iterate individuals for tmax timesteps
// to plot the stress response curve for
// different individuals
void write_simple_iter(std::ofstream &IterFile)
{
    // number of individuals
    int nrep = 20;
    int tmax = 500;
    int tstress = 100;

    // use an individual's genotype to compute its
    // hormone response to a single stressor event at time step tstress

    double hormone, hormone_tplus1;

    Individual ind;

    IterFile << "time;individual;hormone;" << std::endl;

    std::uniform_int_distribution<int> rint(0, Npop - 1);

    // sample individuals
    for (int ind_i = 0; ind_i < nrep; ++ind_i)
    {
        ind = pop[rint(rng_r)];

        // get clearance, stress_influx, influx and hstart from genotype
        double clearance = 0.5*(ind.clearance[0] + ind.clearance[1]);
        double stress_influx = 0.5*(ind.stress_influx[0] +
                               ind.stress_influx[1]);
        double influx = 0.5*(ind.influx[0] + ind.influx[1]);
        double hstart = 0.5*(ind.hstart[0] + ind.hstart[1]);

        // hormone baseline level
        hormone = (clearance > 0) ? influx/clearance : 1.0;
        hormone = hstart;
        hormone = std::clamp(hormone, 0.0, zmax);
        hormone_tplus1 = 0.0;

        IterFile << 0 << ";" << ind_i << ";" << hormone << ";" << std::endl;

        // iterate the stress response for this individual
        for (int timestep = 0; timestep < tmax; ++timestep)
        {
            hormone_tplus1 = (1.0 - clearance)*hormone + influx;
            if (timestep == tstress)
            {
                hormone_tplus1 += stress_influx;
            }
            hormone_tplus1 = std::clamp(hormone_tplus1, 0.0, zmax);
            hormone = hormone_tplus1;

            IterFile << timestep << ";" << (ind_i + 1) << ";"
                     << hormone << ";" << std::endl;
        }
    }
}

// code to read population from tab-separated text file
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
        double clearance; // dummy variable
        double stress_influx; // dummy variable
        double influx; // dummy variable
        double hstart; // dummy variable
        int envtP = 0;
        int alv = 0;
        while (infile && ind_i < Npop) {
            infile >> indi.clearance[0];
            infile >> indi.clearance[1];
            infile >> indi.stress_influx[0];
            infile >> indi.stress_influx[1];
            infile >> indi.influx[0];
            infile >> indi.influx[1];
            infile >> indi.hstart[0];
            infile >> indi.hstart[1];
            infile >> clearance;
            infile >> stress_influx;
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
void write_pop_to_file(std::ofstream &PopFile)
{
    // first write headers
    PopFile << "clearance1" << "\t" << "clearance2" << "\t"
        << "stress_influx1" << "\t" << "stress_influx2" << "\t"
        << "influx1" << "\t" << "influx2" << "\t"
        << "hstart1" << "\t" << "hstart2" << "\t"
        << "clearance" << "\t"
        << "stress_influx" << "\t"
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
        double influx = 0.5*(pop[ind_i].influx[0] +
                        pop[ind_i].influx[1]);
        double hstart = 0.5*(pop[ind_i].hstart[0] +
                        pop[ind_i].hstart[1]);
        PopFile << pop[ind_i].clearance[0] << "\t"
                << pop[ind_i].clearance[1] << "\t"
                << pop[ind_i].stress_influx[0] << "\t"
                << pop[ind_i].stress_influx[1] << "\t"
                << pop[ind_i].influx[0] << "\t"
                << pop[ind_i].influx[1] << "\t"
                << pop[ind_i].hstart[0] << "\t"
                << pop[ind_i].hstart[1] << "\t"
                << clearance << "\t"
                << stress_influx << "\t"
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
    std::ofstream DataFile(base_name.c_str());
    std::ofstream IterFile((base_name + "iters.csv").c_str());

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
        std::ofstream PopFile((base_name + "pop.txt").c_str());
        write_pop_to_file(PopFile);
    }

    // iterate the stress response curves for a subset of individuals
    write_simple_iter(IterFile);
//    // finally write some params
//	write_parameters(DataFile);
}
