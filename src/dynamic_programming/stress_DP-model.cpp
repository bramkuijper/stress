// **********************************************************************************
// Dynamic programming model of stress response, described in:
// "Towards an evolutionary theory of stress responses" by Taborsky et al.
//
// March 2020, Exeter
// **********************************************************************************


//HEADER FILES

#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <string>



// constants, type definitions, etc.

using namespace std;

const int seed        = time(0); // pseudo-random seed
const double lambdaA  = 0.005;   // probability that predator arrives
const double lambdaL  = 0.095;   // probability that predator leaves
const double pAtt     = 0.5;     // probability that predator attacks if present
const double alpha    = 1.0;     // parameter controlling effect of hormone level on pKill
const double beta     = 1.5;     // parameter controlling effect of hormone level on reproductive rate
const double mu       = 0.002;   // background mortality (independent of hormone level and predation risk)
const int maxI        = 1000000; // maximum number of iterations
const int maxT        = 200;     // maximum number of time steps since last saw predator
const int maxH        = 1000;    // maximum hormone level
const int skip        = 10;      // interval between print-outs

ofstream outputfile; // output file
stringstream outfile; // for naming output file

int hormone[maxT];          // hormone level (strategy)
double pKill[maxH];         // probability of being killed by an attacking predator
double repro[maxH];         // reproductive output
double Wopt[maxT];          // fitness immediately after predator has/hasn't attacked, under optimal decision h
double W[maxT][maxH];       // expected fitness at start of time step, before predator does/doesn't attack
double Wnext[maxT][maxH];   // expected fitness at start of next time step
double pPred[maxT];         // probability that predator is present
double totfitdiff;          // fitness difference between optimal strategy in successive iterations

int i;     // iteration




/* SPECIFY FINAL FITNESS */
void FinalFit()
{
  int t,h;

  for (t=1;t<maxT;t++) // note that Wnext is undefined for t=0 because t=1 if predator has just attacked
  {
    for (h=0;h<maxH;h++)
    {
      Wnext[t][h] = 1.0;
    }
  }
}


/* CALCULATE PROBABILITY THAT PREDATOR IS PRESENT */
void Predator()
{
  int t;

  pPred[1] = 1.0-lambdaL; // if predator attacked in last time step
  for (t=2;t<maxT;t++) // if predator did NOT attack in last time step
  {
    pPred[t] = (pPred[t-1]*(1.0-pAtt)*(1.0-lambdaL)+(1.0-pPred[t-1])*lambdaA) / (1.0 - pPred[t-1]*pAtt);
  }

}



/* CALCULATE PROBABILITY OF BEING KILLED BY AN ATTACKING PREDATOR */
void Death()
{
  int h;

  for (h=0;h<maxH;h++)
  {
    pKill[h] = 1.0 - pow(double(h)/double(maxH),alpha);
  }
}


/* CALCULATE PROBABILITY OF REPRODUCING */
void Reproduction()
{
  int h;

  for (h=0;h<maxH;h++)
  {
    repro[h] = 1.0 - pow(double(h)/double(maxH),beta);
  }
}



/* CALCULATE OPTIMAL DECISION FOR EACH t */
void OptDec()
{
  int t,h;
  double fitness;

  // calculate optimal decision h given current t (N.B. t=0 if survived attack)
  for (t=0;t<maxT;t++)
  {
    Wopt[t] = 0.0;
    hormone[t] = 0;
    for (h=0;h<maxH;h++)
    {
      fitness = Wnext[min(maxT-1,t+1)][h]; // fitness as a function of h
      if (fitness>Wopt[t])
      {
        Wopt[t] = fitness; // overwrite Wopt
        hormone[t] = h; // overwrite hormone
      }
    }
  }

  // calculate expected fitness as a function of t and h, before predator does/doesn't attack
  for (t=1;t<maxT;t++) // note that W is undefined for t=0 because t=1 if predator has just attacked
  {
    for (h=0;h<maxH;h++)
    {
      W[t][h] = pPred[t]*pAtt*(1.0-pKill[h])*(1-mu)*(repro[h]+Wopt[0]) // survive attack
                + (1.0-pPred[t]*pAtt)*(1-mu)*(repro[h]+Wopt[t]); // no attack
    }
  }

}



/* OVERWRITE FITNESS ARRAY FROM PREVIOUS ITERATION */
void ReplaceFit()
{
  int t,h;
  double fitdiff;

  fitdiff = 0.0;

  for (t=1;t<maxT;t++)
  {
    for (h=0;h<maxH;h++)
    {
      fitdiff = fitdiff + abs(Wnext[t][h]-W[t][h]);
      Wnext[t][h] = W[t][h];
    }
  }

  totfitdiff = fitdiff;

}



/* PRINT OUT OPTIMAL STRATEGY */
void PrintStrat()
{
  int t;

  outputfile << "t" << "\t" << "hormone" << endl;

  for (t=0;t<maxT;t++)
  {
    outputfile << t << "\t" << hormone[t] << endl;
  }
  outputfile << endl;
  outputfile << "nIterations" << "\t" << i << endl;
  outputfile << endl;
}




/* WRITE PARAMETER SETTINGS TO OUTPUT FILE */
void PrintParams()
{
  outputfile << endl << "PARAMETER VALUES" << endl
       << "lambdaL: " << "\t" << lambdaL << endl
       << "lambdaA: " << "\t" << lambdaA << endl
       << "pAtt: " << "\t" << pAtt << endl
       << "alpha: " << "\t" << alpha << endl
       << "beta: " << "\t" << beta << endl
       << "mu: " << "\t" << mu << endl
       << "maxI: " << "\t" << maxI << endl
       << "maxT: " << "\t" << maxT << endl
       << "maxH: " << "\t" << maxH << endl;
}



/* MAIN PROGRAM */
int main()
{

		///////////////////////////////////////////////////////
		outfile.str("");
		outfile << "stress.txt";
		string outputfilename = outfile.str();
		outputfile.open(outputfilename.c_str());
		///////////////////////////////////////////////////////

        outputfile << "Random seed: " << seed << endl; // write seed to output file

        FinalFit();
        Predator();
        Death();
        Reproduction();

        cout << "i" << "\t" << "totfitdiff" << endl;
        for (i=1;i<=maxI;i++)
          {
          OptDec();
          ReplaceFit();

          if (totfitdiff < 0.000001) break; // strategy has converged on optimal solution, so exit loop
          if (i==maxI) { outputfile << "*** DID NOT CONVERGE WITHIN " << i << " ITERATIONS ***" << endl;}

 		  if (i%skip==0)
            {
              cout << i << "\t" << totfitdiff << endl; // show fitness difference every 'skip' generations
            }
          }

        cout << endl;
        outputfile << endl;

        PrintStrat();
        PrintParams();
        outputfile.close();


  return 0;
}
