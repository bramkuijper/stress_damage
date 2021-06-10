// **********************************************************************************
// Dynamic programming model of stress response with somatic damage.
//
// Survival only.
//
// May 2021, Exeter
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
//const int seed      = << enter seed here >>;

double pLeave;   // probability that predator leaves
double pArrive;  // probability that predator arrives
const double pAttack  = 0.5;     // probability that predator attacks if present
const double alpha    = 1.0;     // parameter controlling effect of hormone level on pKilled
//const double beta     = 1.5;     // parameter controlling effect of hormone level on reproductive rate
const double mu0      = 0.002;   // background mortality (independent of hormone level and predation risk)
const double phi_inv  = 1.0/((sqrt(5.0)+1.0)/2.0); // inverse of golden ratio (for golden section search)
const int maxD        = 20;      // maximum damage level
const double Kmort        = 0.01;    // parameter Kmort controlling increase in mortality with damage level
const double Kfec        = 0.0;    // parameter Kmort controlling increase in mortality with damage level
const int maxI        = 1000000; // maximum number of iterations
const int maxT        = 100;     // maximum number of time steps since last saw predator
const int maxH        = 500;     // maximum hormone level
const int skip        = 10;      // interval between print-outs

ofstream outputfile;  // output file
ofstream fwdCalcfile; // forward calculation output file
ofstream attsimfile;  // simulated attacks output file
stringstream outfile; // for naming output file

// random numbers
mt19937 mt(seed); // random number generator
uniform_real_distribution<double> Uniform(0, 1); // real number between 0 and 1 (uniform)

int hormone[maxT][maxD+1];        // hormone level (strategy)
double pKilled[maxH];             // probability of being killed by an attacking predator
double mu[maxD+1];                // probability of background mortality, as a function of damage
double dnew[maxD+1][maxH];        // new damage level, as a function of previous damage and hormone
double repro[maxD+1];       // reproductive output
double Wopt[maxT][maxD+1];        // fitness immediately after predator has/hasn't attacked, under optimal decision h
double W[maxT][maxD+1][maxH];     // expected fitness at start of time step, before predator does/doesn't attack
double Wnext[maxT][maxD+1][maxH]; // expected fitness at start of next time step
double F[maxT][maxD+1][maxH];     // frequency of individuals at start of time step, before predator does/doesn't attack
double Fnext[maxT][maxD+1][maxH]; // frequency of individuals at start of next time step
double pPred[maxT];               // probability that predator is present
double totfitdiff;                // fitness difference between optimal strategy in successive iterations

int i;     // iteration




/* SPECIFY FINAL FITNESS */
void FinalFit()
{
  int t,d,h;

  for (t=1;t<maxT;t++) // note that Wnext is undefined for t=0 because t=1 if predator has just attacked
  {
    for (d=0;d<=maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        Wnext[t][d][h] = 1.0;
      }
    }
  }
}


/* CALCULATE PROBABILITY THAT PREDATOR IS PRESENT */
void PredProb()
{
  int t;

  pPred[1] = 1.0-pLeave; // if predator attacked in last time step

  for (t=2;t<maxT;t++) // if predator did NOT attack in last time step
  {
//      Pr(predator present at time t | predator did not attack at time t-1)
//      = Pr (predator did not attack at time t-1 | predator present at time t) * 
//              Pr( predator present at time t)  / Pr(predator did not attack at time t-1)
    pPred[t] = (pPred[t-1]*(1.0-pAttack)*(1.0-pLeave)+(1.0-pPred[t-1])*pArrive) / (1.0 - pPred[t-1]*pAttack);
  }
}



/* CALCULATE PROBABILITY OF BEING KILLED BY AN ATTACKING PREDATOR */
void Predation()
{
  int h;

  for (h=0;h<maxH;h++)
  {
    pKilled[h] = 1.0 - pow(double(h)/double(maxH),alpha);
  }
}


/* CALCULATE BACKGROUND MORTALITY */
void Mortality()
{
  int d;

  for (d=0;d<=maxD;d++)
  {
    mu[d] = min(1.0,mu0 + Kmort*double(d));
  }
}



/* CALCULATE DAMAGE */
void Damage()
{
  int d,h;

  for (d=0;d<=maxD;d++)
  {
    for (h=0;h<maxH;h++)
    {
      dnew[d][h] = max(0.0,min(double(maxD),double(d) + 4.0*(double(h)/double(maxH))*(double(h)/double(maxH))-1.0));
    }
  }
}



/* CALCULATE PROBABILITY OF REPRODUCING */
void Reproduction()
{
  int d;

  for (d=0;d<=maxD;d++)
  {
    repro[d] = max(0.0, 1.0 - Kfec*double(d));
  }
} // end Reproduction()



/* CALCULATE OPTIMAL DECISION FOR EACH t */
void OptDec()
{
  int t,h,d,d1,d2,
    LHS,RHS,x1,x2,cal_x1,cal_x2;
  double fitness,fitness_x1,fitness_x2,ddec;

  // calculate optimal decision h given current t and d (N.B. t=0 if survived attack)
  for (t=0;t<maxT;t++)
  {
    for (d=0;d<=maxD;d++)
    {
      // GOLDEN SECTION SEARCH
      // following https://medium.datadriveninvestor.com/golden-section-search-method-peak-index-in-a-mountain-array-leetcode-852-a00f53ed4076
      LHS = 0;
      RHS = maxH;
      x1 = RHS - (round((double(RHS)-double(LHS))*phi_inv));
      x2 = LHS + (round((double(RHS)-double(LHS))*phi_inv));

      while (x1<x2)
      {
        fitness_x1 = Wnext[min(maxT-1,t+1)][d][x1]; // fitness as a function of h=x1
        fitness_x2 = Wnext[min(maxT-1,t+1)][d][x2]; // fitness as a function of h=x2

        if (fitness_x1<fitness_x2)
        {
            LHS = x1;
            x1 = x2;
            x2 = RHS - (round((double(RHS)-double(x1))*phi_inv));
        }
        else
        {
            RHS = x2;
            x2 = x1;
            x1 = LHS + (round((double(x2)-double(LHS))*phi_inv));
        }
      }
      hormone[t][d] = x1; // optimal hormone level
      Wopt[t][d] = fitness_x1; // fitness of optimal decision

    }
  }

  // calculate expected fitness as a function of t, h and d, before predator does/doesn't attack
  for (t=1;t<maxT;t++) // note that W is undefined for t=0 because t=1 if predator has just attacked
  {
    for (d=0;d<=maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        d1=floor(dnew[d][h]); // for linear interpolation
        d2=ceil(dnew[d][h]); // for linear interpolation
        ddec=dnew[d][h]-double(d1); // for linear interpolation
        W[t][d][h] = pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*(repro[d] + (1.0-ddec)*Wopt[0][d1]+ddec*Wopt[0][d2]) // survive attack
                    + (1.0-pPred[t]*pAttack)*(1.0-mu[d])*(repro[d] +(1.0-ddec)*Wopt[t][d1]+ddec*Wopt[t][d2]); // no attack
      }
    }
  }

}



/* OVERWRITE FITNESS ARRAY FROM PREVIOUS ITERATION */
void ReplaceFit()
{
  int t,h,d;
  double fitdiff;

  fitdiff = 0.0;

  for (t=1;t<maxT;t++)
  {
    for (d=0;d<=maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        fitdiff = fitdiff + abs(Wnext[t][d][h]-W[t][d][h]);
        Wnext[t][d][h] = W[t][d][h];
      }
    }
  }

  totfitdiff = fitdiff;
}



/* PRINT OUT OPTIMAL STRATEGY */
void PrintStrat()
{
  int t,d;

  outputfile << "t" << "\t" << "d" << "\t" << "hormone" << endl;

  for (t=0;t<maxT;t++)
  {
    for (d=0;d<=maxD;d++)
    {
      outputfile << t << "\t" << d << "\t" << hormone[t][d] << endl;
    }
  }
  outputfile << endl;
  outputfile << "nIterations" << "\t" << i << endl;
  outputfile << endl;
}




/* WRITE PARAMETER SETTINGS TO OUTPUT FILE */
void PrintParams()
{
  outputfile << endl << "PARAMETER VALUES" << endl
       << "pLeave: " << "\t" << pLeave << endl
       << "pArrive: " << "\t" << pArrive << endl
       << "pAttack: " << "\t" << pAttack << endl
       << "alpha: " << "\t" << alpha << endl
//       << "beta: " << "\t" << beta << endl
       << "mu0: " << "\t" << mu0 << endl
       << "Kmort: " << "\t" << Kmort << endl
       << "Kfec: " << "\t" << Kfec << endl
       << "maxI: " << "\t" << maxI << endl
       << "maxT: " << "\t" << maxT << endl
       << "maxD: " << "\t" << maxD << endl
       << "maxH: " << "\t" << maxH << endl;
}



/* FORWARD CALCULATION TO OBTAIN PER-TIME-STEP MORTALITY FROM STRESSOR VS. DAMAGE */
void fwdCalc()
{
  int t,d,h,d1,d2,h1,h2,i;
  double ddec,predDeaths,damageDeaths,maxfreqdiff;

  for (t=1;t<maxT;t++) // note that F is undefined for t=0 because t=1 if predator has just attacked
  {
    for (d=0;d<=maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        F[t][d][h] = 0.0;
        Fnext[t][d][h] = 0.0;
      }
    }
  }
  F[50][0][0] = 1.0; // initialise all individuals with zero damage, zero hormone and 50 time steps since last attack

  i = 0;
  maxfreqdiff = 1.0;
  while (maxfreqdiff > 0.000001)
  {
      i++;
      predDeaths = 0.0;
      damageDeaths = 0.0;
      for (t=1;t<maxT;t++) // note that F is undefined for t=0 because t=1 if predator has just attacked
      {
        for (d=0;d<=maxD;d++)
        {
          for (h=0;h<maxH;h++)
          {
            d1=floor(dnew[d][h]); // for linear interpolation
            d2=ceil(dnew[d][h]); // for linear interpolation
            ddec=dnew[d][h]-double(d1); // for linear interpolation
            // attack
            h1=hormone[0][d1];
            Fnext[1][d1][h1] += F[t][d][h]*pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*(1.0-ddec);
            h2=hormone[0][d2];
            Fnext[1][d2][h2] += F[t][d][h]*pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*ddec;
            // no attack
            h1=hormone[min(maxT-1,t+1)][d1];
            Fnext[min(maxT-1,t+1)][d1][h1] += F[t][d][h]*(1.0-pPred[t]*pAttack)*(1.0-mu[d])*(1.0-ddec);
            h2=hormone[min(maxT-1,t+1)][d2];
            Fnext[min(maxT-1,t+1)][d2][h2] += F[t][d][h]*(1.0-pPred[t]*pAttack)*(1.0-mu[d])*ddec;
            // deaths from predation and damage
            predDeaths += F[t][d][h]*pPred[t]*pAttack*pKilled[h];
            damageDeaths += F[t][d][h]*(1.0-pPred[t]*pAttack*pKilled[h])*mu[d];
          }
        }
      }

      // NORMALISE AND OVERWRITE FREQUENCIES
      Fnext[1][0][0] = Fnext[1][0][0]/(1.0-predDeaths-damageDeaths); // normalise
      maxfreqdiff = abs(F[1][0][0] - Fnext[1][0][0]);
      for (t=1;t<maxT;t++)
      {
        for (d=0;d<=maxD;d++)
        {
          for (h=0;h<maxH;h++)
          {
            Fnext[t][d][h] = Fnext[t][d][h]/(1.0-predDeaths-damageDeaths); // normalise
            maxfreqdiff = max(maxfreqdiff,abs(F[t][d][h]-Fnext[t][d][h])); // stores largest frequency difference so far
            F[t][d][h] = Fnext[t][d][h]; // next time step becomes this time step
            Fnext[t][d][h] = 0.0; // wipe next time step
          }
        }
      }
      if (i%skip==0)
      {
        cout << i << "\t" << maxfreqdiff << endl; // show fitness difference every 'skip' generations
      }

  }

  ///////////////////////////////////////////////////////
  outfile.str("");
  outfile << "fwdCalcL";
  outfile << std::fixed << pLeave;
  outfile << "A";
  outfile << std::fixed << pArrive;
  outfile << "Kmort";
  outfile << std::fixed << Kmort;
  outfile << "Kfec";
  outfile << std::fixed << Kfec;
  outfile << ".txt";
  string fwdCalcfilename = outfile.str();
  fwdCalcfile.open(fwdCalcfilename.c_str());
  ///////////////////////////////////////////////////////

  fwdCalcfile << "SUMMARY STATS" << endl
    << "predDeaths: " << "\t" << predDeaths << endl
    << "damageDeaths: " << "\t" << damageDeaths << endl
    << endl;

  fwdCalcfile << "\t" << "t" << "\t" << "damage" << "\t" << "hormone" << "\t" << //"repro" << "\t" <<
    "freq" << endl; // column headings in output file

  for (t=1;t<maxT;t++)
  {
    for (d=0;d<=maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        fwdCalcfile << "\t" << t << "\t" << d << "\t" << h << "\t" << setprecision(4) << F[t][d][h] << "\t" << endl; // print data
      }
    }
  }

  fwdCalcfile.close();
}



/* Simulated series of attacks */
void SimAttacks()
{
  int i,time,t,d,h,d1,d2;
  double //r,
    ddec;
  bool attack;

  ///////////////////////////////////////////////////////
  outfile.str("");
  outfile << "simAttacksL";
  outfile << std::fixed << pLeave;
  outfile << "A";
  outfile << std::fixed << pArrive;
  outfile << "Kmort";
  outfile << std::fixed << Kmort;
  outfile << "Kfec";
  outfile << std::fixed << Kfec;
  outfile << ".txt";
  string attsimfilename = outfile.str();
  attsimfile.open(attsimfilename.c_str());
  ///////////////////////////////////////////////////////

  attsimfile << "time" << "\t" << "t" << "\t" << "damage" << "\t" << "hormone" << "\t" << "attack" << endl; // column headings in output file

    // initialise individual (alive, no damage, no offspring, baseline hormone level) and starting environment (predator)
    attack = false;
    time = 0;
    t = 50;
    d = 0;
    //    r = 0.0;
    h = hormone[t][d];

    while (time <= 60)
    {
      if (time > 16 && time < 33) // predator attacks
      {
        t = 0;
        attack = true;
      }
      else
      {
        t++;
        attack = false;
      }
      h = hormone[t][d];
      d1 = floor(dnew[d][h]);
      d2 = ceil(dnew[d][h]);
      ddec = dnew[d][h]-d1;
      if (Uniform(mt)<ddec) d = d2; else d = d1;
      attsimfile << time << "\t" << t << "\t" << d << "\t" << h << "\t" << attack << "\t" << endl; // print data
      time++;
    }

  attsimfile.close();

}



/* MAIN PROGRAM */
int main()
{

    int xGrid,yGrid;
    double risk,autocorr;

    for(xGrid=1;xGrid<=3;xGrid++)
      {
      if (xGrid == 1) risk = 0.05;
      if (xGrid == 2) risk = 0.1;
      if (xGrid == 3) risk = 0.2;
      for(yGrid=1;yGrid<=6;yGrid++)
        {
        if (yGrid == 1) autocorr = 0.0;
        if (yGrid == 2) autocorr = 0.1;
        if (yGrid == 3) autocorr = 0.3;
        if (yGrid == 4) autocorr = 0.5;
        if (yGrid == 5) autocorr = 0.7;
        if (yGrid == 6) autocorr = 0.9;

        pLeave = (1.0 - autocorr)/(1.0+(risk/(1.0-risk)));
        pArrive = 1.0 - pLeave - autocorr;
//    pLeave= 0.095;
//    pArrive = 0.005;

		///////////////////////////////////////////////////////
		outfile.str("");
		outfile << "stressL";
		outfile << std::fixed << pLeave;
		outfile << "A";
		outfile << std::fixed << pArrive;
		outfile << "Kmort";
		outfile << std::fixed << Kmort;
		outfile << "Kfec";
		outfile << std::fixed << Kfec;
		outfile << ".txt";
		string outputfilename = outfile.str();
		outputfile.open(outputfilename.c_str());
		///////////////////////////////////////////////////////

        outputfile << "Random seed: " << seed << endl; // write seed to output file

        // initialize arrays
        FinalFit();
        PredProb();
        Predation();
        Mortality();
        Damage();
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

        fwdCalc();
        SimAttacks();

        }
      }

  return 0;
}
