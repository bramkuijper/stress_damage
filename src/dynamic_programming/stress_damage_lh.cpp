// **********************************************************************************
// Dynamic programming model of stress response with somatic damage.
//
// fecundity and seasonality
//
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
#include <cassert>

// constants, type definitions, etc.
const int seed        = std::time(0); // pseudo-random seed

double pLeave;   // probability that predator leaves
double pArrive;  // probability that predator arrives
double pAttack  = 0.5;     // probability that predator attacks if present
double alpha    = 1.0;     // parameter controlling effect of hormone level on pKilled
//const double beta     = 1.5;     // parameter controlling effect of hormone level on reproductive rate
const double mu0      = 0.002;   // background mortality (independent of hormone level and predation risk)
const double phi_inv  = 1.0/((sqrt(5.0)+1.0)/2.0); // inverse of golden ratio (for golden section search)
const double hmin     = 0.3;     // level of hormone that minimises damage
const double hslope   = 20.0;     // slope parameter controlling increase in damage with deviation from hmin
const int maxD        = 1000;      // maximum damage level
int repair      = 1;      // damage units removed per time step
double Kmort        = 0.0;    // parameter Kmort controlling increase in mortality with damage level
double Kfec        = 0.05;    // parameter Kmort controlling increase in mortality with damage level
const int maxI        = 1000000; // maximum number of iterations
const int maxT        = 50;     // maximum number of time steps since last saw predator
const int maxH        = 100;     // maximum hormone level
const int skip        = 1;      // interval between print-outs
const int maxTs = 20; // duration of a season

int n_repro_bout_length = 1;

std::ofstream outputfile;  // output file
std::ofstream fwdCalcfile; // forward calculation output file
std::ofstream attsimfile;  // simulated attacks output file
std::stringstream outfile; // for naming output file
std::string base_name = "";

// random numbers
std::mt19937 mt(seed); // random number generator
std::uniform_real_distribution<double> Uniform(0, 1); // real number between 0 and 1 (uniform)

///int hormone[maxT+1][maxTs+1][maxD+1];        // hormone level (strategy)
std::vector < std::vector < std::vector<int> > > 
    hormone(maxT+1, std::vector < std::vector <int> >(
                maxTs+1, std::vector<int>(maxD + 1, 0)));

double pKilled[maxH+1];             // probability of being killed by an attacking predator
double mu[maxD+1];                // probability of background mortality, as a function of damage
double dnew[maxD+1][maxH+1];        // new damage level, as a function of previous damage and hormone
double repro[maxTs+1][maxD+1];       // reproductive output

// for the dimensions used here 
// there is not enough memory
// to allocate global C-style arrays 
// (they are allocated nor on the
// heap nor on the stack)
// yet, allocation through vectors works
// I guess because these are allocated on the stack (??)
// anyway, it works.
//
// Wopt[maxT+1][maxTs+1][maxD+1]
std::vector < std::vector < std::vector<double> > > 
    Wopt(maxT+1, std::vector < std::vector <double> >(
                maxTs+1, std::vector<double>(maxD+1, 0.0)));

// W[maxT+1][maxTs+1][maxD+1]
std::vector < std::vector < std::vector < std::vector <double> > > > 
    W(maxT+1, std::vector < std::vector < std::vector <double> > >(
                maxTs+1, std::vector < std::vector <double> >(
                    maxD+1, std::vector<double>(maxH+1, 0.0))));

// reproductive values V    
// V[maxT+1][maxD+1][maxH+1]
std::vector < std::vector < std::vector<double> > > 
    V(maxT+1, std::vector < std::vector <double> >(
                maxD+1, std::vector<double>(maxH+1, 0.0)));

// Wnext[maxT+1][maxTs+1][maxD+1][maxH+1]
std::vector < std::vector < std::vector < std::vector <double> > > > 
    Wnext(maxT+1, std::vector < std::vector < std::vector <double> > >(
                maxTs+1, std::vector < std::vector <double> >(
                    maxD+1, std::vector<double>(maxH+1, 0.0))));

// F[maxT+1][maxTs+1][maxD+1][maxH+1]
std::vector < std::vector < std::vector < std::vector <double> > > > 
    F(maxT+1, std::vector < std::vector < std::vector <double> > >(
                maxTs+1, std::vector < std::vector <double> >(
                    maxD+1, std::vector<double>(maxH+1, 0.0))));


    
double pPred[maxT+1];               // probability that predator is present
double totfitdiff;                // fitness difference between optimal strategy in successive iterations

int i;     // iteration


/* SPECIFY FINAL FITNESS */
void FinalFit()
{
  int t,d,h;

    for (t=1;t<=maxT;++t) // note that Wnext is undefined for t=0 because t=1 if predator has just attacked
    {
        for (d=0;d<=maxD;++d)
        {
            for (h=0;h<=maxH;++h)
            {
               V[t][d][h] =  Wnext[t][maxTs][d][h] = repro[maxTs][d];
            }
        }
    }
} // end FinalFit()


/* CALCULATE PROBABILITY THAT PREDATOR IS PRESENT */
void PredProb()
{
  int t;

  pPred[1] = 1.0 - pLeave; // if predator attacked in last time step

  for (t=2;t<=maxT;++t) // if predator did NOT attack in last time step
  {
//      Pr(predator present at time t | predator did not attack at time t-1)
//      = Pr (predator did not attack at time t-1 | predator present at time t) * 
//              Pr( predator present at time t)  / Pr(predator did not attack at time t-1)
    pPred[t] = (pPred[t-1]*(1.0-pAttack)*(1.0-pLeave)+(1.0-pPred[t-1])*pArrive) / (1.0 - pPred[t-1]*pAttack);
  }
} // end PredProb



/* CALCULATE PROBABILITY OF BEING KILLED BY AN ATTACKING PREDATOR */
void Predation()
{
  int h;

  for (h=0;h<=maxH;++h)
  {
    pKilled[h] = std::max(0.0,1.0 - pow(double(h)/double(maxH),alpha));
  }
} // end Predation()


/* CALCULATE BACKGROUND MORTALITY */
void Mortality()
{
  int d;

  for (d=0;d<=maxD;++d)
  {
    mu[d] = std::min(1.0,mu0 + Kmort*double(d));
  }

  mu[maxD]=1.0;
} // end Mortality()



/* CALCULATE DAMAGE */
void Damage()
{
  int d,h;

  for (d=0;d<=maxD;++d)
  {
    for (h=0;h<=maxH;++h)
    {
      dnew[d][h] = 
          std::max(0.0,
              std::min(double(maxD),
                  double(d) + hslope*pow(h/((double)maxH) - hmin,2) - double(repair)));
    }
  }
} // void Damage()


/* CALCULATE PROBABILITY OF REPRODUCING */
void Reproduction()
{
  int d, ts;

  for (ts = 0; ts <= maxTs; ++ts)
  { 
    for (d=0;d<=maxD;++d)
    {
        // reproductive event
        // in the middle between 0 and maxTs
        // lasting one timestep 
        repro[ts][d] = abs(ts % maxTs - floor(maxTs/2)) < n_repro_bout_length ? 
            std::max(0.0, 1.0 - Kfec*double(d)) 
            : 
            0.0;

        repro[ts][maxD] = 0.0;
    }
  }
} // end Reproduction()



/* CALCULATE OPTIMAL DECISION FOR EACH t */
void OptDec(int const iteration)
{
    // auxiliary variables
    int t,ts,h,d,d1,d2,LHS,RHS,x1,x2,cal_x1,cal_x2;
    double fitness,fitness_x1,fitness_x2,ddec;

    // go from maxTs down to 0
    // start from maxTs - 1, as we need to reach back
    // to array positions given by ts + 1
    for (ts = maxTs - 1; ts >= 0; --ts)
    {
      // calculate optimal decision h given current t, 
      // ts and d (N.B. t=0 if survived attack)
        for (t=0;t<=maxT;++t)
        {
            for (d=0;d<=maxD;++d)
            {
                double fit_max = 0;
                double wx;
                int hval = 0.0;

                for (int h_idx = 0; h_idx <= maxH; ++h_idx)
                {
                    wx = Wnext[std::min(maxT,t+1)][ts + 1][d][h_idx];

                    if (wx > fit_max)
                    {
                        fit_max = wx;
                        hval = h_idx;
                    }
                }
                
                // ts ranges here from MaxTs - 1 to 0
                // i.e., there are no hormone, Wopt values here for MaxTs
                hormone[t][ts][d] = hval; // optimal hormone level
                Wopt[t][ts][d] = fit_max; // fitness of optimal decision
            } // end for d
        } // end for t

          // calculate expected fitness W as a function of t, 
          // h and d, before predator does/doesn't attack
          // later on we will then set Wnext = W and see for 
          // which hormone level fitness is max


        for (t=1;t<=maxT;++t) // note that W is undefined for t=0 because t=1 if predator has just attacked
        {
            for (d=0;d<=maxD;++d)
            {
                for (h=0;h<=maxH;++h)
                {
                    d1=floor(dnew[d][h]); // for linear interpolation
                    d2=ceil(dnew[d][h]); // for linear interpolation
                    ddec=dnew[d][h]-double(d1); // for linear interpolation

                    W[t][ts][d][h] = pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*(
                            repro[ts][d] + 
                            (1.0-ddec)*Wopt[0][ts][d1]+ddec*Wopt[0][ts][d2]
                            ) // survive attack
                        + (1.0-pPred[t]*pAttack)*(1.0-mu[d])*(
                                    repro[ts][d] +
                                    (1.0-ddec)*Wopt[t][ts][d1]+ddec*Wopt[t][ts][d2]
                                ); // no attack

                    Wnext[t][ts][d][h] = W[t][ts][d][h];
                    
                } // end for h
            } // end for d
        } // end for t
    } // end for ts
} // end void OptDec()



/* OVERWRITE FITNESS ARRAY FROM PREVIOUS ITERATION */
void ReplaceFit()
{
    int t,h,d,ts;
    double fitdiff;

    fitdiff = 0.0;

    for (t=1;t<=maxT;t++)
    {
        for (d=0;d<=maxD;++d)
        {
            for (h=0;h<=maxH;++h)
            {
                fitdiff = fitdiff + fabs(V[t][d][h]-W[t][0][d][h]);

                Wnext[t][maxTs][d][h] = W[t][0][d][h];
                hormone[t][maxTs][d] = hormone[t][0][d];

                V[t][d][h] = W[t][0][d][h];
            }
        }
    }

    totfitdiff = fitdiff;
} // void ReplaceFit()



/* PRINT OUT OPTIMAL STRATEGY */
void PrintStrat()
{
  int t,d,ts;

  outputfile << "t" << "\t" << "d" << "\t" << "ts" << "\t" << "hormone" << std::endl;

  for (t=0;t<=maxT;++t)
  {
      for (ts = 0; ts <= maxTs; ++ts)
      {
        for (d=0;d<=maxD;++d)
            {
              outputfile << t << "\t" << d << "\t" << ts << "\t" << hormone[t][ts][d] << std::endl;
            }
      }
  }

  outputfile << std::endl;
  outputfile << "nIterations" << "\t" << i << std::endl;
  outputfile << std::endl;
}




/* WRITE PARAMETER SETTINGS TO OUTPUT FILE */
void PrintParams()
{
  outputfile << std::endl << "PARAMETER VALUES" << std::endl
       << "pLeave: " << "\t" << pLeave << std::endl
       << "pArrive: " << "\t" << pArrive << std::endl
       << "pAttack: " << "\t" << pAttack << std::endl
       << "alpha: " << "\t" << alpha << std::endl
//       << "beta: " << "\t" << beta << std::endl
       << "mu0: " << "\t" << mu0 << std::endl
       << "Kmort: " << "\t" << Kmort << std::endl
       << "Kfec: " << "\t" << Kfec << std::endl
       << "maxI: " << "\t" << maxI << std::endl
       << "maxT: " << "\t" << maxT << std::endl
       << "maxTs: " << "\t" << maxTs << std::endl
       << "maxD: " << "\t" << maxD << std::endl
       << "maxH: " << "\t" << maxH << std::endl
       << "hmin: " << "\t" << hmin << std::endl
       << "n_repro_bout_length: " << "\t" << n_repro_bout_length << std::endl
       << "hslope: " << "\t" << hslope << std::endl
       << "repair: " << "\t" << repair << std::endl;
}



/* FORWARD CALCULATION TO OBTAIN PER-TIME-STEP MORTALITY FROM STRESSOR VS. DAMAGE */
void fwdCalc()
{
  int t,ts,d,h,d1,d2,h1,h2,i;
  double ddec,predDeaths,damageDeaths,bkgrndDeaths,maxfreqdiff;

  for (t=1;t<=maxT;++t) // note that F is undefined for t=0 because t=1 if predator has just attacked
  {
      for (ts = 0; ts<=maxTs; ++ts)
      {
        for (d=0;d<=maxD;++d)
        {
          for (h=0;h<=maxH;++h)
          {
            F[t][ts][d][h] = 0.0;
          }
        }
      }
  }

  F[maxT][0][0][0] = 1.0; // initialise all individuals with zero damage, zero hormone and maxT time steps since last attack, during the first reproductive bout

  i = 0;

  maxfreqdiff = 1.0;

  while (maxfreqdiff > 0.000001)
  {
      i++;
      predDeaths = 0.0;
      damageDeaths = 0.0;
      bkgrndDeaths = 0.0;
      for (t=1;t<=maxT;++t) // note that F is undefined for t=0 because t=1 if predator has just attacked
      {
          for(ts=0; ts<maxTs;++ts)
          { 
            for (d=0;d<=maxD;++d)
            {
              for (h=0;h<=maxH;++h)
              {
                d1=floor(dnew[d][h]); // for linear interpolation
                d2=ceil(dnew[d][h]); // for linear interpolation
                ddec=dnew[d][h]-double(d1); // for linear interpolation

                // attack
                h1=hormone[0][ts][d1];
                F[1][ts+1][d1][h1] += F[t][ts][d][h]*pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*(1.0-ddec);
                h2=hormone[0][ts][d2];
                F[1][ts+1][d2][h2] += F[t][ts][d][h]*pPred[t]*pAttack*(1.0-pKilled[h])*(1.0-mu[d])*ddec;
                // no attack
                h1=hormone[std::min(maxT-1,t+1)][ts][d1];
                F[std::min(maxT-1,t+1)][ts+1][d1][h1] += F[t][ts][d][h]*(1.0-pPred[t]*pAttack)*(1.0-mu[d])*(1.0-ddec);
                h2=hormone[std::min(maxT-1,t+1)][ts][d2];
                F[std::min(maxT-1,t+1)][ts+1][d2][h2] += F[t][ts][d][h]*(1.0-pPred[t]*pAttack)*(1.0-mu[d])*ddec;
                // deaths from predation, damage and background
                predDeaths += F[t][ts][d][h]*pPred[t]*pAttack*pKilled[h];
                damageDeaths += F[t][ts][d][h]*(1.0-pPred[t]*pAttack*pKilled[h])*(mu[d]-mu[0]);
                bkgrndDeaths += F[t][ts][d][h]*(1.0-pPred[t]*pAttack*pKilled[h])*mu[0];
              } // end for h
            } // end for d
          } // end for ts
      } // end for t

      // NORMALISE AND COMPARE FREQUENCIES AT BREEDING POINT (ts = 0 = maxTs)
      double maxdiffcalc = 0.0;

      for (t=1;t<=maxT;t++)
      {
            for (d=0;d<=maxD;d++)
            {
              for (h=0;h<=maxH;h++)
              {
                // normalise
                F[t][maxTs][d][h] = F[t][maxTs][d][h]/(1.0-predDeaths-damageDeaths-bkgrndDeaths); 
                
                // stores largest frequency difference so far
                maxdiffcalc = std::max(maxdiffcalc,
                        fabs(F[t][maxTs][d][h]-F[t][0][d][h])); 

                // frequencies at time ts = maxTs become 
                // frequencies at time ts = 0 (next breeding cycle)
                F[t][0][d][h] = F[t][maxTs][d][h]; 

                for (ts = 1;ts<=maxTs;ts++)
                {
                    F[t][ts][d][h] = 0.0; // wipe frequencies for ts > 0 in next breeding cycle
                } // end for ts
              } // end for h
            } // end for d
      }
      if (i%skip==0)
      {
        std::cout << i << "\t" << maxfreqdiff << std::endl; // show fitness difference every 'skip' generations
      }

      maxfreqdiff = maxdiffcalc;
  } // end while 

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
  std::string fwdCalcfilename = outfile.str();
  fwdCalcfile.open(fwdCalcfilename.c_str());
  ///////////////////////////////////////////////////////

  fwdCalcfile << "SUMMARY STATS" << std::endl
    << "predDeaths: " << "\t" << predDeaths << std::endl
    << "damageDeaths: " << "\t" << damageDeaths << std::endl
    << std::endl;

  fwdCalcfile << "\t" << "t" << "\t" << "ts" << "\t" << "damage" << "\t" << "hormone" << "\t" << //"repro" << "\t" <<
    "freq" << std::endl; // column headings in output file

  for (t=1;t<=maxT;++t)
  {
      for (ts = 0; ts <= maxTs; ++ts)
      {
        for (d=0;d<=maxD;++d)
        {
          for (h=0;h<=maxH;++h)
          {
            fwdCalcfile << "\t" << t << "\t" << ts << "\t" << d << "\t" << h << "\t" << std::setprecision(4) << F[t][ts][d][h] << "\t" << std::endl; // print data
          }
        }
      }
  }

  fwdCalcfile.close();
}



/* Simulated series of attacks */
void SimAttacks(std::string const &base_name)
{
  int i,time_i,t,d,h,d1,d2;
  double //r,
    ddec;
  bool attack;

  std::string attsimfilename = "simAttacks" + base_name;
    attsimfilename += ".txt";
  attsimfile.open(attsimfilename.c_str());
  ///////////////////////////////////////////////////////

  attsimfile << "time" << "\t" << "t" << "\t" << "ts" << "\t" << "damage" << "\t" << "hormone" << "\t" << "attack" << "\t" << "reproduce" << "\t" << "d1" << "\t" << "d2" << "\t" << "ddec" << "\t" << std::endl; // column headings in output file

    // initialise individual (alive, no damage, no offspring, baseline hormone level) and starting environment (predator)
    //
    int time_sim_max = 300;

    attack = false;
    time_i = 0; // time overall
    t = maxT; // time since attack
   
    // time since last reproductive event
    // 1 timestep before treproduce and we start to count
    // time from 0
    int ts = 0;

    d = 0;

    int reproduce;

    for(time_i = 0; time_i < time_sim_max; ++time_i)
    {
      if (time_i > 16 && time_i < 33) // predator attacks repeatedly and after that does not attack
      {
        t = 0;
        attack = true;
      }
      else
      {
        t++;
        attack = false;
      }

      if (t >= maxT)
      {
          t = maxT;
      }

      h = hormone[t][ts % maxTs][d];

      d1 = floor(dnew[d][h]);
      d2 = ceil(dnew[d][h]);
      ddec = dnew[d][h]-d1;

      reproduce = abs(ts % maxTs - floor(maxTs/2)) < n_repro_bout_length;

      if (Uniform(mt)<ddec)
      {
          d = d2; 
      }
      else 
      {
          d = d1;
      }

      attsimfile << time_i << "\t" << t << "\t" << ts << "\t" << d << "\t" << h << "\t" << attack << "\t" << reproduce << "\t" << d1 << "\t" << d2 << "\t" << ddec << "\t" << std::endl; // print data


      ++ts;
    } // for time_i time sim max

  attsimfile.close();

} // end void SimAttacks()


void init_params(int argc, char** argv)
{
    pLeave = atof(argv[1]);
    pArrive = atof(argv[2]);
    pAttack = atof(argv[3]);
    alpha = atof(argv[4]);
    Kmort = atof(argv[5]);
    Kfec = atof(argv[6]);
    repair = atoi(argv[7]);
    n_repro_bout_length = atoi(argv[8]);
    base_name = argv[9];
} // end init_params() 

/* MAIN PROGRAM */
int main(int argc, char** argv)
{
    init_params(argc, argv);

		///////////////////////////////////////////////////////
//		outfile.str("");
//		outfile << "L";
//		outfile << std::fixed << pLeave;
//		outfile << "A";
//		outfile << std::fixed << pArrive;
//		outfile << "Kfec";
//		outfile << std::fixed << Kfec;
//		outfile << "Kmort";
//		outfile << std::fixed << Kmort;
//		outfile << "r";
//		outfile << std::fixed << repair;
//		outfile << "pAtt";
//		outfile << std::fixed << pAttack;
//		outfile << "bout";
//		outfile << std::fixed << n_repro_bout_length;

        std::string outputfilename = "stress_strategy" + base_name;
        outputfilename += ".txt";

		outputfile.open(outputfilename.c_str());
		///////////////////////////////////////////////////////

        outputfile << "Random seed: " << seed << std::endl; // write seed to output file

        Reproduction();
        FinalFit();
        PredProb();
        Predation();
        Mortality();
        Damage();

        std::cout << "i" << "\t" << "totfitdiff" << "\t" << std::endl;

        for (i=1;i<=maxI;++i)
        {
            OptDec(i);
            ReplaceFit();


 		  if (i%skip==0 || totfitdiff < 0.000001)
            {
              std::cout << i << "\t" << totfitdiff << std::endl; // show fitness difference every 'skip' generations
            }
          
          if (totfitdiff < 0.000001) break; // strategy has converged on optimal solution, so exit loop

          if (i==maxI) { outputfile << "*** DID NOT CONVERGE WITHIN " << i << " ITERATIONS ***" << std::endl;}
          }

        std::cout << std::endl;
        outputfile << std::endl;

        PrintStrat();
        PrintParams();
        outputfile.close();

//        fwdCalc();
        SimAttacks(base_name);

  return 0;
}
