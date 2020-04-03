/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include <list>
#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "trimatrixpair.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"

#include "../EATester/Log.h"

#include <ctime>

class DSMGA2 {
public:
    DSMGA2 (int n_ell, int n_nInitial, bool n_twoEdge, bool n_useComparativeMixing, int n_maxGen, int n_maxFe, int fffff, time_t n_maxTime, CLog* n_logger);

    ~DSMGA2 ();

    void selection ();
    /** tournament selection without replacement*/
    void tournamentSelection();

    void oneRun (int iteration, bool output = true);
    int doIt (bool output = true);

    void buildGraph ();
    void mixing ();
	void findMask(Chromosome& ch, list<int>& mask, int startNode);
	void findMask_size(Chromosome& ch, list<int>& mask, int startNode, int bound);
	void buildGraph_sizecheck();
    void restrictedMixing(Chromosome&);
    bool restrictedMixing(Chromosome& ch, list<int>& mask);
    void backMixing(Chromosome& source, list<int>& mask, Chromosome& des);
    void backMixingE(Chromosome& source, list<int>& mask, Chromosome& des);

    bool shouldTerminate ();

    bool foundOptima ();

    int getGeneration () const {
        return generation;
    }

    bool isInP(const Chromosome& ) const;
    void genOrderN();
    void genOrderELL();

    void showStatistics ();

    bool isSteadyState ();

	void getResult(int* bestFFE, time_t* bestTime, int* bestSolution);

//protected:
public:

    int ell;                                  // chromosome length
    int nCurrent;                             // population size
    bool EQ;
    unordered_map<unsigned long, double> pHash; // to check if a chromosome is in the population

	bool twoEdge;
	bool useComparativeMixing;

    list<int> *masks;
    int *selectionIndex;
    int *orderN;                             // for random order
    int *orderELL;                             // for random order
    int selectionPressure;
    int maxGen;
    int maxFe;
	time_t maxTime;
    int repeat;
    int generation;
	time_t startTime;
    int bestIndex;

    Chromosome* population;
    FastCounting* fastCounting;

    TriMatrix<double> graph;

	TriMatrixPair graphTwoEdge;
	TriMatrixPair graphSizeTwoEdge;

    double previousFitnessMean;
    Statistics stFitness;

    // methods
    double computeMI(double, double, double, double) const;


    void findClique(int startNode, list<int>& result);

    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;

    size_t findSize(Chromosome&, list<int>&) const;
    size_t findSize(Chromosome&, list<int>&, Chromosome&) const;

	bool converged();
	double lastMax, lastMean, lastMin;
	int convergeCount;

	double bestFitness;
	int bestFFE;
	time_t bestTime;
	int* bestSolution;

private:
	CLog* logger;
};


#endif /* _DSMGA2_H_ */
