#ifndef LTGA_ORIGINAL_H
#define LTGA_ORIGINAL_H

#define LTGA_ORIGINAL_ARGUMENT_DO_LOCAL_SEARCH "do_local_search"
#define LTGA_ORIGINAL_ARGUMENT_LOCAL_SEARCH_ONE_ITERATION "local_search_one_iteration"
#define LTGA_ORIGINAL_ARGUMENT_WITHOUT_TOURNAMENT_SELECTION "without_tournament_selection"
#define LTGA_ORIGINAL_ARGUMENT_LINKAGE_TREE_RANDOM_ORDER "linkage_tree_random_order"

#include "BinaryCoding.h"
#include "Error.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"

#include "../LTGA/LTGA.h"

#include <ctime>
#include <cstdint>
#include <istream>
#include <vector>

using namespace std;

class CLTGAOriginal : public COptimizer<CBinaryCoding, CBinaryCoding>
{
public:
	CLTGAOriginal(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CLTGAOriginal(CLTGAOriginal *pcOther);

	virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CLTGAOriginal(this); };

	virtual CError eConfigure(istream *psSettings);

	virtual void vInitialize(time_t tStartTime);
	virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);
	bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime, CIndividual<CBinaryCoding, CBinaryCoding> *pcBestIndividual);

	virtual void vRun();

	bool bAreAllIndividualsTheSame() { return c_ltga.areAllSolutionsTheSame(); }
	double dComputeAverageFitnessValue() { return c_ltga.computeAverageObjectiveValue(); }

	uint32_t iGetPopulationSize() { return (uint32_t)c_ltga.getPopulationSize(); }
	void vSetPopulationSize(uint32_t iPopulationSize) { c_ltga.setPopulationSize((int)iPopulationSize); }

private:
	using COptimizer<CBinaryCoding, CBinaryCoding>::b_update_best_individual;
	bool b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime);

	NLTGA::LTGA c_ltga;
	vector<char> v_ltga_individual;
};//class CLTGAOriginal : public COptimizer<CBinaryCoding, CBinaryCoding>

#endif//LTGA_ORIGINAL_H