#ifndef BINARY_OPTIMIZER
#define BINARY_OPTIMIZER

#include "BinaryCoding.h"
#include "Individual.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"

#include <ctime>
#include <cstdint>

using namespace std;

class CBinaryOptimizer : public COptimizer<CBinaryCoding, CBinaryCoding>
{
public:
	CBinaryOptimizer(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CBinaryOptimizer(CBinaryOptimizer *pcOther);

	virtual ~CBinaryOptimizer();

protected:
	using COptimizer<CBinaryCoding, CBinaryCoding>::b_update_best_individual;
	bool b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, int32_t *piBits, double dFitnessValue);

	double d_compute_fitness_value(int32_t *piBits);

private:
	CIndividual<CBinaryCoding, CBinaryCoding> *pc_evaluation_individual;
};//class CBinaryOptimizer : COptimizer<CBinaryCoding, CBinaryCoding>

#endif//BINARY_OPTIMIZER