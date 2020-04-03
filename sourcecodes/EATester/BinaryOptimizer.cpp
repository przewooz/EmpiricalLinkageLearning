#include "BinaryOptimizer.h"

CBinaryOptimizer::CBinaryOptimizer(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcProblem, pcLog, iRandomSeed)
{
	pc_evaluation_individual = pc_create_individual(new CBinaryCoding(pc_problem->pcGetEvaluation()->iGetNumberOfElements(), nullptr));
}//CBinaryOptimizer::CBinaryOptimizer(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

CBinaryOptimizer::CBinaryOptimizer(CBinaryOptimizer *pcOther)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcOther)
{
	pc_evaluation_individual = pcOther->pc_evaluation_individual;
}//CBinaryOptimizer::CBinaryOptimizer(CBinaryOptimizer *pcOther)

CBinaryOptimizer::~CBinaryOptimizer()
{
	if (b_own_params)
	{
		delete pc_evaluation_individual;
	}//if (b_own_params)
}//CBinaryOptimizer::~CBinaryOptimizer()

bool CBinaryOptimizer::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, int32_t *piBits, double dFitnessValue)
{
	return b_update_best_individual(iIterationNumber, tStartTime, dFitnessValue, [&](CBinaryCoding *pcBestGenotype)
	{
		for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
		{
			*(pcBestGenotype->piGetBits() + i) = *(piBits + i);
		}//for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
	});
}//bool CBinaryOptimizer::b_update_best_individual(int32_t *piBits, double dFitnessValue)

double CBinaryOptimizer::d_compute_fitness_value(int32_t *piBits)
{
	pc_evaluation_individual->pcGetGenotype()->vSetBits(piBits);
	pc_evaluation_individual->vIsEvaluated(false);
	pc_evaluation_individual->vEvaluate();
	
	return pc_evaluation_individual->dGetFitnessValue();
}//double CBinaryOptimizer::d_compute_fitness_value(int32_t *piBits)