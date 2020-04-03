#include "LTGAOriginal.h"

#include "CommandParam.h"
#include "PopulationOptimizer.h"
#include "RandUtils.h"
#include "UIntCommandParam.h"

#include <algorithm>

CLTGAOriginal::CLTGAOriginal(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcProblem, pcLog, iRandomSeed), c_ltga(pcProblem->pcGetEvaluation(), pcLog, (int64_t)RandUtils::iRandNumber((uint32_t)0, UINT32_MAX - 1))
{
	v_ltga_individual.resize((size_t)pcProblem->pcGetEvaluation()->iGetNumberOfElements());
}//CLTGAOriginal::CLTGAOriginal(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

CLTGAOriginal::CLTGAOriginal(CLTGAOriginal *pcOther)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcOther), c_ltga(pcOther->pc_problem->pcGetEvaluation(), pcOther->pc_log, (int64_t)RandUtils::iRandNumber((uint32_t)0, UINT32_MAX - 1))
{
	c_ltga.setPopulationSize(pcOther->c_ltga.getPopulationSize());
	c_ltga.setDoLocalSearch(pcOther->c_ltga.getDoLocalSearch());
	c_ltga.setLocalSearchOneIteration(pcOther->c_ltga.getLocalSearchOneIteration());
	c_ltga.setWithoutTournamentSelection(pcOther->c_ltga.getWithoutTournamentSelection());
	c_ltga.setLinkageTreeRandomOrder(pcOther->c_ltga.getLinkageTreeRandomOrder());

	v_ltga_individual.resize(pcOther->v_ltga_individual.size());
}//CLTGAOriginal::CLTGAOriginal(CLTGAOriginal *pcOther)

CError CLTGAOriginal::eConfigure(istream *psSettings)
{
	CError c_error = COptimizer<CBinaryCoding, CBinaryCoding>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_population_size(POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE);
		vSetPopulationSize(p_population_size.iGetValue(psSettings, &c_error));
	}//if (!c_error)

	if (!c_error)
	{
		CBoolCommandParam p_do_local_search(LTGA_ORIGINAL_ARGUMENT_DO_LOCAL_SEARCH);
		c_ltga.setDoLocalSearch(p_do_local_search.bGetValue(psSettings, &c_error));
	}//if (!c_error)

	if (!c_error)
	{
		if (c_ltga.getDoLocalSearch())
		{
			CBoolCommandParam p_local_search_one_iteration(LTGA_ORIGINAL_ARGUMENT_LOCAL_SEARCH_ONE_ITERATION);
			c_ltga.setLocalSearchOneIteration(p_local_search_one_iteration.bGetValue(psSettings, &c_error));
		}//if (c_ltga.getDoLocalSearch())
	}//if (!c_error)

	if (!c_error)
	{
		CBoolCommandParam p_without_tournament_selection(LTGA_ORIGINAL_ARGUMENT_WITHOUT_TOURNAMENT_SELECTION, false, false);
		c_ltga.setWithoutTournamentSelection(p_without_tournament_selection.bGetValue(psSettings, &c_error));
	}//if (!c_error)

	if (!c_error)
	{
		CBoolCommandParam p_linkage_tree_random_order(LTGA_ORIGINAL_ARGUMENT_LINKAGE_TREE_RANDOM_ORDER, false, false);
		c_ltga.setLinkageTreeRandomOrder(p_linkage_tree_random_order.bGetValue(psSettings, &c_error));
	}//if (!c_error)

	return c_error;
}//CError CLTGAOriginal::eConfigure(istream *psSettings)

void CLTGAOriginal::vInitialize(time_t tStartTime)
{
	COptimizer<CBinaryCoding, CBinaryCoding>::vInitialize(tStartTime);

	c_ltga.initialize();
	c_ltga.updateBestPrevGenSolution();

	b_update_best_individual(0, tStartTime);
}//void CLTGAOriginal::vInitialize(time_t tStartTime)

bool CLTGAOriginal::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	return bRunIteration(iIterationNumber, tStartTime, nullptr);
}//bool CLTGAOriginal::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)

bool CLTGAOriginal::bRunIteration(uint32_t iIterationNumber, time_t tStartTime, CIndividual<CBinaryCoding, CBinaryCoding> *pcBestIndividual)
{
	if (pcBestIndividual)
	{
		for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
		{
			v_ltga_individual.at(i) = *(pcBestIndividual->pcGetGenotype()->piGetBits() + i) == 1 ? TRUE : FALSE;
		}//for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)

		c_ltga.runGeneration(v_ltga_individual.data());
	}//if (pcBestIndividual)
	else
	{
		c_ltga.runGeneration(nullptr);
	}//else if (pcBestIndividual)


	bool b_updated = b_update_best_individual(iIterationNumber, tStartTime);

	CString s_log_message;
	
	/*
	s_log_message.Format("best fitness: %f; best unitation: %f; ffe: %u; time: %u; population size: %u; avg fitness: %f; all the same: %d",
		pc_best_individual->dGetFitnessValue(), pc_best_individual->pcGetGenotype()->dGetUnitation(), pc_problem->pcGetEvaluation()->iGetFFE(), 
		(uint32_t)(time(nullptr) - tStartTime), iGetPopulationSize(), dComputeAverageFitnessValue(), bAreAllIndividualsTheSame());
	*/

	s_log_message.Format("best fitness: %f; best unitation: %f; ffe: %u; time: %u; population size: %u",
		pc_best_individual->dGetFitnessValue(), pc_best_individual->pcGetGenotype()->dGetUnitation(), pc_problem->pcGetEvaluation()->iGetFFE(),
		(uint32_t)(time(nullptr) - tStartTime), iGetPopulationSize());

	pc_log->vPrintLine(s_log_message, true);


	return b_updated;
}//bool CLTGAOriginal::bRunIteration(uint32_t iIterationNumber, time_t tStartTime, CIndividual<CBinaryCoding, CBinaryCoding> *pcBestIndividual)

void CLTGAOriginal::vRun()
{
	COptimizer<CBinaryCoding, CBinaryCoding>::vRun();
	c_ltga.ezilaitiniMemory();
}//void CLTGAOriginal::vRun()

bool CLTGAOriginal::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_updated = b_update_best_individual(iIterationNumber, tStartTime, c_ltga.getBestEverEvaluatedObjectiveValue(), [&](CBinaryCoding *pcBestGenotype)
	{
		for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
		{
			*(pcBestGenotype->piGetBits() + i) = *(c_ltga.getBestEverEvaluatedSolution() + i) == TRUE ? 1 : 0;
		}//for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
	});//bool b_updated = b_update_best_individual(iIterationNumber, tStartTime, c_ltga.getBestEverEvaluatedObjectiveValue(), [&](CBinaryCoding *pcBestGenotype)

	return b_updated;
}//bool CLTGAOriginal::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)