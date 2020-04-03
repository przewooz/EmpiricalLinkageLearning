#include "LTGAOriginalPopulationSizing.h"

#include "UIntCommandParam.h"

CLTGAOriginalPopulationSizing::CLTGAOriginalPopulationSizing(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcProblem, pcLog, iRandomSeed), c_params_ltga(pcProblem, pcLog, iRandomSeed)
{

}//CLTGAOriginalPopulationSizing::CLTGAOriginalPopulationSizing(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

CLTGAOriginalPopulationSizing::CLTGAOriginalPopulationSizing(CLTGAOriginalPopulationSizing *pcOther)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcOther), c_params_ltga(pcOther->c_params_ltga)
{

}//CLTGAOriginalPopulationSizing::CLTGAOriginalPopulationSizing(CLTGAOriginalPopulationSizing *pcOther)

CLTGAOriginalPopulationSizing::~CLTGAOriginalPopulationSizing()
{
	v_clear_ltgas();
}//CLTGAOriginalPopulationSizing::~CLTGAOriginalPopulationSizing()

CError CLTGAOriginalPopulationSizing::eConfigure(istream *psSettings)
{
	CError c_error = COptimizer<CBinaryCoding, CBinaryCoding>::eConfigure(psSettings);

	if (!c_error)
	{
		c_error = c_params_ltga.eConfigure(psSettings);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_population_creation_frequency(LTGA_ORIGINAL_POPULATION_SIZING_ARGUMENT_POPULATION_CREATION_FREQUENCY);
		i_population_creation_frequency = p_population_creation_frequency.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_population_size_multiplier(LTGA_ORIGINAL_POPULATION_SIZING_ARGUMENT_POPULATION_SIZE_MULTIPLIER);
		i_population_size_multiplier = p_population_size_multiplier.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CLTGAOriginalPopulationSizing::eConfigure(istream *psSettings)

void CLTGAOriginalPopulationSizing::vInitialize(time_t tStartTime)
{
	COptimizer<CBinaryCoding, CBinaryCoding>::vInitialize(tStartTime);

	v_clear_ltgas();

	i_next_population_size = c_params_ltga.iGetPopulationSize();
	v_add_new_ltga(0, tStartTime);

	b_update_best_individual(0, tStartTime);
}//void CLTGAOriginalPopulationSizing::vInitialize(time_t tStartTime)

bool CLTGAOriginalPopulationSizing::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	v_handle_ltgas(iIterationNumber, tStartTime);

	bool b_updated = false;

	if (b_run_ltgas_iteration(iIterationNumber, tStartTime))
	{
		b_updated = b_update_best_individual(iIterationNumber, tStartTime);
	}//if (b_run_ltgas_iteration(iIterationNumber, tStartTime))


	CString s_log;

	s_log.Format("iteration: %d; time: %u; number of ltgas: %d; best: %f", iIterationNumber, (uint32_t)(time(nullptr) - tStartTime),
		v_ltgas.size(), pc_best_individual->dGetFitnessValue());

	pc_log->vPrintLine(s_log, true);
	pc_log->vPrintEmptyLine(true);


	return b_updated;
}//bool CLTGAOriginalPopulationSizing::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)

bool CLTGAOriginalPopulationSizing::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)
{
	CLTGAOriginal *pc_ltga;
	CIndividual<CBinaryCoding, CBinaryCoding> *pc_ltga_best_individual;

	bool b_updated = false;

	for (size_t i = 0; i < v_ltgas.size(); i++)
	{
		pc_ltga = v_ltgas.at(i);
		pc_ltga_best_individual = pc_ltga->pcGetBestIndividual();

		if (pc_ltga_best_individual)
		{
			if (b_update_best_individual(iIterationNumber, tStartTime, pc_ltga_best_individual))
			{
				b_updated = true;
			}//if (b_update_best_individual(iIterationNumber, tStartTime, pc_ltga_best_individual))
		}//if (pc_ltga_best_individual)
	}//for (size_t i = 0; i < v_ltgas.size(); i++)

	return b_updated;
}//bool CLTGAOriginalPopulationSizing::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)

void CLTGAOriginalPopulationSizing::v_clear_ltgas()
{
	for (size_t i = 0; i < v_ltgas.size(); i++)
	{
		delete v_ltgas.at(i);
	}//for (size_t i = 0; i < v_ltgas.size(); i++)

	v_ltgas.clear();
}//void CLTGAOriginalPopulationSizing::v_clear_ltgas()

bool CLTGAOriginalPopulationSizing::b_run_ltgas_iteration(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_at_least_one_improvement = false;

	for (size_t i = 0; i < v_ltgas.size(); i++)
	{
		if (v_ltgas.at(i)->bRunIteration(iIterationNumber, tStartTime, pc_best_individual))
		{
			b_at_least_one_improvement = true;
		}//if (v_ltgas.at(i)->bRunIteration(iIterationNumber, tStartTime, pc_best_individual))
	}//for (size_t i = 0; i < v_ltgas.size(); i++)

	return b_at_least_one_improvement;
}//bool CLTGAOriginalPopulationSizing::b_run_ltgas_iteration(uint32_t iIterationNumber, time_t tStartTime)

void CLTGAOriginalPopulationSizing::v_handle_ltgas(uint32_t iIterationNumber, time_t tStartTime)
{
	v_ltgas_deletion();
	v_ltgas_creation(iIterationNumber, tStartTime);
}//void CLTGAOriginalPopulationSizing::v_handle_ltgas(uint32_t iIterationNumber, time_t tStartTime)

void CLTGAOriginalPopulationSizing::v_ltgas_deletion()
{
	unordered_set<size_t> s_ltgas_to_delete_indexes;
	s_ltgas_to_delete_indexes.reserve(v_ltgas.size());

	v_get_ltgas_to_delete(&s_ltgas_to_delete_indexes);

	if (!s_ltgas_to_delete_indexes.empty())
	{
		for (size_t i = v_ltgas.size(); i > 0; i--)
		{
			if (s_ltgas_to_delete_indexes.count(i - 1) > 0)
			{
				delete v_ltgas.at(i - 1);
				v_ltgas.erase(v_ltgas.begin() + i - 1);
			}//if (s_ltgas_to_delete_indexes.count(i - 1) > 0)
		}//for (size_t i = v_ltgas.size(); i > 0; i--)
	}//if (!s_ltgas_to_delete_indexes.empty())
}//void CLTGAOriginalPopulationSizing::v_ltgas_deletion()

void CLTGAOriginalPopulationSizing::v_ltgas_creation(uint32_t iIterationNumber, time_t tStartTime)
{
	if (v_ltgas.empty() || iIterationNumber - i_last_population_creation_iteration_number == i_population_creation_frequency)
	{
		v_add_new_ltga(iIterationNumber, tStartTime);
	}//if (v_ltgas.empty() || iIterationNumber - i_last_population_creation_iteration_number == i_population_creation_frequency)
}//void CLTGAOriginalPopulationSizing::v_ltgas_creation(uint32_t iIterationNumber, time_t tStartTime)

void CLTGAOriginalPopulationSizing::v_get_ltgas_to_delete(unordered_set<size_t> *psIndexes)
{
	double *pd_average_fitnesses = new double[v_ltgas.size()];

	for (size_t i = 0; i < v_ltgas.size(); i++)
	{
		*(pd_average_fitnesses + i) = v_ltgas.at(i)->dComputeAverageFitnessValue();

		if (v_ltgas.at(i)->bAreAllIndividualsTheSame())
		{
			psIndexes->insert(i);
		}//if (v_ltgas.at(i)->bAreAllIndividualsTheSame())
	}//for (size_t i = 0; i < v_ltgas.size(); i++)

	bool b_to_delete;

	for (size_t i = 0; i < v_ltgas.size(); i++)
	{
		if (psIndexes->count(i) == 0)
		{
			b_to_delete = false;

			for (size_t j = i + 1; j < v_ltgas.size() && !b_to_delete; j++)
			{
				b_to_delete = *(pd_average_fitnesses + i) < *(pd_average_fitnesses + j);
			}//for (size_t j = i + 1; j < v_ltgas.size() && !b_to_delete; j++)

			if (b_to_delete)
			{
				for (size_t j = 0; j <= i; j++)
				{
					psIndexes->insert(j);
				}//for (size_t j = 0; j <= i; j++)
			}//if (b_to_delete)
		}//if (psIndexes->count(i) == 0)
	}//for (size_t i = 0; i < v_ltgas.size(); i++)

	delete pd_average_fitnesses;
}//void CLTGAOriginalPopulationSizing::v_get_ltgas_to_delete(unordered_set<size_t> *psIndexes)

void CLTGAOriginalPopulationSizing::v_add_new_ltga(uint32_t iIterationNumber, time_t tStartTime)
{
	CLTGAOriginal *pc_ltga = new CLTGAOriginal(&c_params_ltga);

	pc_ltga->vSetPopulationSize(i_next_population_size);
	pc_ltga->vInitialize(tStartTime);
	
	v_ltgas.push_back(pc_ltga);

	i_next_population_size *= i_population_size_multiplier;
	i_last_population_creation_iteration_number = iIterationNumber;
}//void CLTGAOriginalPopulationSizing::v_add_new_ltga(uint32_t iIterationNumber, time_t tStartTime)