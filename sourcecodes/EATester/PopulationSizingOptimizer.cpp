#include "PopulationSizingOptimizer.h"

#include "BinaryCoding.h"
#include "OptimizerUtils.h"
#include "UIntCommandParam.h"


CPopulationSizingCounter::CPopulationSizingCounter() : CPopulationSizingCounter(POPULATION_SIZING_COUNTER_DEFAULT_BASE)
{

}//CPopulationSizingCounter::CPopulationSizingCounter()

CPopulationSizingCounter::CPopulationSizingCounter(uint32_t iBase)
{
	vReset();
	bSetBase(iBase);
}//CPopulationSizingCounter::CPopulationSizingCounter(uint32_t iBase)

uint32_t CPopulationSizingCounter::iIncrement()
{
	uint32_t i_most_significant_change_index = 0;

	if (bIsReset())
	{
		v_numbers.push_back(0);
	}//if (bIsReset())
	else
	{
		for (uint32_t i = 0; i < (uint32_t)v_numbers.size() && i == i_most_significant_change_index; i++)
		{
			v_numbers.at(i)++;

			if (v_numbers.at(i) % i_base == 0)
			{
				v_numbers.at(i) = 0;
				i_most_significant_change_index++;
			}//if (v_numbers.at(i) % i_base == 0)
		}//for (uint32_t i = 0; i < (uint32_t)v_numbers.size() && i == i_most_significant_change_index; i++)

		if (i_most_significant_change_index == (uint32_t)v_numbers.size())
		{
			v_numbers.push_back(1);
		}//if (i_most_significant_change_index == (uint32_t)v_numbers.size())
	}//if (bIsReset())

	i_last_most_significant_change_index = i_most_significant_change_index;

	return i_most_significant_change_index;
}//uint32_t CPopulationSizingCounter::iIncrement()

void CPopulationSizingCounter::vReset()
{
	i_last_most_significant_change_index = 0;
	v_numbers.clear();
}//void CPopulationSizingCounter::vReset()

bool CPopulationSizingCounter::bSetBase(int iBase)
{
	bool b_is_set = false;

	if (bIsReset())
	{
		i_base = iBase;
		b_is_set = true;
	}//if (bIsReset())

	return b_is_set;
}//bool CPopulationSizingCounter::bSetBase(int iBase)


template <class TGenotype, class TFenotype>
CPopulationSizingSingleOptimizer<TGenotype, TFenotype>::CPopulationSizingSingleOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<TGenotype, TFenotype>(pcProblem, pcLog, iRandomSeed)
{

}//CPopulationSizingSingleOptimizer<TGenotype, TFenotype>::CPopulationSizingSingleOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

template <class TGenotype, class TFenotype>
CPopulationSizingSingleOptimizer<TGenotype, TFenotype>::CPopulationSizingSingleOptimizer(CPopulationSizingSingleOptimizer<TGenotype, TFenotype> *pcOther)
	: COptimizer<TGenotype, TFenotype>(pcOther)
{

}//CPopulationSizingSingleOptimizer<TGenotype, TFenotype>::CPopulationSizingSingleOptimizer(CPopulationSizingSingleOptimizer<TGenotype, TFenotype> *pcOther)


template <class TGenotype, class TFenotype>
CPopulationSizingOptimizer<TGenotype, TFenotype>::CPopulationSizingOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<TGenotype, TFenotype>(pcProblem, pcLog, iRandomSeed)
{
	pc_params_optimizer = nullptr;
}//CPopulationSizingOptimizer<TGenotype, TFenotype>::CPopulationSizingOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

template <class TGenotype, class TFenotype>
CPopulationSizingOptimizer<TGenotype, TFenotype>::CPopulationSizingOptimizer(CPopulationSizingOptimizer<TGenotype, TFenotype> *pcOther)
	: COptimizer<TGenotype, TFenotype>(pcOther)
{
	pc_params_optimizer = pcOther->pc_params_optimizer;
}//CPopulationSizingOptimizer<TGenotype, TFenotype>::CPopulationSizingOptimizer(CPopulationSizingOptimizer<TGenotype, TFenotype> *pcOther)

template <class TGenotype, class TFenotype>
CPopulationSizingOptimizer<TGenotype, TFenotype>::~CPopulationSizingOptimizer()
{
	v_clear_optimizers();
	v_clear_params();
}//CPopulationSizingOptimizer<TGenotype, TFenotype>::~CPopulationSizingOptimizer()

template <class TGenotype, class TFenotype>
CError CPopulationSizingOptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	v_clear_params();

	CError c_error = COptimizer<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		c_error = e_configure_params_optimizer(psSettings);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_population_creation_frequency(POPULATION_SIZING_OPTIMIZER_ARGUMENT_POPULATION_CREATION_FREQUENCY);
		i_population_creation_frequency = p_population_creation_frequency.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_population_size_multiplier(POPULATION_SIZING_OPTIMIZER_ARGUMENT_POPULATION_SIZE_MULTIPLIER);
		i_population_size_multiplier = p_population_size_multiplier.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CPopulationSizingOptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)
{
	COptimizer<TGenotype, TFenotype>::vInitialize(tStartTime);

	v_clear_optimizers();

	i_next_population_size = pc_params_optimizer->iGetPopulationSize();
	v_add_new_optimizer(0, tStartTime);

	b_update_best_individual(0, tStartTime);
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)

template <class TGenotype, class TFenotype>
bool CPopulationSizingOptimizer<TGenotype, TFenotype>::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	v_handle_optimizers(iIterationNumber, tStartTime);

	bool b_updated = false;

	if (b_run_optimizers_iteration(iIterationNumber, tStartTime))
	{
		b_updated = b_update_best_individual(iIterationNumber, tStartTime);
	}//if (b_run_optimizers_iteration(iIterationNumber, tStartTime))


	CString s_log;

	s_log.Format("iteration: %d; time: %u; number of optimizers: %d; best: %f", iIterationNumber, (uint32_t)(time(nullptr) - tStartTime),
		v_optimizers.size(), pc_best_individual->dGetFitnessValue());

	pc_log->vPrintLine(s_log, true);
	pc_log->vPrintEmptyLine(true);


	return b_updated;
}//bool CPopulationSizingOptimizer<TGenotype, TFenotype>::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
CError CPopulationSizingOptimizer<TGenotype, TFenotype>::e_configure_params_optimizer(istream *psSettings)
{
	CError c_error;
	
	pc_params_optimizer = OptimizerUtils::pcGetPopulationSizingSingleOptimizer(pc_problem, pc_log, i_random_seed, psSettings, &c_error);

	return c_error;
}//CError CPopulationSizingOptimizer<TGenotype, TFenotype>::e_configure_params_optimizer(istream *psSettings)

template <class TGenotype, class TFenotype>
bool CPopulationSizingOptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)
{
	CPopulationSizingSingleOptimizer<TGenotype, TFenotype> *pc_optimizer;
	CIndividual<TGenotype, TFenotype> *pc_optimizer_best_individual;

	bool b_updated = false;

	for (size_t i = 0; i < v_optimizers.size(); i++)
	{
		pc_optimizer = v_optimizers.at(i);
		pc_optimizer_best_individual = pc_optimizer->pcGetBestIndividual();

		if (pc_optimizer_best_individual)
		{
			if (b_update_best_individual(iIterationNumber, tStartTime, pc_optimizer_best_individual))
			{
				b_updated = true;
			}//if (b_update_best_individual(iIterationNumber, tStartTime, pc_ltga_best_individual))
		}//if (pc_optimizer_best_individual)
	}//for (size_t i = 0; i < v_optimizers.size(); i++)

	return b_updated;
}//bool CPopulationSizingOptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_clear_optimizers()
{
	for (size_t i = 0; i < v_optimizers.size(); i++)
	{
		delete v_optimizers.at(i);
	}//for (size_t i = 0; i < v_optimizers.size(); i++)

	v_optimizers.clear();
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_clear_optimizers()

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_clear_params()
{
	if (b_own_params)
	{
		delete pc_params_optimizer;
		pc_params_optimizer = nullptr;
	}//if (b_own_params)
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_clear_params()

template <class TGenotype, class TFenotype>
bool CPopulationSizingOptimizer<TGenotype, TFenotype>::b_run_optimizers_iteration(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_at_least_one_improvement = false;

	for (size_t i = 0; i < v_optimizers.size(); i++)
	{
		if (v_optimizers.at(i)->bRunIteration(iIterationNumber, tStartTime))
		{
			b_at_least_one_improvement = true;
		}//if (v_optimizers.at(i)->bRunIteration(iIterationNumber, tStartTime))
	}//for (size_t i = 0; i < v_optimizers.size(); i++)

	return b_at_least_one_improvement;
}//bool CPopulationSizingOptimizer<TGenotype, TFenotype>::b_run_optimizers_iteration(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_handle_optimizers(uint32_t iIterationNumber, time_t tStartTime)
{
	v_optimizers_deletion();
	v_optimizers_creation(iIterationNumber, tStartTime);
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_handle_optimizers(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_optimizers_deletion()
{
	unordered_set<size_t> s_optimizers_to_delete_indexes;
	s_optimizers_to_delete_indexes.reserve(v_optimizers.size());

	v_get_optimizers_to_delete(&s_optimizers_to_delete_indexes);

	if (!s_optimizers_to_delete_indexes.empty())
	{
		for (size_t i = v_optimizers.size(); i > 0; i--)
		{
			if (s_optimizers_to_delete_indexes.count(i - 1) > 0)
			{
				delete v_optimizers.at(i - 1);
				v_optimizers.erase(v_optimizers.begin() + i - 1);
			}//if (s_optimizers_to_delete_indexes.count(i - 1) > 0)
		}//for (size_t i = v_optimizers.size(); i > 0; i--)
	}//if (!s_optimizers_to_delete_indexes.empty())
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_optimizers_deletion()

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_optimizers_creation(uint32_t iIterationNumber, time_t tStartTime)
{
	if (v_optimizers.empty() || iIterationNumber - i_last_population_creation_iteration_number == i_population_creation_frequency)
	{
		v_add_new_optimizer(iIterationNumber, tStartTime);
	}//if (v_optimizers.empty() || iIterationNumber - i_last_population_creation_iteration_number == i_population_creation_frequency)
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_optimizers_creation(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_get_optimizers_to_delete(unordered_set<size_t> *psIndexes)
{
	double *pd_average_fitnesses = new double[v_optimizers.size()];

	for (size_t i = 0; i < v_optimizers.size(); i++)
	{
		*(pd_average_fitnesses + i) = v_optimizers.at(i)->dComputeAverageFitnessValue();

		if (v_optimizers.at(i)->bIsSteadyState())
		{
			psIndexes->insert(i);
		}//if (v_optimizers.at(i)->bIsSteadyState())
	}//for (size_t i = 0; i < v_optimizers.size(); i++)

	bool b_to_delete;

	for (size_t i = 0; i < v_optimizers.size(); i++)
	{
		if (psIndexes->count(i) == 0)
		{
			b_to_delete = false;

			for (size_t j = i + 1; j < v_optimizers.size() && !b_to_delete; j++)
			{
				b_to_delete = *(pd_average_fitnesses + i) < *(pd_average_fitnesses + j);
			}//for (size_t j = i + 1; j < v_optimizers.size() && !b_to_delete; j++)

			if (b_to_delete)
			{
				for (size_t j = 0; j <= i; j++)
				{
					psIndexes->insert(j);
				}//for (size_t j = 0; j <= i; j++)
			}//if (b_to_delete)
		}//if (psIndexes->count(i) == 0)
	}//for (size_t i = 0; i < v_optimizers.size(); i++)

	delete pd_average_fitnesses;
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_get_optimizers_to_delete(unordered_set<size_t> *psIndexes)

template <class TGenotype, class TFenotype>
void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_add_new_optimizer(uint32_t iIterationNumber, time_t tStartTime)
{
	CPopulationSizingSingleOptimizer<TGenotype, TFenotype> *pc_optimizer = (CPopulationSizingSingleOptimizer<TGenotype, TFenotype>*)pc_params_optimizer->pcCopy();

	pc_optimizer->vSetPopulationSize(i_next_population_size);
	pc_optimizer->vInitialize(tStartTime);

	v_optimizers.push_back(pc_optimizer);

	i_next_population_size *= i_population_size_multiplier;
	i_last_population_creation_iteration_number = iIterationNumber;
}//void CPopulationSizingOptimizer<TGenotype, TFenotype>::v_add_new_optimizer(uint32_t iIterationNumber, time_t tStartTime)



template class CPopulationSizingSingleOptimizer<CBinaryCoding, CBinaryCoding>;

template class CPopulationSizingOptimizer<CBinaryCoding, CBinaryCoding>;