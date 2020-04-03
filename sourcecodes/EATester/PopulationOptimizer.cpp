#include "PopulationOptimizer.h"

#include "BinaryCoding.h"
#include "CommandParam.h"
#include "GenerationUtils.h"
#include "UIntCommandParam.h"

template <class TGenotype, class TFenotype>
CPopulationOptimizer<TGenotype, TFenotype>::CPopulationOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<TGenotype, TFenotype>(pcProblem, pcLog, iRandomSeed)
{
	pc_generation = nullptr;
	ppc_population = nullptr;

	b_own_injected_genotypes = false;
	ppc_injected_genotypes = nullptr;
}//CPopulationOptimizer<TGenotype, TFenotype>::CPopulationOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

template <class TGenotype, class TFenotype>
CPopulationOptimizer<TGenotype, TFenotype>::CPopulationOptimizer(CPopulationOptimizer<TGenotype, TFenotype>* pcOther)
	: COptimizer<TGenotype, TFenotype>(pcOther)
{
	i_population_size = pcOther->i_population_size;
	pc_generation = pcOther->pc_generation;
	ppc_population = nullptr;

	b_own_injected_genotypes = false;
	i_number_of_injected_genotypes = pcOther->i_number_of_injected_genotypes;
	ppc_injected_genotypes = pcOther->ppc_injected_genotypes;
}//CPopulationOptimizer<TGenotype, TFenotype>::CPopulationOptimizer(CPopulationOptimizer<TGenotype, TFenotype>* pcOther)

template <class TGenotype, class TFenotype>
CPopulationOptimizer<TGenotype, TFenotype>::~CPopulationOptimizer()
{
	v_replace_population(nullptr);
	v_clear_params();
	vInjectGenotypes(0, nullptr);
}//CPopulationOptimizer<TGenotype, TFenotype>::~CPopulationOptimizer()

template <class TGenotype, class TFenotype>
CError CPopulationOptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	v_clear_params();

	CError c_error = COptimizer<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CBoolCommandParam p_injected_genotypes(POPULATION_OPTIMIZER_ARGUMENT_INJECTED_GENOTYPES);
		bool b_injected_genotypes = p_injected_genotypes.bGetValue(psSettings, &c_error);

		if (!c_error)
		{
			if (b_injected_genotypes)
			{
				if (ppc_injected_genotypes != nullptr)
				{
					if (i_number_of_injected_genotypes >= 2)
					{
						pc_generation = GenerationUtils::pcGetEmptyGeneration(pc_problem);

						TGenotype *pc_sample_genotype = pc_generation->pcGenerateEmpty();

						for (uint32_t i = 0; i < i_number_of_injected_genotypes && !c_error; i++)
						{
							if (!pc_sample_genotype->bHasSameParams(*(ppc_injected_genotypes + i)))
							{
								c_error.vSetError(CError::iERROR_CODE_SYSTEM_ARGUMENT_WRONG_VALUE, "injected genotypes are different than expected by provided evaluation");
							}//if (!pc_sample_genotype->bHasSameParams(*(ppc_injected_genotypes + i)))
						}//for (uint32_t i = 0; i < i_number_of_injected_genotypes && !c_error; i++)

						delete pc_sample_genotype;
					}//if (i_number_of_injected_genotypes >= 2)
					else
					{
						c_error.vSetError(CError::iERROR_CODE_SYSTEM_OUT_OF_BOUND_ARGUMENT, "the number of injected genotypes must be greater or equal to 2");
					}//else if (i_number_of_injected_genotypes >= 2)
				}//if (ppc_injected_genotypes != nullptr)
				else
				{
					c_error.vSetError(CError::iERROR_CODE_OPERATOR_NOT_FOUND, "genotypes are not injected");
				}//else if (ppc_injected_genotypes != nullptr)
			}//if (b_injected_genotypes)
			else
			{
				if (ppc_injected_genotypes != nullptr)
				{
					c_error.vSetError(CError::iERROR_CODE_SYSTEM_ARGUMENT_WRONG_VALUE, "genotypes cannot be injected - the corresponding argument is set to 0");
				}//if (ppc_injected_genotypes != nullptr)

				if (!c_error)
				{
					CUIntCommandParam p_population_size(POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE, 2, UINT32_MAX);
					i_population_size = p_population_size.iGetValue(psSettings, &c_error);
				}//if (!c_error)

				if (!c_error)
				{
					pc_generation = GenerationUtils::pcGetGeneration(pc_problem, psSettings, &c_error);
				}//if (!c_error)
			}//else if (b_injected_genotypes)
		}//if (!c_error)
	}//if (!c_error)

	return c_error;
}//CError CPopulationOptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)

template <class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)
{
	delete pc_best_individual;
	pc_best_individual = nullptr;

	v_generation();
	v_evaluation();

	b_update_best_individual(0, tStartTime);
}//void CPopulationOptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)

template <class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::vInjectGenotypes(uint32_t iSize, TGenotype **ppcGenotypes, bool bOwn = false)
{
	if (b_own_injected_genotypes)
	{
		if (ppc_injected_genotypes != nullptr)
		{
			for (uint32_t i = 0; i < i_number_of_injected_genotypes; i++)
			{
				delete *(ppc_injected_genotypes + i);
			}//for (uint32_t i = 0; i < i_number_of_injected_genotypes; i++)
		}//if (ppc_injected_genotypes != nullptr)

		delete ppc_injected_genotypes;
	}//if (b_own_injected_genotypes)

	v_replace_population(nullptr);

	i_population_size = iSize;

	b_own_injected_genotypes = bOwn;
	i_number_of_injected_genotypes = iSize;
	ppc_injected_genotypes = ppcGenotypes;
}//void CPopulationOptimizer<TGenotype, TFenotype>::vInjectGenotypes(uint32_t iSize, TGenotype **ppcGenotypes, bool bOwn = false)

template <class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::v_generation()
{
	CIndividual<TGenotype, TFenotype> **ppc_generated_population = new CIndividual<TGenotype, TFenotype>*[i_population_size];

	for (uint32_t i = 0; i < i_population_size; i++)
	{
		TGenotype *pc_genotype = ppc_injected_genotypes == nullptr ? pc_generation->pcGenerate() : new TGenotype(*(ppc_injected_genotypes + i));
		(*(ppc_generated_population + i)) = pc_create_individual(pc_genotype);
	}//for (uint32_t i = 0; i < i_population_size; i++)

	v_replace_population(ppc_generated_population);
}//void CPopulationOptimizer<TGenotype, TFenotype>::v_generation()

template <class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::v_evaluation()
{
	for (uint32_t i = 0; i < i_population_size; i++)
	{
		(*(ppc_population + i))->vEvaluate();
	}//for (uint32_t i = 0; i < i_population_size; i++)
}//void CPopulationOptimizer<TGenotype, TFenotype>::v_evaluation()

template <class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::v_replace_population(CIndividual<TGenotype, TFenotype> **ppcNewPopulation)
{
	if (ppc_population)
	{
		for (uint32_t i = 0; i < i_population_size; i++)
		{
			delete *(ppc_population + i);
		}//for (uint32_t i = 0; i < i_population_size; i++)
	}//if (ppc_population)

	delete ppc_population;

	ppc_population = ppcNewPopulation;
}//void CPopulationOptimizer<TGenotype, TFenotype>::v_replace_population(CIndividual<TGenotype, TFenotype> **ppcNewPopulation)

template <class TGenotype, class TFenotype>
bool CPopulationOptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_updated = false;

	CIndividual<TGenotype, TFenotype> *pc_individual;

	for (uint32_t i = 0; i < i_population_size; i++)
	{
		pc_individual = *(ppc_population + i);
		b_updated = b_update_best_individual(iIterationNumber, tStartTime, pc_individual) || b_updated;
	}//for (uint32_t i = 0; i < i_population_size; i++)

	return b_updated;
}//bool CPopulationOptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)

template <class TGenotype, class TFenotype>
TGenotype ** CPopulationOptimizer<TGenotype, TFenotype>::pcGetPopulationGenotypes()
{
	TGenotype **ppc_population_genotypes = new TGenotype*[i_population_size];

	for (uint32_t i = 0; i < i_population_size; i++)
	{
		*(ppc_population_genotypes + i) = (*(ppc_population+ i))->pcGetGenotype();
	}//for (uint32_t i = 0; i < i_population_size; i++)

	return ppc_population_genotypes;
}//TGenotype ** CPopulationOptimizer<TGenotype, TFenotype>::pcGetPopulationGenotypes()

template <class TGenotype, class TFenotype>
TFenotype ** CPopulationOptimizer<TGenotype, TFenotype>::pcGetPopulationFenotypes()
{
	TFenotype **ppc_population_fenotypes = new TFenotype*[i_population_size];

	for (uint32_t i = 0; i < i_population_size; i++)
	{
		*(ppc_population_fenotypes + i) = (*(ppc_population + i))->pcGetFenotype();
	}//for (uint32_t i = 0; i < i_population_size; i++)

	return ppc_population_fenotypes;
}//TFenotype ** CPopulationOptimizer<TGenotype, TFenotype>::pcGetPopulationFenotypes()

template<class TGenotype, class TFenotype>
void CPopulationOptimizer<TGenotype, TFenotype>::v_clear_params()
{
	if (b_own_params)
	{
		delete pc_generation;
		pc_generation = nullptr;
	}//if (b_own_params)
}//void CPopulationOptimizer<TGenotype, TFenotype>::v_clear_params()

template class CPopulationOptimizer<CBinaryCoding, CBinaryCoding>;