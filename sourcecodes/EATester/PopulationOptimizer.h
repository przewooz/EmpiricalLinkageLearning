#ifndef POPULATION_OPTIMIZER_H
#define POPULATION_OPTIMIZER_H

#define POPULATION_OPTIMIZER_ARGUMENT_INJECTED_GENOTYPES "injected_genotypes"

#define POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE "population_size"

#define POPULATION_OPTIMIZER_ARGUMENT_MAX_NUMBER_OF_GENE_PATTERNS "max_number_of_gene_patterns"

#include "Error.h"
#include "Generation.h"
#include "Individual.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"

#include <cstdint>
#include <ctime>
#include <istream>

using namespace std;

namespace ChangeDetector
{
	class CDeterministicOptimizerChangeDetector;
}//namespace ChangeDetector

template <class TGenotype, class TFenotype>
class CPopulationOptimizer : public COptimizer<TGenotype, TFenotype>
{
friend class ChangeDetector::CDeterministicOptimizerChangeDetector;

public:
	CPopulationOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CPopulationOptimizer(CPopulationOptimizer<TGenotype, TFenotype> *pcOther);

	virtual ~CPopulationOptimizer();

	virtual CError eConfigure(istream *psSettings);

	virtual void vInitialize(time_t tStartTime);

	void vInjectGenotypes(uint32_t iSize, TGenotype **ppcGenotypes, bool bOwn = false);

	CIndividual<TGenotype, TFenotype> **ppcGetPopulation() { return ppc_population; };
	uint32_t iGetPopulationSize() { return i_population_size; };

protected:
	virtual void v_generation();
	virtual void v_evaluation();

	virtual void v_replace_population(CIndividual<TGenotype, TFenotype> **ppcNewPopulation);

	using COptimizer<TGenotype, TFenotype>::b_update_best_individual;
	virtual bool b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime);

	TGenotype **pcGetPopulationGenotypes();
	TFenotype **pcGetPopulationFenotypes();

	void v_clear_params();

	uint32_t i_population_size;

	CGeneration<TGenotype> *pc_generation;

	CIndividual<TGenotype, TFenotype> **ppc_population;

private:
	bool b_own_injected_genotypes;
	uint32_t i_number_of_injected_genotypes;
	TGenotype **ppc_injected_genotypes;
};//class CPopulationOptimizer

#endif//POPULATION_OPTIMIZER_H