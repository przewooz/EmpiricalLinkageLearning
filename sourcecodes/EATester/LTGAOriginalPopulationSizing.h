#ifndef LTGA_ORIGINAL_POPULATION_SIZING_H
#define LTGA_ORIGINAL_POPULATION_SIZING_H

#define LTGA_ORIGINAL_POPULATION_SIZING_ARGUMENT_POPULATION_CREATION_FREQUENCY "population_creation_frequency"
#define LTGA_ORIGINAL_POPULATION_SIZING_ARGUMENT_POPULATION_SIZE_MULTIPLIER "population_size_multiplier"

#include "BinaryCoding.h"
#include "Error.h"
#include "Log.h"
#include "LTGAOriginal.h"
#include "Optimizer.h"
#include "Problem.h"

#include <ctime>
#include <cstdint>
#include <istream>
#include <unordered_set>
#include <vector>

using namespace std;

class CLTGAOriginalPopulationSizing : public COptimizer<CBinaryCoding, CBinaryCoding>
{
public:
	CLTGAOriginalPopulationSizing(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CLTGAOriginalPopulationSizing(CLTGAOriginalPopulationSizing *pcOther);

	virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CLTGAOriginalPopulationSizing(this); };

	virtual ~CLTGAOriginalPopulationSizing();

	virtual CError eConfigure(istream *psSettings);

	virtual void vInitialize(time_t tStartTime);
	virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

private:
	using COptimizer<CBinaryCoding, CBinaryCoding>::b_update_best_individual;
	virtual bool b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime);

	void v_clear_ltgas();

	bool b_run_ltgas_iteration(uint32_t iIterationNumber, time_t tStartTime);
	
	void v_handle_ltgas(uint32_t iIterationNumber, time_t tStartTime);
	void v_ltgas_deletion();
	void v_ltgas_creation(uint32_t iIterationNumber, time_t tStartTime);
	void v_get_ltgas_to_delete(unordered_set<size_t> *psIndexes);
	void v_add_new_ltga(uint32_t iIterationNumber, time_t tStartTime);

	vector<CLTGAOriginal*> v_ltgas;

	uint32_t i_population_creation_frequency;
	uint32_t i_population_size_multiplier;

	uint32_t i_next_population_size;
	uint32_t i_last_population_creation_iteration_number;

	CLTGAOriginal c_params_ltga;
};//class CLTGAOriginalPopulationSizing : public COptimizer<CBinaryCoding, CBinaryCoding>

#endif//LTGA_ORIGINAL_POPULATION_SIZING_H