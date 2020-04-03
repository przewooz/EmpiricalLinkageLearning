#include "P3BasedOptimizer.h"

#include "StringCommandParam.h"

#include "../P3/Configuration.h"

#include <tuple>
#include <vector>

CP3BasedOptimizer::CP3BasedOptimizer(string sOptimizerName, CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcProblem, pcLog, iRandomSeed), s_optimizer_name(sOptimizerName)
{
	pc_recorder = nullptr;

	pc_config = new Configuration();
	pc_rand = new Random(iRandomSeed);
}//CP3BasedOptimizer::CP3BasedOptimizer(string sOptimizerName, CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

CP3BasedOptimizer::CP3BasedOptimizer(CP3BasedOptimizer *pcOther)
	: COptimizer<CBinaryCoding, CBinaryCoding>(pcOther), s_optimizer_name(pcOther->s_optimizer_name)
{
	pc_recorder = nullptr;

	pc_config = pcOther->pc_config;
	pc_rand = pcOther->pc_rand;
}//CP3BasedOptimizer::CP3BasedOptimizer(CP3BasedOptimizer *pcOther)

CP3BasedOptimizer::~CP3BasedOptimizer()
{
	delete pc_recorder;

	if (b_own_params)
	{
		delete pc_config;
		delete pc_rand;
	}//if (b_own_params)
}//CP3BasedOptimizer::~CP3BasedOptimizer()

CError CP3BasedOptimizer::eConfigure(istream *psSettings)
{
	CError c_error = COptimizer::eConfigure(psSettings);

	if (!c_error)
	{
		CFilePathCommandParam p_config_file_path(OPTIMIZER_ARGUMENT_CONFIG_FILE_PATH);
		CString s_config_path = p_config_file_path.sGetValue(psSettings, &c_error);

		if (!c_error)
		{
			uint32_t i_config_file_path_length = s_config_path.GetLength() + 1;
			char *pc_config_file_path = new char[i_config_file_path_length];
			strcpy_s(pc_config_file_path, i_config_file_path_length, s_config_path);

			pc_config->parse(pc_config_file_path);
			pc_config->set_injected_evaluation(pc_problem->pcGetEvaluation());
			pc_config->set("optimizer", s_optimizer_name);
			pc_config->set("problem", string("Injected"));
			pc_config->set("verbosity", 0);

			delete pc_config_file_path;
		}//if (!c_error)
	}//if (!c_error)

	return c_error;
}//CError CP3BasedOptimizer::eConfigure(istream *psSettings)

void CP3BasedOptimizer::vInitialize(time_t tStartTime)
{
	COptimizer<CBinaryCoding, CBinaryCoding>::vInitialize(tStartTime);

	evaluation::pointer c_evaluation_ptr = pc_config->get<evaluation::pointer>("problem");
	optimize::pointer c_optimizer_ptr = pc_config->get<optimize::pointer>("optimizer");

	delete pc_recorder;
	pc_recorder = new Middle_Layer(*pc_config, c_evaluation_ptr(*pc_config, 0), tStartTime, pc_log);

	pc_optimizer = c_optimizer_ptr(*pc_rand, *pc_recorder, *pc_config);
}//void CP3BasedOptimizer::vInitialize(time_t tStartTime)

bool CP3BasedOptimizer::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	pc_optimizer->iterate();
	
	return b_update_best_individual(iIterationNumber, tStartTime);
}//bool CP3BasedOptimizer::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)

bool CP3BasedOptimizer::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)
{
	tuple<float, float, int, time_t, vector<bool>> t_summary = pc_recorder->results.best();

	double d_best_individual_candidate_fitness = (double)get<0>(t_summary);

	bool b_updated = COptimizer<CBinaryCoding, CBinaryCoding>::b_update_best_individual(iIterationNumber, tStartTime, d_best_individual_candidate_fitness, [&](CBinaryCoding *pcBestGenotype)
	{
		vector<bool> v_best_bits = get<4>(t_summary);

		for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
		{
			*(pcBestGenotype->piGetBits() + i) = (int32_t)v_best_bits.at(i);
		}//for (uint16_t i = 0; i < pc_problem->pcGetEvaluation()->iGetNumberOfElements(); i++)
	});

	return b_updated;
}//bool CP3BasedOptimizer::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime)