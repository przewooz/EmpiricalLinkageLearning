#ifndef SYSTEM_H
#define SYSTEM_H

#define SYSTEM_ARGUMENT_DESTINATION_FILE "destination_file"
#define SYSTEM_ARGUMENT_MODE "mode"
#define SYSTEM_ARGUMENT_MODE_OPTIMIZATION "optimization"
#define SYSTEM_ARGUMENT_MODE_DYNAMIC_OPTIMIZATION "dynamic_optimization"
#define SYSTEM_ARGUMENT_MODE_GENERATE_RANDOM_MAX_SAT "generate_random_max_sat"
#define SYSTEM_ARGUMENT_MODE_PLANTED_MAX_SAT_2_PROBLEM_GENERATION "planted_max_sat_2_problem_generation"
#define SYSTEM_ARGUMENT_MODE_AVIS_SUBSET_SUM_KNAPSACK_PROBLEM_GENERATOR "avis_subset_sum_knapsack_problem_generator"
#define SYSTEM_ARGUMENT_MODE_LOCAL_SEARCH_REPEATABILITY_CHECKING "local_search_repeatability_checking"
#define SYSTEM_ARGUMENT_MODE_INJECTED_POPULATION_OPTIMIZER_REPEATABILITY_CHECKING "injected_population_optimizer_repeatability_checking"
#define SYSTEM_ARGUMENT_MODE_CHANGE_DETECTION "change_detection"

#define SYSTEM_ARGUMENT_RANDOM_SEED "random_seed"
#define SYSTEM_ARGUMENT_LOG_SYSTEM_OFFSET "log_system"


#define SYSTEM_ARGUMENT_NUMBER_OF_VARIABLES "number_of_variables"

#define SYSTEM_ARGUMENT_NUMBER_OF_CLAUSES "number_of_clauses"
#define SYSTEM_ARGUMENT_CLAUSE_LENGTH "clause_length"

#define SYSTEM_ARGUMENT_P "p"
#define SYSTEM_ARGUMENT_R "r"

#define SYSTEM_ARGUMENT_OUTPUT_FILE_PATH "output_file_path" 

#define SYSTEM_ARGUMENT_NUMBER_OF_REPEATS "number_of_repeats" 

#define SYSTEM_ARGUMENT_NUMBER_OF_INJECTED_GENOTYPES "number_of_injected_genotypes" 

#include "CommandExecution.h"
#include "Error.h"
#include "Log.h"
#include "StringUtils.h"
#include "util\tools.h"

#include <atlstr.h>
#include <cstdint>
#include <fstream>
#include <unordered_map>





class CSystem
{
static uint32_t iERROR_CODE_SYSTEM_IN_FILE_NOT_FOUND;
static uint32_t iERROR_CODE_SYSTEM_SETTNGS_FILE_NOT_FOUND;

public:
	CSystem(int iArgCount, char **ppcArgValues, CError *pcError);
	
	~CSystem();

	CError eRun();

private:
	void v_register_command(CCommandExecution *pcCommand) { m_commands.insert(pair<const CString, CCommandExecution*>(pcCommand->sGetName(), pcCommand)); };
	void v_register_commands();

	CError e_create_log_file(CString  sFileName, CString  *psFinalFilename);
	CError e_read_arguments_from_file(int iArgCount, char **ppcArgValues);
	CError e_load_settings(CString sSettingsFile);

	ofstream *pf_log;

	CString  s_result_summary_file;
	ofstream f_result_summary;

	CString s_log_file;
	CLog c_log;

	uint32_t i_random_seed;

	ifstream f_settings;

	unordered_map<CString, CString> m_arguments;
	unordered_map<CString, CCommandExecution*> m_commands;
};//class CSystem

#endif//SYSTEM_H