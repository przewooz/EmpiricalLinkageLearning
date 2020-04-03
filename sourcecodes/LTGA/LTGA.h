/*
* LTGA.c
*
* Copyright (c) Peter A.N. Bosman
*
* The software in this file is the proprietary information of
* Peter A.N. Bosman.
*
* IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
* DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
* INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
* TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
* OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
* REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
* AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
* USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
*
* Linkage Tree Genetic Algorithm
*
* In this implementation, maximization is assumed.
*
* The software in this file is the result of (ongoing) scientific research.
* The following people have been actively involved in this research over
* the years:
* - Peter A.N. Bosman
* - Dirk Thierens
*
* The most recent publication that this code was used in, is:
*
* P.A.N. Bosman and D. Thierens. More Concise and Robust Linkage Learning by
* Filtering and Combining Linkage Hierarchies. In C. Blum and E. Alba, editors,
* Proceedings of the Genetic and Evolutionary Computation Conference -
* GECCO-2013, pages 359-366, ACM Press, New York, New York, 2013.
*/

#ifndef _LTGA_H_
#define _LTGA_H_

#define _CRT_SECURE_NO_WARNINGS

#define FALSE 0
#define TRUE 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "../EATester/BinaryCoding.h"
#include "../EATester/Evaluation.h"
#include "../EATester/Log.h"

namespace NLTGA
{
	class LTGA
	{
	public:
		LTGA(CEvaluation<CBinaryCoding> *evaluation, CLog *log, int64_t random_seed);

		void initialize();
		void ezilaitiniMemory(void);

		void runGeneration(char *best_solution);

		void updateBestPrevGenSolution();

		bool areAllSolutionsTheSame();
		double computeAverageObjectiveValue();

		double getBestEverEvaluatedObjectiveValue() { return best_ever_evaluated_objective_value; }
		char *getBestEverEvaluatedSolution() { return best_ever_evaluated_solution; }

		int getPopulationSize() { return population_size; }
		void setPopulationSize(int population_size) { this->population_size = population_size; }
		
		bool getDoLocalSearch() { return do_local_search; }
		void setDoLocalSearch(bool do_local_search) { this->do_local_search = do_local_search; }

		bool getLocalSearchOneIteration() { return local_search_one_iteration; }
		void setLocalSearchOneIteration(bool local_search_one_iteration) { this->local_search_one_iteration = local_search_one_iteration; }

		bool getWithoutTournamentSelection() { return without_tournament_selection; }
		void setWithoutTournamentSelection(bool without_torunament_selection) { this->without_tournament_selection = without_torunament_selection; }

		bool getLinkageTreeRandomOrder() { return linkage_tree_random_order; }
		void setLinkageTreeRandomOrder(bool linkage_tree_random_order) { this->linkage_tree_random_order = linkage_tree_random_order; }

		//int main(int argc, char **argv);
		//void run();

	private:
		void *Malloc(long size);
		void interpretCommandLine(int argc, char **argv);
		void parseCommandLine(int argc, char **argv);
		void parseOptions(int argc, char **argv, int *index);
		void printAllInstalledProblems(void);
		void optionError(char **argv, int index);
		void parseParameters(int argc, char **argv, int *index);
		void printUsage(void);
		void checkOptions(void);
		void printVerboseOverview(void);
		double randomRealUniform01(void);
		int randomInt(int maximum);
		int *randomPermutation(int n);
		void randomPermutation(int *permutation, int n);
		char *installedProblemName(int index);
		int numberOfInstalledProblems(void);
		void installedProblemEvaluation(int index, char *parameters, double *objective_value, double *constraint_value);
		void onemaxFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value);
		void deceptiveTrap4TightEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value);
		void deceptiveTrap4LooseEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value);
		void deceptiveTrap5TightEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value);
		void deceptiveTrap5LooseEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value);
		void deceptiveTrapKTightEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value, int k);
		void deceptiveTrapKLooseEncodingFunctionProblemEvaluation(char *parameters, double *objective_value, double *constraint_value, int k);
		void initializeMemory();
		//void initializeRandomNumberGenerator();
		void initializePopulationAndFitnessValues();
		void writeGenerationalStatistics();
		void writeGenerationalSolutions(char is_final_generation);
		void writeGenerationalSolutionsBest(char is_final_generation);
		void writeGenerationalSolutionsBestEver();
		//void writeRunningTime(char *filename);
		char checkTerminationCondition();
		char checkNumberOfEvaluationsTerminationCondition();
		char checkVTRTerminationCondition();
		void determineBestSolutionInCurrentPopulation(int *index_of_best);
		char checkFitnessVarianceTermination();
		void selectFinalSurvivors();
		char betterFitness(double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y);
		char equalFitness(double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y);
		void makeOffspring(char *best_solution);
		void selectForLearningLT();
		void learnLT();
		double *estimateParametersForSingleBinaryMarginal(int *indices, int number_of_indices, int *factor_size);
		int determineNearestNeighbour(int index, double **S_matrix, int mpm_length);
		void printLT();
		double log2(double x);
		void generateAndEvaluateNewSolutionsToFillOffspring(char *best_solution);
		char *generateNewSolution(int which, double *obj, double *con, char *best_solution);
		//long getMilliSecondsRunning();
		void localSearch(char *parameters, double *objective_value, double *constraint_value);

		char	**population,                           /* Population solutions. */
				**selection,                            /* Selected solutions from which to learn the LT. */
				**offspring,                            /* Offspring solutions. */
				 *best_prevgen_solution,                /* The best solution found in all previous generations. */
				 *best_ever_evaluated_solution;         /* The best ever evaluated solution. */
		short     write_generational_statistics,        /* Whether to compute and write statistics every generation (0 = no). */
				  write_generational_solutions,         /* Whether to write the population every generation (0 = no). */
				  print_verbose_overview,               /* Whether to print a overview of settings (0 = no). */
				  use_vtr,                              /* Whether to terminate at the value-to-reach (VTR) (0 = no). */
				  vtr_hit_has_happened,                 /* Whether the vtr has been hit yet. */
				  print_lt_contents;                    /* Whether to print the LT structure every generation. */
		int       problem_index,                        /* The index of the optimization problem. */
				  number_of_parameters,                 /* The number of parameters to be optimized. */
				  number_of_generations,                /* The current generation count. */
				  population_size,                      /* The size of the population. */
				  selection_size,                       /* The size of the selection. */
			      offspring_size,                       /* The size of the offspring. */
			      maximum_number_of_evaluations,        /* The maximum number of evaluations. */
			      no_improvement_stretch,               /* The number of subsequent generations without an improvement. */
			    **mpm,                                  /* The marginal product model (mpm). */
				 *mpm_number_of_indices,                /* The number of variables in each factor in the mpm. */
				  mpm_length,                           /* The number of factors in the mpm. */
				**lt,                                   /* The linkage tree (lt). */
				 *lt_number_of_indices,                 /* The number of variables in each factor in the lt. */
				  lt_length;                            /* The number of factors in the lt. */
		long      number_of_evaluations,                /* The current number of times a function evaluation was performed. */
				  timestamp_start;                      /* The time stamp in milliseconds for when the algorithm was started. */
		double   *objective_values,                     /* Objective values for population members. */
				 *constraint_values,                    /* Sum of all constraint violations for population members. */
				 *objective_values_offspring,           /* Objective values of offspring solutions. */
				 *constraint_values_offspring,          /* Sum of all constraint violations of offspring solutions. */
				  best_prevgen_objective_value,         /* Objective value of best solution in all previous generations. */
				  best_prevgen_constraint_value,        /* Constraint value of best solution in all previous generations. */
				  fitness_variance_tolerance,           /* The minimum fitness variance level that is allowed. */
				  vtr,                                  /* The value-to-reach (function value of best solution that is feasible). */
			    **MI_matrix,                            /* Mutual information between any two variables. */
			      best_ever_evaluated_objective_value,  /* The best ever evaluated objective value. */
			      best_ever_evaluated_constraint_value; /* The best ever evaluated constraint value. */
		int64_t   random_seed,                          /* The seed used for the random-number generator. */
				  random_seed_changing;                 /* Internally used variable for randomly setting a random seed. */

		bool	  do_local_search;
		bool	  local_search_one_iteration;
		bool	  without_tournament_selection;
		bool      linkage_tree_random_order;

		CEvaluation<CBinaryCoding> *evaluation;
		CLog *log;

		CBinaryCoding evaluation_individual;
	};
}

#endif