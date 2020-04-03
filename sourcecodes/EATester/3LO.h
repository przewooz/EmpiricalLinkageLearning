#ifndef C3LO_OPTIMIZER_H
#define C3LO_OPTIMIZER_H


#include "BinaryCoding.h"
#include "BinaryOptimizer.h"
#include "CommandParam.h"
#include "Error.h"
#include "Log.h"
#include "MathUtils.h"
#include "Optimizer.h"
#include "Problem.h"
#include  "util\timer.h"
#include  "util\tools.h"


#include <istream>
#include <algorithm>



namespace ThreeLO
{
	#define  _3LO_ARGUMENT_USE_DSM									"use_dsm" 

	#define  _3LO_TRIBE_MAX_ROUNDS_NUM								5
	#define  _3LO_HIGH_LVL_GENE_VALUE_GENERATION_SUBGENE_MAX		10


	class  C3LOIndividual;
	class  C3LOPattern;
	class  CMessyGene;
	class  C3LOPatternPair;
	class  C3LOHighLvlGeneFreq;
	class  CGenotypeInfoNode;


	class C3LO : public CBinaryOptimizer
	{
		friend class C3LOSingle;
		friend class C3LOIndividual;
	public:
		static uint32_t iERROR_PARENT_C3LOOptimizer;
		static uint32_t iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0;



		C3LO(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		C3LO(C3LO *pcOther);
		~C3LO();

		virtual CError eConfigure(istream *psSettings);

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new C3LO(this); };

		virtual void vInitialize(time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);


		//for individuals
		double dComputeFitness(int32_t *piBits);

		CString  sAdditionalSummaryInfo();

	private:
		void  v_execute_super_pop(uint32_t iIterationNumber, time_t tStartTime);
		void  v_add_new_3lo_pop();

		bool  b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int iPatternLevel);
		void  v_report_scraps_crossings_trees(C3LOSingle *pc3LOSingle, CString  *psScraps, CString  *psCrossings, CString  *psTrees);
		void  v_report_high_level_genes(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel);
		void  v_report_high_level_genes_numbers(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel);
		void  v_report_groupped_trees(CString  *psGrouppedTrees);

		void  v_join_children_trees_if_necessary(vector<C3LOPattern  *> *pvTreesAtLevel, int iPatternLevel);
		void  v_join_children_refresh_grouped_trees(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees);
		void  v_join_children_assign_to_tree_groups(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees);


		int  i_templ_length;

		vector <C3LOSingle *>  v_3lo_pops;
		vector<vector<C3LOPattern *>>  v_grouped_children_trees;

		TimeCounters::CTimeCounter  c_time_counter;
		time_t t_start;
		int32_t  *pi_best_genotype;//for COptimizer best individual update


		double  **pd_dsm_matrix;
		int  **pi_gene_realtions;
		int  i_relation_length_max;

		int  i_dsm_reporting_enumeration_tool;
		C3LOSingle  *pc_3lo_single_tool;
		C3LOSingle  *pc_3lo_additional_pop;

		bool b_use_dsm;
	};//class C3LO : public CBinaryOptimizer


	class C3LOSingle : public CBinaryOptimizer
	{
		friend class C3LO;
		friend class C3LOIndividual;				
	public:
		static uint32_t iERROR_PARENT_C3LOSingleOptimizer;
		static uint32_t iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0;



		C3LOSingle(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		C3LOSingle(C3LOSingle *pcOther);
		~C3LOSingle();

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new C3LOSingle(this); };

		virtual void vInitialize(time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);
		virtual bool bRunIteration_3lo_berlin(uint32_t iIterationNumber, time_t tStartTime);//just for history

		void  vSortAllSubpopsByFitness();

		vector<vector<C3LOPattern *>>  *pvGetHigherLevelTreesGenes() { return(&v_higher_level_trees_genes); }

		//for individuals
		double dComputeFitness(int32_t *piBits);
		bool  bCanIndBeAddedAtAnyLevel(C3LOIndividual *pcInd, bool  bComparePhenotypes = true, bool  bCompareGenotypes = false);

		void  vShufflePatternsToCrossSingleLevel(int iLevel);
		void  vShufflePatternsToCross();
		void  vShufflePatternsAndFreqs();
		bool  bLinkageGenerated() { return(b_linkage_generated); };

		

	private:
		void  v_initialize_orders();
		void  v_add_order_by_genotype(bool  bReversed);
		void  v_add_order_random(vector <int>  *pvNewOrder, vector <int>  *pvOppositeNewOrder, bool  bWithOpposite);

		void  v_insert_gene_pattern_part_level(int  iPatternPartLevel);
		void  v_add_gene_pattern_tree(C3LOPattern  *pcTree, int  iPatternLevel);
		double  d_get_most_similar_pattern_pair_offset(int  *piMostSimlarPairOffset, vector<C3LOPatternPair>  *pvPatternPairs);

		void  v_remove_pattern_from_buffer(vector<C3LOPattern  *>  *pvPartsBuffer, C3LOPattern  *pcPatternToRemove);
		void  v_remove_pattern_from_pairs_buffer(vector<C3LOPatternPair>  *pvPatternPairsBuffer, C3LOPattern  *pcPatternToRemove);


		bool  b_add_higher_level_gene(C3LOPattern  *pcHigherLevelGene, int  iGeneLevel);
		//bool  b_ind_on_lvl(C3LOIndividual *pcInd, int  iEvolutionLvl);
		bool  b_run_flat_pop_for_ind(C3LOIndividual  *pcIndMain, int iPatternLevel, vector<C3LOIndividual  *>  *pvLevelsOffspringsDifferentToParent);
		void  v_tiny_restricted_evolution(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentDonator, int  iPatternLevel, int  iPopDestLvl, int iRepeatation = -1);

		void  v_process_berlin(C3LOIndividual  *pcParentMain, int  iPopLevel);
		bool  b_process_berlin(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int  *piGenotypeToShuffle, int  *piShuffle, int iPatternLevel);
		void  v_process(C3LOIndividual  *pcParentMain, int  iPopLevel);
		bool  b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int  *piGenotypeToShuffle, int  *piShuffle, int iPatternLevel);


		int i_get_parrent_offset_respect_inner_level(vector<C3LOIndividual *>  *pvIndToChoose);
		int  i_get_next_gene_offset(int iBaseGeneOffset, int  **piGeneRelations, int  *piIncludedGenesMask);
		int  i_get_pattern_in_touch_if_possible(int  *piGeneGroup, vector<C3LOPattern *>  *pvPatternPool);
		int  i_get_pattern_in_touch_in_differences(int  *piGeneGroup, int  *piDifferenceMask, vector<C3LOPattern *>  *pvPatternPool);
		void  v_get_all_original_pattern_parts_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel);
		void  v_get_all_intercepting_patterns_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel);
		C3LOPattern*  pc_get_any_original_pattern_part_extending_gene_group(int  *piGeneGroup, int iPatternLevel);
		bool  b_shuffle_the_genotype(int  *piGenotypeToShuffle, C3LOIndividual  *pcCandidate, int iPatternLevel);
		void  v_add_ind_to_process_in_next_generation(C3LOIndividual  *pcIndividualToProcessInNextGeneration);
		
		void  v_update_high_lvl_genes_and_trees();
		void  v_create_trees_from_ind(C3LOIndividual  *pcIndiv, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel);
		void  v_create_trees_from_ind_2(C3LOIndividual  *pcIndSource, C3LOIndividual  *pcIndDest, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel);
		void  v_create_trees_from_dsm(vector<C3LOPattern *>  *pvTreesDest);
		void  v_build_create_intercept_patterns(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies);
		void  v_add_intercept_pattern(int  *piPatternDefinition, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies);
		void  v_build_pattern_tree_dsm_like(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel);
		void  v_build_trees_from_dsm(double  **pdDSM, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel);
		void  v_build_dsm(double  **pdDSM, vector<C3LOPattern  *>  *pvGenePatternsParts);
		void  v_build_gene_relations_from_dsm(double  **pdDSM, int  **piGeneRelations);

		void  v_build_pattern_tree(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel);
		bool  b_cross_individuals_by_tree_branch(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOPattern *pcCrossingTree, int  iPatternLevel);
		bool  b_cross_individuals_by_tree_branch_or_part(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, int  iPatternLevel);
		int   i_cross_individuals_by_crossing_tree(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, bool bMainGroup, int  iPatternLevel);
		int   i_cross_individuals_by_crossing_tree_berlin(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, bool bMainGroup, int  iPatternLevel);
		bool  b_cross_individuals_by_tree_branch_longest_reasonable_branch(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOIndividual  **pcOffspring, C3LOIndividual  **pcReversedOffspring, C3LOPattern *pcCrossingTree, int  iPatternLevel);
		bool  b_cross_individuals_by_tree(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOPattern *pcCrossingTree, int iPatternLevel);

		int  i_cross_dsm_glueing_similarities_from_different_start_genes(int  **piGeneRelations, int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel);
		int  i_cross_dsm(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel);
		int  i_cross_intercept_scraps(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel);
		int   i_cross_fluent_scraps(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel);
		int   i_cross_fluent_scraps_berlin(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel);
		bool  b_cross_infect_individual_by_tree(C3LOIndividual  *pcInd, C3LOIndividual  *pcDonator, C3LOPattern *pcCrossingTree, int iPatternLevel);
		void  v_cross_infect_individual_by_tree(int  *piGenotype, int  *piGenotypeDonator, C3LOPattern *pcCrossingTree);

		int  i_check_number_of_building_block_occurences(int  *piGenotype0, C3LOPattern *pcCrossingTree, int iPatternLevel);


		void  v_remove_from_levels(C3LOIndividual *pcInd);
		double  d_get_covering_at_level(int  iLevel);
		void  v_add_linkage_covering_at_level(int  iLevel);

		int  i_get_ind_number();

		C3LOPattern*  pc_get_tree_root_for_pattern(C3LOPattern  *pcCrossingTree, int  iPatternLevel);



		void  v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd);
		void  v_update_higher_genes_freq_by_ind_this_pop(C3LOIndividual  *pcInd);
		void  v_count_higher_genes_freq();

		void  v_filter_gene_freq_by_ind(int iBestGeneValueNumber, C3LOIndividual *pcInd, int  iPatternLevel, C3LOPattern *pcHighLvlGene, int *piGenotypeBuffer);
		void  v_count_higher_genes_freq_for_ind(C3LOIndividual  *pcInd, int iPatternLevel);
		void  v_update_higher_genes_freq_by_ind_for_best(C3LOIndividual  *pcInd);
		void  v_update_higher_genes_freq_by_ind_this_pop_for_best(C3LOIndividual  *pcInd);


		C3LOIndividual* pc_get_random_individual(int iLevel);
		C3LOIndividual *pc_get_individual_with_lvl_tournament(C3LOIndividual *pcDifferentToIndiv = NULL);
		C3LOIndividual *pc_get_individual_from_lvl(int  iEvolutionLvl, C3LOIndividual *pcDifferentToIndiv = NULL);
		bool  b_add_at_level(C3LOIndividual *pcInd, int  iLevel);
		bool  b_add_start_individual_copy(C3LOIndividual *pcInd, int  iLevel);
		C3LOIndividual* pc_insert_into_pop(vector  <C3LOIndividual  *>  *pvPopDest, C3LOIndividual  *pcIndToInsert, bool bCheckPhenotype = true, bool  bCheckGenotype = false);
		bool  b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes = true, bool  bCompareGenotypes = false);
		bool  b_can_ind_be_added_at_any_level(int  *piGeneString, bool  bComparePhenotypes = true, bool  bCompareGenotypes = false);
		bool  b_can_ind_be_added_at_pattern_level(int iPatternLevel, int  *piGeneString, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/);
		bool  b_produce_new_linkage(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, int  iPatternLevel);

		void  v_save_trees_and_patterns(CString  sPatternFile);
		void  v_load_isis_dsm(double  ***pdDSM, CString  sSource);
		void  v_report_dsm(double  **pdDSM, CString  sDest);
		void  v_report_dsm(int  **piDSM, CString  sDest);
		void  v_report_dsm_comparison(double  **pdDSM, double  **pdDSMisis, CString  sDest);
		int  i_dsm_get_miss_number(double  **pdDSM, double  **pdDSMisis);

		void  v_update_genes_values_occurrences(C3LOIndividual *pcInd);

		
	private:
		int  i_templ_length;
		vector<vector<int>> v_optimization_orders;//add levels and inner levels
		vector  <C3LOIndividual  *>  *pv_population;
		vector<C3LOIndividual  *>  v_individuals_to_process_next_gen;

		//vector<vector<C3LOPattern *>>  v_gene_patterns_parts_original;
		vector<vector<C3LOPattern *>>  v_gene_patterns_parts_interceptions_and_differencies;
		vector<vector<C3LOPattern *>>  v_gene_patterns_parts;
		vector<vector<C3LOPattern *>>  v_gene_patterns_trees;
		vector<vector<C3LOPattern *>>  v_higher_level_trees_genes;

		vector  <vector<C3LOIndividual  *>>  v_population_levels;//add inner levels
		vector  <vector<C3LOIndividual  *>>  v_population_levels_start_individuals;//add inner levels

		vector  <int  *>  v_genes_marked_by_linkage;//add inner levels

		CGenotypeInfoNode  *pc_genotype_root;

		double  **pd_dsm_matrix;
		int  **pi_gene_realtions;
		int  i_relation_length_max;


		//stats generation
		int  i_random_indiv_adds;
		time_t t_start;


		//for single iteration:
		C3LOIndividual  *pc_best;

		int  *pi_genotype_tool;
		int32_t  *pi_best_genotype;//for COptimizer best individual update

		CString  s_debug;

		int i_pop_to_add;
		int i_upgraded_new_ind_at_last_insert;
		int i_upgraded_new_ind;
		int i_worse_but_different_used;
		int i_worse_but_different_used_at_last_insert;
		int  i_crossings_with_effect_in_the_pop;
		int  i_all_crossings;
		int  i_all_effective_crossings;
		int  i_better_but_cannot_be_added_at_last_insert;
		int  i_better_but_cannot_be_added;
		int  i_tribe_effective_in_this_iteration;
		int  i_brutal_random_search;
		int  i_brutal_random_search_effective;
		double  d_unit_max_cur;
		int  i_cur_gen;
		bool  b_tribe_effective;
		int  i_the_same;
		int  i_improved;
		int  i_improved_other;
		int  i_shuffled_individuals;
		bool  b_linkage_generated;

		TimeCounters::
		CTimeCounter  c_time_counter;

		int  i_dsm_reporting_enumeration_tool;
		int  i_last_m_from_dsm;

		double d_best_fitness_prev;
		double d_best_fitness_cur;
		//C3LOIndividual  *pc_best_ind_0;
		
		C3LOIndividual  *pc_best_ind;
		C3LOIndividual  *pc_last_ind;

		C3LO  *pc_parent;
		vector<C3LOSingle *>  *pv_other_pops;

		int  i_linkage_generations;
		double  d_linkage_ffe;
		double  d_linkage_time;

		int  *pi_last_start_genotype;
		vector<int> v_last_optimization_order;
		bool  b_use_last_start_genotype;

		int  ***pppi_genes_values_occurrences;
	};//class CDummyOptimizer : public CBinaryOptimizer




	class  C3LOIndividual 
	{
	#define  C3LO_INDIVIDUAL_MAX_GENES_FOR_BS					12
	#define  C3LO_INDIVIDUAL_MAX_HIGH_LEVEL_GENES_VALUES					20100

		friend class C3LO;
		friend class C3LOSingle;
		friend class CGenotypeInfoNode;
		friend class C3LOHighLvlGeneFreq;
		friend class C3LOPattern;

	public:
		C3LOIndividual(int iLevel, int  iTemplLength, CProblem<CBinaryCoding, CBinaryCoding > *pcProblem, vector<vector<int>> *pvOptimizationOrders, C3LOSingle  *pcParent);
		~C3LOIndividual();

		void vRandomInit();
		double  dGetBestComputedFitness(int  **piBestPhenotype);
		double  dComputeFitness(int  iPatternLevel);
		//double  dComputeFitnessOfCurrentPhenotype(int  iPhenotype);
		void  vSetOptimized(bool  bOptimized);
		void  vGenerateOrder();


		void  vIncrementalBrutalUpdate(vector<C3LOPattern *> *pvGenePatterns, vector<C3LOPattern *> *pvGenePatternPartsOriginal, int iPatternLevel);
		bool  bLinkedGenesBrutalSearch(C3LOPattern *pcLinkedGenes, int iPatternLevel, vector<C3LOPattern *> *pvGenePatternPartsOriginal, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent);
		bool  bLinkedGenesBrutalSearchZero(C3LOPattern *pcLinkedGenes, vector<C3LOPattern *> *pvGenePatternPartsOriginal, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent);
		bool  bLinkedGenesBrutalSearchHigher(C3LOPattern *pcLinkedGenes, int iPatternLevel, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent);


		bool  bGetTreesGenerated(int  iPatternLevel);
		void  vGeneratePatternSingle(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel);
		void  vGeneratePatternSingle(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel);
		void  vGetPatternsParts(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel);
		void  vGetPatternsParts(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel);
		

		void  vCopyTo(C3LOIndividual  *pcOther);
		//void  vCopyPhenotypeTo(C3LOIndividual  *pcOther);

		bool  bComparePhenotype(int  iPatternLevel, int iOrder, C3LOIndividual  *pcOther);
		bool  bComparePhenotype(int  iPatternLevel, int iOrder, int  *piGeneString);
		bool  bCompareGenotype(C3LOIndividual  *pcOther);
		bool  bCompareGenotype(int  *piGeneString);
		int*  piGetGenotype() { return(pi_genotype); }
		int*  piGetPhenotype(int  iPatternLevel, int iOrder);
		int*  piGetBestPhenotype();

		void  vSetGenotype(int  *piNewGenotype);

		double  dUnitation(int  iPatternLevel = 0, int iOrder = 0);
		int  iGetOptDistance(int  iPatternLevel, int iOrder);


		void  vSaveGenePatternParts(CString  sDebugName);
		void  vSave(FILE  *pfReport, int  *piProblemBlockStructure, bool bDoNotForcePhenotype = false);
		//double  dComputeFitness(int  iPatternLevel, bool  bOptimize = true);
		//double  dComputeFitnessOfCurrentPhenotype(int  iPhenotype);*/

	private:
		double  d_compute_fitness_at_level_zero(int  iOrder, int  *piGenotypeToOptimize, int *piPhenotype);
		double  d_compute_fitness_at_level(int  iPatternLevel, int  iOrder, int  *piGenotypeToOptimize, int *piPhenotype);

		bool  b_brutal_check(vector<int>  *pvGenesShuffled, C3LOIndividual  *pcOptimizationBuf, C3LOIndividual  *pcBestThatFar, int  *piBestThatFarInitialized, int  iPos);
		bool  b_brutal_check_high_level(vector<C3LOPattern  *>  *pvHighLevelGenesShuffled, C3LOIndividual  *pcOptimizationBuf, C3LOIndividual  *pcBestThatFar, int  *piBestThatFarInitialized, int iPatternLevel, int  iPos);
		C3LOPattern*  pc_get_pattern_part(int  iGenePosition, int *piFilledGenes, vector<C3LOPattern *> *pvGenePatternPartsOriginal);

		void  vGetBestGenes
			(
				int  iBestGenesNumber, vector<C3LOHighLvlGeneFreq *> *pvResult, vector<C3LOHighLvlGeneFreq  *>  *pvGenesPool,
				int  *piOriginalGenotype, int  *piGenotypeBuffer,
				int  iPatternLevel
			);
		
		void  v_get_patterns_parts_zero_level(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts);
		void  v_get_patterns_parts_zero_level(vector<C3LOPattern  *>  *pvGenePatternsParts);
		void  v_get_patterns_parts_high_level(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel);
		void  v_get_patterns_parts_high_level(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel);
		//bool  b_phenotypes_the_same(int  iPosition, C3LOIndividual *pcOther);
		bool  b_phenotypes_the_same(int iPatternLevel, int iOrder, int  iPosition, C3LOIndividual *pcOther);
		void  v_add_phenotype_level(int iPhenotypeLevel);
		void  v_rem_phenotype_level(int iPhenotypeLevel);

		//double  d_compute_fitness_at_level_zero(bool  bOptimize);
		//double  d_compute_fitness_at_level(int  iPatternLevel);

		void  v_print_genotype(int  *piGenotypeToPrint, FILE  *pfDest);

		C3LOSingle  *pc_parent;
		int  i_level;
		int  i_level_inner;
		vector<int>  v_optimizations_on_levels;

		bool  b_trees_generated;
		vector<bool>  v_pattern_generated;
		vector<vector<C3LOPattern *>>  v_gene_patterns_parts;
		vector<vector<C3LOPattern *>>  v_gene_patterns_trees;

				

		CProblem<CBinaryCoding, CBinaryCoding > *pc_problem;
		vector<vector<int>> *pv_optimization_orders;
		vector<vector<int>> v_optimization_order;



		//data
		int  i_templ_length;
		int  *pi_genotype;
		vector<int  **> vpi_optimized_genotypes;
		vector<int  **> vpi_optimized_genotypes_for_linkage;
		vector<double  *> vpd_fit_values;
		vector<int *> v_optimized;

		//double  d_fit_value;

	};//class  C3LOIndividual



	class  C3LOPattern
	{
	#define  C3LO_PATTERN_TAB_INCREASE					"  "
		friend class C3LO;
		friend class C3LOSingle;
		friend class C3LOIndividual;
		friend class C3LOPatternPair;
	public:
		C3LOPattern(int iTemplLength);
		~C3LOPattern();

		vector<C3LOHighLvlGeneFreq  *>  *pvGetHighLevelGeneFreq() { return(&v_high_level_gene_freq); }
		vector<CMessyGene>  *pvGetPattern() { return(&v_pattern); };

		vector<C3LOPattern *>  *pvGetChildren() { return(&v_children); };

		int  iGetPatternLevel() { return(i_pattern_level); }

		void  vGetCommonAndUNcommonPositions(C3LOPattern  *pcOther, int *piCommonPos, int *piUncommonPos, int *piUnmarkedPos);
		void  vGetCommonAndUNcommonPositions(int  *piDifferences, int *piCommonPos, int *piUncommonPos, int *piUnmarkedPos);
		int  iCommonPositions(C3LOPattern  *pcOther, bool  bCheckValues = false, bool  bAllMustBeTheSame = true);
		int  iDifferentPositions(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int iPatternLevel);

		void  vZeroHighLvlGeneFreq();
		void  vBrutalValuesGenerationZeroLvl(C3LOIndividual *pcIndForContext, int iPatternLevel);
		void  vBrutalValuesGenerationHighLvl(C3LOIndividual *pcIndForContext, int iPatternLevel, vector<C3LOPattern *> *pvHighLvlGenesFroLowerLevel);
		void  vCountHighLvlGeneFreq(vector<C3LOIndividual *>  *pvPopoulation, int iPatternLevel);
		void  vShuffleFreqs();
		double  dCheckHighLevelGeneValues();

		C3LOPattern*  pcGetBestCrossingNode(int iPatternLevel, C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int *piDiffrentPositions = NULL, double  *pdDifrPosPerc = NULL);
		void  vGetRandomTreeBranch(int  *piPhenotype0, int *piPhenotype1, vector<C3LOPattern  *>  *pvTreeBranch, int iPatternLevel);
		C3LOPattern*  pcGetLongestReasonableBranch(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int *piDifferenceSearchTool, int iPatternLevel);
		void  vGetBestCovering(vector<C3LOPattern*> *pvDiffrCoveringPatterns, int  *piDifferences, int iAcceptableLevel);
		bool  bMarkedPhenotypeDifferences(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int iPatternLevel);
		bool  bMarkedPhenotypeDifferences(int  *piPhenotype0, int *piPhenotype1);
		bool  bUnMarkedPhenotypeDifferences(int  *piPhenotype0, int *piPhenotype1, int *piDifferenceSearchTool);


		void  vJoinSinglePatternAdd(C3LOPattern  *pcPatOther);
		void  vJoinTwoPatterns(C3LOPattern  *pcPat0, C3LOPattern  *pcPat1);


		bool  bDoesIndividualContainMarkedBlock(int  *piGenotype0, C3LOIndividual  *pcInd);
		bool  bDoesIndividualContainMarkedBlock(int  *piGenotype0, int  *piControlledGenotype);

		bool  bAmIInsideGeneGroup(int  *piGeneGroup);
		bool  bDoITouchGeneGroup(int  *piGeneGroup);
		bool  bDoITouchGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask);
		bool  bDoIExtendGeneGroup(int  *piGeneGroup);
		bool  bDoIExtendGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask);
		bool  bAmIThere(vector<C3LOPattern  *>  *pvGenePatternsParts, bool  bCheckValues = false);
		bool  bGeneThere(int  iGenePosition);

		void  vCopyPatternAndFreq(C3LOPattern  *pcOther);

		void  vSavePattern(FILE  *pfDest, int  *piProblemBlockStructure, CString  sTab = "", bool bOnlyRoot = false);
		void  vSaveHighLevelGene(FILE  *pfDest, int  *piProblemBlockStructure);


		void  vGeneratePatternTable();
		int  iGetSubtraction(C3LOPattern *pcOther, int  *piResult);
		int  iGetAND(C3LOPattern *pcOther, int  *piResult);
		bool  bEqualToDefinition(int  *piPatternDefinition);
		int  iGetHitNumber() { return(i_number_of_hits); }


		void  vRemoveFromDSM_SimilarityList(C3LOPattern  *pcOther);
		void  vGetAll_DSM_PairsWithSimilarity(double  dDSM_Similarity, vector<C3LOPatternPair> *pvDest);
		double  dGetDSM_Similarity(C3LOPattern  *pcOther, double  **pdDSM);
		double  dGetMaxDSM();


	private:
		void  v_check_high_level_gene_for_ind(C3LOIndividual *pcInd, int iPatternLevel);
		void  v_brutal_values_generate_zero_lvl(vector<CMessyGene>  *pvGenesChosen, int  iGeneOffset, int  *piGenotypeBuf);
		void  v_brutal_values_generate_high_lvl(vector<C3LOPattern *>  *pvHighLvlGenes, int  iGeneOffset, int  *piGenotypeBuf);

		int  i_pattern_level;
		int  i_templ_length;
		vector<CMessyGene>  v_pattern;

		vector<C3LOPattern *>  v_parents;
		vector<C3LOPattern *>  v_children;

		vector<C3LOPatternPair>  v_pattern_distances;


		vector<C3LOHighLvlGeneFreq  *>  v_high_level_gene_freq;

		int  *pi_pattern_table;
		int  i_number_of_hits;

		double  d_dsm_value;
	};//class  C3LOPattern


	class  C3LOPatternPair
	{
	public:
		C3LOPattern  *pc_first;
		C3LOPattern  *pc_second;
		double  d_similarity;

		double  dGetDSM_Similarity();

	};//class  C3LOPatternPair




	class  C3LOHighLvlGeneFreq
	{
	public:
		C3LOHighLvlGeneFreq(int  iTemplLength);
		~C3LOHighLvlGeneFreq();

		void  vCopy(C3LOHighLvlGeneFreq  *pcOther);

		void  vSetTemplLen(int  iTemplLength);
		void  vCreateTest(C3LOPattern  *pcPattern, int  iVal, int iPatternLevel);
		void  vCreate(C3LOPattern  *pcPattern, C3LOIndividual *pcInd, int iPatternLevel);
		void  vCreate(C3LOPattern  *pcPattern, int  *piGenotype);
		bool  bCheckHighLevelGene(C3LOIndividual *pcInd, int iPatternLevel);
		bool  bCheckHighLevelGene(int *piGenotypeToCheck);
		bool  bCheckAll(int iValue);

		void  vInfect(C3LOIndividual  *pcInd);
		void  vInfectPhenotype(C3LOIndividual  *pcInd, int iPatternLevel, int iOrder);
		void  vInfectPhenotype(int  *piGenotypeToInfcet);


		void  vSave(FILE  *pfDest, CString  sTab = " ");

		int  iGetGeneNumber() { return(i_gene_number); }

	private:
		int  *pi_gene_def;
		int  i_gene_number;

		int  i_templ_length;

	};//class  C3LOHighLvlGeneFreq




	class  CGenotypeInfoNode
	{
	public:
		CGenotypeInfoNode();
		~CGenotypeInfoNode();

		int  iGetSuccOpt() { return(i_succ_opt); }
		C3LOIndividual *pcCheckIndividual(C3LOIndividual *pcIndToCheck);

	private:
		int  i_succ_opt;
		int  i_value;
		vector<CGenotypeInfoNode *>  v_children;
		vector<C3LOIndividual *>  v_individuals;

		CGenotypeInfoNode  *pc_get_child_node(int  iVal);
		C3LOIndividual *pc_check_individual(C3LOIndividual *pcIndToCheck, int  iPos);


	};//class  CGenotypeInfoNode



	class  CMessyGene
	{
	public:
		CMessyGene();
		CMessyGene(int  iGeneVal, int  iGenePos);
		CMessyGene(CMessyGene  *pcOther);
		~CMessyGene();

		int  iGeneVal() { return(i_gene_val); }
		int  iGenePos() { return(i_gene_pos); }

	private:

		int  i_gene_val;
		int  i_gene_pos;
	};//class  CMessyGene
	
}//namespace 3LO


#endif//C3LO_OPTIMIZER_H
