#include "3LO.h"

#include "RandUtils.h"

using  namespace ThreeLO;


uint32_t C3LO::iERROR_PARENT_C3LOOptimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_C3LOOptimizer");

uint32_t C3LOSingle::iERROR_PARENT_C3LOSingleOptimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_C3LOSingleOptimizer");

uint32_t C3LOSingle::iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0 = CError::iADD_ERROR("iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0");




//---------------------------------------------C3LO-------------------------------------------------------
C3LO::C3LO(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	pd_dsm_matrix = NULL;	
	pi_gene_realtions = NULL;
	pc_3lo_single_tool = NULL;
	pi_best_genotype = NULL;
	pc_3lo_additional_pop = NULL;
	b_use_dsm = false;
}//C3LO::C3LO(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


C3LO::C3LO(C3LO *pcOther)	: CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "Brak implementacji: C3LO::C3LO(C3LO *pcOther) : CBinaryOptimizer(pcOther)", "BRAK", MB_OK);
}//C3LO::C3LO(C3LO *pcOther)

C3LO::~C3LO()
{
	if (pi_best_genotype != NULL)  delete  pi_best_genotype;


	for (int ii = 0; ii < v_3lo_pops.size(); ii++)
		delete  v_3lo_pops.at(ii);


	if (pd_dsm_matrix != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			delete  pd_dsm_matrix[ii];
		delete  pd_dsm_matrix;
	}//if (pd_dsm_matrix != NULL)


	if (pi_gene_realtions != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			delete  pi_gene_realtions[ii];
		delete  pi_gene_realtions;
	}//if (pi_gene_realtions != NULL)

	if (pc_3lo_single_tool != NULL)  delete  pc_3lo_single_tool;
	if (pc_3lo_additional_pop != NULL)  delete  pc_3lo_additional_pop;
	
}//C3LO::~C3LO()


CError C3LO::eConfigure(istream *psSettings)
{
	CError c_error = CBinaryOptimizer::eConfigure(psSettings);

	if (c_error == false)
	{
		CBoolCommandParam p_use_dsm(_3LO_ARGUMENT_USE_DSM, false, false);
		b_use_dsm = p_use_dsm.bGetValue(psSettings, &c_error);
	}//if (c_error == false)

	return c_error;
}//CError C3LO::eConfigure(istream *psSettings)


double C3LO::dComputeFitness(int32_t *piBits)
{
	double d_fitness_value = d_compute_fitness_value(piBits);
	return(d_fitness_value);
}//double C3LO::dComputeFitness(int32_t *piBits)



void C3LO::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);
	t_start = tStartTime;

	CError  c_err(iERROR_PARENT_C3LOOptimizer);
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	pc_3lo_single_tool = new C3LOSingle(pc_problem, pc_log, i_random_seed);
	pc_3lo_single_tool->pc_parent = this;
	pc_3lo_single_tool->vInitialize(tStartTime);
	pc_3lo_single_tool->pv_other_pops = &v_3lo_pops;

	if (pc_3lo_additional_pop != NULL)  delete  pc_3lo_additional_pop;
	pc_3lo_additional_pop = NULL;


	if (i_templ_length <= 0)
	{
		c_err.vSetError(C3LOSingle::iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)


	pc_log->vPrintLine("Initializing...", true);


	c_time_counter.vSetStartNow();
	pc_log->vPrintLine("DONE...", true);



	pd_dsm_matrix = new double*[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pd_dsm_matrix[ii] = new double[i_templ_length];


	pi_gene_realtions = new int*[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_realtions[ii] = new int[i_templ_length];

	i_relation_length_max = 0;
	i_dsm_reporting_enumeration_tool = 0;
	
	pi_best_genotype = new int32_t[i_templ_length];
}//void C3LO::vInitialize(time_t tStartTime)


void  C3LO::v_add_new_3lo_pop()
{
	C3LOSingle *pc_new_pop;
	pc_new_pop = new C3LOSingle(pc_problem, pc_log, i_random_seed);
	pc_new_pop->pc_parent = this;

	pc_new_pop->vInitialize(t_start);
	
	v_3lo_pops.push_back(pc_new_pop);

}//void  C3LO::v_add_new_3lo_pop()








void  C3LO::v_join_children_assign_to_tree_groups(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees)
{
	CString  s_buf;
	vector  <int>  v_fitting_groupped_trees_offsets;
	C3LOPattern  *pc_tree_group;


	for (int i_tree_group = 0; i_tree_group < pvGroupedChildrenTrees->size(); i_tree_group++)
	{
		pvGroupedChildrenTrees->at(i_tree_group)->vGeneratePatternTable();
		if (pcChildTree->bDoITouchGeneGroup(pvGroupedChildrenTrees->at(i_tree_group)->pi_pattern_table) == true)
		{
			v_fitting_groupped_trees_offsets.push_back(i_tree_group);
		}//if (pcChildTree->bDoITouchGeneGroup(pvGroupedChildrenTrees->at(ii)->pi_pattern_table) == true)
	}//for (int i_tree_group = 0; i_tree_group < pvGroupedChildrenTrees->size(); i_tree_group++)


	if (v_fitting_groupped_trees_offsets.size() == 1)
	{
		pc_tree_group = pvGroupedChildrenTrees->at(v_fitting_groupped_trees_offsets.at(0));

		if (pcChildTree->bDoIExtendGeneGroup(pc_tree_group->pi_pattern_table) == true)
			::MessageBox(NULL, "extend tree group - error", "extend tree group - error", MB_OK);

		pc_tree_group->v_children.push_back(pcChildTree);
	}//if (v_fitting_groupped_trees_offsets.size() == 1)
	else
	{
		s_buf.Format("wrong number of tree groups: %d", v_fitting_groupped_trees_offsets.size());
		::MessageBox(NULL, s_buf, s_buf, MB_OK);
	}//else  if (v_fitting_groupped_trees_offsets.size() == 1)

}//void  C3LO::v_join_children_assign_to_tree_groups(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees)




void  C3LO::v_join_children_refresh_grouped_trees(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees)
{
	vector  <int>  v_fitting_groupped_trees_offsets;


	for (int i_tree_group = 0; i_tree_group < pvGroupedChildrenTrees->size(); i_tree_group++)
	{
		pvGroupedChildrenTrees->at(i_tree_group)->vGeneratePatternTable();
		if (pcChildTree->bDoITouchGeneGroup(pvGroupedChildrenTrees->at(i_tree_group)->pi_pattern_table) == true)
		{
			v_fitting_groupped_trees_offsets.push_back(i_tree_group);
		}//if (pcChildTree->bDoITouchGeneGroup(pvGroupedChildrenTrees->at(ii)->pi_pattern_table) == true)
	}//for (int i_tree_group = 0; i_tree_group < pvGroupedChildrenTrees->size(); i_tree_group++)


	C3LOPattern  *pc_tree_group;


	if (v_fitting_groupped_trees_offsets.size() == 0)
	{
		pc_tree_group = new C3LOPattern(i_templ_length);
		pc_tree_group->v_pattern = pcChildTree->v_pattern;
		pc_tree_group->vGeneratePatternTable();

		pvGroupedChildrenTrees->push_back(pc_tree_group);
	}//if (v_fitting_groupped_trees_offsets.size() == 0)



	if (v_fitting_groupped_trees_offsets.size() == 1)
	{
		pc_tree_group = pvGroupedChildrenTrees->at(v_fitting_groupped_trees_offsets.at(0));

		if (pcChildTree->bDoIExtendGeneGroup(pc_tree_group->pi_pattern_table) == true)
		{
			for (int ii = 0; ii < pcChildTree->v_pattern.size(); ii++)
				pc_tree_group->pi_pattern_table[pcChildTree->v_pattern.at(ii).iGenePos()] = 1;

			pc_tree_group->v_pattern.clear();

			for (int ii = 0; ii < i_templ_length; ii++)
			{
				if  (pc_tree_group->pi_pattern_table[pcChildTree->v_pattern.at(ii).iGenePos()] == 1)
					pc_tree_group->v_pattern.push_back(CMessyGene(1, ii));
			}//for (int ii = 0; ii < i_templ_length; ii++)
		}//if (pcChildTree->bDoIExtendGeneGroup(pc_tree_group->pi_pattern_table) == true)
	}//if (v_fitting_groupped_trees_offsets.size() == 1)



	if (v_fitting_groupped_trees_offsets.size() > 1)
	{
		pc_tree_group = pvGroupedChildrenTrees->at(v_fitting_groupped_trees_offsets.at(0));

		for (int ii = v_fitting_groupped_trees_offsets.size() - 1; ii >= 1; ii--)//thanks this the trees may be ereased immedietly
		{
			pc_tree_group->vJoinSinglePatternAdd(pvGroupedChildrenTrees->at(v_fitting_groupped_trees_offsets.at(ii)));
			//delete  pvGroupedChildrenTrees->at(v_fitting_groupped_trees_offsets.at(ii));
			pvGroupedChildrenTrees->erase(pvGroupedChildrenTrees->begin() + v_fitting_groupped_trees_offsets.at(ii));
		}//for (int ii = 1; ii < v_fitting_groupped_trees_offsets.size(); ii++)
	}//if (v_fitting_groupped_trees_offsets.size() > 1)


}//void  C3LO::v_join_children_refresh_grouped_trees(C3LOPattern  *pcChildTree, vector<C3LOPattern *> *pvGroupedChildrenTrees)



void  C3LO::v_join_children_trees_if_necessary(vector<C3LOPattern  *> *pvTreesAtLevel, int iPatternLevel)
{
	while (iPatternLevel >= v_grouped_children_trees.size())  v_grouped_children_trees.push_back(vector<C3LOPattern  *>());

	for (int i_child_tree = 0; i_child_tree < pvTreesAtLevel->size(); i_child_tree++)
		v_join_children_refresh_grouped_trees(pvTreesAtLevel->at(i_child_tree), &(v_grouped_children_trees.at(iPatternLevel)));


	for (int ii = 0; ii < v_grouped_children_trees.at(iPatternLevel).size(); ii++)
		v_grouped_children_trees.at(iPatternLevel).at(ii)->v_children.clear();


	for (int i_child_tree = 0; i_child_tree < pvTreesAtLevel->size(); i_child_tree++)
		v_join_children_assign_to_tree_groups(pvTreesAtLevel->at(i_child_tree), &(v_grouped_children_trees.at(iPatternLevel)));


	C3LOPattern  *pc_joined_pattern, *pc_other_pattern;
	for (int i_groupped_trees = 0; i_groupped_trees < v_grouped_children_trees.at(iPatternLevel).size(); i_groupped_trees++)
	{
		if (v_grouped_children_trees.at(iPatternLevel).at(i_groupped_trees)->v_children.size() > 1)
		{

			pc_joined_pattern = v_grouped_children_trees.at(iPatternLevel).at(i_groupped_trees)->v_children.at(0);

			for (int ii = 1; ii < v_grouped_children_trees.at(iPatternLevel).at(i_groupped_trees)->v_children.size(); ii++)
			{
				pc_other_pattern = v_grouped_children_trees.at(iPatternLevel).at(i_groupped_trees)->v_children.at(ii);
				pc_joined_pattern->vJoinSinglePatternAdd(pc_other_pattern);

				//delete  pc_other_pattern;

				for (int i_tree = 0; i_tree < pvTreesAtLevel->size(); i_tree++)
				{
					if (pvTreesAtLevel->at(i_tree) == pc_other_pattern)
					{
						pvTreesAtLevel->erase(pvTreesAtLevel->begin() + i_tree);
						i_tree--;
					}//if (pvTreesAtLevel->at(i_tree) == pc_other_pattern)
				}//for (int i_tree = 0; i_tree < pvTreesAtLevel->size(); i_tree++)
			}//for (int ii = 1; ii < v_grouped_children_trees.at(iPatternLevel).at(i_groupped_trees)->v_children.size(); ii++)
		}//if (v_grouped_children_trees.at(iPatternLevel).at(ii)->v_children.size() > 1)
	}//for (int i_groupped_trees = 0; i_groupped_trees < v_grouped_children_trees.at(iPatternLevel).size(); i_groupped_trees++)


}//void  C3LO::v_join_children_trees_if_necessary(vector<C3LOPattern  *> *pvTreesAtLevel, int iPatternLevel)






bool bPatternLarger(C3LOPattern *elem1, C3LOPattern *elem2)
{
	return elem1->pvGetPattern()->size() < elem2->pvGetPattern()->size();
}//bool bPatternLarger (C3LOHighLvlGeneFreq *elem1, C3LOHighLvlGeneFreq *elem2)


void C3LO::v_execute_super_pop(uint32_t iIterationNumber, time_t tStartTime)
{
	CString  s_buf;

	s_buf.Format("Executing SUPERpop");
	pc_log->vPrintLine(s_buf, true);

	//v_build_consolidated_dsm();
	//v_build_gene_relations();

	bool  b_generate_new_trees = false;
	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		if (v_3lo_pops.at(i_3lo_pop)->bLinkageGenerated() == true)  b_generate_new_trees = true;
		v_3lo_pops.at(i_3lo_pop)->b_linkage_generated = false;
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)



	for  (int  i_pat_lvl = 0; i_pat_lvl < pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size(); i_pat_lvl++)
		pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).clear();
	

	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		for (int i_pat_lvl = 0; i_pat_lvl < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.size(); i_pat_lvl++)
		{
			while (pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size() < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.size())
				pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.push_back(vector<C3LOPattern *> ());

			for (int i_part = 0; i_part < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).size(); i_part++)
			{
				pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).push_back
				(
					v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).at(i_part)
				);
			}//for (int i_part = 0; i_part < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.at(0).size(); i_part++)
		}//for (int i_pat_lvl = 0; i_pat_lvl < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.size(); i_pat_lvl++)
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)


	std::sort(pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.at(0).begin(), pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.at(0).end(), bPatternLarger);


	if (b_generate_new_trees == true)
	{
		pc_log->vPrintLine("join linkage for superpop", true);

		for (int i_pat_lvl = 0; i_pat_lvl < pc_3lo_single_tool->v_gene_patterns_trees.size(); i_pat_lvl++)
		{
			for (int i_tree = 0; i_tree < pc_3lo_single_tool->v_gene_patterns_trees.at(i_pat_lvl).size(); i_tree++)
				delete  pc_3lo_single_tool->v_gene_patterns_trees.at(i_pat_lvl).at(i_tree);
		}//for (int i_pat_lvl = 0; i_pat_lvl < pc_3lo_single_tool->v_gene_patterns_trees.size(); i_pat_lvl++)
		pc_3lo_single_tool->v_gene_patterns_trees.clear();


		C3LOPattern  *pc_pattern;
		for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
		{
			for (int i_trees_lvl = 0; i_trees_lvl < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_trees.size(); i_trees_lvl++)
			{
				while (pc_3lo_single_tool->v_gene_patterns_trees.size() < i_trees_lvl + 1)  pc_3lo_single_tool->v_gene_patterns_trees.push_back(vector<C3LOPattern*>());

				for (int i_tree = 0; i_tree < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_trees.at(i_trees_lvl).size(); i_tree++)
				{

					if (
						v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_trees.at(i_trees_lvl).at(i_tree)->bAmIThere(&(pc_3lo_single_tool->v_gene_patterns_trees.at(i_trees_lvl)))
						==
						false
						)
					{
						pc_pattern = new C3LOPattern(i_templ_length);
						pc_pattern->vCopyPatternAndFreq(v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_trees.at(i_trees_lvl).at(i_tree));
						
						pc_3lo_single_tool->v_gene_patterns_trees.at(i_trees_lvl).push_back(pc_pattern);
					}//if  (				 
				}//for (int i_trees = 0; i_trees < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_trees.at(i_trees_lvl).size(); i_trees++)			
			}//for (int i_part = 0; i_part < v_3lo_pops.at(i_3lo_pop)->v_gene_patterns_parts_interceptions_and_differencies.at(0).size(); i_part++)
		}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)

		//pc_3lo_single_tool->v_count_higher_genes_freq_for_ind(pc_3lo_single_tool->pc_best_ind);

	}//if  (b_generate_new_trees == true)*/


	pc_3lo_single_tool->bRunIteration(iIterationNumber, tStartTime);
}//void C3LO::v_execute_super_pop(uint32_t iIterationNumber, time_t tStartTime)


void  C3LO::v_report_high_level_genes_numbers(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel)
{
	if (iPatternLevel >= pc3LOSingle->v_gene_patterns_trees.size())  return;


	CString  s_info;
	CString  s_buf;
	int i_gene_number;

	s_info = "";
	for (int i_high_lvl_gene = 0; i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size(); i_high_lvl_gene++)
	{		
		i_gene_number = pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).at(i_high_lvl_gene)->pvGetHighLevelGeneFreq()->size();
		s_buf.Format("%d", i_gene_number);
		if (i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size() - 1) s_buf += "/";

		s_info += s_buf;
	}//for (int i_high_lvl_gene = 0; i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size(); i_high_lvl_gene++)

	psHighLevelGenes->Format("[%d] HLGenes:", iPatternLevel);
	(*psHighLevelGenes) += s_info;
}//void  C3LO::v_report_high_level_genes_numbers(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel)



void  C3LO::v_report_high_level_genes(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel)
{
	if (iPatternLevel >= pc3LOSingle->v_gene_patterns_trees.size())  return;

	
	CString  s_info;
	CString  s_buf;
	int  i_ok_genes, i_halved_genes;
	i_ok_genes = 0;
	i_halved_genes = 0;
	double  d_gene_perc;
	s_info = "";
	for (int i_high_lvl_gene = 0; i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size(); i_high_lvl_gene++)
	{
		d_gene_perc = pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).at(i_high_lvl_gene)->dCheckHighLevelGeneValues();
		if (d_gene_perc == 1) i_ok_genes++;
		if (d_gene_perc == 0.5) i_halved_genes++;

		if (d_gene_perc == 1)
			s_buf = "OK";
		else
			s_buf.Format("%.1lf", d_gene_perc);
		if (i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size() - 1) s_buf += "/";

		s_info += s_buf;
	}//for (int i_high_lvl_gene = 0; i_high_lvl_gene < pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size(); i_high_lvl_gene++)

	psHighLevelGenes->Format("[%d] HLGenes (%d/%d/%d):", iPatternLevel, i_ok_genes, i_halved_genes, pc3LOSingle->v_gene_patterns_trees.at(iPatternLevel).size());

	(*psHighLevelGenes) += s_info;

}//void  C3LO::v_report_high_level_genes(C3LOSingle *pc3LOSingle, CString  *psHighLevelGenes, int iPatternLevel)



void  C3LO::v_report_scraps_crossings_trees(C3LOSingle *pc3LOSingle, CString  *psScraps, CString  *psCrossings, CString  *psTrees)
{
	CString  s_buf;

	*psScraps = "";
	for (int i_scrap = 0; i_scrap < pc3LOSingle->v_gene_patterns_parts.size(); i_scrap++)
	{
		if (i_scrap + 1 == pc3LOSingle->v_gene_patterns_parts.size())
			s_buf.Format("%d", pc3LOSingle->v_gene_patterns_parts.at(i_scrap).size());
		else
			s_buf.Format("%d/", pc3LOSingle->v_gene_patterns_parts.at(i_scrap).size());

		*psScraps += s_buf;
	}//for (int i_cross = 0; i_cross < pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size(); i_cross++)

	*psCrossings = "";
	for (int i_cross = 0; i_cross < pc3LOSingle->v_gene_patterns_parts_interceptions_and_differencies.size(); i_cross++)
	{
		if (i_cross + 1 == pc3LOSingle->v_gene_patterns_parts_interceptions_and_differencies.size())
			s_buf.Format("%d", pc3LOSingle->v_gene_patterns_parts_interceptions_and_differencies.at(i_cross).size());
		else
			s_buf.Format("%d/", pc3LOSingle->v_gene_patterns_parts_interceptions_and_differencies.at(i_cross).size());

		*psCrossings += s_buf;
	}//for (int i_cross = 0; i_cross < pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size(); i_cross++)

	*psTrees = "";
	for (int i_tree = 0; i_tree < pc3LOSingle->v_gene_patterns_trees.size(); i_tree++)
	{
		if (i_tree + 1 == pc3LOSingle->v_gene_patterns_trees.size())
			s_buf.Format("%d", pc3LOSingle->v_gene_patterns_trees.at(i_tree).size());
		else
			s_buf.Format("%d/", pc3LOSingle->v_gene_patterns_trees.at(i_tree).size());

		*psTrees += s_buf;
	}//for (int i_cross = 0; i_cross < pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size(); i_cross++)

}//void  C3LO::v_report_scraps_crossings_trees(C3LOSingle *pc3LOSingle, CString  *psScraps, CString  *psCrossings, CString  *psTrees)




void  C3LO::v_report_groupped_trees(CString  *psGrouppedTrees)
{
	CString  s_buf;

	*psGrouppedTrees = "";
	for (int i_tree = 0; i_tree < v_grouped_children_trees.size(); i_tree++)
	{
		if (i_tree + 1 == v_grouped_children_trees.size())
			s_buf.Format("%d", v_grouped_children_trees.at(i_tree).size());
		else
			s_buf.Format("%d/", v_grouped_children_trees.at(i_tree).size());

		*psGrouppedTrees += s_buf;
	}//for (int i_cross = 0; i_cross < pc_3lo_single_tool->v_gene_patterns_parts_interceptions_and_differencies.size(); i_cross++)

}//void  C3LO::v_report_groupped_trees(CString  *psGrouppedTrees)



CString  C3LO::sAdditionalSummaryInfo()
{
	CString  s_result;

	int  i_total_individuals;
	int  i_total_linkage_generations;
	double  d_total_linkage_ffe;
	double  d_total_linkage_time;

	i_total_linkage_generations = 0;
	d_total_linkage_ffe = 0;
	d_total_linkage_time = 0;
	i_total_individuals = 0;


	i_total_individuals += (int)pc_3lo_single_tool->v_population_levels.at(0).size();
	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		i_total_linkage_generations += v_3lo_pops.at(i_3lo_pop)->i_linkage_generations;
		i_total_individuals += (int)v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0).size();
		
		d_total_linkage_ffe += v_3lo_pops.at(i_3lo_pop)->d_linkage_ffe;
		d_total_linkage_time += v_3lo_pops.at(i_3lo_pop)->d_linkage_time;
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)

	s_result.Format("\tllinkGens:\t%d\tlinkTime:\t%.2lf\tlinkFFE:\t%.0lf\ttotalInd:\t%d\t", i_total_linkage_generations, d_total_linkage_time, d_total_linkage_ffe, i_total_individuals);

	return(s_result);
}//CString  C3LO::sAdditionalSummaryInfo()


bool C3LO::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	CString  s_buf;
	double  d_time_passed;
	bool  b_add_new_pop;

	if (v_3lo_pops.size() == 0)  
		v_add_new_3lo_pop();
	else
	{
		b_add_new_pop = true;
		for (int ii = 0; ii < v_3lo_pops.size(); ii++)
		{
			if (v_3lo_pops.at(ii)->d_best_fitness_prev < v_3lo_pops.at(ii)->d_best_fitness_cur)  b_add_new_pop = false;
		}//for (int ii = 0; ii < v_3lo_pops.size(); ii++)

		if (pc_3lo_single_tool->d_best_fitness_prev < pc_3lo_single_tool->d_best_fitness_cur)  b_add_new_pop = false;


		if (b_add_new_pop == true)
		{
			//try superPOP with additional linkage

			//if additional po is not ready - prepare it
			if (pc_3lo_additional_pop == NULL)
			{
				s_buf.Format("Executing new additional pop");
				pc_log->vPrintLine(s_buf, true);

				pc_3lo_additional_pop = new C3LOSingle(pc_problem, pc_log, i_random_seed);
				pc_3lo_additional_pop->pc_parent = this;
				pc_3lo_additional_pop->vInitialize(tStartTime);
				pc_3lo_additional_pop->bRunIteration(iIterationNumber, tStartTime);
			}//if (pc_3lo_additional_pop == NULL)

			/*C3LOIndividual  *pc_last_individual;
			pc_last_individual = pc_3lo_single_tool->pc_last_ind;
			pc_3lo_single_tool->pc_last_ind = NULL;

			for (int ii = 0; ii < pc_3lo_single_tool->v_population_levels.at(0).size(); ii++)
			{
				if (pc_3lo_single_tool->v_population_levels.at(0).at(ii) == pc_last_individual)
					pc_3lo_single_tool->v_population_levels.at(0).erase(pc_3lo_single_tool->v_population_levels.at(0).begin() + ii);
			}//for (int ii = 0; ii < pc_3lo_single_tool->v_population_levels.at(0).size(); ii++)*/

			v_3lo_pops.push_back(pc_3lo_additional_pop);


			/*MessageBox(NULL, "test", "test", MB_OK);

			pc_3lo_single_tool->b_use_last_start_genotype = true;*/
			v_execute_super_pop(iIterationNumber, tStartTime);

			/*if (pc_last_individual->bComparePhenotype(0, 0, pc_3lo_single_tool->pc_last_ind->piGetPhenotype(0, 0)) == true)
				MessageBox(NULL, "the same", "the same", MB_OK);
			else
				MessageBox(NULL, "different", "different", MB_OK);*/



			//if (pc_3lo_single_tool->b_can_ind_be_added_at_pattern_level(pc_last_individual->i_level, pc_last_individual->piGetPhenotype(0, 0), true, true) == true)
				//pc_3lo_single_tool->v_population_levels.at(0).push_back(pc_last_individual);

			if (pc_3lo_single_tool->d_best_fitness_prev < pc_3lo_single_tool->d_best_fitness_cur)
			{
				//if adding the linkage from pc_3lo_additional_pop has improved anything - leave the already added pop
				pc_3lo_additional_pop = NULL;//mark that the new additional po is to be created
				//MessageBox(NULL, "adding new pop", "adding new pop", MB_OK);
			}//if (pc_3lo_single_tool->d_best_fitness_prev < pc_3lo_single_tool->d_best_fitness_cur)
			else
			{
				//remove the additional pop from evolving pops
				v_3lo_pops.erase(v_3lo_pops.begin() + v_3lo_pops.size() - 1);
			}//else  if (pc_3lo_single_tool->d_best_fitness_prev < pc_3lo_single_tool->d_best_fitness_cur)

		}//if (b_add_new_pop == true)
	}//else  if (v_3lo_pops.size() == 0)  


	 //if (b_add_new_pop == true)  v_add_new_3lo_pop();



	if (v_3lo_pops.size() > 0)
	{
		v_execute_super_pop(iIterationNumber, tStartTime);
	}//if ( (b_add_new_pop == true)&&(v_3lo_pops.size() > 0) )



	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		s_buf.Format("Executing population %d:", i_3lo_pop);
		pc_log->vPrintLine(s_buf, true);
		v_3lo_pops.at(i_3lo_pop)->bRunIteration(iIterationNumber, tStartTime);
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)


	pc_log->vPrintLine("", true);
	pc_log->vPrintLine("", true);
	s_buf.Format("Population state report (%d):", v_3lo_pops.size());
	pc_log->vPrintLine(s_buf, true);


	double  d_overall_best_fitness;
	int  *pi_overall_best_phenotype, *pi_phenotype_buf;
	int  i_total_linkage_generations;
	double  d_total_linkage_ffe;
	double  d_total_linkage_time;

	i_total_linkage_generations = 0;
	d_total_linkage_ffe = 0;
	d_total_linkage_time = 0;
	CString  s_info, s_scraps, s_crossings, s_trees;
	CString  s_high_level_genes;

	int  i_total_individuals;
	i_total_individuals = 0;
	pi_overall_best_phenotype = NULL;
	if (pc_3lo_single_tool->pc_best_ind != NULL)
	{
		v_report_scraps_crossings_trees(pc_3lo_single_tool, &s_scraps, &s_crossings, &s_trees);

		
		s_info.Format
		(
			"[Population SUPERpop] best: %.4lf last:%.4lf ZERO: %.4lf/%.4lf size:%d  scraps: %s  crossings: %s  trees: %s  LatsM: %d  linkageGens: %d",
			pc_3lo_single_tool->pc_best_ind->dGetBestComputedFitness(&pi_phenotype_buf), pc_3lo_single_tool->pc_last_ind->dGetBestComputedFitness(&pi_phenotype_buf),
			pc_3lo_single_tool->d_best_fitness_cur, pc_3lo_single_tool->d_best_fitness_prev, 
			(int)pc_3lo_single_tool->v_population_levels.at(0).size(), 
			
			s_scraps, s_crossings, s_trees,
			pc_3lo_single_tool->i_last_m_from_dsm, pc_3lo_single_tool->i_linkage_generations
		);
		pc_log->vPrintLine(s_info, true);

		d_overall_best_fitness = pc_3lo_single_tool->pc_best_ind->dGetBestComputedFitness(&pi_overall_best_phenotype);
		//d_overall_best_fitness = pc_3lo_single_tool->pc_best_ind->dComputeFitness(-1);
		//pi_overall_best_phenotype = pc_3lo_single_tool->pc_best_ind->piGetPhenotype(0, 0);
		i_total_linkage_generations += pc_3lo_single_tool->i_linkage_generations;

		i_total_individuals += (int)pc_3lo_single_tool->v_population_levels.at(0).size();
	}//if (pc_3lo_single_tool->pc_best_ind != NULL)
	else
	{
		d_overall_best_fitness = v_3lo_pops.at(0)->pc_best_ind->dGetBestComputedFitness(&pi_overall_best_phenotype);
		//d_overall_best_fitness = v_3lo_pops.at(0)->pc_best_ind->dComputeFitness(-1);
		//pi_overall_best_phenotype = v_3lo_pops.at(0)->pc_best_ind->piGetPhenotype(0, 0);
	}//else  if (pc_3lo_single_tool->pc_best_ind != NULL)


	
	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		v_report_scraps_crossings_trees(v_3lo_pops.at(i_3lo_pop), &s_scraps, &s_crossings, &s_trees);

		s_info.Format
		(
			"[Population %d] best: %.4lf last:%.4lf ZERO: %.4lf/%.4lf  size:%d  scraps: %s  crossings: %s  trees: %s  LatsM: %d  linkageGens: %d",
			i_3lo_pop,
			v_3lo_pops.at(i_3lo_pop)->pc_best_ind->dGetBestComputedFitness(&pi_phenotype_buf), v_3lo_pops.at(i_3lo_pop)->pc_last_ind->dGetBestComputedFitness(&pi_phenotype_buf),
			v_3lo_pops.at(i_3lo_pop)->d_best_fitness_cur, v_3lo_pops.at(i_3lo_pop)->d_best_fitness_prev,
			(int)v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0).size(), 

			s_scraps, s_crossings, s_trees,
			v_3lo_pops.at(i_3lo_pop)->i_last_m_from_dsm, v_3lo_pops.at(i_3lo_pop)->i_linkage_generations
		);
		i_total_individuals += (int)v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0).size();
		pc_log->vPrintLine(s_info, true);

		if (v_3lo_pops.at(i_3lo_pop)->pc_best_ind->dGetBestComputedFitness(&pi_phenotype_buf) > d_overall_best_fitness)
		{
			d_overall_best_fitness = v_3lo_pops.at(i_3lo_pop)->pc_best_ind->dGetBestComputedFitness(&pi_overall_best_phenotype);
			//d_overall_best_fitness = v_3lo_pops.at(i_3lo_pop)->pc_best_ind->dComputeFitness(-1);
			//pi_overall_best_phenotype = v_3lo_pops.at(i_3lo_pop)->pc_best_ind->piGetPhenotype(0, 0);
		}//if (v_3lo_pops.at(i_3lo_pop)->pc_best_ind->dComputeFitness(-1) > d_overall_best)

		i_total_linkage_generations += v_3lo_pops.at(i_3lo_pop)->i_linkage_generations;
		d_total_linkage_ffe += v_3lo_pops.at(i_3lo_pop)->d_linkage_ffe;
		d_total_linkage_time += v_3lo_pops.at(i_3lo_pop)->d_linkage_time;
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)


	double  d_ffe_total = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	c_time_counter.bGetTimePassed(&d_time_passed);
	s_buf.Format("SUMARY best: %.4lf linkGens: %d   pops: %d (inds:%d) [time:%.2lf]", d_overall_best_fitness, i_total_linkage_generations, v_3lo_pops.size(), i_total_individuals, d_time_passed);
	pc_log->vPrintLine(s_buf, true);
	s_buf.Format
		(
		"SUMARY LINKAGE: time %.2lf/%.2lf (%.4lf)    FFE: %.2lf/%.2lf (%.4lf)", 
		d_total_linkage_time, d_time_passed, d_total_linkage_time / d_time_passed,
		d_total_linkage_ffe, d_ffe_total, d_total_linkage_ffe / d_ffe_total
		);
	pc_log->vPrintLine(s_buf, true);

	v_report_groupped_trees(&s_buf);
	s_buf = "SUMMARY LINKAGE QUALITY/ groupped trees:" + s_buf;
	pc_log->vPrintLine(s_buf, true);

	/*v_report_high_level_genes(pc_3lo_single_tool, &s_high_level_genes, 0);
	pc_log->vPrintLine(s_high_level_genes, true);

	for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		v_report_high_level_genes(v_3lo_pops.at(i_3lo_pop), &s_high_level_genes, 0);
		pc_log->vPrintLine(s_high_level_genes, true);
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)*/

	/*for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)
	{
		v_report_high_level_genes_numbers(v_3lo_pops.at(i_3lo_pop), &s_high_level_genes, 0);
		pc_log->vPrintLine(s_high_level_genes, true);
	}//for (int i_3lo_pop = 0; i_3lo_pop < v_3lo_pops.size(); i_3lo_pop++)*/

	
	pc_log->vPrintLine("", true);
	pc_log->vPrintLine("", true);



	if (pi_overall_best_phenotype != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_best_genotype[ii] = (int32_t)(pi_overall_best_phenotype[ii]);

		b_update_best_individual(iIterationNumber, tStartTime, pi_best_genotype, d_overall_best_fitness);
	}//if (pi_best_phenotype != NULL)


	return(true);
}//bool C3LO::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)











bool  C3LO::b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int iPatternLevel)
{
	CString  s_buf;

	pcParentMain->dComputeFitness(iPatternLevel);
	pcParentOther->dComputeFitness(iPatternLevel);

	int  i_max_ind_level;


	vector<C3LOPattern *>  v_allowed_trees;


	int  i_result = 0;

	//s_buf.Format("PROCESSING... : ");
	//pc_log->vPrintLine(s_buf, true);
	int  i_preserve_inner_level;
	i_preserve_inner_level = pcParentMain->i_level_inner;
	i_result = pc_3lo_single_tool->i_cross_dsm_glueing_similarities_from_different_start_genes(pi_gene_realtions, piCurrentGenotype, pcParentMain, pcParentOther, iPatternLevel);

	pcParentMain->i_level_inner = i_preserve_inner_level;
	//i_result = i_cross_dsm(piCurrentGenotype, pcParentMain, pcParentOther, i_pattern_level);
	//i_result = i_cross_fluent_scraps(piCurrentGenotype, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, i_pattern_level);
	//i_result = i_cross_intercept_scraps(piCurrentGenotype, pcParentMain, pcParentOther, i_pattern_level);


	if (i_result == 2)  return(true);

	return(false);
}//bool  C3LO::b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int iPatternLevel)



//---------------------------------------------C3LOSingle-------------------------------------------------------
C3LOSingle::C3LOSingle(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	pv_population = new vector  <C3LOIndividual  *>();

	i_random_indiv_adds = 0;

	pc_genotype_root = new  CGenotypeInfoNode();

	pi_best_genotype = NULL;
	pi_genotype_tool = NULL;
	pd_dsm_matrix = NULL;
	pc_last_ind = NULL;

	pc_parent = NULL;
	pv_other_pops = NULL;

	pi_last_start_genotype = NULL;
}//C3LO::C3LO(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, , uint32_t iRandomSeed)


C3LOSingle::C3LOSingle(C3LOSingle *pcOther)
	: CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "Brak implementacji: C3LO::C3LO(C3LO *pcOther) : CBinaryOptimizer(pcOther)", "BRAK", MB_OK);
}//C3LO::C3LO(C3LO *pcOther)


C3LOSingle::~C3LOSingle()
{
	if (pi_best_genotype != NULL)  delete  pi_best_genotype;
	if (pi_genotype_tool != NULL)  delete  pi_genotype_tool;
	if (pi_last_start_genotype != NULL)  delete  pi_last_start_genotype;
		
	
	
	//v_delete_population(pv_population);
	if (pv_population != NULL)
	{
		for (int ii = 0; ii < (int)pv_population->size(); ii++)
		{
			delete  pv_population->at(ii);
		}//for  (int  ii = 0; ii < (int) pv_population->size(); ii++)
		delete  pv_population;
	}//if  (pvPopulationToDelete  !=  NULL)


	//for (int ii = 0; ii < v_new_individuals.size(); ii++)
		//delete  v_new_individuals.at(ii);


	//if  (pi_genes_marked_by_linkage != NULL)  delete  pi_genes_marked_by_linkage;

	delete  pc_genotype_root;


	for (int ii = 0; ii < v_genes_marked_by_linkage.size(); ii++)
		delete  v_genes_marked_by_linkage.at(ii);



	if (pd_dsm_matrix != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			delete  pd_dsm_matrix[ii];
		delete  pd_dsm_matrix;
	}//if (pd_dsm_matrix != NULL)


	if (pi_gene_realtions != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			delete  pi_gene_realtions[ii];
		delete  pi_gene_realtions;
	}//if (pi_gene_realtions != NULL)


	if (pc_parent->b_use_dsm == true)
	{
		if (pppi_genes_values_occurrences != NULL)
		{
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				for (int jj = 0; jj < i_templ_length; jj++)
					delete pppi_genes_values_occurrences[ii][jj];

				delete pppi_genes_values_occurrences[ii];
			}//for (int ii = 0; ii < i_templ_length; ii++)

			delete pppi_genes_values_occurrences;
		}//if (ppi_genes_values_occurrences != NULL)
	}//if (pc_parent->b_use_dsm == true)


	//not an owner...
	/*for  (int  ii = 0; ii < v_gene_patterns.size(); ii++)
	delete  v_gene_patterns.at(ii);

	for  (int  ii = 0; ii < v_gene_patterns_parts.size(); ii++)
	delete  v_gene_patterns_parts.at(ii);*/
}//C3LO::~C3LO()




void  C3LOSingle::v_add_order_random(vector <int>  *pvNewOrder, vector <int>  *pvOppositeNewOrder, bool  bWithOpposite)
{
	pvNewOrder->clear();

	vector <int>  v_genotype_order;
	//vector <int>  v_new_order;

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		v_genotype_order.push_back(ii);
	}//for  (int ii = 0; ii < i_templ_length; ii++)

	int  i_gene_pos_offset;
	while (v_genotype_order.size() > 0)
	{
		//i_gene_pos_offset = pc_random_gen->Next(0, v_genotype_order.size());
		i_gene_pos_offset = RandUtils::iRandNumber(0, v_genotype_order.size() - 1);
		pvNewOrder->push_back(v_genotype_order.at(i_gene_pos_offset));

		v_genotype_order.erase(v_genotype_order.begin() + i_gene_pos_offset);
	}//while  (v_genotype_order.size() > 0)

	//v_optimization_orders.push_back(v_new_order);


	if ( (bWithOpposite == true)&&(pvOppositeNewOrder !=  NULL) )
	{
		for (int ii = 0; ii < i_templ_length; ii++)
		{
			pvOppositeNewOrder->push_back(pvNewOrder->at(i_templ_length - ii - 1));
		}//for  (int  ii = 0; ii < i_templ_length; ii++)	

		//v_optimization_orders.push_back(v_opposite_order);
	}//if ( (bWithOpposite == true)&&(pvOppositeNewOrder !=  NULL) )
}//void  C3LO::v_add_order_random()



void  C3LOSingle::v_add_order_by_genotype(bool  bReversed)
{
	vector <int>  v_new_order;

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (bReversed == false)
			v_new_order.push_back(ii);
		else
			v_new_order.push_back(i_templ_length - ii - 1);
	}//for  (int ii = 0; ii < i_templ_length; ii++)

	v_optimization_orders.push_back(v_new_order);
}//void  C3LO::v_add_order_by_genotype(bool  bReversed)




void  C3LOSingle::v_initialize_orders()
{
	v_optimization_orders.clear();
	v_add_order_by_genotype(false);
	//v_add_order_by_genotype(true);

	//for  (int  ii = 0; ii < 50; ii++)
	//	v_add_order_random(true);

}//void  C3LO::v_initialize_orders()



void  C3LOSingle::v_insert_gene_pattern_part_level(int  iPatternPartLevel)
{
	if (iPatternPartLevel  <  0)  return;

	while (v_gene_patterns_parts.size() <= iPatternPartLevel)  v_gene_patterns_parts.push_back(vector<C3LOPattern *>());
	while (v_gene_patterns_parts_interceptions_and_differencies.size() <= iPatternPartLevel)  v_gene_patterns_parts_interceptions_and_differencies.push_back(vector<C3LOPattern *>());
	while (v_gene_patterns_trees.size() <= iPatternPartLevel)  v_gene_patterns_trees.push_back(vector<C3LOPattern *>());

	
//	while (v_gene_patterns_parts_original.size() <= iPatternPartLevel)  v_gene_patterns_parts_original.push_back(vector<C3LOPattern *>());
	
}//void  C3LO::v_insert_gene_pattern_part_level(int  iPatternPartLevel)



void  C3LOSingle::v_add_gene_pattern_tree(C3LOPattern  *pcTree, int  iPatternLevel)
{
	if (iPatternLevel  <  0)  return;

	while (v_gene_patterns_trees.size() <= iPatternLevel)  v_gene_patterns_trees.push_back(vector<C3LOPattern *>());

	if (pcTree == NULL)  return;
	v_gene_patterns_trees.at(iPatternLevel).push_back(pcTree);
}//void  C3LO::v_add_gene_pattern_tree(C3LOPattern  *pcPattern, int  iPatternLevel)




double  C3LOSingle::d_get_most_similar_pattern_pair_offset(int  *piMostSimlarPairOffset, vector<C3LOPatternPair>  *pvPatternPairs)
{
	double  d_most_similar = 0;
	int  i_most_similar_offset = -1;
	for (int ii = 0; ii < pvPatternPairs->size(); ii++)
	{
		if (pvPatternPairs->at(ii).d_similarity > d_most_similar)
		{
			d_most_similar = pvPatternPairs->at(ii).d_similarity;
			i_most_similar_offset = ii;
		}//if  (v_pattern_pairs.at(ii).d_similarity > d_most_similar)
	}//for  (int ii = 0; ii < v_pattern_pairs.size(); ii++)

	*piMostSimlarPairOffset = i_most_similar_offset;
	return(d_most_similar);
}//double  C3LO::d_get_most_similar_pattern_pair_offset()





void  C3LOSingle::v_remove_pattern_from_buffer(vector<C3LOPattern  *>  *pvPartsBuffer, C3LOPattern  *pcPatternToRemove)
{
	//PRW FINISHED HERE

	for (int ii = 0; ii < pvPartsBuffer->size(); ii++)
	{
		if (pvPartsBuffer->at(ii) == pcPatternToRemove)
		{
			pvPartsBuffer->erase(pvPartsBuffer->begin() + ii);
			return;
		}//if  (pvPartsBuffer->at(ii) == pcPatternToRemove)
	}//for  (int  ii = 0; ii < pvPartsBuffer->size() ii++)*/
}//void  C3LO::v_remove_pattern_from_buffer(vector<C3LOPattern  *>  *pvPartsBuffer, C3LOPattern  *pcPatternToRemove)



void  C3LOSingle::v_remove_pattern_from_pairs_buffer(vector<C3LOPatternPair>  *pvPatternPairsBuffer, C3LOPattern  *pcPatternToRemove)
{
	for (int ii = 0; ii < pvPatternPairsBuffer->size(); ii++)
	{
		if (
			(pvPatternPairsBuffer->at(ii).pc_first == pcPatternToRemove) ||
			(pvPatternPairsBuffer->at(ii).pc_second == pcPatternToRemove)
			)
		{
			pvPatternPairsBuffer->erase(pvPatternPairsBuffer->begin() + ii);
			ii--;
		}//if  (pvPartsBuffer->at(ii) == pcPatternToRemove)
	}//for  (int  ii = 0; ii < pvPartsBuffer->size() ii++)*/	


}//void  C3LO::v_remove_pattern_from_pairs_buffer(vector<C3LOPatternPair>  *pvPatternPairsBuffer, C3LOPattern  *pcPatternToRemove)






bool  C3LOSingle::b_add_higher_level_gene(C3LOPattern  *pcHigherLevelGene, int  iGeneLevel)
{
	if (iGeneLevel  <  0)  return(false);

	while (v_higher_level_trees_genes.size() <= iGeneLevel)  v_higher_level_trees_genes.push_back(vector<C3LOPattern *>());

	if (pcHigherLevelGene == NULL)  return(false);
	if (pcHigherLevelGene->v_pattern.size() == 0)  return(false);

	//check ifit is not contained 
	int  i_common, i_uncommon, i_unmarked;
	for (int ii = 0; ii < v_higher_level_trees_genes.at(iGeneLevel).size(); ii++)
	{
		v_higher_level_trees_genes.at(iGeneLevel).at(ii)->vGetCommonAndUNcommonPositions(pcHigherLevelGene, &i_common, &i_uncommon, &i_unmarked);

		if (i_unmarked == 0)  return(false);
	}//for  (int  ii = 0; ii < this->v_gene_higher_level_trees_genes.at(iGeneLevel).size(); ii++)


	 //now check if it contains...
	for (int ii = 0; ii < v_higher_level_trees_genes.at(iGeneLevel).size(); ii++)
	{
		pcHigherLevelGene->vGetCommonAndUNcommonPositions(v_higher_level_trees_genes.at(iGeneLevel).at(ii), &i_common, &i_uncommon, &i_unmarked);

		if (i_unmarked == 0)
		{
			delete  v_higher_level_trees_genes.at(iGeneLevel).at(ii);
			v_higher_level_trees_genes.at(iGeneLevel).erase(v_higher_level_trees_genes.at(iGeneLevel).begin() + ii);
			ii--;
		}//if  ( (i_uncommon == 0)&&(i_unmarked == 0)  )
	}//for  (int  ii = 0; ii < this->v_gene_higher_level_trees_genes.at(iGeneLevel).size(); ii++)

	v_higher_level_trees_genes.at(iGeneLevel).push_back(pcHigherLevelGene);



	return(true);
}//void  C3LO::v_insert_gene_pattern_part_level(int  iPatternPartLevel)



int  C3LOSingle::i_get_ind_number()
{
	int  i_cur_pop_size;
	i_cur_pop_size = 0;
	for (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)
	{
		i_cur_pop_size += v_population_levels.at(i_lvl).size();
	}//for  (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)

	return(i_cur_pop_size);
}//int  C3LO::i_get_ind_number()




void  C3LOSingle::v_update_higher_genes_freq_by_ind_for_best(C3LOIndividual  *pcInd)
{
	for (int i_3lo_pop = 0; i_3lo_pop < pc_parent->v_3lo_pops.size(); i_3lo_pop++)
		pc_parent->v_3lo_pops.at(i_3lo_pop)->v_update_higher_genes_freq_by_ind_this_pop_for_best(pcInd);

	pc_parent->pc_3lo_single_tool->v_update_higher_genes_freq_by_ind_this_pop_for_best(pcInd);

}//void  C3LOSingle::v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd)



void  C3LOSingle::v_update_higher_genes_freq_by_ind_this_pop_for_best(C3LOIndividual  *pcInd)
{
	if (pc_best_ind == NULL)  return;

	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		if (v_gene_patterns_trees.at(i_pattern_level).size() > 1)
		{
			//if (v_population_levels.size() > iPatternLevel)
			{
				for (int i_hl_gene = 0; i_hl_gene < v_gene_patterns_trees.at(i_pattern_level).size(); i_hl_gene++)
				{
					if (pc_parent != NULL)
					{
						v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->v_check_high_level_gene_for_ind
						(
							pcInd, i_pattern_level + 1
						);

						v_filter_gene_freq_by_ind(2, pc_best_ind, i_pattern_level, v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene), pi_genotype_tool);
					}//if (pc_parent != NULL)
				}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
			}//if  (v_population_levels.size() > iPatternLevel)
		}//if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
	}//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	return;
}//void  C3LOSingle::v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd)







void  C3LOSingle::v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd)
{
	for (int i_3lo_pop = 0; i_3lo_pop < pc_parent->v_3lo_pops.size(); i_3lo_pop++)
		pc_parent->v_3lo_pops.at(i_3lo_pop)->v_update_higher_genes_freq_by_ind_this_pop(pcInd);

	pc_parent->pc_3lo_single_tool->v_update_higher_genes_freq_by_ind_this_pop(pcInd);
	
}//void  C3LOSingle::v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd)



void  C3LOSingle::v_update_higher_genes_freq_by_ind_this_pop(C3LOIndividual  *pcInd)
{
	return;
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		if (v_gene_patterns_trees.at(i_pattern_level).size() > 1)
		{
			//if (v_population_levels.size() > iPatternLevel)
			{
				for (int i_hl_gene = 0; i_hl_gene < v_gene_patterns_trees.at(i_pattern_level).size(); i_hl_gene++)
				{
					if (pc_parent != NULL)
					{
						v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->v_check_high_level_gene_for_ind
							(
								pcInd, i_pattern_level + 1
							);
					}//if (pc_parent != NULL)
				}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
			}//if  (v_population_levels.size() > iPatternLevel)
		}//if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
	}//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	return;
}//void  C3LOSingle::v_update_higher_genes_freq_by_ind(C3LOIndividual  *pcInd)



void  C3LOSingle::v_filter_gene_freq_by_ind(int iBestGeneValueNumber, C3LOIndividual *pcInd, int  iPatternLevel, C3LOPattern *pcHighLvlGene, int *piGenotypeBuffer)
{
	vector<C3LOHighLvlGeneFreq *>  *pv_gene_freqs;
	vector<C3LOHighLvlGeneFreq *>  v_chosen_genes;
	int  i_old_level;

	pv_gene_freqs = pcHighLvlGene->pvGetHighLevelGeneFreq();
	i_old_level = pcInd->i_level;
	pcInd->i_level = iPatternLevel + 1;
	pcInd->dComputeFitness(iPatternLevel + 1);
	pcInd->
		vGetBestGenes
		(
			2, &v_chosen_genes, pv_gene_freqs,
			pcInd->vpi_optimized_genotypes.at(iPatternLevel + 1)[0], piGenotypeBuffer,
			iPatternLevel + 2
		);
	pcInd->i_level = i_old_level;

	//preserve chosen genes
	for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)
	{
		for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
		{
			if (v_chosen_genes.at(i_chosen) == pv_gene_freqs->at(i_gene_freq))
			{
				pv_gene_freqs->erase(pv_gene_freqs->begin() + i_gene_freq);
				i_gene_freq--;
			}//if (v_chosen_genes.at(i_chosen) == pv_gene_freqs->at(i_gene_freq))
		}//for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
	}//for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)

	for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
		delete  pv_gene_freqs->at(i_gene_freq);
	pv_gene_freqs->clear();

	for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)
		pv_gene_freqs->push_back(v_chosen_genes.at(i_chosen));

	v_chosen_genes.clear();
}//void  C3LOSingle::v_filter_gene_freq_by_ind(int iBestGeneValueNumber, C3LOIndividual *pcInd, int  iPatternLevel, C3LOPattern *pcHighLvlGene)






void  C3LOSingle::v_count_higher_genes_freq_for_ind(C3LOIndividual  *pcInd, int iPatternLevel)
{
	double  d_ffe_before, d_ffe_after;
	double  d_time_before, d_time_after;
	int  i_old_level;

	d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	c_time_counter.bGetTimePassed(&d_time_before);

	vector<C3LOHighLvlGeneFreq *>  *pv_gene_freqs;
	vector<C3LOHighLvlGeneFreq *>  v_chosen_genes;
	int  *pi_genotype_buffer;
	pi_genotype_buffer = new int[i_templ_length];


	//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
		{
			//if (v_population_levels.size() > iPatternLevel)
			{
				for (int i_hl_gene = 0; i_hl_gene < v_gene_patterns_trees.at(iPatternLevel).size(); i_hl_gene++)
				{
					v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vZeroHighLvlGeneFreq();

					if (pc_parent != NULL)
					{
						//if (pc_parent->pc_3lo_single_tool == this)//i am superpop...
						{
							if (iPatternLevel == 0)
							{
								v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vBrutalValuesGenerationZeroLvl(pcInd, iPatternLevel);
							}//if (i_pattern_level == 0)
							else
							{
								v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vBrutalValuesGenerationHighLvl(pcInd, iPatternLevel, &(v_gene_patterns_trees.at(iPatternLevel - 1)) );
							}//else  if (i_pattern_level == 0)

							for (int i_3lo_pop = 0; i_3lo_pop < pc_parent->v_3lo_pops.size(); i_3lo_pop++)
							{
								if (pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.size() > 0)
								{
									v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
									(
										&(pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0)), iPatternLevel + 1
									);

									v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
									(
										&(pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels_start_individuals.at(0)), iPatternLevel + 1
									);
								}//if  (pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.size() > 0)
							}//for (int ii = 0; pc_parent->v_3lo_pops.size(); ii++)*/


							if (pc_parent->pc_3lo_single_tool->v_population_levels.size() > 0)
							{
								v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
								(
									&(pc_parent->pc_3lo_single_tool->v_population_levels.at(0)), iPatternLevel + 1
								);

								v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
								(
									&(pc_parent->pc_3lo_single_tool->v_population_levels_start_individuals.at(0)), iPatternLevel + 1
								);
							}//if (pc_parent->pc_3lo_single_tool->v_population_levels.size() > 0)
						}//if (pc_parent->pc_3lo_single_tool == this)//i am superpop...


						v_filter_gene_freq_by_ind(2, pcInd, iPatternLevel, v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene), pi_genotype_buffer);

						/*pv_gene_freqs = v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->pvGetHighLevelGeneFreq();
						i_old_level = pcInd->i_level;
						pcInd->i_level = i_pattern_level + 1;
						pcInd->dComputeFitness(i_pattern_level + 1);
						pcInd->
							vGetBestGenes
							(
								2, &v_chosen_genes, pv_gene_freqs,
								pcInd->vpi_optimized_genotypes.at(i_pattern_level + 1)[0], pi_genotype_buffer,
								i_pattern_level + 2
							);
						pcInd->i_level = i_old_level;

						//preserve chosen genes
						for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)
						{
							for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
							{
								if (v_chosen_genes.at(i_chosen) == pv_gene_freqs->at(i_gene_freq))
								{
									pv_gene_freqs->erase(pv_gene_freqs->begin() + i_gene_freq);
									i_gene_freq--;
								}//if (v_chosen_genes.at(i_chosen) == pv_gene_freqs->at(i_gene_freq))
							}//for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
						}//for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)

						for (int i_gene_freq = 0; i_gene_freq < pv_gene_freqs->size(); i_gene_freq++)
							delete  pv_gene_freqs->at(i_gene_freq);
						pv_gene_freqs->clear();

						for (int i_chosen = 0; i_chosen < v_chosen_genes.size(); i_chosen++)
							pv_gene_freqs->push_back(v_chosen_genes.at(i_chosen));

						v_chosen_genes.clear();//*/


						 /*else
						 {
						 if (v_population_levels.size() > 0)
						 {
						 v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
						 (
						 &(v_population_levels.at(0)), i_pattern_level + 1
						 );
						 }//if (v_population_levels.size() > 0)
						 }//else  if (pc_parent->pc_3lo_single_tool == this)//i am superpop...*/
					}//if (pc_parent != NULL)


				}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
			}//if  (v_population_levels.size() > iPatternLevel)
		}//if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
	}//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)

	delete  pi_genotype_buffer;


	d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	pc_parent->c_time_counter.bGetTimePassed(&d_time_after);

	CString  s_buf;
	s_buf.Format("COUNT GENE FREQ COST (%.2lf sec  %.0lf ffe)\n", d_time_after - d_time_before, d_ffe_after - d_ffe_before);
	pc_parent->pc_log->vPrintLine(s_buf, true);
	//Tools::vShow(0);

	return;


	for (int i_lvl = 0; i_lvl < v_higher_level_trees_genes.size(); i_lvl++)
	{
		//first compute individual fitness at particular level
		for (int i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)
		{
			v_population_levels.at(0).at(i_ind)->dComputeFitness(i_lvl);
		}//for  (int  i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)


		for (int i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
		{
			v_higher_level_trees_genes.at(i_lvl).at(i_hl_gene)->vCountHighLvlGeneFreq(&(v_population_levels.at(0)), i_lvl);
		}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
	}//for  (int  i_lvl = 0; i_lvl < v_higher_level_trees_genes.size(); i_lvl++)
}//void  C3LOSingle::v_count_higher_genes_freq_for_ind(C3LOIndividual  *pcInd)


void  C3LOSingle::v_count_higher_genes_freq()
{
	return;
	double  d_ffe_before, d_ffe_after;
	double  d_time_before, d_time_after;

	d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	c_time_counter.bGetTimePassed(&d_time_before);
	
	
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		if (v_gene_patterns_trees.at(i_pattern_level).size() > 1)
		{
			//if (v_population_levels.size() > iPatternLevel)
			{
				for (int i_hl_gene = 0; i_hl_gene < v_gene_patterns_trees.at(i_pattern_level).size(); i_hl_gene++)
				{
					v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vZeroHighLvlGeneFreq();

					if (pc_parent != NULL)
					{
						//if (pc_parent->pc_3lo_single_tool == this)//i am superpop...
						{
							for (int i_3lo_pop = 0; i_3lo_pop < pc_parent->v_3lo_pops.size(); i_3lo_pop++)
							{
								if (pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.size() > 0)
								{
									v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
									(
										&(pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0)), i_pattern_level + 1
									);

									v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
									(
										&(pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels_start_individuals.at(0)), i_pattern_level + 1
									);
								}//if  (pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.size() > 0)
							}//for (int ii = 0; pc_parent->v_3lo_pops.size(); ii++)*/


							if (pc_parent->pc_3lo_single_tool->v_population_levels.size() > 0)
							{
								v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
								(
									&(pc_parent->pc_3lo_single_tool->v_population_levels.at(0)), i_pattern_level + 1
								);

								v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
								(
									&(pc_parent->pc_3lo_single_tool->v_population_levels_start_individuals.at(0)), i_pattern_level + 1
								);
							}//if (pc_parent->pc_3lo_single_tool->v_population_levels.size() > 0)
						}//if (pc_parent->pc_3lo_single_tool == this)//i am superpop...
						/*else
						{
							if (v_population_levels.size() > 0)
							{
								v_gene_patterns_trees.at(i_pattern_level).at(i_hl_gene)->vCountHighLvlGeneFreq
								(
									&(v_population_levels.at(0)), i_pattern_level + 1
								);
							}//if (v_population_levels.size() > 0)
						}//else  if (pc_parent->pc_3lo_single_tool == this)//i am superpop...*/
					}//if (pc_parent != NULL)

					
				}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
			}//if  (v_population_levels.size() > iPatternLevel)
		}//if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
	}//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)


	 d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	 pc_parent->c_time_counter.bGetTimePassed(&d_time_after);

	 CString  s_buf;
	 s_buf.Format("COUNT GENE FREQ COSTa (%.2lf sec  %.0lf ffe)\n", d_time_after - d_time_before, d_ffe_after - d_ffe_before);
	 pc_parent->pc_log->vPrintLine(s_buf, true);
	 //Tools::vShow(0);

	return;


	for (int i_lvl = 0; i_lvl < v_higher_level_trees_genes.size(); i_lvl++)
	{
		//first compute individual fitness at particular level
		for (int i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)
		{
			v_population_levels.at(0).at(i_ind)->dComputeFitness(i_lvl);
		}//for  (int  i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)


		for (int i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
		{
			v_higher_level_trees_genes.at(i_lvl).at(i_hl_gene)->vCountHighLvlGeneFreq(&(v_population_levels.at(0)), i_lvl);
		}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
	}//for  (int  i_lvl = 0; i_lvl < v_higher_level_trees_genes.size(); i_lvl++)
}//void  C3LO::v_count_higher_genes_freq()




void C3LOSingle::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);

	t_start = tStartTime;

	

	CError  c_err(iERROR_PARENT_C3LOSingleOptimizer);
	CString   s_buf;
	CString   s_levels_info;
	CString   s_tribe_eff_info;


	
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();
	v_initialize_orders();//i_templ_length needed


	if (i_templ_length <= 0)
	{
		c_err.vSetError(C3LOSingle::iERROR_CODE_3LO_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)


	pi_best_genotype = new int32_t[i_templ_length];
	pi_genotype_tool = new int32_t[i_templ_length];

	pi_last_start_genotype = new int[i_templ_length];
	b_use_last_start_genotype = false;

	
	pc_log->vPrintLine("Initializing...", true);


	double  d_fit_max_cur;
	C3LOIndividual  *pc_best;


	double  d_buf;
	int  i_max_lvl;



	d_fit_max_cur = 0;
	pc_best = NULL;


	v_insert_gene_pattern_part_level(0);
	v_add_gene_pattern_tree(NULL, 0);

	s_debug = "";

	i_pop_to_add = 1;
	i_upgraded_new_ind_at_last_insert = 0;
	i_upgraded_new_ind = 0;
	i_worse_but_different_used = 0;
	i_worse_but_different_used_at_last_insert = 0;
	b_add_higher_level_gene(NULL, 0);
	i_crossings_with_effect_in_the_pop = 0;
	i_all_crossings = 0;
	i_all_effective_crossings = 0;

	i_better_but_cannot_be_added_at_last_insert = 0;
	i_better_but_cannot_be_added = 0;
	i_tribe_effective_in_this_iteration = 0;

	i_brutal_random_search = 0;
	i_brutal_random_search_effective = 0;

	d_fit_max_cur = 0;
	d_unit_max_cur = 0;
	i_cur_gen = 0;

	b_tribe_effective = false;


	c_time_counter.vSetStartNow();
	pc_log->vPrintLine("DONE...", true);

	i_dsm_reporting_enumeration_tool = 0;
	i_last_m_from_dsm = 1;


	pd_dsm_matrix = new double*[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pd_dsm_matrix[ii] = new double[i_templ_length];


	pi_gene_realtions = new int*[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_realtions[ii] = new int[i_templ_length];

	i_relation_length_max = 0;

	d_best_fitness_prev = 0;
	d_best_fitness_cur = 0;
	pc_best_ind = NULL;
	pc_last_ind = NULL;

	i_linkage_generations = 0;
	d_linkage_ffe = 0;
	d_linkage_time = 0;

	b_linkage_generated = false;

	if (pc_parent->b_use_dsm == true)
	{
		pppi_genes_values_occurrences = new int**[i_templ_length];

		for (int ii = 0; ii < i_templ_length; ii++)
		{
			pppi_genes_values_occurrences[ii] = new int*[i_templ_length];

			for (int jj = 0; jj < i_templ_length; jj++)
			{
				pppi_genes_values_occurrences[ii][jj] = new int[4];

				for (int kk = 0; kk < 4; kk++)
				{
					pppi_genes_values_occurrences[ii][jj][kk] = 0;
				}//for (int kk = 0; kk < 4; kk++)
			}//for (int jj = 0; jj < i_templ_length; jj++)
		}//for (int ii = 0; ii < i_templ_length; ii++)
	}//if (pc_parent->b_use_dsm == true)
}//void CDummyOptimizer::vInitialize(time_t tStartTime)




C3LOIndividual* C3LOSingle::pc_get_random_individual(int iLevel)
{
	C3LOIndividual  *pc_indiv_new;

	pc_indiv_new = new C3LOIndividual(iLevel, i_templ_length, pc_problem, &v_optimization_orders, this);
	pc_indiv_new->vRandomInit();
	
	if (b_use_last_start_genotype == true)
	{
		pc_indiv_new->vSetGenotype(pi_last_start_genotype);
		pc_indiv_new->v_optimization_order.at(0) = v_last_optimization_order;
	}//if (b_use_last_start_genotype == true)

	//store last start genotype
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_last_start_genotype[ii] = pc_indiv_new->pi_genotype[ii];

	//and opitimization order
	v_last_optimization_order = pc_indiv_new->v_optimization_order.at(0);

	pc_indiv_new->dComputeFitness(-1);//we want the genotype to be optimized//we do not want it any more...

	i_random_indiv_adds++;

	return(pc_indiv_new);
}//C3LOIndividual  *pc_get_random_individual()



C3LOIndividual *C3LOSingle::pc_get_individual_with_lvl_tournament(C3LOIndividual *pcDifferentToIndiv /*= NULL*/)
{
	C3LOIndividual *pc_parent_0, *pc_parent_1;

	pc_parent_0 = pc_get_individual_from_lvl(0, pcDifferentToIndiv);
	pc_parent_1 = pc_get_individual_from_lvl(0, pcDifferentToIndiv);

	if (pc_parent_0 == NULL)  return(pc_parent_1);
	if (pc_parent_1 == NULL)  return(pc_parent_0);

	if (pc_parent_0->i_level > pc_parent_1->i_level)
		return(pc_parent_0);
	else
		return(pc_parent_1);

	return(NULL);
}//C3LOIndividual *C3LO::pc_get_individual_from_lvl(int  iEvolutionLvl, C3LOIndividual *pcDifferentToIndiv /*= NULL*/)


C3LOIndividual *C3LOSingle::pc_get_individual_from_lvl(int  iEvolutionLvl, C3LOIndividual *pcDifferentToIndiv /*= NULL*/)
{
	if (iEvolutionLvl < 0)  return(NULL);
	if (v_population_levels.size() <= iEvolutionLvl)  return(NULL);

	if (v_population_levels.at(iEvolutionLvl).size() <  1)  return(NULL);

	int  i_rand_ind;

	int  i_rand_rounds = 0;
	bool  b_rand_again = true;
	while (b_rand_again == true)
	{
		//i_rand_ind = pc_random_gen->Next(0, v_population_levels.at(iEvolutionLvl).size());
		i_rand_ind = RandUtils::iRandNumber(0, v_population_levels.at(iEvolutionLvl).size() -1);
		if (v_population_levels.at(iEvolutionLvl).at(i_rand_ind) != pcDifferentToIndiv)  return(v_population_levels.at(iEvolutionLvl).at(i_rand_ind));

		i_rand_rounds++;

		if (i_rand_rounds >= _3LO_TRIBE_MAX_ROUNDS_NUM)  b_rand_again = false;
	}//while  (b_rand_again == true)


	return(NULL);
}//C3LOIndividual *C3LO::pc_get_individual_from_lvl(int  iEvolutionLvl, C3LOIndividual *pcDifferentToIndiv /*= NULL*/)





C3LOIndividual* C3LOSingle::pc_insert_into_pop(vector  <C3LOIndividual  *>  *pvPopDest, C3LOIndividual  *pcIndToInsert, bool bCheckPhenotype /*= true*/, bool  bCheckGenotype /*= false*/)
{
	double  d_fit;
	d_fit = pcIndToInsert->dComputeFitness(0);//also forces the creation of the genotype


	if (pvPopDest != pv_population)
	{
		for (int ii = 0; ii < pv_population->size(); ii++)
		{
			if (pv_population->at(ii)->dComputeFitness(0) == pcIndToInsert->dComputeFitness(0))
			{
				if (pv_population->at(ii)->bComparePhenotype(0, 0, pcIndToInsert) == true)
				{
					return(pv_population->at(ii));
				}//if  (pv_population->at(ii)->bComparePhenotype(pcIndToInsert) ==  true)
			}//if  (pv_population->at(ii)->dComputeFitness() == pcIndToInsert->dComputeFitness())
		}//for  (int  ii = 0; ii < pv_population->size(); ii++)
	}//if  (pvPopDest != pv_population)


	for (int ii = 0; ii < pvPopDest->size(); ii++)
	{
		if (pcIndToInsert == pvPopDest->at(ii))  return(pvPopDest->at(ii));
	}//for  (int ii = 0;  ii < pvPopDest->size(); ii++)


	 //find individuals of the same fitness
	for (int ii = pvPopDest->size() - 1; ii >= 0; ii--)
	{
		if (pvPopDest->at(ii)->dComputeFitness(0) > pcIndToInsert->dComputeFitness(0))
		{
			if (ii  <  pvPopDest->size() - 1)
			{
				pvPopDest->insert(pvPopDest->begin() + ii + 1, pcIndToInsert);
				return(pcIndToInsert);
			}//if  (ii  <  pv_population->size() - 1)
			else
			{
				pvPopDest->push_back(pcIndToInsert);
				return(pcIndToInsert);
			}//else  if  (ii  <  pv_population->size() - 1)
		}//if  (pv_population->at(ii)->dComputeFitness() > pc_indiv->dComputeFitness())

		if (pvPopDest->at(ii)->dComputeFitness(0) == pcIndToInsert->dComputeFitness(0))
		{
			if (bCheckPhenotype == true)
			{
				if (pvPopDest->at(ii)->bComparePhenotype(0,0, pcIndToInsert) == true)
				{
					return(pvPopDest->at(ii));
				}//if  (pv_population->at(ii)->bComparePhenotype(pcIndToInsert) ==  true)

				if (bCheckGenotype == true)
				{
					if (pvPopDest->at(ii)->bCompareGenotype(pcIndToInsert) == true)
					{
						return(pvPopDest->at(ii));
					}//if  (pv_population->at(ii)->bComparePhenotype(pcIndToInsert) ==  true)
				}//if  (bCheckGenotype == true)
			}//if  (bCheckPhenotype == true)
		}//if  (pv_population->at(ii)->dComputeFitness() == pc_indiv->dComputeFitness())
	}//for  (int  ii = 0; ii < pv_population->size(); ii++)



	 //new best individual...
	pvPopDest->insert(pvPopDest->begin(), pcIndToInsert);
	return(pcIndToInsert);
}//C3LOIndividual* C3LO::pc_insert_into_pop(vector  <C3LOIndividual  *>  *pvPopDest, C3LOIndividual  *pcIndToInsert, bool bCheckPhenotype /*= true*/)




bool  C3LOSingle::b_add_start_individual_copy(C3LOIndividual *pcInd, int  iLevel)
{
	if (iLevel < 0)  return(false);

	while (v_population_levels_start_individuals.size() <= iLevel)
		v_population_levels_start_individuals.push_back(vector<C3LOIndividual *>());


	for (int ii = 0; ii < v_population_levels_start_individuals.at(iLevel).size(); ii++)
	{
		if (iLevel == 0)
		{
			if (v_population_levels_start_individuals.at(iLevel).at(ii)->bCompareGenotype(pcInd) == true)  return(false);
		}//if  (bComparePhenotypes == true)

		if (iLevel > 0)
		{
			if (v_population_levels_start_individuals.at(iLevel).at(ii)->bComparePhenotype(iLevel + 1, 0, pcInd) == true)  return(false);
		}//if  (bComparePhenotypes == true)	
	}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)

	C3LOIndividual  *pc_ind_copy;
	pc_ind_copy = new C3LOIndividual(iLevel, i_templ_length, pc_problem, &v_optimization_orders, this);
	pc_ind_copy->vRandomInit();
	pc_ind_copy->vSetGenotype(pcInd->piGetPhenotype(iLevel, 0));
	pc_ind_copy->dComputeFitness(-1);
	
	v_population_levels_start_individuals.at(iLevel).push_back(pc_ind_copy);
	

	return(true);
}//bool  C3LO::b_add_at_level(C3LOIndividual *pcInd, int  iLevel, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)



bool  C3LOSingle::b_add_at_level(C3LOIndividual *pcInd, int  iLevel)
{
	if (iLevel < 0)  return(false);

	while (v_population_levels.size() <= iLevel)
		v_population_levels.push_back(vector<C3LOIndividual *>());


	for (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
	{
		if (iLevel == 0)
		{
			if (v_population_levels.at(iLevel).at(ii)->bCompareGenotype(pcInd) == true)  return(false);
		}//if  (bComparePhenotypes == true)

		if (iLevel > 0)
		{
			if (v_population_levels.at(iLevel).at(ii)->bComparePhenotype(iLevel+1, 0, pcInd) == true)  return(false);			
		}//if  (bComparePhenotypes == true)	
	}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)

	pcInd->i_level = iLevel;
	v_population_levels.at(iLevel).push_back(pcInd);

	return(true);
}//bool  C3LO::b_add_at_level(C3LOIndividual *pcInd, int  iLevel, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)


/*
bool  C3LO::b_ind_on_lvl(C3LOIndividual *pcInd, int  iEvolutionLvl)
{
	if (iEvolutionLvl  <  0)  return(false);
	if (iEvolutionLvl >= v_population_levels.size())  return(false);

	for (int i_ind = 0; i_ind < v_population_levels.at(iEvolutionLvl).size(); i_ind++)
	{
		if (v_population_levels.at(iEvolutionLvl).at(i_ind) == pcInd)  return(true);
	}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)

	return(false);
}//bool  C3LO::b_ind_on_lvl(C3LOIndividual *pcInd, , int  iEvolutionLvl)*/



bool  C3LOSingle::b_can_ind_be_added_at_any_level(int  *piGeneString, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)
{
	for (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)
	{
		if (b_can_ind_be_added_at_pattern_level(i_lvl, piGeneString, bComparePhenotypes, bCompareGenotypes) == false)  return(false);
	}//for  (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)

	return(true);
}//bool  C3LO::b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)



bool  C3LOSingle::b_can_ind_be_added_at_pattern_level(int iPatternLevel, int  *piGeneString, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)
{
	if (iPatternLevel >= v_population_levels.size())  return(true);


	for (int i_ind = 0; i_ind < v_population_levels.at(iPatternLevel).size(); i_ind++)
	{
		if (bComparePhenotypes == true)
		{
			if (v_population_levels.at(iPatternLevel).at(i_ind)->bComparePhenotype(0, 0, piGeneString) == true)  return(false);
		}//if  (bComparePhenotypes == true)

		if (bCompareGenotypes == true)
		{
			if (v_population_levels.at(iPatternLevel).at(i_ind)->bCompareGenotype(piGeneString) == true)  return(false);
		}//if  (bComparePhenotypes == true)	
	}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)


	if (pv_other_pops != NULL)
	{
		for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
		{
			if (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPatternLevel)
			{
				for (int i_ind = 0; i_ind < pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).size(); i_ind++)
				{
					if (bComparePhenotypes == true)
					{
						if (pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).at(i_ind)->bComparePhenotype(0, 0, piGeneString) == true)  return(false);
					}//if  (bComparePhenotypes == true)

					if (bCompareGenotypes == true)
					{
						if (pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).at(i_ind)->bCompareGenotype(piGeneString) == true)  return(false);
					}//if  (bComparePhenotypes == true)	
				}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
			}//if  (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPatternLevel)
		}//for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
	}//if (pv_other_pops != NULL)


	return(true);
}//bool  C3LO::b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)



bool  C3LOSingle::b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)
{
	for (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)
	{
		for (int i_ind = 0; i_ind < v_population_levels.at(i_lvl).size(); i_ind++)
		{
			if (bComparePhenotypes == true)
			{
				if (v_population_levels.at(i_lvl).at(i_ind)->bComparePhenotype(0,0, pcInd) == true)  return(false);
			}//if  (bComparePhenotypes == true)

			if (bCompareGenotypes == true)
			{
				if (v_population_levels.at(i_lvl).at(i_ind)->bCompareGenotype(pcInd) == true)  return(false);
			}//if  (bComparePhenotypes == true)	
		}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
	}//for  (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)

	return(true);
}//bool  C3LO::b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)



double  C3LOSingle::d_get_covering_at_level(int  iLevel)
{
	if (iLevel >= v_genes_marked_by_linkage.size())  return(0);


	double  d_pat_coverage;
	d_pat_coverage = 0;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (v_genes_marked_by_linkage.at(iLevel)[ii] > 0)  d_pat_coverage++;
	}//for  (int  ii = 0; ii < i_templ_length; ii++)
	d_pat_coverage = d_pat_coverage / i_templ_length;

	return(d_pat_coverage);
}//double  C3LO::d_get_covering_at_level(int  iLevel)


void  C3LOSingle::v_add_linkage_covering_at_level(int  iLevel)
{
	while (v_genes_marked_by_linkage.size() <= iLevel)
	{
		int  *pi_linkage_covering;
		pi_linkage_covering = new int[i_templ_length];

		for (int ii = 0; ii < i_templ_length; ii++)
			pi_linkage_covering[ii] = 0;

		v_genes_marked_by_linkage.push_back(pi_linkage_covering);
	}//while  (v_genes_marked_by_linkage.size() <= iLevel)
}//void  C3LO::v_add_linkage_covering_at_level(int  iLevel)



bool  C3LOSingle::b_produce_new_linkage(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, int  iPatternLevel)
{
	int  *pi_phenotype_main, *pi_phenotype_donator;

	//pcInd0->dComputeFitness(iPatternLevel); //force phenotype creation
	//pcInd1->dComputeFitness(iPatternLevel); //force phenotype creation

	//pi_phenotype_main = pcInd0->piGetPhenotype(iPatternLevel, 0);
	//pi_phenotype_donator = pcInd1->piGetPhenotype(iPatternLevel, 0);

	if (pcInd0->i_level == 0)  return(false);
	if (pcInd1->i_level == 0)  return(false);

	pi_phenotype_main = pcInd0->vpi_optimized_genotypes_for_linkage.at(iPatternLevel)[0];
	pi_phenotype_donator = pcInd1->vpi_optimized_genotypes_for_linkage.at(iPatternLevel)[0];

	

	v_add_linkage_covering_at_level(iPatternLevel);

	for (int i_gen_off = 0; i_gen_off < i_templ_length; i_gen_off++)
	{
		if (pi_phenotype_main[i_gen_off] != pi_phenotype_donator[i_gen_off])
		{
			if (v_genes_marked_by_linkage.at(iPatternLevel)[i_gen_off] == 0)  return(true);
			//if  (pi_genes_marked_by_linkage[i_gen_off] == 0)  return(true);
		}//if  (pi_phenotype_0[i_gen_off]  ==  pi_phenotype_1[i_gen_off])
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(false);
}//bool  C3LO::b_produce_new_linkage_high_level(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, int  iLevel)





bool  C3LOSingle::bCanIndBeAddedAtAnyLevel(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)
{
	for (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)
	{
		for (int i_ind = 0; i_ind < v_population_levels.at(i_lvl).size(); i_ind++)
		{
			if (bComparePhenotypes == true)
			{
				if (v_population_levels.at(i_lvl).at(i_ind)->bComparePhenotype(0,0, pcInd) == true)  return(false);
			}//if  (bComparePhenotypes == true)

			if (bCompareGenotypes == true)
			{
				if (v_population_levels.at(i_lvl).at(i_ind)->bCompareGenotype(pcInd) == true)  return(false);
			}//if  (bComparePhenotypes == true)	
		}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
	}//for  (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)

	return(true);
}//bool  C3LO::b_can_ind_be_added_at_any_level(C3LOIndividual *pcInd, bool  bComparePhenotypes /*= true*/, bool  bCompareGenotypes /*= false*/)







bool bIndividualOrder(C3LOIndividual *elem1, C3LOIndividual *elem2)
{
	return elem1->dComputeFitness(-1) > elem2->dComputeFitness(-1);
}//bool bIndividualOrder(C3LOIndividual *elem1, C3LOIndividual *elem2)

void  C3LOSingle::vSortAllSubpopsByFitness()
{
	CString  s_buf;

	for (int i_level = 0; i_level < v_population_levels.size(); i_level++)
	{
		std::sort(v_population_levels.at(i_level).begin(), v_population_levels.at(i_level).end(), bIndividualOrder);

		/*for (int ii = 0; ii < v_population_levels.at(i_level).size(); ii++)
		{
			s_buf.Format("%.4lf", v_population_levels.at(i_level).at(ii)->dComputeFitness(-1));
			MessageBox(NULL, s_buf, s_buf, MB_OK);
		}//for (int  ii = 0; ii < v_population_levels.at(i_level).size(); ii++)*/

	}//for (int i_level = 0; i_level < v_population_levels.size() i_level++)
}//void  C3LOSingle::vSortAllSubpopsByFitness()



bool  C3LOSingle::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	double  d_time_passed;

	CString  s_buf;
	CString  s_levels_info;

	b_linkage_generated = false;

	i_the_same = 0;
	i_improved = 0;
	i_improved_other = 0;
	i_shuffled_individuals = 0;


	pc_log->vPrintLine("", true);
	

	if (pc_parent->b_use_dsm == false)
	{
		//while (d_get_covering_at_level(0) < 1)
		//while (i_last_m_from_dsm > 0)
		if (pv_other_pops == NULL)//we do not detect linkage for the super pop
		{
			if (d_get_covering_at_level(0) == 0)
			{
				C3LOIndividual  *pc_ind_for_linkage;

				v_add_linkage_covering_at_level(0);

				pc_ind_for_linkage = pc_get_random_individual(0);
				pc_ind_for_linkage->i_level = 1;

				v_create_trees_from_ind(pc_ind_for_linkage, &(v_gene_patterns_trees.at(0)), 0);

				delete pc_ind_for_linkage;

				s_buf.Format("(%d:%.4lf) ", 0, d_get_covering_at_level(0));
				pc_log->vPrintLine(s_buf, true);
			}//while (d_get_covering_at_level(0))
		}//if  (pv_other_pops  ==  NULL)//we do not detect linkage for the super pop
	}//if (pc_parent->b_use_dsm == false)

	std::sort(v_gene_patterns_parts_interceptions_and_differencies.at(0).begin(), v_gene_patterns_parts_interceptions_and_differencies.at(0).end(), bPatternLarger);
	//s_buf.Format("first pat size: %d", v_gene_patterns_parts_interceptions_and_differencies.at(0).at(0)->v_pattern.size());
	//MessageBox(NULL, s_buf, s_buf, MB_OK);


	
	/*FILE  *pf_test;
	s_debug.Format("%d_", i_cur_gen);
	s_buf.Format("zzz_pats_super_pop.txt", i_cur_gen);
	v_save_trees_and_patterns(s_buf);//*/
	


	C3LOIndividual  *pc_indiv_new;

	pc_indiv_new = pc_get_random_individual(0);
	v_update_higher_genes_freq_by_ind_for_best(pc_indiv_new);
	
	if (pc_indiv_new == NULL)  return(false);
	//optimize new individual
	//pc_indiv_new->i_level = 1;
	//v_count_higher_genes_freq();
	vShufflePatternsAndFreqs();//so the individuals is optimized with different orders...
	//vShufflePatternsToCross();//so the individuals is processed with different orders...
	
	pc_indiv_new->i_level = v_gene_patterns_trees.size() + 1;
	pc_indiv_new->i_level = 1;
	pc_indiv_new->dComputeFitness(-1);


	/*v_update_higher_genes_freq_by_ind(pc_indiv_new);

	for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size() + 2; i_pattern_lvl++)
	{
		if (b_can_ind_be_added_at_pattern_level(0, pc_indiv_new->piGetPhenotype(i_pattern_lvl, 0), true, true) == true)
		{
			pc_indiv_new->vSetGenotype(pc_indiv_new->piGetPhenotype(i_pattern_lvl, 0));
			pc_indiv_new->vSetOptimized(false);
		}//if (b_can_ind_be_added_at_pattern_level(0, pc_ind_upgraded->piGetPhenotype(1, 0), true, true) == true)
	}//for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size() + 2; i_pattern_lvl++)*/
	
	if (b_can_ind_be_added_at_pattern_level(0, pc_indiv_new->piGetPhenotype(1, 0), true, true) == true)
		pc_indiv_new->vSetGenotype(pc_indiv_new->piGetPhenotype(1,0));

	pc_indiv_new->vSetOptimized(false);
	pc_indiv_new->i_level = 0;//*/
	pc_indiv_new->dComputeFitness(-1);
	
	b_add_start_individual_copy(pc_indiv_new, 0);

	
	/*if (pv_other_pops != NULL)
	{
		FILE  *pf_test;
		pf_test = fopen("zz_ind_test.txt", "w+");
		pc_indiv_new->vSave(pf_test, NULL);
		fclose(pf_test);
	}//if  (pv_other_pops != NULL)*/
	

	if (b_add_at_level(pc_indiv_new, 0) == false) return(false);
	

	v_process(pc_indiv_new, 0);
	
	if (pc_parent->b_use_dsm == true)
	{
		v_update_genes_values_occurrences(pc_indiv_new);

		for (int i_tre_lev = 0; i_tre_lev < v_gene_patterns_trees.size(); i_tre_lev++)
		{
			for (int i_tree = 0; i_tree < v_gene_patterns_trees.at(i_tre_lev).size(); i_tree++)
				delete v_gene_patterns_trees.at(i_tre_lev).at(i_tree);
			v_gene_patterns_trees.at(i_tre_lev).clear();

			v_gene_patterns_parts_interceptions_and_differencies.at(i_tre_lev).clear();
		}//for (int i_tre_lev = 0; i_tre_lev < v_gene_patterns_trees.size(); i_tre_lev++)

		v_create_trees_from_dsm(&(v_gene_patterns_trees.at(0)));
		pc_parent->v_join_children_trees_if_necessary(&(v_gene_patterns_trees.at(0)), 0);
	}//if (pc_parent->b_use_dsm == true)

	C3LOIndividual  *pc_ind_upgraded;

	pc_ind_upgraded = new C3LOIndividual(v_gene_patterns_trees.size() + 1, i_templ_length, pc_problem, &v_optimization_orders, this);
	//pc_ind_upgraded = new C3LOIndividual(1, i_templ_length, pc_problem, &v_optimization_orders, this);
	pc_ind_upgraded->vRandomInit();
	pc_ind_upgraded->vSetGenotype(pc_indiv_new->piGetPhenotype(0, 0));
	pc_ind_upgraded->vSetOptimized(false);
	pc_ind_upgraded->dComputeFitness(-1);
	if (b_can_ind_be_added_at_pattern_level(0, pc_ind_upgraded->piGetPhenotype(1, 0), true, true) == true)
	{
		pc_indiv_new->vSetGenotype(pc_ind_upgraded->piGetPhenotype(1, 0));
		pc_indiv_new->vSetOptimized(false);
		pc_indiv_new->dComputeFitness(-1);
	}//if (b_can_ind_be_added_at_pattern_level(0, pc_ind_upgraded->piGetPhenotype(1, 0), true, true) == true)*/

	/*for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size() + 2; i_pattern_lvl++)
	{
		if (b_can_ind_be_added_at_pattern_level(0, pc_ind_upgraded->piGetPhenotype(i_pattern_lvl, 0), true, true) == true)
		{
			pc_indiv_new->vSetGenotype(pc_ind_upgraded->piGetPhenotype(i_pattern_lvl, 0));
			pc_indiv_new->vSetOptimized(false);
		}//if (b_can_ind_be_added_at_pattern_level(0, pc_ind_upgraded->piGetPhenotype(1, 0), true, true) == true)
	}//for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size() + 2; i_pattern_lvl++)
	pc_indiv_new->dComputeFitness(-1);//*/


	v_update_higher_genes_freq_by_ind_for_best(pc_indiv_new);
	//v_update_higher_genes_freq_by_ind(pc_indiv_new);

	pc_last_ind = pc_indiv_new;

	delete  pc_ind_upgraded;

	

	double  d_best_fitness_0;
	double  d_best_fitness, d_fit_buf;
	int  *pi_curr_phenotype, *pi_best_phenotype;
	C3LOIndividual  *pc_best, *pc_best_0;

	//int  i_best_opt_distance, i_dist_buf;
	//double  d_best_opt_distance_fit;

	int  i_old_level;
	pi_best_phenotype = NULL;
	d_best_fitness = 0;
	d_best_fitness_0 = 0;
	pc_best = NULL;
	pc_best_0 = NULL;
	//i_best_opt_distance = i_templ_length;
	for (int i_pop_lvl = 0; i_pop_lvl < v_population_levels.size(); i_pop_lvl++)
	{
		for (int i_ind = 0; i_ind < v_population_levels.at(i_pop_lvl).size(); i_ind++)
		{
			//int i_pattern_lvl = -1;
			for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size(); i_pattern_lvl++)
			{
				i_old_level = v_population_levels.at(i_pop_lvl).at(i_ind)->i_level;
				v_population_levels.at(i_pop_lvl).at(i_ind)->i_level = i_pattern_lvl + 2;
				d_fit_buf = v_population_levels.at(i_pop_lvl).at(i_ind)->dComputeFitness(i_pattern_lvl + 2);
				v_population_levels.at(i_pop_lvl).at(i_ind)->i_level = i_old_level;

				pi_curr_phenotype = v_population_levels.at(i_pop_lvl).at(i_ind)->piGetPhenotype(i_pattern_lvl + 2, 0);
				//d_fit_buf = dComputeFitness(pi_curr_phenotype);
				if (d_best_fitness < d_fit_buf)
				{
					d_best_fitness = d_fit_buf;
					pi_best_phenotype = pi_curr_phenotype;

					pc_best = v_population_levels.at(i_pop_lvl).at(i_ind);
				}//if (d_best_fitness < d_fit_buf)

				if (i_pattern_lvl == 0)
				{

					if (d_best_fitness_0 < d_fit_buf)
					{
						d_best_fitness_0 = d_fit_buf;
						pc_best_0 = v_population_levels.at(i_pop_lvl).at(i_ind);
					}//if (d_best_fitness < d_fit_buf)
				}//if (i_pop_lvl == 0)


			}//for (int i_pattern_lvl = 0; v_gene_patterns_trees.size(); i_pattern_lvl++)		

		}//for  (int  i_ind = 0; i_ind < v_population_levels.at(i_pop_lvl).size(); i_pop_lvl++)

	}//for (int i_pop_lvl = 0; i_pop_lvl < v_population_levels.size(); i_pop_lvl++)

	//pc_best_ind_0 = pc_best_0;
	pc_best_ind = pc_best;

	d_best_fitness_prev = d_best_fitness_cur;
	d_best_fitness_cur = d_best_fitness_0;


	if (d_best_fitness_prev < d_best_fitness_cur)
	{
		s_buf.Format("best fitness increase: %.4lf -> %.4lf", d_best_fitness_prev,  d_best_fitness_cur);
		pc_log->vPrintLine(s_buf, true);

		if (pv_other_pops == NULL)//we do not detect linkage for the super pop
		{
			//clear the trees first...
			for (int i_tre_lev = 0; i_tre_lev < v_gene_patterns_trees.size(); i_tre_lev++)
			{
				for (int i_tree = 0; i_tree < v_gene_patterns_trees.at(i_tre_lev).size(); i_tree++)
					delete v_gene_patterns_trees.at(i_tre_lev).at(i_tree);
				v_gene_patterns_trees.at(i_tre_lev).clear();

				v_gene_patterns_parts_interceptions_and_differencies.at(i_tre_lev).clear();
			}//for (int i_tre_lev = 0; i_tre_lev < v_gene_patterns_trees.size(); i_tre_lev++)


			C3LOIndividual  *pc_ind_for_linkage;

			v_add_linkage_covering_at_level(0);

			pc_ind_for_linkage = pc_get_random_individual(0);
			pc_ind_for_linkage->i_level = 1;
			pc_ind_for_linkage->vSetGenotype(pc_best_0->piGetPhenotype(0, 0));
			pc_ind_for_linkage->vSetOptimized(false);

			int i_gene_trees_at_previous_level;
			i_gene_trees_at_previous_level = i_templ_length + 1;
			v_create_trees_from_ind(pc_ind_for_linkage, &(v_gene_patterns_trees.at(0)), 0);
			pc_parent->v_join_children_trees_if_necessary(&(v_gene_patterns_trees.at(0)), 0);
			v_count_higher_genes_freq_for_ind(pc_ind_for_linkage, 0);
			

			int  i_current_level = 1;
			
			while ( (i_gene_trees_at_previous_level > v_gene_patterns_trees.at(i_current_level-1).size())&&(v_gene_patterns_trees.at(i_current_level-1).size() > 1) )
			{
				pc_ind_for_linkage->i_level = i_current_level + 1;
				i_gene_trees_at_previous_level = v_gene_patterns_trees.at(i_current_level-1).size();

				v_add_linkage_covering_at_level(i_current_level);
				v_insert_gene_pattern_part_level(i_current_level);
				v_create_trees_from_ind(pc_ind_for_linkage, &(v_gene_patterns_trees.at(i_current_level)), i_current_level);


				pc_parent->v_join_children_trees_if_necessary(&(v_gene_patterns_trees.at(i_current_level)), i_current_level);
				v_count_higher_genes_freq_for_ind(pc_ind_for_linkage, i_current_level);

				i_current_level++;

				//::MessageBox(NULL, "", "", MB_OK);
			}//while ( (i_gene_trees_at_previous_level < v_gene_patterns_trees.at(0).size())&&(v_gene_patterns_trees.at(0).size() > 1) )*/
					
			
			delete pc_ind_for_linkage;

			s_buf.Format("(%d:%.4lf) ", 0, d_get_covering_at_level(0));
			pc_log->vPrintLine(s_buf, true);


			CString  s_pop_id;

			s_pop_id = "";
			if (pc_parent != NULL)
			{
				for (int ii = 0; ii < pc_parent->v_3lo_pops.size(); ii++)
				{
					if (this == pc_parent->v_3lo_pops.at(ii))
					{
						s_pop_id.Format("%.2d", ii);
					}//if (this == pc_parent->v_3lo_pops.at(ii))
				}//for (int ii = 0; ii < pc_parent->v_3lo_pops.size(); ii++)

				if (pc_parent->pc_3lo_single_tool == this)  s_pop_id = "superpop";
			}//if (pc_parent != NULL)


			/*if (s_pop_id != "superpop")
			{
				//FILE  *pf_test;
				s_debug.Format("%d_", i_cur_gen);
				//s_buf.Format("zzz_pats_%s_%.2d.txt", s_pop_id, i_cur_gen);
				s_buf.Format("zzz_pats_%s_%.3d_.txt", s_pop_id, i_cur_gen);
				v_save_trees_and_patterns(s_buf);

				s_buf.Format("zzz_DSM_%s_%.3d_.txt", s_pop_id, i_cur_gen);
				v_report_dsm(pd_dsm_matrix, s_buf);
			}//if  (s_pop_id != "")*/


			/*if (s_pop_id != "")
			{
				FILE  *pf_test;
				s_debug.Format("%d_", i_cur_gen);
				s_buf.Format("zzz_pats_%s_%.2d.txt", s_pop_id, i_cur_gen);
				v_save_trees_and_patterns(s_buf);
			}//if  (s_pop_id != "")*/


			/*C3LOIndividual  *pc_ind_for_test;
			pc_ind_for_test = pc_get_random_individual(0);
			pc_ind_for_test->i_level = 2;
			pc_ind_for_test->dComputeFitness(-1);

			FILE  *pf_test;
			pf_test = fopen("zzz__test.txt", "w+");
			pc_ind_for_test->vSave(pf_test, NULL);
			fclose(pf_test);


			v_add_gene_pattern_tree(NULL, 1);
			v_add_linkage_covering_at_level(1);
			v_create_trees_from_ind(pc_ind_for_test, &(v_gene_patterns_trees.at(1)), 1);

			::MessageBox(NULL, "a", "a", MB_OK);

			pc_log->vPrintLine(s_buf, true);*/
		}//if  (pv_other_pops  ==  NULL)//we do not detect linkage for the super pop
	}//if (d_best_fitness_prev < d_best_fitness_cur)



	if (pi_best_phenotype != NULL)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_best_genotype[ii] = (int32_t)(pi_best_phenotype[ii]);

		b_update_best_individual(iIterationNumber, tStartTime, pi_best_genotype, d_best_fitness);
	}//if (pi_best_phenotype != NULL)



	//pattern coverage
	CString  s_coverage;
	s_coverage = "";

	for (int ii = 0; ii < v_gene_patterns_trees.size(); ii++)
	{
		s_buf.Format("(%d:%.4lf) ", ii, d_get_covering_at_level(ii));
		s_coverage += s_buf;
	}//for  (int  ii = 0; ii < v_gene_patterns_trees.size(); ii++)


	CString  s_block_info;

	CString  s_phenotype_part;
	int  i_ones, i_zeros, i_other;
	vector <CString>  v_other;

	s_block_info = "";

	/*for  (int  i_gene = 0; i_gene < i_templ_length; i_gene+=10)
	{
	i_ones = 0;
	i_zeros = 0;
	i_other = 0;
	v_other.clear();
	for  (int  ii = 0; ii < v_population_levels.at(0).size(); ii++)
	{
	s_phenotype_part = v_population_levels.at(0).at(ii)->sGetPhenotype(0, i_gene, i_gene + 10);

	if  (s_phenotype_part == "1111111111")  i_ones++;
	if  (s_phenotype_part == "0000000000")  i_zeros++;

	if  (
	(s_phenotype_part != "1111111111")&&
	(s_phenotype_part != "0000000000")
	)
	{
	i_other++;

	if (!( std::find(v_other.begin(), v_other.end(), s_phenotype_part) != v_other.end() ))
	v_other.push_back(s_phenotype_part);

	}//if  (

	}//for  (int  ii = 0; ii < v_population_levels.at(0).size(); ii++)

	s_buf.Format("[%d: %d/%d %d/%d] ", i_gene, i_ones, i_zeros, i_other, (int) v_other.size());
	s_block_info += s_buf;
	}//for  (int  i_gene = 0; i_gene < i_templ_length; i_gene++)*/


	CString  s_pop_id;

	s_pop_id = "";
	if (pc_parent != NULL)
	{
		if (pc_parent->pc_3lo_single_tool == this)  s_pop_id = "superpop";
	}//if (pc_parent != NULL)


	/*if (s_pop_id != "superpop")
	{
		//FILE  *pf_test;
		s_debug.Format("%d_", i_cur_gen);
		//s_buf.Format("zzz_pats_%s_%.2d.txt", s_pop_id, i_cur_gen);
		s_buf.Format("zzz_pats_%s_%.3d_.txt", s_pop_id, i_cur_gen);
		v_save_trees_and_patterns(s_buf);

		s_buf.Format("zzz_DSM_%s_%.3d_.txt", s_pop_id, i_cur_gen);
		v_report_dsm(pd_dsm_matrix, s_buf);
	}//if  (s_pop_id != "")*/




	//d_time_passed = (uint32_t)(time(nullptr) - t_start);
	c_time_counter.bGetTimePassed(&d_time_passed);

	int  i_highest_hit = 0;

	for (int ii = 0; ii < v_gene_patterns_parts_interceptions_and_differencies.at(0).size(); ii++)
	{
		if (i_highest_hit < v_gene_patterns_parts_interceptions_and_differencies.at(0).at(ii)->i_number_of_hits)
			i_highest_hit = v_gene_patterns_parts_interceptions_and_differencies.at(0).at(ii)->i_number_of_hits;
	}//for (int ii = 0; ii < v_gene_patterns_parts_interceptions_and_differencies.size(); ii++)


	s_buf.Format
	(
		"blocked improvements: %d %s Iteration: %d (next gen ind: %d new shuffled ind: %d) (upgraded/upgraded other/the same: %d/%d/%d) BestFit: %.4lf (zero: %.4lf) last: %.4lf [u:%.4lf] PopSize:%d  Patterns[coverage: %s] FFE:%.0lf FFE_OPTs: %d time passed:%.4lf  %s   original parts: %d parts: %d   interceptParts:%d hh:%d ",
		i_better_but_cannot_be_added, s_block_info,
		i_cur_gen, v_individuals_to_process_next_gen.size(), i_shuffled_individuals,
		i_improved, i_improved_other, i_the_same,

		d_best_fitness, d_best_fitness_0, pc_last_ind->dComputeFitness(-1),  d_unit_max_cur,
		i_get_ind_number(),

		s_coverage,
		//pc_fitness->dGetFuncEval(),
		(double)pc_problem->pcGetEvaluation()->iGetFFE(),
		pc_genotype_root->iGetSuccOpt(),
		d_time_passed,

		s_levels_info,
		0,//v_gene_patterns_parts_original.at(0).size(),
		v_gene_patterns_parts.at(0).size(),
		v_gene_patterns_parts_interceptions_and_differencies.at(0).size(),
		i_highest_hit
	);


	pc_log->vPrintLine(s_buf, true);


	b_use_last_start_genotype = false;
	i_cur_gen++;
	return(true);
}//bool C3LO::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)






bool  C3LOSingle::bRunIteration_3lo_berlin(uint32_t iIterationNumber, time_t tStartTime)
{
	double  d_time_passed;

	CString  s_buf;
	CString  s_levels_info;

	i_the_same = 0;
	i_improved = 0;
	i_improved_other = 0;
	i_shuffled_individuals = 0;
		
	vector<C3LOIndividual  *>  v_individuals_to_process;
	v_individuals_to_process = v_individuals_to_process_next_gen;
	v_individuals_to_process_next_gen.clear();


	C3LOIndividual  *pc_indiv_main;

	while (d_get_covering_at_level(0) < 1)
	{
		v_add_linkage_covering_at_level(0);

		pc_indiv_main = pc_get_random_individual(0);
		pc_indiv_main->i_level = 1;

		v_create_trees_from_ind(pc_indiv_main, &(v_gene_patterns_trees.at(0)), 0);

		delete pc_indiv_main;

		s_buf.Format("(%d:%.4lf) ", 0, d_get_covering_at_level(0));
		pc_log->vPrintLine(s_buf, true);
	}//while (d_get_covering_at_level(0))


	
	if (v_individuals_to_process.size() == 0)
	{
		pc_indiv_main = pc_get_random_individual(0);
		if (pc_indiv_main != NULL)  b_add_at_level(pc_indiv_main, 0);

		v_add_ind_to_process_in_next_generation(pc_indiv_main);
		//v_individuals_to_process_next_gen.push_back(pc_indiv_main);


		for (int ii = 0; ii < i_templ_length; ii++)
			pi_best_genotype[ii] = (int32_t)(pc_indiv_main->piGetPhenotype(0,0)[ii]);

		b_update_best_individual(iIterationNumber, tStartTime, pi_best_genotype, d_best_fitness_cur);
	}//if  (pv_population->size() == 0)
	else
	{
		C3LOIndividual  *pc_ind_current;
		for (int ind_to_proc = 0; ind_to_proc < v_individuals_to_process.size(); ind_to_proc++)
		{
			pc_ind_current = v_individuals_to_process.at(ind_to_proc);

			double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
			for (int i_level = 0; i_level < v_population_levels.size(); i_level++)
			{
				v_process(pc_ind_current, i_level);
			}//for (int i_level = 0; i_level < v_population_levels.size(); i_level++)
			double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();


			if (ind_to_proc % 1 == 0)
			{
				c_time_counter.bGetTimePassed(&d_time_passed);
				s_buf.Format("Ind processing: %d/%d  time:%.4lf  FFEchange: %.2lf =(%.2lf-%.2lf) Levels:%d", ind_to_proc + 1, v_individuals_to_process.size(), d_time_passed, d_ffe_after - d_ffe_before, d_ffe_after, d_ffe_before, v_population_levels.size());
				pc_log->vPrintLine(s_buf, true);
			}//if (ind_to_proc % 50 == 0)

		}//for (int ind_to_proc = 0; ind_to_proc < v_individuals_to_process.size(); ind_to_proc++)


		v_count_higher_genes_freq();



		double  d_best_fitness, d_fit_buf;
		int  *pi_curr_phenotype, *pi_best_phenotype;

		pi_best_phenotype = NULL;
		d_best_fitness = 0;
		for (int i_pop_lvl = 0; i_pop_lvl < v_population_levels.size(); i_pop_lvl++)
		{
			for (int i_ind = 0; i_ind < v_population_levels.at(i_pop_lvl).size(); i_ind++)
			{
				for (int i_pattern_lvl = 0; i_pattern_lvl < v_gene_patterns_trees.size(); i_pattern_lvl++)
				{
					pi_curr_phenotype = v_population_levels.at(i_pop_lvl).at(i_ind)->piGetPhenotype(i_pattern_lvl, 0);
					d_fit_buf = dComputeFitness(pi_curr_phenotype);
					if (d_best_fitness < d_fit_buf)
					{
						d_best_fitness = d_fit_buf;
						pi_best_phenotype = pi_curr_phenotype;
					}//if (d_best_fitness < d_fit_buf)
				}//for (int i_pattern_lvl = 0; v_gene_patterns_trees.size(); i_pattern_lvl++)		
			}//for  (int  i_ind = 0; i_ind < v_population_levels.at(i_pop_lvl).size(); i_pop_lvl++)
		}//for (int i_pop_lvl = 0; i_pop_lvl < v_population_levels.size(); i_pop_lvl++)



		if (pi_best_phenotype != NULL)
		{
			for (int ii = 0; ii < i_templ_length; ii++)
				pi_best_genotype[ii] = (int32_t)(pi_best_phenotype[ii]);

			b_update_best_individual(iIterationNumber, tStartTime, pi_best_genotype, d_best_fitness);
		}//if (pi_best_phenotype != NULL)



		 //pattern coverage
		CString  s_coverage;
		s_coverage = "";

		for (int ii = 0; ii < v_gene_patterns_trees.size(); ii++)
		{
			s_buf.Format("(%d:%.4lf) ", ii, d_get_covering_at_level(ii));
			s_coverage += s_buf;
		}//for  (int  ii = 0; ii < v_gene_patterns_trees.size(); ii++)


		CString  s_block_info;

		CString  s_phenotype_part;
		int  i_ones, i_zeros, i_other;
		vector <CString>  v_other;

		s_block_info = "";

		/*for  (int  i_gene = 0; i_gene < i_templ_length; i_gene+=10)
		{
		i_ones = 0;
		i_zeros = 0;
		i_other = 0;
		v_other.clear();
		for  (int  ii = 0; ii < v_population_levels.at(0).size(); ii++)
		{
		s_phenotype_part = v_population_levels.at(0).at(ii)->sGetPhenotype(0, i_gene, i_gene + 10);

		if  (s_phenotype_part == "1111111111")  i_ones++;
		if  (s_phenotype_part == "0000000000")  i_zeros++;

		if  (
		(s_phenotype_part != "1111111111")&&
		(s_phenotype_part != "0000000000")
		)
		{
		i_other++;

		if (!( std::find(v_other.begin(), v_other.end(), s_phenotype_part) != v_other.end() ))
		v_other.push_back(s_phenotype_part);

		}//if  (

		}//for  (int  ii = 0; ii < v_population_levels.at(0).size(); ii++)

		s_buf.Format("[%d: %d/%d %d/%d] ", i_gene, i_ones, i_zeros, i_other, (int) v_other.size());
		s_block_info += s_buf;
		}//for  (int  i_gene = 0; i_gene < i_templ_length; i_gene++)*/


		/*FILE  *pf_test;
		s_debug.Format("%d_", i_cur_gen);
		s_buf.Format("zzz_pats_%.2d.txt", i_cur_gen);
		v_save_trees_and_patterns(s_buf);//*/



		
		//d_time_passed = (uint32_t)(time(nullptr) - t_start);
		c_time_counter.bGetTimePassed(&d_time_passed);

		s_buf.Format
		(
			"blocked improvements: %d %s Iteration: %d (next gen ind: %d new shuffled ind: %d) (upgraded/upgraded other/the same: %d/%d/%d) BestFit: %.4lf [u:%.4lf] PopSize:%d  Patterns[coverage: %s] FFE:%.0lf FFE_OPTs: %d time passed:%.4lf  %s   original parts: %d parts: %d",
			i_better_but_cannot_be_added, s_block_info,
			i_cur_gen, v_individuals_to_process_next_gen.size(),  i_shuffled_individuals,
			i_improved, i_improved_other, i_the_same,

			d_best_fitness, d_unit_max_cur,
			i_get_ind_number(),

			s_coverage,
			//pc_fitness->dGetFuncEval(),
			(double)pc_problem->pcGetEvaluation()->iGetFFE(),
			pc_genotype_root->iGetSuccOpt(),
			d_time_passed,
			
			s_levels_info,
			0,//v_gene_patterns_parts_original.at(0).size(),
			v_gene_patterns_parts.at(0).size()
		);


		pc_log->vPrintLine(s_buf, true);


		i_cur_gen++;
	}//if (v_individuals_to_process.size() == 0)
	return(true);


	uint16_t i_number_of_bits = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	//osobnik w Twojej konwencji
	//int *pi_bits = nullptr;

	//wyliczenie fitnessu
	//double d_fitness_value = d_compute_fitness_value(pi_bits);

	//bool b_updated = false;

	//aktualizacja najlepszego znalezionego osobnika - w przypadku gdy metoda produkuje tylko jednego osobnika podczas pojedynczej iteracji, 
	//metoda b_update_best_individual zostanie wywolana tylko raz z bitami i fitnessem tego konkretnego osobnika
	//w ponizszym przykladzie for jest tylko po to by pokazac, ze metode b_update_best_individual trzeba wywolac dla kazdego "nowego" osobnika
	//for (int i = 0; i < 10; i++)
	//{
	//	b_updated = b_update_best_individual(iIterationNumber, tStartTime, pi_bits, d_fitness_value) || b_updated;
	//}//for (int i = 0; i < 10; i++)

	//bity najlepszego znalezionego do teraz osobnika
	//pc_best_individual->pcGetGenotype()->piGetBits();

	int *pi_bits = new int[i_number_of_bits];

	for (uint16_t i = 0; i < i_number_of_bits; i++)
	{
		//wywolanie utilsow do liczb pseudolosowych (metoda powinna korzystac jedynie z RandUtils do tego typu rzeczy 
		//- w przeciwnym wypadku podany jako parametr seed nie bedzie uwzgledniony)
		*(pi_bits + i) = RandUtils::iRandNumber(0, 1);
	}//for (uint16_t i = 0; i < i_number_of_bits; i++)

	double d_fitness_value = d_compute_fitness_value(pi_bits);
	bool b_updated = b_update_best_individual(iIterationNumber, tStartTime, pi_bits, d_fitness_value);



	//::MessageBox(NULL, "", "", MB_OK);

	//if ((time(nullptr) - tStartTime) % t_log_frequency == 0)
	{
		CString s_log_message;
		s_log_message.Format("best fitness: %f; best unitation: %f; ffe: %u; time: %u; iteration: %u", pc_best_individual->dGetFitnessValue(),
			pc_best_individual->pcGetGenotype()->dGetUnitation(), pc_problem->pcGetEvaluation()->iGetFFE(), (uint32_t)(time(nullptr) - tStartTime),
			iIterationNumber);

		pc_log->vPrintLine(s_log_message, true);
	}//if ((time(nullptr) - tStartTime) % t_log_frequency == 0)


	return b_updated;
}//bool C3LO::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)






bool  C3LOSingle::b_run_flat_pop_for_ind(C3LOIndividual  *pcIndMain, int iPatternLevel, vector<C3LOIndividual  *>  *pvLevelsOffspringsDifferentToParent)
{
	C3LOIndividual  *pc_indiv_donator;
	C3LOIndividual  *pc_indiv_new;

	double  d_ind_fitness, d_donator_fitness;

	int  i_max_lvl_number;
	i_max_lvl_number = v_population_levels.size() - 1;


	//for  (int  i_donator_lvl = i_max_lvl_number; i_donator_lvl >= 0; i_donator_lvl--)
	for (int i_donator_lvl = 0; i_donator_lvl <= i_max_lvl_number; i_donator_lvl++)
	{
		pc_indiv_donator = pc_get_individual_with_lvl_tournament(pcIndMain);

		if (pc_indiv_donator != NULL)
		{
			d_ind_fitness = pcIndMain->dComputeFitness(iPatternLevel);
			d_donator_fitness = pc_indiv_donator->dComputeFitness(iPatternLevel);

			v_tiny_restricted_evolution(pcIndMain, pc_indiv_donator, iPatternLevel, 0, 1);


			if (d_ind_fitness > pcIndMain->dComputeFitness(iPatternLevel))
			{
				if (v_gene_patterns_trees.size() > iPatternLevel + 1)
				{
					b_run_flat_pop_for_ind(pcIndMain, iPatternLevel + 1, pvLevelsOffspringsDifferentToParent);
				}//if  (v_gene_patterns_trees.size() > iPatternLevel + 1)

			}//if (d_ind_fitness  >  pcIndMain->dComputeFitness(iPatternLevel))


		}//if  (pc_indiv_donator  != NULL)
	}//for  (int  i_donator_lvl = i_max_lvl_number; i_max_lvl_number >= 0; i_max_lvl_number++)


	 /*C3LOIndividual  *pc_indiv_new;
	 pc_indiv_new = pc_get_random_individual();
	 if  (pc_indiv_new != NULL)  b_add_at_level(pc_indiv_new, 0);*/
	return(false);
}//bool  C3LO::b_run_flat_pop_for_ind(C3LOIndividual  *pcIndMain, int  iEvolutionLvl)



 
int C3LOSingle::i_get_parrent_offset_respect_inner_level(vector<C3LOIndividual *>  *pvIndToChoose)
{
	if (pvIndToChoose->size() < 1) return(-1);

	int i_lowest_inner_level;
	
	i_lowest_inner_level = pvIndToChoose->at(0)->i_level_inner;

	for (int ii = 0; ii < pvIndToChoose->size(); ii++)
	{
		if (i_lowest_inner_level > pvIndToChoose->at(ii)->i_level_inner)  i_lowest_inner_level = pvIndToChoose->at(ii)->i_level_inner;
	}//for (int ii = 0; ii < pvIndToChoose->size(); ii++)


	vector<int>  v_par_candidates;
	for (int ii = 0; ii < pvIndToChoose->size(); ii++)
	{
		if (i_lowest_inner_level == pvIndToChoose->at(ii)->i_level_inner)  v_par_candidates.push_back(ii);
	}//for (int ii = 0; ii < pvIndToChoose->size(); ii++)

	if (v_par_candidates.size() == 0)  ::MessageBox(NULL, "no ind found on the chosen inner level", "ERROR", MB_OK);

	int  i_result;
	i_result = v_par_candidates.at(RandUtils::iRandNumber(0, v_par_candidates.size() - 1));

	return(i_result);
}//int C3LO::i_get_parrent_offset_respect_inner_level(vector<C3LOIndividual *>  *pvIndToChoose)



void  C3LOSingle::v_process(C3LOIndividual  *pcParentMain, int  iPopLevel)
{
	CString  s_buf;
	vector<C3LOIndividual *> v_pop_level_copy;
	//bool b_crossing_effective;
	int  *pi_genotype, *pi_original_genotype;
	int  *pi_genotype_to_shuffle;
	int  i_shuffle;

	double  d_genotype_fitness_on_start;
	double  d_genotype_fitness;

	int  i_other_parent_offset;
	C3LOIndividual *pc_other_parent;
	pi_genotype = new int[i_templ_length];
	pi_genotype_to_shuffle = new int[i_templ_length];

	int  i_shuffled_ind_created = 0;
	i_shuffle = 0;

	double  d_ffe_sum;
	d_ffe_sum = 0;

	int i_pattern_level = 0;
	//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		v_pop_level_copy = v_population_levels.at(iPopLevel);

		if (pv_other_pops != NULL)
		{
			for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
			{
				if (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPopLevel)
				{
					for (int i_ind = 0; i_ind < pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPopLevel).size(); i_ind++)
					{
						v_pop_level_copy.push_back(pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPopLevel).at(i_ind));
					}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
				}//if  (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPatternLevel)
			}//for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
		}//if (pv_other_pops != NULL)


		//b_crossing_effective = false;
		pi_original_genotype = pcParentMain->piGetPhenotype(i_pattern_level, 0);

		if (pi_original_genotype != NULL)
		{

			for (int ii = 0; ii < i_templ_length; ii++)
				pi_genotype[ii] = pi_original_genotype[ii];
			d_genotype_fitness_on_start = dComputeFitness(pi_genotype);

			while (v_pop_level_copy.size() > 0)
			{
				//i_other_parent_offset = RandUtils::iRandNumber(0, v_pop_level_copy.size() - 1);
				i_other_parent_offset = i_get_parrent_offset_respect_inner_level(&v_pop_level_copy);
				
				pc_other_parent = v_pop_level_copy.at(i_other_parent_offset);
				v_pop_level_copy.erase(v_pop_level_copy.begin() + i_other_parent_offset);


				double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();

				if (pc_other_parent != pcParentMain)
				{
					int  i_level_before = pcParentMain->i_level_inner;
					if (b_process(pcParentMain, pc_other_parent, pi_genotype, pi_genotype_to_shuffle, &i_shuffle, i_pattern_level) == true) //b_crossing_effective = true;
					{
						i_improved++;

						double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
						s_buf.Format("CROSSED  FFEchange: %.2lf OTHER ind left: %d  INNER LEVELS: %d/%d/%d", d_ffe_after - d_ffe_before, v_pop_level_copy.size(), i_level_before, pcParentMain->i_level_inner, pc_other_parent->i_level_inner);
						//if (pv_other_pops != NULL)  pc_log->vPrintLine(s_buf, true);

						d_ffe_sum += d_ffe_after - d_ffe_before;
					}//if (b_process(pcParentMain, pc_other_parent, pi_genotype, pi_genotype_to_shuffle, &i_shuffle, i_pattern_level) == true) //b_crossing_effective = true;
					else
					{
						double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
						s_buf.Format("NOT CROSSED  FFEchange: %.2lf OTHER ind left: %d  INNER LEVELS: %d/%d/%d", d_ffe_after - d_ffe_before, v_pop_level_copy.size(), i_level_before, pcParentMain->i_level_inner, pc_other_parent->i_level_inner);
						//if (pv_other_pops != NULL)  pc_log->vPrintLine(s_buf, true);

						d_ffe_sum += d_ffe_after - d_ffe_before;
					}//else  if (b_process(pcParentMain, pc_other_parent, pi_genotype, pi_genotype_to_shuffle, &i_shuffle, i_pattern_level) == true) //b_crossing_effective = true;

				}//if (pc_other_parent != pcParentMain)

				
			}//while ((v_pop_level.size() > 0) && (b_crossing_effective == false))

					


			d_genotype_fitness = dComputeFitness(pi_genotype);
			
		}//if (pi_original_genotype !=  NULL)
	}//for  (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)

	s_buf.Format("V_PROCES FINISHED FFE COST: %.2lf   fintess:(%.4lf->%.4lf)", d_ffe_sum, d_genotype_fitness_on_start, d_genotype_fitness);
	pc_log->vPrintLine(s_buf, true);

	delete pi_genotype;
	delete pi_genotype_to_shuffle;

}//void  C3LO::v_process(C3LOIndividual  *pcParentMain, int  iPopLevel)


void  C3LOSingle::v_process_berlin(C3LOIndividual  *pcParentMain, int  iPopLevel)
{
	vector<C3LOIndividual *> v_pop_level_copy;
	bool b_crossing_effective;
	int  *pi_genotype, *pi_original_genotype;
	int  *pi_genotype_to_shuffle;
	int  i_shuffle;

	double  d_genotype_fitness_on_start;
	double  d_genotype_fitness;

	int  i_other_parent_offset;
	C3LOIndividual *pc_other_parent;
	pi_genotype = new int[i_templ_length];
	pi_genotype_to_shuffle = new int[i_templ_length];

	int  i_shuffled_ind_created = 0;
	i_shuffle = 0;
	
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		v_pop_level_copy = v_population_levels.at(iPopLevel);
		b_crossing_effective = false;
		pi_original_genotype = pcParentMain->piGetPhenotype(i_pattern_level, 0);

		if (pi_original_genotype != NULL)
		{

			for (int ii = 0; ii < i_templ_length; ii++)
				pi_genotype[ii] = pi_original_genotype[ii];
			d_genotype_fitness_on_start = dComputeFitness(pi_genotype);

			while (v_pop_level_copy.size() > 0)
			{
				i_other_parent_offset = RandUtils::iRandNumber(0, v_pop_level_copy.size() - 1);
				pc_other_parent = v_pop_level_copy.at(i_other_parent_offset);
				v_pop_level_copy.erase(v_pop_level_copy.begin() + i_other_parent_offset);


				if (pc_other_parent != pcParentMain)
				{
					if (b_process(pcParentMain, pc_other_parent, pi_genotype, pi_genotype_to_shuffle, &i_shuffle, i_pattern_level) == true)  b_crossing_effective = true;
				}//if (pc_other_parent != pcParentMain)
			}//while ((v_pop_level.size() > 0) && (b_crossing_effective == false))


			if ((b_crossing_effective == false) && (i_shuffle == 1))
			{
				//the we shuffle the individual with the block from ParentOther
				int i_offspring_level = 0;
				if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
				if (i_offspring_level < iPopLevel)  i_offspring_level = iPopLevel;
				if (i_offspring_level < i_pattern_level)  i_offspring_level = i_pattern_level;

				//if ofsspring level is higher than parents - make a new individual
				C3LOIndividual *pc_indiv_new;
				pc_indiv_new = new C3LOIndividual(i_offspring_level, i_templ_length, pc_problem, &v_optimization_orders, this);
				pc_indiv_new->vRandomInit();

				if (b_shuffle_the_genotype(pi_genotype_to_shuffle, pc_indiv_new, i_pattern_level) == true)
				{
					//SHUFFLING!
					i_shuffled_individuals++;

					pc_indiv_new->vSetGenotype(pi_genotype_to_shuffle);
					pc_indiv_new->vSetOptimized(false);
					pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);

					b_add_at_level(pc_indiv_new, i_offspring_level);

					v_add_ind_to_process_in_next_generation(pc_indiv_new);
				}//if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new) == true)
				else
					delete pc_indiv_new;
			}//if ( (b_crossing_effective == false)&&(i_shuffle == 1) )


			d_genotype_fitness = dComputeFitness(pi_genotype);
			//if ((b_crossing_effective == true) && (d_genotype_fitness <= d_genotype_fitness_on_start))  ::MessageBox(NULL, "IMPOSSIBLE", "IMPOSSIBLE", MB_OK);

			if (d_genotype_fitness > d_genotype_fitness_on_start)
			{
				if (b_can_ind_be_added_at_any_level(pi_genotype) == true)
				{
					i_improved++;

					int i_offspring_level = 0;
					if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
					if (i_offspring_level < iPopLevel)  i_offspring_level = iPopLevel;
					if (i_offspring_level < i_pattern_level)  i_offspring_level = i_pattern_level;

					//if ofsspring level is the same as parents - replace parent
					if (i_offspring_level == pcParentMain->i_level)
					{
						pcParentMain->vSetGenotype(pi_genotype);
						pcParentMain->vSetOptimized(false);
						pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);

						v_add_ind_to_process_in_next_generation(pcParentMain);
						//v_individuals_to_process_next_gen.push_back(pcParentMain);
					}//if (i_offspring_level == pcParentOther->i_level)
					else
					{
						//if ofsspring level is higher than parents - make a new individual
						C3LOIndividual *pc_indiv_new;
						pc_indiv_new = new C3LOIndividual(i_offspring_level, i_templ_length, pc_problem, &v_optimization_orders, this);
						pc_indiv_new->vRandomInit();
						pc_indiv_new->vSetGenotype(pi_genotype);
						pc_indiv_new->vSetOptimized(false);
						pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);


						b_add_at_level(pc_indiv_new, i_offspring_level);

						v_add_ind_to_process_in_next_generation(pc_indiv_new);
						//v_individuals_to_process_next_gen.push_back(pc_indiv_new);
					}//else  if ofsspring level is the same as parents - replace parent

				}//if (b_can_ind_be_added_at_any_level(pi_genotype) == true)
				else
				{
					//WHAT TO DO HERE???
					//nothing for now...
					i_the_same++;
				}//else  if (b_can_ind_be_added_at_any_level(pi_genotype) == true)
			}//if (d_genotype_fitness > pcParentMain->dComputeFitness(i_pattern_level))		
		}//if (pi_original_genotype !=  NULL)
	}//for  (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)

	delete pi_genotype;
	delete pi_genotype_to_shuffle;

}//void  C3LO::v_process(C3LOIndividual  *pcParentMain, int  iPopLevel)



bool  C3LOSingle::b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int  *piGenotypeToShuffle, int  *piShuffle, int iPatternLevel)
{
	CString  s_buf;

	//int  i_max_pattern_level;
	//i_max_pattern_level = v_gene_patterns_trees.size() - 1;

	//s_buf.Format("i_max_pattern_level = %d", i_max_pattern_level);
	//::MessageBox(NULL, s_buf, s_buf, MB_OK);
	pcParentMain->dComputeFitness(iPatternLevel);
	pcParentOther->dComputeFitness(iPatternLevel);

	int  i_max_ind_level;


	if (pcParentMain->i_level  >  pcParentOther->i_level)
		i_max_ind_level = pcParentMain->i_level;
	else
		i_max_ind_level = pcParentOther->i_level;

	vector<C3LOPattern *>  v_allowed_trees;


	int  i_result = 0;

	//first produce linkage...
	bool  b_at_least_one_update = false;
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		/*if (b_produce_new_linkage(pcParentMain, pcParentOther, i_pattern_level) == true)
		{
			v_update_high_lvl_genes_and_trees();
			v_create_trees_from_ind(pcParentMain, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);
			v_update_high_lvl_genes_and_trees();
		}//if (b_produce_new_linkage(pcParentMain, pcParentOther, i_level) == true)


		if (b_produce_new_linkage(pcParentMain, pcParentOther, i_pattern_level) == true)
		{
			v_update_high_lvl_genes_and_trees();
			v_create_trees_from_ind(pcParentOther, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);
			v_update_high_lvl_genes_and_trees();
		}//if (b_produce_new_linkage(pcParentMain, pcParentOther, i_level) == true)*/

		//i_result = i_cross_dsm_glueing_similarities_from_different_start_genes(pi_gene_realtions, piCurrentGenotype, pcParentMain, pcParentOther, i_pattern_level);
		i_result = i_cross_dsm(piCurrentGenotype, pcParentMain, pcParentOther, i_pattern_level);
		//i_result = i_cross_fluent_scraps(piCurrentGenotype, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, i_pattern_level);
		//i_result = i_cross_intercept_scraps(piCurrentGenotype, pcParentMain, pcParentOther, i_pattern_level);
		
	}//for (int i_level = 0; i_level < i_max_pattern_level; i_level++)


	if (i_result == 2)  return(true);

	return(false);
}//bool  C3LO::b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int iPatternLevel)



bool  C3LOSingle::b_process_berlin(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int  *piGenotypeToShuffle, int  *piShuffle, int iPatternLevel)
{
	int  i_max_pattern_level;
	i_max_pattern_level = v_gene_patterns_trees.size() - 1;

	pcParentMain->dComputeFitness(i_max_pattern_level);
	pcParentOther->dComputeFitness(i_max_pattern_level);

	int  i_max_ind_level;

	
	if  (pcParentMain->i_level  >  pcParentOther->i_level)
		i_max_ind_level = pcParentMain->i_level;
	else
		i_max_ind_level = pcParentOther->i_level;

	vector<C3LOPattern *>  v_allowed_trees;
	

	int  i_result = 0;

	//first produce linkage...
	bool  b_at_least_one_update = false;
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		if (b_produce_new_linkage(pcParentMain, pcParentOther, i_pattern_level) == true)
		{
			v_update_high_lvl_genes_and_trees();
			//v_create_trees_from_ind_2(pcParentMain, pcParentOther, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);
			v_create_trees_from_ind(pcParentMain, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);

			
			v_update_high_lvl_genes_and_trees();
		}//if (b_produce_new_linkage(pcParentMain, pcParentOther, i_level) == true)

		
		if (b_produce_new_linkage(pcParentMain, pcParentOther, i_pattern_level) == true)
		{
			v_update_high_lvl_genes_and_trees();
			v_create_trees_from_ind(pcParentOther, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);

			//vector<C3LOPattern *>  v_gene_patterns_trees_buf;
			//v_gene_patterns_trees_buf = v_gene_patterns_trees.at(i_pattern_level);
			//v_gene_patterns_trees.at(i_pattern_level).clear();

			//v_build_pattern_tree(&v_gene_patterns_trees_buf, &(v_gene_patterns_trees.at(i_pattern_level)), i_pattern_level);

			v_update_high_lvl_genes_and_trees();
		}//if (b_produce_new_linkage(pcParentMain, pcParentOther, i_level) == true)


		//double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		
		i_result = i_cross_fluent_scraps(piCurrentGenotype, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, i_pattern_level);
		//if (b_cross_fluent_scraps(piCurrentGenotype, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, i_pattern_level) == true)  b_at_least_one_update = true;

		//double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();

		//CString  s_buf;
		//s_buf.Format("Single i_cross_fluent_scraps FFEchange: %.2lf =(%.2lf-%.2lf)", d_ffe_after - d_ffe_before, d_ffe_after, d_ffe_before);
		//pc_log->vPrintLine(s_buf, true);

		


		/*//OLD ONE
		v_allowed_trees.clear();

		for (int ii = 0; ii < v_gene_patterns_trees.at(i_pattern_level).size(); ii++)
		{
			if (v_gene_patterns_trees.at(i_pattern_level).at(ii)->bMarkedPhenotypeDifferences(piCurrentGenotype, pcParentOther->piGetPhenotype(i_pattern_level, 0)) == true)  v_allowed_trees.push_back(v_gene_patterns_trees.at(i_pattern_level).at(ii));
		}//for  (int ii = 0; ii < v_gene_patterns_trees.size(); ii++)


		if (v_allowed_trees.size() > 0)
		{
			int i_chosen_tree;
			i_chosen_tree = RandUtils::iRandNumber(0, v_allowed_trees.size() - 1);
			b_cross_individuals_by_tree_branch_or_part(piCurrentGenotype, pcParentMain, pcParentOther, v_allowed_trees.at(i_chosen_tree), i_pattern_level);
		}//if  (v_allowed_trees.size() > 0)*/
		
	}//for (int i_level = 0; i_level < i_max_pattern_level; i_level++)


	if (i_result > 0)  return(true);
	
	return(false);
}//bool  C3LO::b_process(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int *piCurrentGenotype, int iPatternLevel)




void  C3LOSingle::v_tiny_restricted_evolution(C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentDonator, int  iPatternLevel, int  iPopDestLvl, int iRepeatation /*= -1*/)
{
	if (v_gene_patterns_trees.size()  <  iPatternLevel)  return;

	C3LOIndividual  *pc_parent_main;

	pc_parent_main = pcParentMain;
	bool  b_continue_tribe = true;
	bool  b_updated_after_linkage_gen;
	vector<C3LOPattern *>  v_allowed_trees;


	/*FILE  *pf_test;
	pf_test  = fopen("test.txt", "w+");
	fprintf(pf_test, "MAIN:\n");
	pc_parent_main->vSave(pf_test);
	fprintf(pf_test, "DONATOR:\n");
	pcParentDonator->vSave(pf_test);*/

	int  i_iteration = 0;
	b_updated_after_linkage_gen = false;
	while (b_continue_tribe == true)
	{
		b_continue_tribe = false;

		//pc_parent_main->b_optimized = false;
		//pcParentDonator->b_optimized = false;
		double  d_main_fitness, d_donator_fitness;
		d_main_fitness = pc_parent_main->dComputeFitness(iPatternLevel); //force phenotype creation
		d_donator_fitness = pcParentDonator->dComputeFitness(iPatternLevel); //force phenotype creation



		if (b_produce_new_linkage(pc_parent_main, pcParentDonator, iPatternLevel) == true)
		{
			if (iPatternLevel > 0)
			{
				int  ig = 0;
				ig++;
			}//if  (iPatternLevel > 0)

			if (s_debug != "")
			{
				CString  s_debug_name;

				s_debug_name.Format("%s_patterns_00_before_level_%d.txt", s_debug, iPatternLevel);
				v_save_trees_and_patterns(s_debug_name);
			}//if  (b_debug  == true)


			v_update_high_lvl_genes_and_trees();
			v_create_trees_from_ind(pc_parent_main, &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);

			if (s_debug != "")
			{
				CString  s_debug_name;

				s_debug_name.Format("%s_patterns_parts_from_ind.txt", s_debug);
				pc_parent_main->vSaveGenePatternParts(s_debug_name);
			}//if  (b_debug  == true)



			if (s_debug != "")
			{
				CString  s_debug_name;

				s_debug_name.Format("%s_patterns_01_after_level_%d.txt", s_debug, iPatternLevel);
				v_save_trees_and_patterns(s_debug_name);
			}//if  (b_debug  == true)


			 /*FILE  *pf_before;
			 pf_before = fopen("zz_before.txt", "w+");
			 fprintf(pf_before, "\n\nPATTERNS:\n");
			 for  (int  ii = 0; ii < v_gene_patterns_trees.size(); ii++)
			 v_gene_patterns_trees.at(0).at(ii)->vSavePattern(pf_before, pc_fitness);
			 fclose(pf_before);//*/



			vector<C3LOPattern *>  v_gene_patterns_trees_buf;
			v_gene_patterns_trees_buf = v_gene_patterns_trees.at(iPatternLevel);
			v_gene_patterns_trees.at(iPatternLevel).clear();

			v_build_pattern_tree(&v_gene_patterns_trees_buf, &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);

			v_update_high_lvl_genes_and_trees();


			if (s_debug != "")
			{
				CString  s_debug_name;

				s_debug_name.Format("%s_patterns_02_built_level_%d.txt", s_debug, iPatternLevel);
				v_save_trees_and_patterns(s_debug_name);
			}//if  (b_debug  == true)


			 /*FILE  *pf_after;
			 pf_after = fopen("zz_after.txt", "w+");
			 fprintf(pf_after, "\n\nPATTERNS:\n");
			 for  (int  ii = 0; ii < v_gene_patterns_trees.at(0).size(); ii++)
			 v_gene_patterns_trees.at(0).at(ii)->vSavePattern(pf_after, pc_fitness);
			 fclose(pf_after);*/


			 //::vShow("patterns!");//*/
		}//if  (b_produce_new_linkage(pc_parent_main, pcParentDonator, i_lvl + 1) == true)


		if (b_produce_new_linkage(pc_parent_main, pcParentDonator, iPatternLevel) == true)
		{
			v_update_high_lvl_genes_and_trees();
			v_create_trees_from_ind(pcParentDonator, &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);

			vector<C3LOPattern *>  v_gene_patterns_trees_buf;
			v_gene_patterns_trees_buf = v_gene_patterns_trees.at(iPatternLevel);
			v_gene_patterns_trees.at(iPatternLevel).clear();

			v_build_pattern_tree(&v_gene_patterns_trees_buf, &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);

			v_update_high_lvl_genes_and_trees();
		}//if  (b_produce_new_linkage(pc_parent_main, pcParentDonator, i_lvl + 1) == true)



		if (v_allowed_trees.size() == 0)
		{
			for (int ii = 0; ii < v_gene_patterns_trees.at(iPatternLevel).size(); ii++)
			{
				if (v_gene_patterns_trees.at(iPatternLevel).at(ii)->bMarkedPhenotypeDifferences(pc_parent_main, pcParentDonator, iPatternLevel) == true)  v_allowed_trees.push_back(v_gene_patterns_trees.at(iPatternLevel).at(ii));
			}//for  (int ii = 0; ii < v_gene_patterns_trees.size(); ii++)
		}//if  (v_allowed_trees.size() == 0)

		if (v_allowed_trees.size() > 0)
		{
			int i_chosen_tree;
			//i_chosen_tree = pc_random_gen->Next(0, v_allowed_trees.size());
			i_chosen_tree = RandUtils::iRandNumber(0, v_allowed_trees.size() - 1);

			//b_cross_individuals_by_tree_branch(pc_parent_main, pcParentDonator, &pc_offspring, &pc_reversed_offspring, v_allowed_trees.at(i_chosen_tree), iPatternLevel);
			//b_cross_individuals_by_tree_branch_longest_reasonable_branch(pc_parent_main, pcParentDonator, &pc_offspring, &pc_reversed_offspring, v_allowed_trees.at(i_chosen_tree), iPatternLevel);
			b_cross_individuals_by_tree_branch_or_part(pc_parent_main->piGetPhenotype(iPatternLevel, 0), pc_parent_main, pcParentDonator, v_allowed_trees.at(i_chosen_tree), iPatternLevel);


			if (pc_parent_main->dComputeFitness(iPatternLevel) <= d_main_fitness)
			{
				C3LOIndividual  *pc_best_update;
				C3LOIndividual  *pc_best_different_ind;

				pc_best_update = NULL;
				pc_best_different_ind = NULL;
				//i_brutal_random_search++;
				//if (pcParentMain->bLinkedGenesBrutalSearch(v_allowed_trees.at(i_chosen_tree), iPatternLevel, &(v_gene_patterns_parts_original.at(iPatternLevel)),  &pc_best_update, &pc_best_different_ind) == true)
				if (pcParentMain->bLinkedGenesBrutalSearch(v_allowed_trees.at(i_chosen_tree), iPatternLevel, &(v_gene_patterns_parts.at(iPatternLevel)), &pc_best_update, &pc_best_different_ind) == true)
				{
					if (pc_best_update != NULL)
					{
						pc_best_update->i_level = pcParentMain->i_level + 1;
						pc_best_update->vCopyTo(pcParentMain);
					}//if (pc_best_update != NULL)
				}//if  (pcParentMain->bLinkedGenesBrutalSearch(v_allowed_trees.at(i_chosen_tree), iPatternLevel, &pc_best_update, &pc_best_different_ind)  ==  true)
				
			}//if (pc_parent_main->dComputeFitness(iPatternLevel) <= d_main_fitness)

			v_allowed_trees.erase(v_allowed_trees.begin() + i_chosen_tree);
			if (v_allowed_trees.size() > 0)  b_continue_tribe = true;
		}//if  (v_allowed_trees.size() > 0)
		else
		{
			if (b_updated_after_linkage_gen == true)
			{
				b_updated_after_linkage_gen = false;
				b_continue_tribe = true;
			}//if  (b_updated_after_linkage_gen  ==  true)
		}//else  if  (v_allowed_trees.size() > 0)

		 //fprintf(pf_test, "\nOFFSPRING:\n");


		if (d_main_fitness < pc_parent_main->dComputeFitness(iPatternLevel))  b_updated_after_linkage_gen = true;

		i_iteration++;
		if (iRepeatation > 0)
		{
			if (i_iteration >= iRepeatation)  b_continue_tribe = false;
		}//if  (iRepeatation > 0) 
	}//while  (b_continue_tribe == true)


}//void  C3LO::v_tiny_restricted_evolution(C3LOIndividual  *pcParent0, C3LOIndividual  *pcParent1, C3LOIndividual  **pcOffspring)




bool  C3LOSingle::b_cross_individuals_by_tree_branch_longest_reasonable_branch(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOIndividual  **pcOffspring, C3LOIndividual  **pcReversedOffspring, C3LOPattern *pcCrossingTree, int  iPatternLevel)
{
	//C3LOIndividual  *pc_best_acceptable_offspring;
	//pc_best_acceptable_offspring = NULL;

	C3LOPattern  *pc_pattern_to_cross;

	pc_pattern_to_cross = pcCrossingTree->pcGetLongestReasonableBranch(pcInd0, pcInd1, pi_genotype_tool, iPatternLevel);
	if (pc_pattern_to_cross == NULL) return(false);


	if (pc_pattern_to_cross->bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == true)
	{
		b_cross_infect_individual_by_tree(pcInd0, pcInd1, pc_pattern_to_cross, iPatternLevel);
		if ((*pcOffspring)->dComputeFitness(iPatternLevel) > pcInd0->dComputeFitness(iPatternLevel))  return(true);



		if (((*pcOffspring)->dComputeFitness(iPatternLevel) < pcInd0->dComputeFitness(iPatternLevel)) && (*pcReversedOffspring == NULL))
		{
			b_cross_infect_individual_by_tree(pcInd1, pcInd0, pc_pattern_to_cross, iPatternLevel);

			if (pc_pattern_to_cross == pcCrossingTree)//if we are at the top of the tree
			{
				if (
					(
					((*pcOffspring)->dComputeFitness(iPatternLevel) > pcInd0->dComputeFitness(iPatternLevel)) &&
						((*pcReversedOffspring)->dComputeFitness(iPatternLevel) > pcInd1->dComputeFitness(iPatternLevel))
						)
					||
					(
					((*pcOffspring)->dComputeFitness(iPatternLevel) < pcInd0->dComputeFitness(iPatternLevel)) &&
						((*pcReversedOffspring)->dComputeFitness(iPatternLevel) < pcInd1->dComputeFitness(iPatternLevel))
						)
					)
				{
					b_add_higher_level_gene(NULL, pcCrossingTree->iGetPatternLevel() + 1);//create level of high level genes
					v_insert_gene_pattern_part_level(pcCrossingTree->iGetPatternLevel() + 1);
					v_add_gene_pattern_tree(NULL, pcCrossingTree->iGetPatternLevel() + 1);
					v_update_high_lvl_genes_and_trees();
				}//if  (
			}//if  (i_pattern == 0)//if we are at the top of the tree



			if ((*pcReversedOffspring)->dComputeFitness(iPatternLevel) < pcInd1->dComputeFitness(iPatternLevel))  *pcReversedOffspring = NULL;

		}//if (((*pcOffspring)->dComputeFitness(iPatternLevel) < pcInd0->dComputeFitness(iPatternLevel)) && (*pcReversedOffspring == NULL))

	}//if  (pcCrossingTree->bMarkedPhenotypeDifferences(pcInd0, pcInd1) == true)


	return(false);
}//bool  C3LO::b_cross_individuals_by_tree_branch_longest_reasonable_branch(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOIndividual  **pcOffspring, C3LOIndividual  **pcReversedOffspring, C3LOPattern *pcCrossingTree, int  iPatternLevel)



 
void  C3LOSingle::v_get_all_intercepting_patterns_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel)
{
	for (int i_pattern = 0; i_pattern < v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).size(); i_pattern++)
	{
		if (pc_parent->b_use_dsm == false || v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).at(i_pattern)->pvGetPattern()->size() < i_templ_length)
		{
			//if (v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern)->bDoITouchGeneGroup(piDifferenceMask) == true)
			if (v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).at(i_pattern)->bDoITouchGeneGroup(piDifferenceMask) == true)
				pvPatternPartsMarkingDifferences->push_back(v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).at(i_pattern));
			//pvPatternPartsMarkingDifferences->push_back(v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern));
		}//if (pc_parent->b_use_dsm == false || v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).at(i_pattern)->pvGetPattern()->size() < i_templ_length)
	}//for (int i_pattern = 0; i_pattern < v_gene_patterns_parts_original.at(iPatternLevel).size(); i_pattern++)
}//void  C3LO::v_get_all_original_pattern_parts_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel)



void  C3LOSingle::v_get_all_original_pattern_parts_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel)
{
	//for (int i_pattern = 0; i_pattern < v_gene_patterns_parts_original.at(iPatternLevel).size(); i_pattern++)
	for (int i_pattern = 0; i_pattern < v_gene_patterns_parts.at(iPatternLevel).size(); i_pattern++)
	{
		//if (v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern)->bDoITouchGeneGroup(piDifferenceMask) == true)
		if (v_gene_patterns_parts.at(iPatternLevel).at(i_pattern)->bDoITouchGeneGroup(piDifferenceMask) == true)
			pvPatternPartsMarkingDifferences->push_back(v_gene_patterns_parts.at(iPatternLevel).at(i_pattern));
			//pvPatternPartsMarkingDifferences->push_back(v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern));
	}//for (int i_pattern = 0; i_pattern < v_gene_patterns_parts_original.at(iPatternLevel).size(); i_pattern++)
}//void  C3LO::v_get_all_original_pattern_parts_marking_the_differences(int  *piDifferenceMask, vector<C3LOPattern *> *pvPatternPartsMarkingDifferences, int iPatternLevel)


C3LOPattern*  C3LOSingle::pc_get_any_original_pattern_part_extending_gene_group(int  *piGeneGroup, int iPatternLevel)
{
	//if (v_gene_patterns_parts_original.size() < iPatternLevel)  return(NULL);
	if (v_gene_patterns_parts.size() < iPatternLevel)  return(NULL);

	C3LOPattern  *pc_result;
	pc_result = NULL;

	
	vector<C3LOPattern *>  v_acceptable_patterns;
	//for (int ii = 0; ii < v_gene_patterns_parts_original.at(iPatternLevel).size(); ii++)
	for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)
	{
		//if (v_gene_patterns_parts_original.at(iPatternLevel).at(ii)->bDoIExtendGeneGroup(piGeneGroup) == true)
			//v_acceptable_patterns.push_back(v_gene_patterns_parts_original.at(iPatternLevel).at(ii));
		if (v_gene_patterns_parts.at(iPatternLevel).at(ii)->bDoIExtendGeneGroup(piGeneGroup) == true)
			v_acceptable_patterns.push_back(v_gene_patterns_parts.at(iPatternLevel).at(ii));
	}//for (int ii = 0; ii < v_gene_patterns_parts_original.size(); ii++)
	
	if (v_acceptable_patterns.size() == 0)  return(NULL);


	pc_result = v_acceptable_patterns.at(RandUtils::iRandNumber(0, v_acceptable_patterns.size() - 1));
	

	return(pc_result);
}//C3LOPattern*  C3LO::pc_get_any_original_pattern_part_extending_gene_group(int  *piGeneGroup, int iPatternLevel)




bool  C3LOSingle::b_shuffle_the_genotype(int  *piGenotypeToShuffle, C3LOIndividual  *pcIndSupportingOptOrder, int iPatternLevel)
{
	return(false);
	int  *pi_shuffling_mask;
	int  *pi_test_phenotype;
	pi_shuffling_mask = new int[i_templ_length];
	pi_test_phenotype = new int[i_templ_length];
	vector<int>  v_available_genes;

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_shuffling_mask[ii] = 0;


	C3LOPattern  *pc_pattern_part;
	int  i_chosen_offset;
	int  i_processed_gene = -1;
	bool  b_done = false;
	pc_pattern_part = NULL;

	while  (b_done ==  false)
	{
		if (pc_pattern_part == NULL)
		{
			v_available_genes.clear();

			for (int ii = 0; ii < i_templ_length; ii++)
			{
				if (pi_shuffling_mask[ii] == 0)
				{
					if (v_genes_marked_by_linkage.at(0)[ii] > 0)  v_available_genes.push_back(ii);
				}//if  (pi_shuffling_mask[ii] == 0)
			}//for (int ii = 0; ii < i_templ_length; ii++)

			if (v_available_genes.size() > 0)
			{
				i_chosen_offset = RandUtils::iRandNumber(0, v_available_genes.size() - 1);
				i_processed_gene = v_available_genes.at(i_chosen_offset);

				pi_shuffling_mask[i_processed_gene] = 1;
				pc_pattern_part = pc_get_any_original_pattern_part_extending_gene_group(pi_shuffling_mask, 0);//it must return something, because from the v_genes_marked_by_linkage.at(0)[ii] we KNOW it was marked by some pattern
			}//if (v_available_genes.size() > 0)
			else
			{
				b_done = true;
			}//else  if (v_available_genes.size() > 0)
		}//if (pc_pattern_part == NULL)
		else
		{
			pc_pattern_part = pc_get_any_original_pattern_part_extending_gene_group(pi_shuffling_mask, 0);
		}//else  if (pc_pattern_part == NULL)

		if (pc_pattern_part != NULL)
		{
			for  (int ii = 0; ii < pc_pattern_part->v_pattern.size(); ii++)
				pi_shuffling_mask[pc_pattern_part->v_pattern.at(ii).iGenePos()] = 1;

			//we shuffle marked genes...
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				if (pi_shuffling_mask[ii] == 1)
					piGenotypeToShuffle[ii] = RandUtils::iRandNumber(0, 1);
			}//for (int ii = 0; ii < i_templ_length; ii++)

			//now check if you can add it to the population...
			pcIndSupportingOptOrder->d_compute_fitness_at_level_zero(0, piGenotypeToShuffle, pi_test_phenotype);

			//check if the phenotype is "addable"
			if (b_can_ind_be_added_at_any_level(pi_test_phenotype, true, true) == true)
			{
				delete  pi_test_phenotype;
				delete  pi_shuffling_mask;
				return(true);
			}//if (b_can_ind_be_added_at_any_level(pi_test_phenotype, true, true) == true)
			
		}//if (pc_pattern_part != NULL)
	}//while  (b_done ==  false)



	delete  pi_test_phenotype;
	delete  pi_shuffling_mask;
	return(false);
}//bool  C3LO::b_shuffle_the_genotype(int  *piGenotypeToShuffle, C3LOIndividual  *pcCandidate, int iPatternLevel)



void  C3LOSingle::v_add_ind_to_process_in_next_generation(C3LOIndividual  *pcIndividualToProcessInNextGeneration)
{
	for (int ii = 0; ii < v_individuals_to_process_next_gen.size(); ii++)
	{
		if (v_individuals_to_process_next_gen.at(ii) == pcIndividualToProcessInNextGeneration)  return;
	}//for (int ii = 0; ii < v_individuals_to_process_next_gen.size(); ii++)
	
	v_individuals_to_process_next_gen.push_back(pcIndividualToProcessInNextGeneration);
}//void  C3LO::v_add_ind_to_process_in_next_generation(C3LOIndividual  *pcParentMain)





int  C3LOSingle::i_get_pattern_in_touch_in_differences(int  *piGeneGroup, int  *piDifferenceMask, vector<C3LOPattern *>  *pvPatternPool)
{
	vector<int>  v_patterns_in_touch_offsets;

	for (int ii = 0; ii < pvPatternPool->size(); ii++)
	{
		//if (pvPatternPool->at(ii)->bDoITouchGeneGroupWithMask(piGeneGroup, piDifferenceMask) == true)
		if (pvPatternPool->at(ii)->bDoIExtendGeneGroupWithMask(piGeneGroup, piDifferenceMask) == true)
		{
			if (pvPatternPool->at(ii)->bDoIExtendGeneGroupWithMask(piGeneGroup, piDifferenceMask) == true)
				v_patterns_in_touch_offsets.push_back(ii);
		}//if (pvPatternPool->at(ii)->bDoITouchGeneGroupWithMask(piGeneGroup, piDifferenceMask) == true)
	}//for (int ii = 0; ii < pvPatternPool->size(); ii++)

	if (v_patterns_in_touch_offsets.size() > 0)  return(v_patterns_in_touch_offsets.at(RandUtils::iRandNumber(0, v_patterns_in_touch_offsets.size() - 1)));

	return(-1);
}//int  C3LO::i_get_pattern_in_touch_in_differences(int  *piGeneGroup, int  *piDifferenceMask, vector<C3LOPattern *>  *pvPatternPool)



int  C3LOSingle::i_get_pattern_in_touch_if_possible(int  *piGeneGroup, vector<C3LOPattern *>  *pvPatternPool)
{
	vector<int>  v_patterns_in_touch_offsets;

	for (int ii = 0; ii < pvPatternPool->size(); ii++)
	{
		if (pvPatternPool->at(ii)->bDoITouchGeneGroup(piGeneGroup) == true)
			v_patterns_in_touch_offsets.push_back(ii);
	}//for (int ii = 0; ii < pvPatternPool->size(); ii++)

	if (v_patterns_in_touch_offsets.size() > 0)  return(RandUtils::iRandNumber(0, v_patterns_in_touch_offsets.size() - 1));

	return(RandUtils::iRandNumber(0, pvPatternPool->size() - 1));
}//int  C3LO::i_get_pattern_in_touch_if_possible(int  *piGeneGroup, vector<C3LOPattern *>  *pvPatternPool)



int  C3LOSingle::i_get_next_gene_offset(int iBaseGeneOffset, int  **piGeneRelations, int  *piIncludedGenesMask)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGeneRelations[iBaseGeneOffset][ii] < 0)  return(-1);

		if (piGeneRelations[iBaseGeneOffset][ii] >= 0)
		{
			if (piIncludedGenesMask[piGeneRelations[iBaseGeneOffset][ii]] == 0)  return(piGeneRelations[iBaseGeneOffset][ii]);
		}//if (piGeneRelations[iBaseGeneOffset][ii] >= 0)

	}//for (int ii = 0; ii < i_templ_length; ii++)

	return(-1);
}//int  C3LO::i_get_next_gene_offset(int iBaseGeneOffset, int  **piGeneRelations, int  *piIncludedGenesMask)




int  C3LOSingle::i_cross_dsm_glueing_similarities_from_different_start_genes(int  **piGeneRelations, int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)
{
	CString  s_buf, s_info;
	double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();


	int  *pi_other_parent_phenotype;
	pi_other_parent_phenotype = pcParentOther->piGetPhenotype(iPatternLevel, 0);

	
	vector<int>  v_differencing_genes;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGenotype0[ii] != pi_other_parent_phenotype[ii])  v_differencing_genes.push_back(ii);
	}//for (int ii = 0; ii < i_templ_length; ii++)


	if (v_differencing_genes.size() == 0)  return(0);

	//shuffle the genes
	int  i_buf;
	vector<int>  v_differencing_genes_shuffled;
	while (v_differencing_genes.size() > 0)
	{
		i_buf = RandUtils::iRandNumber(0, v_differencing_genes.size() - 1);
		v_differencing_genes_shuffled.push_back(v_differencing_genes.at(i_buf));
		v_differencing_genes.erase(v_differencing_genes.begin() + i_buf);
	}//while (v_differencing_genes.size() > 0)

	int  i_differences_at_start;
	i_differences_at_start = v_differencing_genes_shuffled.size();


	int  *pi_included_genes_mask;
	pi_included_genes_mask = new int[i_templ_length];

	C3LOPattern c_current_pattern(i_templ_length);

	int  i_crossing_result;
	int  i_start_gene_offset, i_next_gene_offset;
	bool  b_no_genes_to_add;
	int  i_current_length;
	i_current_length = 0;

	//FILE  *pf_test;
	
	while (v_differencing_genes_shuffled.size() > 0)
	{
		i_current_length++;

		/*s_buf.Format("zzz_patterns_%.3d.txt", i_current_length);
		pf_test = fopen(s_buf, "w+");
		fprintf(pf_test, "CURRENT LENGTH: %d \n\n", i_current_length);
		fprintf(pf_test, "STARTING GENES SHUFFLED:");

		s_info = "";
		for (int ii = 0; ii < v_differencing_genes_shuffled.size(); ii++)
		{
			s_buf.Format("%d ", v_differencing_genes_shuffled.at(ii));
			s_info += s_buf;
		}//for (int ii = 0; ii < v_differencing_genes_shuffled.size(); ii++)
		s_info += "\n\n\n";
		fprintf(pf_test, s_info);*/


		if ( (i_current_length >= i_templ_length)||(i_current_length > i_relation_length_max + 1))  v_differencing_genes_shuffled.clear();



		for (int i_start = 0; i_start < v_differencing_genes_shuffled.size(); i_start++)
		{
			i_start_gene_offset = v_differencing_genes_shuffled.at(i_start);

			c_current_pattern.v_pattern.clear();
			for (int ii = 0; ii < i_templ_length; ii++)
				pi_included_genes_mask[ii] = 0;

			pi_included_genes_mask[i_start_gene_offset] = 1;
			c_current_pattern.v_pattern.push_back(CMessyGene(1, i_start_gene_offset));

			b_no_genes_to_add = false;

			while ((c_current_pattern.v_pattern.size() < i_current_length) && (b_no_genes_to_add == false))
			{
				i_next_gene_offset = i_get_next_gene_offset(c_current_pattern.v_pattern.at(c_current_pattern.v_pattern.size() - 1).iGenePos(), piGeneRelations, pi_included_genes_mask);

				if (i_next_gene_offset <= 0)
				{
					b_no_genes_to_add = true;
					v_differencing_genes_shuffled.erase(v_differencing_genes_shuffled.begin() + i_start);
					i_start--;
				}//if (i_next_gene_offset <= 0)
				else
				{
					pi_included_genes_mask[i_next_gene_offset] = 1;
					c_current_pattern.v_pattern.push_back(CMessyGene(1, i_next_gene_offset));
				}//else  if (i_next_gene_offset <= 0)

			}//while ((c_current_pattern.v_pattern.size() < i_current_length) && (b_no_genes_to_add == false))


			if (c_current_pattern.v_pattern.size() == i_current_length)
			{
				/*s_info = "";
				for (int ii = 0; ii < c_current_pattern.v_pattern.size(); ii++)
				{
					s_buf.Format("%d ", c_current_pattern.v_pattern.at(ii).iGenePos());
					s_info += s_buf;
				}//for (int ii = 0; ii < v_differencing_genes_shuffled.size(); ii++)
				s_info += "\n";

				fprintf(pf_test, s_info);*/

				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, &c_current_pattern, true, iPatternLevel);
				if (i_crossing_result == 2)
				{
					double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
					s_buf.Format("Size:%d FFEchange: %.2lf Differences available (now/atstart): %d/%d max length: %d", c_current_pattern.v_pattern.size(), d_ffe_after - d_ffe_before, v_differencing_genes_shuffled.size(), i_differences_at_start, i_relation_length_max);
					//pc_log->vPrintLine(s_buf, true);

					delete  pi_included_genes_mask;
					return(2);
				}//if (b_crossing_result == true)
			}//if (c_current_pattern.v_pattern.size() == i_current_length)


		}//for (int i_start = 0; i_start < v_differencing_genes_shuffled.size(); i_start++)

		//fclose(pf_test);
	}//while (v_differencing_genes_shuffled.size() > 0)
	
	delete  pi_included_genes_mask;
	return(0);
	 
	 //::MessageBox(NULL, "", "", MB_OK);
		

	/*	



		pc_pattern_to_add = v_allowed_intercepting_patterns_scraps.at(i_pattern_offset);
		v_allowed_intercepting_patterns_scraps.erase(v_allowed_intercepting_patterns_scraps.begin() + i_pattern_offset);
		b_at_least_one_new_gene_marked = false;
		i_patterns_used++;


		if (v_allowed_intercepting_patterns_scraps.size() == 0)
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, true, iPatternLevel);
		else
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, false, iPatternLevel);

		if (i_crossing_result == 2)
		{
			//fclose(pf_test);
			//::MessageBox(NULL, "end", "end", MB_OK);
			double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
			//s_buf.Format("PAT_TYPE:pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
			s_buf.Format("Size:%d FFEchange: %.2lf GLued scraps: %d/%d", pc_pattern_to_add->v_pattern.size(), d_ffe_after - d_ffe_before, i_patterns_used, v_allowed_intercepting_patterns_scraps.size());
			pc_log->vPrintLine(s_buf, true);

			delete  pi_difference_mask;
			return(2);
		}//if (b_crossing_result == true)

	}//while (v_allowed_based_patterns.size() > 0)

	 //fclose(pf_test);
	 //::MessageBox(NULL, "end", "end", MB_OK);


	delete  pi_difference_mask;*/
	return(0);
}//int  C3LO::i_cross_dsm_glueing_similarities_from_different_start_genes(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)




int  C3LOSingle::i_cross_dsm(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)
{
	if (v_gene_patterns_parts_interceptions_and_differencies.size() <= iPatternLevel)  return(0);

	CString  s_buf;
	double  d_ffe_before, d_time_before;
	double  d_ffe_after, d_time_after;
	d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	c_time_counter.bGetTimePassed(&d_time_before);


	int  *pi_difference_mask;
	int  *pi_other_parent_phenotype;
	pi_difference_mask = new int[i_templ_length];

	bool  b_any_differences;
	b_any_differences = false;
	pi_other_parent_phenotype = pcParentOther->piGetPhenotype(0, 0);
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGenotype0[ii] == pi_other_parent_phenotype[ii])
			pi_difference_mask[ii] = 0;
		else
		{
			pi_difference_mask[ii] = 1;
			b_any_differences = true;
		}//else  if (pi_difference_mask[ii] == pi_other_parent_phenotype[ii])
	}//for (int ii = 0; ii < i_templ_length; ii++)


	if (b_any_differences == false)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (b_any_differences == false)


	vector<C3LOPattern *>  v_allowed_intercepting_patterns_scraps;
	v_allowed_intercepting_patterns_scraps.reserve(v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel).size());
	v_get_all_intercepting_patterns_marking_the_differences(pi_difference_mask, &v_allowed_intercepting_patterns_scraps, iPatternLevel);


	C3LOPattern *pc_pattern_to_add;
	int  i_pattern_offset;
	int  i_crossing_result;
	bool  b_at_least_one_new_gene_marked;


	//d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	//c_time_counter.bGetTimePassed(&d_time_before);

	int  i_patterns_used = 0;
	while (v_allowed_intercepting_patterns_scraps.size() > 0)
	{

		//i_pattern_offset = RandUtils::iRandNumber(0, v_allowed_intercepting_patterns_scraps.size() - 1);
		i_pattern_offset = 0;


		pc_pattern_to_add = v_allowed_intercepting_patterns_scraps.at(i_pattern_offset);
		v_allowed_intercepting_patterns_scraps.erase(v_allowed_intercepting_patterns_scraps.begin() + i_pattern_offset);
		b_at_least_one_new_gene_marked = false;
		i_patterns_used++;


		if (v_allowed_intercepting_patterns_scraps.size() == 0)
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, true, iPatternLevel);
		else
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, false, iPatternLevel);

		if (i_crossing_result == 2)
		{
			//fclose(pf_test);
			//::MessageBox(NULL, "end", "end", MB_OK);
			d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
			c_time_counter.bGetTimePassed(&d_time_after);
			//s_buf.Format("PAT_TYPE:pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());

			double  d_ffe_per_sec;
			if (d_time_after - d_time_before > 0)
				d_ffe_per_sec = (d_ffe_after - d_ffe_before) / (d_time_after - d_time_before);
			else
				d_ffe_per_sec = 0;

			s_buf.Format
				(
				"Size:%d FFEchange: %.2lf Time change: %.2lf FFE/sec: %.4lf Scraps (used/left): %d/%d Fitness: %.4lf", 
				pc_pattern_to_add->v_pattern.size(), 
					d_ffe_after - d_ffe_before, d_time_after - d_time_before, d_ffe_per_sec,
					i_patterns_used, v_allowed_intercepting_patterns_scraps.size(),
					pcParentMain->dComputeFitness(-1)
				);
			pc_log->vPrintLine(s_buf, true);

			delete  pi_difference_mask;
			return(2);
		}//if (b_crossing_result == true)

	}//while (v_allowed_based_patterns.size() > 0)

	 //fclose(pf_test);
	 //::MessageBox(NULL, "end", "end", MB_OK);

	d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
	c_time_counter.bGetTimePassed(&d_time_after);
	//s_buf.Format("PAT_TYPE:pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());

	double  d_ffe_per_sec;
	if (d_time_after - d_time_before > 0)
		d_ffe_per_sec = (d_ffe_after - d_ffe_before) / (d_time_after - d_time_before);
	else
		d_ffe_per_sec = 0;

	s_buf.Format
	(
		"NO UPDATE FFEchange: %.2lf Time change: %.2lf FFE/sec: %.4lf Scraps (used/left): %d/%d",
		d_ffe_after - d_ffe_before, d_time_after - d_time_before, d_ffe_per_sec,
		i_patterns_used, v_allowed_intercepting_patterns_scraps.size()
	);
	//pc_log->vPrintLine(s_buf, true);


	delete  pi_difference_mask;
	return(0);
}//int  C3LO::i_cross_dsm(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)




int  C3LOSingle::i_cross_intercept_scraps(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)
{
	CString  s_buf;
	double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();


	int  *pi_difference_mask;
	int  *pi_other_parent_phenotype;
	pi_difference_mask = new int[i_templ_length];

	bool  b_any_differences;
	b_any_differences = false;
	pi_other_parent_phenotype = pcParentOther->piGetPhenotype(iPatternLevel, 0);
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGenotype0[ii] == pi_other_parent_phenotype[ii])
			pi_difference_mask[ii] = 0;
		else
		{
			pi_difference_mask[ii] = 1;
			b_any_differences = true;
		}//else  if (pi_difference_mask[ii] == pi_other_parent_phenotype[ii])
	}//for (int ii = 0; ii < i_templ_length; ii++)


	if (b_any_differences == false)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (b_any_differences == false)


	vector<C3LOPattern *>  v_allowed_intercepting_patterns_scraps;
	v_get_all_intercepting_patterns_marking_the_differences(pi_difference_mask, &v_allowed_intercepting_patterns_scraps, iPatternLevel);


	C3LOPattern *pc_pattern_to_add;
	int  i_pattern_offset;
	int  i_crossing_result;
	bool  b_at_least_one_new_gene_marked;

	int  i_patterns_used = 0;
	while (v_allowed_intercepting_patterns_scraps.size() > 0)
	{

		/*i_pattern_offset = i_get_pattern_in_touch_in_differences(pi_pattern_sum, pi_difference_mask, &v_allowed_original_patterns_scraps);
		if (i_pattern_offset < 0)
		{
			i_pattern_offset = RandUtils::iRandNumber(0, v_allowed_original_patterns_scraps.size() - 1);

			for (int ii = 0; ii < i_templ_length; ii++)
				pi_pattern_sum[ii] = 0;
			c_pattern_current.v_pattern.clear();
			i_glued_patterns = 1;
		}//if (i_pattern_offset < 0)

		if (i_glued_patterns != i_used_patterns)
			s_patterns_divided = "DIVIDED";
		else
			s_patterns_divided = "";*/

		i_pattern_offset = RandUtils::iRandNumber(0, v_allowed_intercepting_patterns_scraps.size() - 1);


		pc_pattern_to_add = v_allowed_intercepting_patterns_scraps.at(i_pattern_offset);
		v_allowed_intercepting_patterns_scraps.erase(v_allowed_intercepting_patterns_scraps.begin() + i_pattern_offset);
		b_at_least_one_new_gene_marked = false;
		i_patterns_used++;


		/*//pattern_sum - pattern_to_add
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = pi_pattern_sum[ii];*/

		/*b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 1)  b_difference_valid = true;

			pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 0;
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)*/

		


		if (v_allowed_intercepting_patterns_scraps.size() == 0)
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, true, iPatternLevel);
		else
			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, NULL, NULL, pcParentMain, pcParentOther, pc_pattern_to_add, false, iPatternLevel);

		if (i_crossing_result == 2)
		{
			//fclose(pf_test);
			//::MessageBox(NULL, "end", "end", MB_OK);
			double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
			//s_buf.Format("PAT_TYPE:pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
			s_buf.Format("Size:%d FFEchange: %.2lf GLued scraps: %d/%d", pc_pattern_to_add->v_pattern.size(), d_ffe_after - d_ffe_before, i_patterns_used, v_allowed_intercepting_patterns_scraps.size());
			pc_log->vPrintLine(s_buf, true);

			delete  pi_difference_mask;
			return(2);
		}//if (b_crossing_result == true)

	}//while (v_allowed_based_patterns.size() > 0)

	 //fclose(pf_test);
	 //::MessageBox(NULL, "end", "end", MB_OK);


	delete  pi_difference_mask;
	return(0);
}//int  C3LO::i_cross_intercept_scraps(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)





int  C3LOSingle::i_cross_fluent_scraps(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)
{
	CString  s_buf;
	double  d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();


	int  *pi_difference_mask;
	int  *pi_other_parent_phenotype;
	pi_difference_mask = new int[i_templ_length];
	
	bool  b_any_differences;
	b_any_differences = false;
	pi_other_parent_phenotype = pcParentOther->piGetPhenotype(iPatternLevel, 0);
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGenotype0[ii] == pi_other_parent_phenotype[ii])
			pi_difference_mask[ii] = 0;
		else
		{
			pi_difference_mask[ii] = 1;
			b_any_differences = true;
		}//else  if (pi_difference_mask[ii] == pi_other_parent_phenotype[ii])
	}//for (int ii = 0; ii < i_templ_length; ii++)


	if (b_any_differences == false)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (b_any_differences == false)


	vector<C3LOPattern *>  v_allowed_original_patterns_scraps;
	v_get_all_original_pattern_parts_marking_the_differences(pi_difference_mask, &v_allowed_original_patterns_scraps, iPatternLevel);

	if (v_allowed_original_patterns_scraps.size() == 0)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (v_allowed_original_patterns_scraps.size() == 0)



	bool  b_difference_valid;
	int  *pi_pattern_sum;
	pi_pattern_sum = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_pattern_sum[ii] = 0;

	int  *pi_pattern_difference_tool;
	pi_pattern_difference_tool = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_pattern_difference_tool[ii] = 0;

	C3LOPattern  c_pattern_difference(i_templ_length);
	C3LOPattern  c_pattern_current(i_templ_length);
	bool  b_at_least_one_new_gene_marked;
	C3LOPattern *pc_pattern_to_add;
	int  i_pattern_offset;
	int  i_crossing_result;

	int  i_used_patterns = 0;
	int  i_glued_patterns = 0;

	//FILE  *pf_test;
	//pf_test = fopen("zzz__test.txt", "w+");
	//::MessageBox(NULL, "start", "start", MB_OK);
	CString  s_patterns_divided;

	while (v_allowed_original_patterns_scraps.size() > 0)
	{
		i_glued_patterns++;
		i_used_patterns++;

		
		i_pattern_offset = i_get_pattern_in_touch_if_possible(pi_pattern_sum, &v_allowed_original_patterns_scraps);
		//i_pattern_offset = i_get_pattern_in_touch_in_differences(pi_pattern_sum, pi_difference_mask, &v_allowed_original_patterns_scraps);
		if (i_pattern_offset < 0)
		{
			i_pattern_offset = RandUtils::iRandNumber(0, v_allowed_original_patterns_scraps.size() - 1);

			for (int ii = 0; ii < i_templ_length; ii++)
				pi_pattern_sum[ii] = 0;
			c_pattern_current.v_pattern.clear();
			i_glued_patterns = 1;
		}//if (i_pattern_offset < 0)

		if  (i_glued_patterns != i_used_patterns)
			s_patterns_divided = "DIVIDED";
		else
			s_patterns_divided = "";
			

		pc_pattern_to_add = v_allowed_original_patterns_scraps.at(i_pattern_offset);
		v_allowed_original_patterns_scraps.erase(v_allowed_original_patterns_scraps.begin() + i_pattern_offset);
		b_at_least_one_new_gene_marked = false;

		
		//pattern_sum - pattern_to_add
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = pi_pattern_sum[ii];

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if  (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 1)  b_difference_valid = true;

			pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 0;
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			b_difference_valid = false;
			
			for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)
			{
				if  (pi_pattern_difference_tool[ii] == 1)  b_difference_valid = true;
			}//for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)

			if (b_difference_valid == true)
			{
				c_pattern_difference.v_pattern.clear();
				for (int ii = 0; ii < i_templ_length; ii++)
				{
					if (pi_pattern_difference_tool[ii] == 1)
						c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
				}//for (int ii = 0; ii < i_templ_length; ii++)

				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

				if (i_crossing_result == 2)
				{

					double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
					//s_buf.Format("pattern_sum - pattern_to_add  FFEchange: %.2lf =(%.2lf-%.2lf) GLued pattern scraps: %d/%d", d_ffe_after - d_ffe_before, d_ffe_after, d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
					s_buf.Format("PAT_TYPE:1 pattern_sum - pattern_to_add  Size:%d FFEchange: %.2lf GLued scraps: %d/%d (used:%d) %s", c_pattern_difference.v_pattern.size(), d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size(), i_used_patterns, s_patterns_divided);
					pc_log->vPrintLine(s_buf, true);

					//fclose(pf_test);
					//::MessageBox(NULL, "end", "end", MB_OK);

					delete  pi_difference_mask;
					delete  pi_pattern_sum;
					delete  pi_pattern_difference_tool;
					return(2);
				}//if (b_crossing_result == true)
			}//if (b_difference_valid == true)
		}//if (b_difference_valid == true)




		//pattern_to_add - pattern_sum
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = 0;

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
			{
				b_difference_valid = true;
				pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
			}//if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			b_difference_valid = false;

			for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)
			{
				if (pi_pattern_difference_tool[ii] == 1)  b_difference_valid = true;
			}//for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)

			if (b_difference_valid == true)
			{
				c_pattern_difference.v_pattern.clear();
				for (int ii = 0; ii < i_templ_length; ii++)
				{
					if (pi_pattern_difference_tool[ii] == 1)
						c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
				}//for (int ii = 0; ii < i_templ_length; ii++)

				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

				if (i_crossing_result == 2)
				{
					double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
					//s_buf.Format("PAT_TYPE:pattern_to_add - pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
					s_buf.Format("PAT_TYPE:2 pattern_to_add - pattern_sum  Size:%d FFEchange: %.2lf GLued scraps: %d/%d (used:%d) %s", c_pattern_difference.v_pattern.size(), d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size(), i_used_patterns, s_patterns_divided);
					pc_log->vPrintLine(s_buf, true);

					//fclose(pf_test);
					//::MessageBox(NULL, "end", "end", MB_OK);

					delete  pi_difference_mask;
					delete  pi_pattern_sum;
					delete  pi_pattern_difference_tool;
					return(2);
				}//if (b_crossing_result == true)
			}//if (b_difference_valid == true)
		}//if (b_difference_valid == true)




		//pattern_to_add AND pattern_sum
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = 0;

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 1)
			{
				b_difference_valid = true;
				pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
			}//if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			c_pattern_difference.v_pattern.clear();
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				if (pi_pattern_difference_tool[ii] == 1)
					c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
			}//for (int ii = 0; ii < i_templ_length; ii++)

			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

			if (i_crossing_result == 2)
			{
				//fclose(pf_test);
				//::MessageBox(NULL, "end", "end", MB_OK);
				double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
				//s_buf.Format("PAT_TYPE:pattern_to_add AND pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
				s_buf.Format("PAT_TYPE:3 pattern_to_add AND pattern_sum  Size:%d FFEchange: %.2lf GLued scraps: %d/%d (used:%d) %s", c_pattern_difference.v_pattern.size(), d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size(), i_used_patterns, s_patterns_divided);
				pc_log->vPrintLine(s_buf, true);

				delete  pi_difference_mask;
				delete  pi_pattern_sum;
				delete  pi_pattern_difference_tool;
				return(2);
			}//if (b_crossing_result == true)
		}//if (b_difference_valid == true)







		

		//finally sum...
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
			{
				b_at_least_one_new_gene_marked = true;
				c_pattern_current.v_pattern.push_back(CMessyGene(1, pc_pattern_to_add->v_pattern.at(ii).iGenePos()));
			}//if (pi_current_pattern[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)

			pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if ( (b_at_least_one_new_gene_marked == true)||(v_allowed_original_patterns_scraps.size() == 0) )
		{
			//c_current_pattern.vSavePattern(pf_test, NULL);

			if (v_allowed_original_patterns_scraps.size() == 0)
				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_current, true, iPatternLevel);
			else
				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_current, false, iPatternLevel);

			if (i_crossing_result == 2)
			{
				//fclose(pf_test);
				//::MessageBox(NULL, "end", "end", MB_OK);
				double  d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
				//s_buf.Format("PAT_TYPE:pattern_sum  FFEchange: %.2lf GLued scraps: %d/%d", d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size());
				s_buf.Format("PAT_TYPE:4 pattern_sum  Size:%d FFEchange: %.2lf GLued scraps: %d/%d (used:%d) %s", c_pattern_current.v_pattern.size(), d_ffe_after - d_ffe_before, i_glued_patterns, v_allowed_original_patterns_scraps.size(), i_used_patterns, s_patterns_divided);
				pc_log->vPrintLine(s_buf, true);

				delete  pi_difference_mask;
				delete  pi_pattern_sum;
				delete  pi_pattern_difference_tool;
				return(2);
			}//if (b_crossing_result == true)
		}//if ( (b_at_least_one_new_gene_marked == true)||(v_allowed_original_patterns_scraps.size() == 0) )
		
	}//while (v_allowed_based_patterns.size() > 0)

	//fclose(pf_test);
	//::MessageBox(NULL, "end", "end", MB_OK);
	

	delete  pi_difference_mask;
	delete  pi_pattern_sum;
	delete  pi_pattern_difference_tool;
	return(0);
}//bool  C3LO::b_fluent_scraps_crossing(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)




int  C3LOSingle::i_cross_fluent_scraps_berlin(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)
{
	int  *pi_difference_mask;
	int  *pi_other_parent_phenotype;
	pi_difference_mask = new int[i_templ_length];

	bool  b_any_differences;
	b_any_differences = false;
	pi_other_parent_phenotype = pcParentOther->piGetPhenotype(iPatternLevel, 0);
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGenotype0[ii] == pi_other_parent_phenotype[ii])
			pi_difference_mask[ii] = 0;
		else
		{
			pi_difference_mask[ii] = 1;
			b_any_differences = true;
		}//else  if (pi_difference_mask[ii] == pi_other_parent_phenotype[ii])
	}//for (int ii = 0; ii < i_templ_length; ii++)


	if (b_any_differences == false)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (b_any_differences == false)


	vector<C3LOPattern *>  v_allowed_original_patterns_scraps;
	v_get_all_original_pattern_parts_marking_the_differences(pi_difference_mask, &v_allowed_original_patterns_scraps, iPatternLevel);

	if (v_allowed_original_patterns_scraps.size() == 0)
	{
		delete  pi_difference_mask;
		return(0);
	}//if (v_allowed_original_patterns_scraps.size() == 0)



	bool  b_difference_valid;
	int  *pi_pattern_sum;
	pi_pattern_sum = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_pattern_sum[ii] = 0;

	int  *pi_pattern_difference_tool;
	pi_pattern_difference_tool = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_pattern_difference_tool[ii] = 0;

	C3LOPattern  c_pattern_difference(i_templ_length);
	C3LOPattern  c_pattern_current(i_templ_length);
	bool  b_at_least_one_new_gene_marked;
	C3LOPattern *pc_pattern_to_add;
	int  i_pattern_offset;
	int  i_crossing_result;
	int  i_best_crossing_result;


	//FILE  *pf_test;
	//pf_test = fopen("zzz__test.txt", "w+");
	//::MessageBox(NULL, "start", "start", MB_OK);
	i_best_crossing_result = 0;
	while (v_allowed_original_patterns_scraps.size() > 0)
	{
		i_pattern_offset = i_get_pattern_in_touch_if_possible(pi_pattern_sum, &v_allowed_original_patterns_scraps);
		pc_pattern_to_add = v_allowed_original_patterns_scraps.at(i_pattern_offset);
		v_allowed_original_patterns_scraps.erase(v_allowed_original_patterns_scraps.begin() + i_pattern_offset);
		b_at_least_one_new_gene_marked = false;



		//pattern_sum - pattern_to_add
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = pi_pattern_sum[ii];

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 1)  b_difference_valid = true;

			pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 0;
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			b_difference_valid = false;

			for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)
			{
				if (pi_pattern_difference_tool[ii] == 1)  b_difference_valid = true;
			}//for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)

			if (b_difference_valid == true)
			{
				c_pattern_difference.v_pattern.clear();
				for (int ii = 0; ii < i_templ_length; ii++)
				{
					if (pi_pattern_difference_tool[ii] == 1)
						c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
				}//for (int ii = 0; ii < i_templ_length; ii++)

				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

				if (i_best_crossing_result < i_crossing_result)  i_best_crossing_result = i_crossing_result;
				if (i_crossing_result == 2)
				{
					//fclose(pf_test);
					//::MessageBox(NULL, "end", "end", MB_OK);

					delete  pi_difference_mask;
					delete  pi_pattern_sum;
					delete  pi_pattern_difference_tool;
					return(2);
				}//if (b_crossing_result == true)
			}//if (b_difference_valid == true)
		}//if (b_difference_valid == true)




		 //pattern_to_add - pattern_sum
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = 0;

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
			{
				b_difference_valid = true;
				pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
			}//if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			b_difference_valid = false;

			for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)
			{
				if (pi_pattern_difference_tool[ii] == 1)  b_difference_valid = true;
			}//for (int ii = 0; (ii < i_templ_length) && (b_difference_valid == false); ii++)

			if (b_difference_valid == true)
			{
				c_pattern_difference.v_pattern.clear();
				for (int ii = 0; ii < i_templ_length; ii++)
				{
					if (pi_pattern_difference_tool[ii] == 1)
						c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
				}//for (int ii = 0; ii < i_templ_length; ii++)

				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

				if (i_best_crossing_result < i_crossing_result)  i_best_crossing_result = i_crossing_result;
				if (i_crossing_result == 2)
				{
					//fclose(pf_test);
					//::MessageBox(NULL, "end", "end", MB_OK);

					delete  pi_difference_mask;
					delete  pi_pattern_sum;
					delete  pi_pattern_difference_tool;
					return(2);
				}//if (b_crossing_result == true)
			}//if (b_difference_valid == true)
		}//if (b_difference_valid == true)




		 //pattern_to_add AND pattern_sum
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_pattern_difference_tool[ii] = 0;

		b_difference_valid = false;
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 1)
			{
				b_difference_valid = true;
				pi_pattern_difference_tool[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
			}//if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if (b_difference_valid == true)
		{
			c_pattern_difference.v_pattern.clear();
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				if (pi_pattern_difference_tool[ii] == 1)
					c_pattern_difference.v_pattern.push_back(CMessyGene(1, ii));
			}//for (int ii = 0; ii < i_templ_length; ii++)

			i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_difference, false, iPatternLevel);

			if (i_best_crossing_result < i_crossing_result)  i_best_crossing_result = i_crossing_result;
			if (i_crossing_result == 2)
			{
				//fclose(pf_test);
				//::MessageBox(NULL, "end", "end", MB_OK);

				delete  pi_difference_mask;
				delete  pi_pattern_sum;
				delete  pi_pattern_difference_tool;
				return(2);
			}//if (b_crossing_result == true)
		}//if (b_difference_valid == true)









		 //finally sum...
		for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)
		{
			if (pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)
			{
				b_at_least_one_new_gene_marked = true;
				c_pattern_current.v_pattern.push_back(CMessyGene(1, pc_pattern_to_add->v_pattern.at(ii).iGenePos()));
			}//if (pi_current_pattern[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] == 0)

			pi_pattern_sum[pc_pattern_to_add->v_pattern.at(ii).iGenePos()] = 1;
		}//for (int ii = 0; ii < pc_pattern_to_add->v_pattern.size(); ii++)

		if ((b_at_least_one_new_gene_marked == true) || (v_allowed_original_patterns_scraps.size() == 0))
		{
			//c_current_pattern.vSavePattern(pf_test, NULL);

			if (v_allowed_original_patterns_scraps.size() == 0)
				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_current, true, iPatternLevel);
			else
				i_crossing_result = i_cross_individuals_by_crossing_tree(piGenotype0, piGenotypeToShuffle, piShuffle, pcParentMain, pcParentOther, &c_pattern_current, false, iPatternLevel);

			if (i_best_crossing_result < i_crossing_result)  i_best_crossing_result = i_crossing_result;

			if (i_crossing_result == 2)
			{
				//fclose(pf_test);
				//::MessageBox(NULL, "end", "end", MB_OK);

				delete  pi_difference_mask;
				delete  pi_pattern_sum;
				delete  pi_pattern_difference_tool;
				return(2);
			}//if (b_crossing_result == true)
		}//if ( (b_at_least_one_new_gene_marked == true)||(v_allowed_original_patterns_scraps.size() == 0) )

	}//while (v_allowed_based_patterns.size() > 0)

	 //fclose(pf_test);
	 //::MessageBox(NULL, "end", "end", MB_OK);


	delete  pi_difference_mask;
	delete  pi_pattern_sum;
	delete  pi_pattern_difference_tool;
	return(i_best_crossing_result);



	//new data
	/*if (pc_chosen_pattern != NULL)
	{
	int  i_effect_already_exists_in_the_pop;

	double  d_ind_0_start_fitness;
	d_ind_0_start_fitness = dComputeFitness(piGenotype0);


	//first make a copy for eventual restoration
	for (int ii = 0; ii < i_templ_length; ii++)
	pi_genotype_tool[ii] = piGenotype0[ii];

	v_cross_infect_individual_by_tree(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pc_chosen_pattern);

	double  d_ind_0_end_fitness;
	d_ind_0_end_fitness = dComputeFitness(piGenotype0);
	if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	{
	if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
	{
	pcParentMain->vSetGenotype(piGenotype0);
	pcParentMain->vSetOptimized(false);
	pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);

	pcParentMain->vGenerateOrder();
	pcParentMain->vSetGenotype(pcParentMain->piGetPhenotype(iPatternLevel, 0));
	pcParentMain->vSetOptimized(false);
	pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);


	for (int ii = 0; ii < i_templ_length; ii++)
	piGenotype0[ii] = pcParentMain->piGetPhenotype(iPatternLevel, 0)[ii];

	return(true);
	}//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
	else
	{
	//the we shuffle the individual with the block from ParentOther
	int i_offspring_level = 0;
	if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
	if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
	if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

	//if ofsspring level is higher than parents - make a new individual
	C3LOIndividual *pc_indiv_new;
	pc_indiv_new = new C3LOIndividual(i_templ_length, pc_problem, &v_optimization_orders, this);
	pc_indiv_new->vRandomInit();

	if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new, iPatternLevel) == true)
	{
	//SHUFFLING!
	i_shuffled_individuals++;

	//the we shuffle the individual with the block from ParentOther
	int i_offspring_level = 0;
	if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
	if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
	if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

	pc_indiv_new->vSetGenotype(piGenotype0);
	pc_indiv_new->vSetOptimized(false);
	pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);

	b_add_at_level(pc_indiv_new, i_offspring_level);

	v_add_ind_to_process_in_next_generation(pc_indiv_new);
	}//if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new) == true)
	else
	delete pc_indiv_new;

	//takeback first
	for (int ii = 0; ii < i_templ_length; ii++)
	piGenotype0[ii] = pi_genotype_tool[ii];
	}//else  if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
	}//if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	else
	{
	//if  (d_ind_0_end_fitness == d_ind_0_start_fitness) DO NOTHING!
	if (d_ind_0_end_fitness < d_ind_0_start_fitness)
	{
	//takeback first
	for (int ii = 0; ii < i_templ_length; ii++)
	piGenotype0[ii] = pi_genotype_tool[ii];
	}//if (d_ind_0_end_fitness < d_ind_0_start_fitness)
	}//else  if (d_ind_0_end_fitness > d_ind_0_start_fitness)



	double d_ind_1_start_fitness, d_ind_1_end_fitness;
	d_ind_1_start_fitness = dComputeFitness(pcParentOther->piGetPhenotype(iPatternLevel, 0));

	//and a copy for eventual restoration
	for (int ii = 0; ii < i_templ_length; ii++)
	pi_genotype_tool[ii] = pcParentOther->piGetPhenotype(iPatternLevel, 0)[ii];

	v_cross_infect_individual_by_tree(pi_genotype_tool, piGenotype0, pc_chosen_pattern);

	d_ind_1_end_fitness = dComputeFitness(pi_genotype_tool);

	if (d_ind_1_end_fitness > d_ind_1_start_fitness)
	{
	//first check if we have something new...
	if (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)
	{
	i_improved_other++;
	//the ParentOther has improved - compute what is the improved parent version level
	int i_offspring_level = 0;
	if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
	if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
	if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

	//if ofsspring level is the same as parents - replace parent
	if (i_offspring_level == pcParentOther->i_level)
	{
	pcParentOther->vSetGenotype(pi_genotype_tool);
	pcParentOther->vSetOptimized(false);
	pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

	pcParentOther->vSetGenotype(pcParentOther->piGetPhenotype(iPatternLevel, 0));
	pcParentOther->vGenerateOrder();
	pcParentOther->vSetOptimized(false);
	pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

	v_add_ind_to_process_in_next_generation(pcParentOther);
	//v_individuals_to_process_next_gen.push_back(pcParentOther);
	}//if (i_offspring_level == pcParentOther->i_level)
	else
	{
	//if ofsspring level is higher than parents - make a new individual
	C3LOIndividual *pc_indiv_new;
	pc_indiv_new = new C3LOIndividual(i_templ_length, pc_problem, &v_optimization_orders, this);
	pc_indiv_new->vRandomInit();
	pc_indiv_new->vSetGenotype(pi_genotype_tool);
	pc_indiv_new->vSetOptimized(false);
	pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);


	b_add_at_level(pc_indiv_new, i_offspring_level);

	v_add_ind_to_process_in_next_generation(pc_indiv_new);
	//v_individuals_to_process_next_gen.push_back(pc_indiv_new);


	}//else  if ofsspring level is the same as parents - replace parent


	}//if  (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)

	if (pc_chosen_pattern == pcCrossingTree)//if we are at the top of the tree
	{
	//TURN IT ON
	//if (d_ind_1_end_fitness < d_ind_1_start_fitness)
	//{
	//b_add_higher_level_gene(NULL, pcCrossingTree->iGetPatternLevel() + 1);//create level of high level genes
	//v_insert_gene_pattern_part_level(pcCrossingTree->iGetPatternLevel() + 1);
	//v_add_gene_pattern_tree(NULL, pcCrossingTree->iGetPatternLevel() + 1);
	//v_update_high_lvl_genes_and_trees();
	//}//if  (
	}//if (pc_chosen_pattern == pcCrossingTree)//if we are at the top of the tree
	}//if  ((*pcOffspring)->dComputeFitness() < pcInd0->dComputeFitness())

	}//if (pc_chosen_pattern != NULL)


	return(false);*/
}//bool  C3LO::b_fluent_scraps_crossing(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, int  iPatternLevel)



int  C3LOSingle::i_check_number_of_building_block_occurences(int  *piGenotype0, C3LOPattern *pcCrossingTree, int iPatternLevel)
{
	/*FILE  *pf_test;
	pf_test = fopen("zzzzz_bb_test.txt", "a");

	fprintf(pf_test, "controlled genotype: \n");
	for (int ii = 0; ii < i_templ_length; ii++)
		fprintf(pf_test, "%d", piGenotype0[ii]);
	fprintf(pf_test, "\n");*/
	

	//pcCrossingTree->vSavePattern(pf_test, NULL, "", true);

	int  i_building_block_occurences_counter;

	i_building_block_occurences_counter = 0;

	if (pv_other_pops != NULL)
	{
		for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
		{
			if (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPatternLevel)
			{
				for (int i_ind = 0; i_ind < pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).size(); i_ind++)
				{
					//pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).at(i_ind)->vSave(pf_test, NULL);


					if (
						pcCrossingTree->bDoesIndividualContainMarkedBlock
						(
							piGenotype0,
							pv_other_pops->at(i_3lo_pop)->v_population_levels.at(iPatternLevel).at(i_ind)
						)
						== true
						)
					{
						//fprintf(pf_test, "true\n");
						i_building_block_occurences_counter++;
					}//if  {
					//else
						//fprintf(pf_test, "false\n");
				}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
			}//if  (pv_other_pops->at(i_3lo_pop)->v_population_levels.size() > iPatternLevel)
		}//for (int i_3lo_pop = 0; i_3lo_pop < pv_other_pops->size(); i_3lo_pop++)
	}//if (pv_other_pops != NULL)


	if (v_population_levels.size() > iPatternLevel)
	{
		for (int i_ind = 0; i_ind < v_population_levels.at(iPatternLevel).size(); i_ind++)
		{
			//v_population_levels.at(iPatternLevel).at(i_ind)->vSave(pf_test, NULL);

			if (
				pcCrossingTree->bDoesIndividualContainMarkedBlock
				(
					piGenotype0,
					v_population_levels.at(iPatternLevel).at(i_ind)
				)
				== true
				)
			{
				//fprintf(pf_test, "true\n");
				i_building_block_occurences_counter++;
			}//if  (
			//else
				//fprintf(pf_test, "false\n");


		}//for (int i_ind = 0; i_ind < v_population_levels->v_population_levels.at(iPatternLevel).size(); i_ind++)
	}//if (v_population_levels.size() > iPatternLevel)


	//fclose(pf_test);
	return(i_building_block_occurences_counter);
}//int  C3LOSingle::i_check_number_of_building_block_occurences(int  *piGenotype0, C3LOPattern *pcCrossingTree)



C3LOPattern*  C3LOSingle::pc_get_tree_root_for_pattern(C3LOPattern  *pcCrossingTree, int  iPatternLevel)
{
	if (iPatternLevel >= v_gene_patterns_trees.size())  return(NULL);

	
	for (int i_tree = 0; i_tree < v_gene_patterns_trees.at(iPatternLevel).size(); i_tree++)
	{
		v_gene_patterns_trees.at(iPatternLevel).at(i_tree)->vGeneratePatternTable();
		if (pcCrossingTree->bDoITouchGeneGroup(v_gene_patterns_trees.at(iPatternLevel).at(i_tree)->pi_pattern_table) == true)  return(v_gene_patterns_trees.at(iPatternLevel).at(i_tree));
	}//for (int i_tree = 0; i_tree < v_gene_patterns_trees.size(); i_tree++)
	

	return(NULL);
}//C3LOPattern*  C3LOSingle::pc_get_tree_root_for_pattern(C3LOPattern  *pcCrossingTree, int  iPatternLevel)


int  C3LOSingle::i_cross_individuals_by_crossing_tree(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle,  C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, bool bMainGroup, int  iPatternLevel)
{
	/*if (
		(pcCrossingTree->bUnMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pi_genotype_tool) == false) ||
		(pcCrossingTree->bMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0)) == false)
		)
	{
		return(0);
	}// if (*/

	int  *pi_genotype_donator;
	


	int  i_effect_already_exists_in_the_pop;

	double  d_ind_0_start_fitness;
	d_ind_0_start_fitness = dComputeFitness(piGenotype0);

	pi_genotype_donator = pcParentOther->piGetPhenotype(0, 0);

	//first make a copy for eventual restoration
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_tool[ii] = piGenotype0[ii];

	v_cross_infect_individual_by_tree(piGenotype0, pi_genotype_donator, pcCrossingTree);

	double  d_ind_0_end_fitness;
	d_ind_0_end_fitness = dComputeFitness(piGenotype0);
	if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	{
		//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
		if (b_can_ind_be_added_at_pattern_level(pcParentMain->i_level, piGenotype0, true, true) == true)
		{
			pcParentMain->vSetGenotype(piGenotype0);
			pcParentMain->vSetOptimized(false);
			pcParentMain->dComputeFitness(-1);
			/*pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);

			pcParentMain->vGenerateOrder();
			pcParentMain->vSetGenotype(pcParentMain->piGetPhenotype(iPatternLevel, 0));
			pcParentMain->vSetOptimized(false);
			pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);*/

			int  i_level_inner_max;
			i_level_inner_max = pcParentMain->i_level_inner;
			if (i_level_inner_max < pcParentOther->i_level_inner)  i_level_inner_max = pcParentOther->i_level_inner;

			pcParentMain->i_level_inner = i_level_inner_max + 1;


			/*for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pcParentMain->piGetPhenotype(iPatternLevel, 0)[ii];*/

			return(2);
		}//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)

		 //takeback first
		for (int ii = 0; ii < i_templ_length; ii++)
			piGenotype0[ii] = pi_genotype_tool[ii];
	}//if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	else
	{
		if (d_ind_0_end_fitness == d_ind_0_start_fitness) //DO NOTHING! - it preserves the change, but the corssing remains unsuccessful
		{
			if (b_can_ind_be_added_at_pattern_level(pcParentMain->i_level, piGenotype0, true, true) == true)
			{
				/*//nothing changes - try random crossing

				//takeback
				for (int ii = 0; ii < i_templ_length; ii++)
					piGenotype0[ii] = pi_genotype_tool[ii];

				vector<CMessyGene> v_differences_marked_by_pattern;
				

				for (int i_gene = 0; i_gene < pcCrossingTree->pvGetPattern()->size(); i_gene++)
				{
					if (piGenotype0[pcCrossingTree->pvGetPattern()->at(i_gene).iGenePos()] != pi_genotype_donator[pcCrossingTree->pvGetPattern()->at(i_gene).iGenePos()])
						v_differences_marked_by_pattern.push_back(CMessyGene(1, pcCrossingTree->pvGetPattern()->at(i_gene).iGenePos()));
				}//for (int i_gene = 0; i_gene < pcCrossingTree->pvGetPattern()->size(); i_gene++)
				

				//now perform uniform crossover on the gene subgroup
				if (v_differences_marked_by_pattern.size() > 0)
				{
					for (int i_diff_gene = 0; i_diff_gene < v_differences_marked_by_pattern.size(); i_diff_gene++)
					{
						if (RandUtils::dRandNumber(0, 1) < 0.5)
						{
							piGenotype0[v_differences_marked_by_pattern.at(i_diff_gene).iGenePos()] = pi_genotype_donator[v_differences_marked_by_pattern.at(i_diff_gene).iGenePos()];
						}//if (RandUtils::dRandNumber(1) < 0.5)
					}//for (int i_diff_gene = 0; i_diff_gene < pcCrossingTree->pvGetPattern()->size(); i_diff_gene++)

					d_ind_0_end_fitness = dComputeFitness(piGenotype0);

					if (d_ind_0_end_fitness < d_ind_0_start_fitness)//it got worse - use previous crossing
						v_cross_infect_individual_by_tree(piGenotype0, pi_genotype_donator, pcCrossingTree);

					if (d_ind_0_end_fitness > d_ind_0_start_fitness)//it has improved! - use it and report success
					{
						pcParentMain->vSetGenotype(piGenotype0);
						pcParentMain->vSetOptimized(false);
						pcParentMain->dComputeFitness(-1);

						return(2);
					}//if (d_ind_0_end_fitness > d_ind_0_start_fitness)//it has improved! - use it and report success

					//if the fitness is the same then just use this genotype version
				}//if (v_differences_marked_by_pattern.size() > 0)
				else
					v_cross_infect_individual_by_tree(piGenotype0, pi_genotype_donator, pcCrossingTree);*/
				



				pcParentMain->vSetGenotype(piGenotype0);
				pcParentMain->vSetOptimized(false);
				pcParentMain->dComputeFitness(-1);
				
				return(0);
			}//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)

			 //takeback first
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pi_genotype_tool[ii];
		}//if  (d_ind_0_end_fitness == d_ind_0_start_fitness) //DO NOTHING! - it preserves the change, but the corssing remains unsuccessful


		if (d_ind_0_end_fitness < d_ind_0_start_fitness)
		{
			//takeback first
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pi_genotype_tool[ii];
		}//if (d_ind_0_end_fitness < d_ind_0_start_fitness)
	}//else  if (d_ind_0_end_fitness > d_ind_0_start_fitness)


	return(0);
}//bool  C3LO::



int  C3LOSingle::i_cross_individuals_by_crossing_tree_berlin(int  *piGenotype0, int  *piGenotypeToShuffle, int  *piShuffle, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, bool bMainGroup, int  iPatternLevel)
{
	if (
		(pcCrossingTree->bUnMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pi_genotype_tool) == false) ||
		(pcCrossingTree->bMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0)) == false)
		)
	{
		return(0);
	}// if (



	int  i_effect_already_exists_in_the_pop;

	double  d_ind_0_start_fitness;
	d_ind_0_start_fitness = dComputeFitness(piGenotype0);


	//first make a copy for eventual restoration
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_tool[ii] = piGenotype0[ii];

	v_cross_infect_individual_by_tree(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pcCrossingTree);

	double  d_ind_0_end_fitness;
	d_ind_0_end_fitness = dComputeFitness(piGenotype0);
	if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	{
		if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
		{
			pcParentMain->vSetGenotype(piGenotype0);
			pcParentMain->vSetOptimized(false);
			pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);

			pcParentMain->vGenerateOrder();
			pcParentMain->vSetGenotype(pcParentMain->piGetPhenotype(iPatternLevel, 0));
			pcParentMain->vSetOptimized(false);
			pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);


			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pcParentMain->piGetPhenotype(iPatternLevel, 0)[ii];

			return(2);
		}//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
		else
		{
			*piShuffle = 1;
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotypeToShuffle[ii] = piGenotype0[ii];

			/*//the we shuffle the individual with the block from ParentOther
			int i_offspring_level = 0;
			if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
			if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
			if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

			//if ofsspring level is higher than parents - make a new individual
			C3LOIndividual *pc_indiv_new;
			pc_indiv_new = new C3LOIndividual(i_templ_length, pc_problem, &v_optimization_orders, this);
			pc_indiv_new->vRandomInit();

			if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new, iPatternLevel) == true)
			{
			//SHUFFLING!
			i_shuffled_individuals++;

			//the we shuffle the individual with the block from ParentOther
			int i_offspring_level = 0;
			if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
			if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
			if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

			pc_indiv_new->vSetGenotype(piGenotype0);
			pc_indiv_new->vSetOptimized(false);
			pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);

			b_add_at_level(pc_indiv_new, i_offspring_level);

			v_add_ind_to_process_in_next_generation(pc_indiv_new);
			}//if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new) == true)
			else
			delete pc_indiv_new;*/

			//takeback first
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pi_genotype_tool[ii];
		}//else  if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
	}//if (d_ind_0_end_fitness > d_ind_0_start_fitness)
	else
	{
		//if  (d_ind_0_end_fitness == d_ind_0_start_fitness) DO NOTHING!
		if (d_ind_0_end_fitness < d_ind_0_start_fitness)
		{
			//takeback first
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotype0[ii] = pi_genotype_tool[ii];
		}//if (d_ind_0_end_fitness < d_ind_0_start_fitness)
	}//else  if (d_ind_0_end_fitness > d_ind_0_start_fitness)



	double d_ind_1_start_fitness, d_ind_1_end_fitness;
	d_ind_1_start_fitness = dComputeFitness(pcParentOther->piGetPhenotype(iPatternLevel, 0));

	//and a copy for eventual restoration
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_tool[ii] = pcParentOther->piGetPhenotype(iPatternLevel, 0)[ii];

	v_cross_infect_individual_by_tree(pi_genotype_tool, piGenotype0, pcCrossingTree);

	d_ind_1_end_fitness = dComputeFitness(pi_genotype_tool);

	if (d_ind_1_end_fitness > d_ind_1_start_fitness)
	{
		//first check if we have something new...
		if (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)
		{
			i_improved_other++;
			//the ParentOther has improved - compute what is the improved parent version level
			int i_offspring_level = 0;
			if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
			if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
			if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

			//if ofsspring level is the same as parents - replace parent
			if (i_offspring_level == pcParentOther->i_level)
			{
				pcParentOther->vSetGenotype(pi_genotype_tool);
				pcParentOther->vSetOptimized(false);
				pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

				pcParentOther->vSetGenotype(pcParentOther->piGetPhenotype(iPatternLevel, 0));
				pcParentOther->vGenerateOrder();
				pcParentOther->vSetOptimized(false);
				pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

				v_add_ind_to_process_in_next_generation(pcParentOther);
				//v_individuals_to_process_next_gen.push_back(pcParentOther);
			}//if (i_offspring_level == pcParentOther->i_level)
			else
			{
				//if ofsspring level is higher than parents - make a new individual
				C3LOIndividual *pc_indiv_new;
				pc_indiv_new = new C3LOIndividual(i_offspring_level, i_templ_length, pc_problem, &v_optimization_orders, this);
				pc_indiv_new->vRandomInit();
				pc_indiv_new->vSetGenotype(pi_genotype_tool);
				pc_indiv_new->vSetOptimized(false);
				pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);


				b_add_at_level(pc_indiv_new, i_offspring_level);

				v_add_ind_to_process_in_next_generation(pc_indiv_new);
				//v_individuals_to_process_next_gen.push_back(pc_indiv_new);
			}//else  if ofsspring level is the same as parents - replace parent

			return(1);
		}//if  (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)
		else
		{
			*piShuffle = 1;
			for (int ii = 0; ii < i_templ_length; ii++)
				piGenotypeToShuffle[ii] = pi_genotype_tool[ii];

			/*//the we shuffle the individual with the block from ParentOther
			int i_offspring_level = 0;
			if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
			if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
			if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

			//if ofsspring level is higher than parents - make a new individual
			C3LOIndividual *pc_indiv_new;
			pc_indiv_new = new C3LOIndividual(i_templ_length, pc_problem, &v_optimization_orders, this);
			pc_indiv_new->vRandomInit();

			if (b_shuffle_the_genotype(pi_genotype_tool, pc_indiv_new, iPatternLevel) == true)
			{
			//SHUFFLING!
			i_shuffled_individuals++;

			//the we shuffle the individual with the block from ParentOther
			int i_offspring_level = 0;
			if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
			if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
			if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

			pc_indiv_new->vSetGenotype(pi_genotype_tool);
			pc_indiv_new->vSetOptimized(false);
			pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);

			b_add_at_level(pc_indiv_new, i_offspring_level);

			v_add_ind_to_process_in_next_generation(pc_indiv_new);
			}//if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new) == true)
			else
			delete pc_indiv_new;*/

		}//else  if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)

		if (bMainGroup == true)//if we are at the top of the tree
		{
			//TURN IT ON
			/*if (d_ind_1_end_fitness < d_ind_1_start_fitness)
			{
			b_add_higher_level_gene(NULL, pcCrossingTree->iGetPatternLevel() + 1);//create level of high level genes
			v_insert_gene_pattern_part_level(pcCrossingTree->iGetPatternLevel() + 1);
			v_add_gene_pattern_tree(NULL, pcCrossingTree->iGetPatternLevel() + 1);
			v_update_high_lvl_genes_and_trees();
			}//if  (*/
		}//if (pc_chosen_pattern == pcCrossingTree)//if we are at the top of the tree
	}//if  ((*pcOffspring)->dComputeFitness() < pcInd0->dComputeFitness())



	return(0);
}//bool  C3LO::


bool  C3LOSingle::b_cross_individuals_by_tree_branch_or_part(int  *piGenotype0, C3LOIndividual  *pcParentMain, C3LOIndividual  *pcParentOther, C3LOPattern *pcCrossingTree, int  iPatternLevel)
{
	C3LOPattern  *pc_chosen_pattern;

	pc_chosen_pattern = NULL;

	if (
		(pcCrossingTree->bUnMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pi_genotype_tool) == true) &&
		(pcCrossingTree->bMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0)) == true)
		)
	{
		pc_chosen_pattern = pcCrossingTree;
	}// if (
	else
	{
		vector<C3LOPattern *>  v_gene_pattern_branch;
		pcCrossingTree->vGetRandomTreeBranch(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), &v_gene_pattern_branch, iPatternLevel);

		int  i_chosen_branch;
		while ((pc_chosen_pattern == NULL) && (v_gene_pattern_branch.size() > 0))
		{
			i_chosen_branch = RandUtils::iRandNumber(0, v_gene_pattern_branch.size() - 1);

			if (
				(v_gene_pattern_branch.at(i_chosen_branch)->bUnMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pi_genotype_tool) == true) &&
				(v_gene_pattern_branch.at(i_chosen_branch)->bMarkedPhenotypeDifferences(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0)) == true)
				)
			{
				pc_chosen_pattern = v_gene_pattern_branch.at(i_chosen_branch);
			}// if (

			v_gene_pattern_branch.erase(v_gene_pattern_branch.begin() + i_chosen_branch);		
		}//while ((pc_chosen_pattern == NULL) && (v_gene_pattern_branch.size() > 0))
	
	}//else  if  (

	
	if (pc_chosen_pattern != NULL)
	{
		int  i_effect_already_exists_in_the_pop;

		double  d_ind_0_start_fitness;
		d_ind_0_start_fitness = dComputeFitness(piGenotype0);
		

		//first make a copy for eventual restoration
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_genotype_tool[ii] = piGenotype0[ii];

		v_cross_infect_individual_by_tree(piGenotype0, pcParentOther->piGetPhenotype(iPatternLevel, 0), pc_chosen_pattern);
	
		double  d_ind_0_end_fitness;
		d_ind_0_end_fitness = dComputeFitness(piGenotype0);
		if (d_ind_0_end_fitness > d_ind_0_start_fitness)
		{
			if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
			{
				pcParentMain->vSetGenotype(piGenotype0);
				pcParentMain->vSetOptimized(false);
				pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);

				pcParentMain->vGenerateOrder();
				pcParentMain->vSetGenotype(pcParentMain->piGetPhenotype(iPatternLevel, 0));
				pcParentMain->vSetOptimized(false);
				pcParentMain->dComputeFitness(v_gene_patterns_trees.size() - 1);


				for (int ii = 0; ii < i_templ_length; ii++)
					piGenotype0[ii] = pcParentMain->piGetPhenotype(iPatternLevel, 0)[ii];

				return(true);
			}//if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
			else
			{
				//the we shuffle the individual with the block from ParentOther
				int i_offspring_level = 0;
				if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
				if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
				if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

				//if ofsspring level is higher than parents - make a new individual
				C3LOIndividual *pc_indiv_new;
				pc_indiv_new = new C3LOIndividual(i_offspring_level, i_templ_length, pc_problem, &v_optimization_orders, this);
				pc_indiv_new->vRandomInit();
				
				if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new, iPatternLevel) == true)
				{
					//SHUFFLING!
					i_shuffled_individuals++;

					//the we shuffle the individual with the block from ParentOther
					int i_offspring_level = 0;
					if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
					if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
					if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

					pc_indiv_new->vSetGenotype(piGenotype0);
					pc_indiv_new->vSetOptimized(false);
					pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);

					b_add_at_level(pc_indiv_new, i_offspring_level);

					v_add_ind_to_process_in_next_generation(pc_indiv_new);
				}//if (b_shuffle_the_genotype(piGenotype0, pc_indiv_new) == true)
				else
					delete pc_indiv_new;

				//takeback first
				for (int ii = 0; ii < i_templ_length; ii++)
					piGenotype0[ii] = pi_genotype_tool[ii];
			}//else  if (b_can_ind_be_added_at_any_level(piGenotype0, true, true) == true)
		}//if (d_ind_0_end_fitness > d_ind_0_start_fitness)
		else
		{
			//if  (d_ind_0_end_fitness == d_ind_0_start_fitness) DO NOTHING!
			if (d_ind_0_end_fitness < d_ind_0_start_fitness)
			{
				//takeback first
				for (int ii = 0; ii < i_templ_length; ii++)
					piGenotype0[ii] = pi_genotype_tool[ii];
			}//if (d_ind_0_end_fitness < d_ind_0_start_fitness)
		}//else  if (d_ind_0_end_fitness > d_ind_0_start_fitness)
		


		double d_ind_1_start_fitness, d_ind_1_end_fitness;
		d_ind_1_start_fitness = dComputeFitness(pcParentOther->piGetPhenotype(iPatternLevel, 0));
				
		//and a copy for eventual restoration
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_genotype_tool[ii] = pcParentOther->piGetPhenotype(iPatternLevel, 0)[ii];

		v_cross_infect_individual_by_tree(pi_genotype_tool, piGenotype0, pc_chosen_pattern);

		d_ind_1_end_fitness = dComputeFitness(pi_genotype_tool);

		if (d_ind_1_end_fitness > d_ind_1_start_fitness)
		{
			//first check if we have something new...
			if (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)
			{
				i_improved_other++;
				//the ParentOther has improved - compute what is the improved parent version level
				int i_offspring_level = 0;
				if (i_offspring_level < pcParentMain->i_level)  i_offspring_level = pcParentMain->i_level;
				if (i_offspring_level < pcParentOther->i_level)  i_offspring_level = pcParentOther->i_level;
				if (i_offspring_level < iPatternLevel)  i_offspring_level = iPatternLevel;

				//if ofsspring level is the same as parents - replace parent
				if (i_offspring_level == pcParentOther->i_level)
				{
					pcParentOther->vSetGenotype(pi_genotype_tool);
					pcParentOther->vSetOptimized(false);
					pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

					pcParentOther->vSetGenotype(pcParentOther->piGetPhenotype(iPatternLevel, 0));
					pcParentOther->vGenerateOrder();
					pcParentOther->vSetOptimized(false);
					pcParentOther->dComputeFitness(v_gene_patterns_trees.size() - 1);

					v_add_ind_to_process_in_next_generation(pcParentOther);
					//v_individuals_to_process_next_gen.push_back(pcParentOther);
				}//if (i_offspring_level == pcParentOther->i_level)
				else
				{
					//if ofsspring level is higher than parents - make a new individual
					C3LOIndividual *pc_indiv_new;
					pc_indiv_new = new C3LOIndividual(i_offspring_level, i_templ_length, pc_problem, &v_optimization_orders, this);
					pc_indiv_new->vRandomInit();
					pc_indiv_new->vSetGenotype(pi_genotype_tool);
					pc_indiv_new->vSetOptimized(false);
					pc_indiv_new->dComputeFitness(v_gene_patterns_trees.size() - 1);


					b_add_at_level(pc_indiv_new, i_offspring_level);

					v_add_ind_to_process_in_next_generation(pc_indiv_new);
					//v_individuals_to_process_next_gen.push_back(pc_indiv_new);
					
					
				}//else  if ofsspring level is the same as parents - replace parent


			}//if  (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == true)

			if (pc_chosen_pattern == pcCrossingTree)//if we are at the top of the tree
			{
				//TURN IT ON
				/*if (d_ind_1_end_fitness < d_ind_1_start_fitness)
				{
					b_add_higher_level_gene(NULL, pcCrossingTree->iGetPatternLevel() + 1);//create level of high level genes
					v_insert_gene_pattern_part_level(pcCrossingTree->iGetPatternLevel() + 1);
					v_add_gene_pattern_tree(NULL, pcCrossingTree->iGetPatternLevel() + 1);
					v_update_high_lvl_genes_and_trees();
				}//if  (*/
			}//if (pc_chosen_pattern == pcCrossingTree)//if we are at the top of the tree
		}//if  ((*pcOffspring)->dComputeFitness() < pcInd0->dComputeFitness())
	
	}//if (pc_chosen_pattern != NULL)


	return(false);
}//bool  C3LO::




bool  C3LOSingle::b_cross_individuals_by_tree_branch(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOPattern *pcCrossingTree, int  iPatternLevel)
{
	//C3LOIndividual  *pc_best_acceptable_offspring;
	//pc_best_acceptable_offspring = NULL;

	//if (bUnMarkedPhenotypeDifferences(pcInd0, pcInd1, piDifferenceSearchTool, iPatternLevel) == true)  return(this);
	if (pcCrossingTree->bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == true)
	{
		C3LOPattern *pc_best_crossing_node;

		vector<C3LOPattern *>  v_gene_pattern_branch;
		pcCrossingTree->vGetRandomTreeBranch(pcInd0->piGetPhenotype(iPatternLevel, 0), pcInd1->piGetPhenotype(iPatternLevel, 0), &v_gene_pattern_branch, iPatternLevel);


		/*int i_chosen_tree;
		i_chosen_tree = pc_random_gen->Next(0, v_allowed_trees.size());*/

		double  d_ind0_fit, d_ind1_fit;
		d_ind0_fit = pcInd0->dComputeFitness(iPatternLevel);
		d_ind1_fit = pcInd1->dComputeFitness(iPatternLevel);


		int  i_effect_already_exists_in_the_pop;
		for (int i_pattern = 0; i_pattern < 1;/*v_gene_pattern_branch.size()*/ i_pattern++)
		{
			b_cross_infect_individual_by_tree(pcInd0, pcInd1, v_gene_pattern_branch.at(i_pattern), iPatternLevel);
			if (pcInd0->dComputeFitness(iPatternLevel) > d_ind0_fit)  return(true);



			b_cross_infect_individual_by_tree(pcInd1, pcInd0, v_gene_pattern_branch.at(i_pattern), iPatternLevel);

			if (i_pattern == 0)//if we are at the top of the tree
			{
				if (
					(
						(pcInd0->dComputeFitness(iPatternLevel) > d_ind0_fit) &&
						(pcInd1->dComputeFitness(iPatternLevel) > d_ind1_fit)
						)
					||
					(
						(d_ind0_fit >= pcInd0->dComputeFitness(iPatternLevel)) &&
						(d_ind1_fit >= pcInd1->dComputeFitness(iPatternLevel))
						)
					)
				{
					b_add_higher_level_gene(NULL, pcCrossingTree->iGetPatternLevel() + 1);//create level of high level genes
					v_insert_gene_pattern_part_level(pcCrossingTree->iGetPatternLevel() + 1);
					v_add_gene_pattern_tree(NULL, pcCrossingTree->iGetPatternLevel() + 1);
					v_update_high_lvl_genes_and_trees();
				}//if  (
			}//if  (i_pattern == 0)//if we are at the top of the tree
	
		}//if  (pc_best_crossing_node !=  NULL)*/
	}//if  (pcCrossingTree->bMarkedPhenotypeDifferences(pcInd0, pcInd1) == true)


	return(false);
}//bool  C3LO::




bool  C3LOSingle::b_cross_individuals_by_tree(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOPattern *pcCrossingTree, int iPatternLevel)
{
	//first do crossing with current tree node



	if (pcCrossingTree->bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == true)
	{
		C3LOPattern *pc_best_crossing_node;

		pc_best_crossing_node = pcCrossingTree->pcGetBestCrossingNode(iPatternLevel, pcInd0, pcInd1);

		if (pc_best_crossing_node != NULL)
		{
			return(b_cross_infect_individual_by_tree(pcInd0, pcInd1, pc_best_crossing_node, iPatternLevel));
		}//if  (pc_best_crossing_node !=  NULL)
	}//if  (pcCrossingTree->bMarkedPhenotypeDifferences(pcInd0, pcInd1) == true)

	return(false);
}//void  C3LO::v_cross_individuals_by_tree(C3LOIndividual  *pcInd0, C3LOIndividual  *pcInd1, C3LOPattern *pcCrossingTree)




bool  C3LOSingle::b_cross_infect_individual_by_tree(C3LOIndividual  *pcInd, C3LOIndividual  *pcDonator, C3LOPattern *pcCrossingTree, int iPatternLevel)
{
	double  d_ind_fit;

	d_ind_fit = pcInd->dComputeFitness(iPatternLevel);
	pcDonator->dComputeFitness(iPatternLevel);


	for (int i_gene = 0; i_gene < i_templ_length; i_gene++)
		pi_genotype_tool[i_gene] = pcInd->vpi_optimized_genotypes.at(iPatternLevel)[0][i_gene];

	for (int ii = 0; ii < pcCrossingTree->pvGetPattern()->size(); ii++)
		pi_genotype_tool[pcCrossingTree->pvGetPattern()->at(ii).iGenePos()] = pcDonator->vpi_optimized_genotypes.at(iPatternLevel)[0][pcCrossingTree->pvGetPattern()->at(ii).iGenePos()];

	i_all_crossings++;

	double  d_crossed_fit;
	d_crossed_fit = dComputeFitness(pi_genotype_tool);

	if (d_crossed_fit > d_ind_fit)
	{
		i_all_effective_crossings++;

		if (b_can_ind_be_added_at_any_level(pi_genotype_tool, true, true) == false)
		{
			i_crossings_with_effect_in_the_pop++;
			return(false);
		}//if (b_can_ind_be_added_at_any_level(pi_genotype_tool) == false)

		pcInd->vSetGenotype(pi_genotype_tool);
		pcInd->i_level++;
	}//if (d_crossed_fit > d_ind_fit)

	
	return(true);
}//void  C3LO::v_cross_infect_individual_by_tree(C3LOIndividual  *pcInd, C3LOIndividual  *pcDonator, C3LOPattern *pcCrossingTree)






void  C3LOSingle::v_cross_infect_individual_by_tree(int  *piGenotype, int  *piGenotypeDonator, C3LOPattern *pcCrossingTree)
{
	for (int ii = 0; ii < pcCrossingTree->pvGetPattern()->size(); ii++)
		piGenotype[pcCrossingTree->pvGetPattern()->at(ii).iGenePos()] = piGenotypeDonator[pcCrossingTree->pvGetPattern()->at(ii).iGenePos()];
}//bool  C3LO::b_cross_infect_individual_by_tree(int  *piGenotype, C3LOIndividual  *pcDonator, C3LOPattern *pcCrossingTree, int iPatternLevel)






void  C3LOSingle::v_remove_from_levels(C3LOIndividual *pcInd)
{
	for (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)
	{
		for (int i_ind = 0; i_ind < v_population_levels.at(i_lvl).size(); i_ind++)
		{
			if (v_population_levels.at(i_lvl).at(i_ind) == pcInd)
				v_population_levels.at(i_lvl).erase
				(
					v_population_levels.at(i_lvl).begin() + i_ind
				);

		}//for  (int ii = 0; ii < v_population_levels.at(iLevel).size(); ii++)
	}//for  (int i_lvl = 0; i_lvl < v_population_levels.size(); i_lvl++)

}//void  C3LO::v_remove_from_levels(C3LOIndividual *pcInd)






void  C3LOSingle::v_add_intercept_pattern(int  *piPatternDefinition, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies)
{
	for (int ii = 0; ii < pvGenePatternsPartsInterceptionsAndDifferencies->size(); ii++)
	{
		if (pvGenePatternsPartsInterceptionsAndDifferencies->at(ii)->bEqualToDefinition(piPatternDefinition))
		{
			pvGenePatternsPartsInterceptionsAndDifferencies->at(ii)->i_number_of_hits++;
			return;
		}//if (pvGenePatternsPartsInterceptionsAndDifferencies->at(ii)->bEqualToDefinition(piPatternDefinition))
	}//for (int ii = 0; ii < pvGenePatternsPartsInterceptionsAndDifferencies->size(); ii++)


	C3LOPattern  *pc_new_pattern;
	pc_new_pattern = new C3LOPattern(i_templ_length);
	
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piPatternDefinition[ii] == 1)  pc_new_pattern->v_pattern.push_back(CMessyGene(1, ii));
	}//for (int ii = 0; ii < i_templ_length; ii++)


	pvGenePatternsPartsInterceptionsAndDifferencies->push_back(pc_new_pattern);
}//void  C3LO::v_add_intercept_pattern(int  *piPatternDefinition, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies)





void  C3LOSingle::v_build_create_intercept_patterns(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies)
{
	int  *pi_pattern_generation_tool;
	pi_pattern_generation_tool = new int[i_templ_length];


	C3LOPattern *pc_part_0, *pc_part_1;
	int  i_resulting_genes;

	for (int i_pat_part_0 = 0; i_pat_part_0 < pvGenePatternsParts->size(); i_pat_part_0++)
	{
		pc_part_0 = pvGenePatternsParts->at(i_pat_part_0);
		
		for (int i_pat_part_1 = 0; i_pat_part_1 < pvGenePatternsParts->size(); i_pat_part_1++)
		{
			pc_part_1 = pvGenePatternsParts->at(i_pat_part_1);


			i_resulting_genes = pc_part_0->iGetSubtraction(pc_part_1, pi_pattern_generation_tool);
			if (i_resulting_genes > 1)
			{
				if (i_resulting_genes < pc_part_0->v_pattern.size())
					v_add_intercept_pattern(pi_pattern_generation_tool, pvGenePatternsPartsInterceptionsAndDifferencies);
			}//if (pc_part_0->iGetSubtraction(pc_part_1, pi_pattern_generation_tool) > 1)

			i_resulting_genes = pc_part_1->iGetSubtraction(pc_part_0, pi_pattern_generation_tool);
			if (i_resulting_genes > 1)
			{
				if (i_resulting_genes < pc_part_1->v_pattern.size())
					v_add_intercept_pattern(pi_pattern_generation_tool, pvGenePatternsPartsInterceptionsAndDifferencies);
			}//if (pc_part_0->iGetSubtraction(pc_part_1, pi_pattern_generation_tool) > 1)
			
			i_resulting_genes = pc_part_0->iGetAND(pc_part_1, pi_pattern_generation_tool);
			if (i_resulting_genes > 1)
			{
				v_add_intercept_pattern(pi_pattern_generation_tool, pvGenePatternsPartsInterceptionsAndDifferencies);
			}//if (pc_part_0->iGetSubtraction(pc_part_1, pi_pattern_generation_tool) > 1)
		}//for (int i_pat_part_1 = 0; i_pat_part_1 < pvGenePatternsParts->size(); i_pat_part_1++)	
	}//for (int i_pat_part_0 = 0; i_pat_part_0 < pvGenePatternsParts->size(); i_pat_part_0++)


	delete  pi_pattern_generation_tool;
}//void  C3LO::v_build_create_intercept_patterns(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsInterceptionsAndDifferencies)



void  C3LOSingle::v_load_isis_dsm(double  ***pdDSM, CString  sSource)
{
	CString  s_line, s_buf;

	FILE  *pf_isis;
	pf_isis = fopen(sSource, "r");

	s_line = Tools::sReadLine(pf_isis);
	s_line = Tools::sReadLine(pf_isis);

	int  i_finish_index;
	int  i_gene_first, i_gene_second;


	//first create dsm matrix...
	int  i_templ_length = 784;
	(*pdDSM) = new double*[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		(*pdDSM)[ii] = new double[i_templ_length];

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int iy = 0; iy < i_templ_length; iy++)
			(*pdDSM)[ix][iy] = 0;
	}//for (int ix = 0; ix < i_templ_length; ix++)


	while (!feof(pf_isis))
	{
		s_line = Tools::sReadLine(pf_isis);
		i_gene_first = Tools::iExtractFromString(s_line, 0, &i_finish_index);
		i_gene_second = Tools::iExtractFromString(s_line, i_finish_index, &i_finish_index);

		(*pdDSM)[i_gene_first][i_gene_second] = 1;
		(*pdDSM)[i_gene_second][i_gene_first] = 1;

		//s_buf.Format("%d/%d", i_gene_first, i_gene_second);
		//::MessageBox(NULL, s_buf, s_buf, MB_OK);
	}//while  (!feof(pf_isis))

	(*pdDSM)[0][0] = 0;

	fclose(pf_isis);
}//void  v_load_isis_dsm(double  **pdDSM, CString  sSource)



void  C3LOSingle::v_report_dsm(double  **pdDSM, CString  sDest)
{
	FILE  *pf_dsm;
	pf_dsm = fopen(sDest, "w+");

	fprintf(pf_dsm, "\t");
	for (int ix = 0; ix < i_templ_length; ix++)
		fprintf(pf_dsm, "%d\t", ix);
	fprintf(pf_dsm, "\n");

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		fprintf(pf_dsm, "%d\t", ix);

		for (int iy = 0; iy < i_templ_length; iy++)
			fprintf(pf_dsm, "%.0lf\t", pdDSM[ix][iy]);

		fprintf(pf_dsm, "\n");
	}//for (int ix = 0; ix < i_templ_length; ix++)


	fclose(pf_dsm);
}//void  v_report_dsm(double  **pdDSM, CString  sDest)



void  C3LOSingle::v_report_dsm(int  **piDSM, CString  sDest)
{
	FILE  *pf_dsm;
	pf_dsm = fopen(sDest, "w+");

	fprintf(pf_dsm, "\t");
	for (int ix = 0; ix < i_templ_length; ix++)
		fprintf(pf_dsm, "%d\t", ix);
	fprintf(pf_dsm, "\n");

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		fprintf(pf_dsm, "%d\t", ix);

		for (int iy = 0; iy < i_templ_length; iy++)
			fprintf(pf_dsm, "%d\t", piDSM[ix][iy]);

		fprintf(pf_dsm, "\n");
	}//for (int ix = 0; ix < i_templ_length; ix++)


	fclose(pf_dsm);
}//void  C3LO::v_report_dsm(int  **pdDSM, CString  sDest)



int  C3LOSingle::i_dsm_get_miss_number(double  **pdDSM, double  **pdDSMisis)
{


	int  i_miss_number, i_overhead_number, i_hit_number;

	i_miss_number = 0;
	i_overhead_number = 0;
	i_hit_number = 0;

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int iy = 0; iy < i_templ_length; iy++)
		{
			if ((pdDSMisis[ix][iy] > 0) && (pdDSM[ix][iy] == 0))  i_miss_number++;
		}//for (int iy = 0; iy < i_templ_length; iy++)
	}//for (int ix = 0; ix < i_templ_length; ix++)


	i_last_m_from_dsm = i_miss_number;

	return(i_last_m_from_dsm);
}//int  C3LO::i_dsm_get_miss_number(double  **pdDSM, double  **pdDSMisis)





void  C3LOSingle::v_report_dsm_comparison(double  **pdDSM, double  **pdDSMisis, CString  sDest)
{
	FILE  *pf_comp;
	pf_comp = fopen(sDest, "w+");


	fprintf(pf_comp, "\t");
	for (int ix = 0; ix < i_templ_length; ix++)
		fprintf(pf_comp, "%d\t", ix);
	fprintf(pf_comp, "\n");

	int  i_miss_number, i_overhead_number, i_hit_number;

	i_miss_number = 0;
	i_overhead_number = 0;
	i_hit_number = 0;

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		fprintf(pf_comp, "%d\t", ix);

		for (int iy = 0; iy < i_templ_length; iy++)
		{
			if ((pdDSMisis[ix][iy] == 0) && (pdDSM[ix][iy] == 0))  fprintf(pf_comp, "%.0lf\t", pdDSM[ix][iy]);
			
			if ((pdDSMisis[ix][iy] == 0) && (pdDSM[ix][iy] > 0))
			{
				fprintf(pf_comp, "x%.0lf\t", pdDSM[ix][iy]);
				i_overhead_number++;
			}//if ((pdDSMisis[ix][iy] == 0) && (pdDSM[ix][iy] > 0))


			if ((pdDSMisis[ix][iy] > 0) && (pdDSM[ix][iy] == 0))
			{
				i_miss_number++;
				fprintf(pf_comp, "M%.0lf\t", pdDSM[ix][iy]);
			}//if ((pdDSMisis[ix][iy] > 0) && (pdDSM[ix][iy] == 0))
			
			if ((pdDSMisis[ix][iy] > 0) && (pdDSM[ix][iy] > 0))
			{
				fprintf(pf_comp, "%.0lf\t", pdDSM[ix][iy]);
				i_hit_number++;
			}//if ((pdDSMisis[ix][iy] > 0) && (pdDSM[ix][iy] > 0))
		}//for (int iy = 0; iy < i_templ_length; iy++)

		fprintf(pf_comp, "\n");
	}//for (int ix = 0; ix < i_templ_length; ix++)


	fprintf(pf_comp, "\n");
	fprintf(pf_comp, "\n");
	fprintf(pf_comp, "\n");

	fprintf(pf_comp, "HIT: %d\n", i_hit_number);
	fprintf(pf_comp, "OVERHEAD: %d\n", i_overhead_number);
	fprintf(pf_comp, "MISS: %d\n", i_miss_number);

	i_last_m_from_dsm = i_miss_number;


	fclose(pf_comp);

}//void  C3LO::v_report_dsm_comparison(double  **pdDSM, double  **pdDSMisis, CString  sDest)




void  C3LOSingle::v_build_trees_from_dsm(double  **pdDSM, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel)
{
	vector  <C3LOPattern *>  v_trees_to_join;

	//create leafs
	C3LOPattern  *pc_pattern;

	if (iLevel == 0)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
		{
			pc_pattern = new C3LOPattern(i_templ_length);
			pc_pattern->i_pattern_level = iLevel;
			pc_pattern->v_pattern.push_back(CMessyGene(1, ii));
			//pc_pattern->vGeneratePatternTable();

			v_trees_to_join.push_back(pc_pattern);
		}//for (int ii = 0; ii < i_templ_length; ii++)
	}//if  (iLevel == 0)
	else
	{
		for (int ii = 0; ii < v_gene_patterns_trees.at(iLevel - 1).size(); ii++)
		{
			pc_pattern = new C3LOPattern(i_templ_length);
			pc_pattern->i_pattern_level = iLevel;
			pc_pattern->v_pattern = v_gene_patterns_trees.at(iLevel - 1).at(ii)->v_pattern;
			//pc_pattern->vGeneratePatternTable();

			v_trees_to_join.push_back(pc_pattern);
		}//for (int ii = 0; ii < v_gene_patterns_trees.at(iLevel - 1).size(); ii++)
	}//else  if  (iLevel == 0)


	for (int i_pat_0 = 0; i_pat_0 < v_trees_to_join.size() - 1; i_pat_0++)
	{
		for (int i_pat_1 = i_pat_0 + 1; i_pat_1 < v_trees_to_join.size(); i_pat_1++)
		{
			v_trees_to_join.at(i_pat_0)->dGetDSM_Similarity(v_trees_to_join.at(i_pat_1), pdDSM);
		}//for (int i_pat_1 = i_pat_0 + 1; i_pat_1 < v_trees_to_join.size(); i_pat_1++)
	}//for (int i_pat_0 = 0; i_pat_0 < v_trees_to_join.size(); i_pat_0++)

	//::MessageBox(NULL, "", "", MB_OK);
	//v_report_dsm(pdDSM, "zz_dsm.txt");
	

	vector<C3LOPatternPair>  v_possible_pairs;
	double  d_dsm_max, d_dsm_buf;
	int  i_tree_pair_to_join;

	d_dsm_max = 1;
	//while (v_trees_to_join.size() > 1)
	while ( (d_dsm_max >= 1)&&(v_trees_to_join.size() > 1) )
	{
		d_dsm_max = 0;
		for (int ii = 0; ii < v_trees_to_join.size(); ii++)
		{
			d_dsm_buf = v_trees_to_join.at(ii)->dGetMaxDSM();
			if (d_dsm_max < d_dsm_buf)  d_dsm_max = d_dsm_buf;
		}//for (int ii = 0; ii < v_trees_to_join.size(); ii++)

		if (d_dsm_max >= 1)//if similarity is below 1 then it is only based on the length, no real similarity was discovered
		{
			//get possible pairs
			v_possible_pairs.clear();
			for (int ii = 0; ii < v_trees_to_join.size(); ii++)
				v_trees_to_join.at(ii)->vGetAll_DSM_PairsWithSimilarity(d_dsm_max, &v_possible_pairs);

			if (v_possible_pairs.size() == 0)  ::MessageBox(NULL, "Can't make a tree!", "Can't make a tree!", MB_OK);

			//create joined pattern...
			i_tree_pair_to_join = RandUtils::iRandNumber(0, v_possible_pairs.size() - 1);
			pc_pattern = new C3LOPattern(i_templ_length);
			pc_pattern->i_pattern_level = iLevel;
			pc_pattern->vJoinTwoPatterns(v_possible_pairs.at(i_tree_pair_to_join).pc_first, v_possible_pairs.at(i_tree_pair_to_join).pc_second);
			pc_pattern->d_dsm_value = d_dsm_max;
			v_trees_to_join.push_back(pc_pattern);

			//remove joined patterns
			for (int ii = 0; ii < v_trees_to_join.size(); ii++)
			{
				if (ii >= 0)
				{
					if (v_trees_to_join.at(ii) == v_possible_pairs.at(i_tree_pair_to_join).pc_first)
					{
						v_trees_to_join.erase(v_trees_to_join.begin() + ii);
						ii--;
					}//if (v_trees_to_join.at(ii) == v_possible_pairs.at(i_tree_pair_to_join).pc_first)
				}//if (ii >= 0)

				if (ii >= 0)
				{
					if (v_trees_to_join.at(ii) == v_possible_pairs.at(i_tree_pair_to_join).pc_second)
					{
						v_trees_to_join.erase(v_trees_to_join.begin() + ii);
						ii--;
					}//if (v_trees_to_join.at(ii) == v_possible_pairs.at(i_tree_pair_to_join).pc_first)
				}//if  (ii  >= 0)
			}//for (int ii = 0; ii < v_trees_to_join.size(); ii++)


			//remove joined trees from cached dsm similarity
			for (int ii = 0; ii < v_trees_to_join.size(); ii++)
			{
				v_trees_to_join.at(ii)->vRemoveFromDSM_SimilarityList(v_possible_pairs.at(i_tree_pair_to_join).pc_first);
				v_trees_to_join.at(ii)->vRemoveFromDSM_SimilarityList(v_possible_pairs.at(i_tree_pair_to_join).pc_second);
			}//for (int ii = 0; ii < v_trees_to_join.size(); ii++)


			//update cached similarities...
			for (int ii = 0; ii < v_trees_to_join.size(); ii++)
				v_trees_to_join.at(ii)->dGetDSM_Similarity(pc_pattern, pdDSM);

			if (pc_pattern->bAmIThere(&(v_gene_patterns_parts_interceptions_and_differencies.at(iLevel))) == false)
				v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).push_back(pc_pattern);
		}//if (d_dsm_max > 0)
		
	}//while (v_trees_to_join.size() > 1)



	for (int ii = 0; ii < v_trees_to_join.size(); ii++)
	{
		/*pc_pattern = new C3LOPattern(i_templ_length);
		pc_pattern->i_pattern_level = iLevel;
		pc_pattern->d_dsm_value = 0;
		pc_pattern->v_pattern = v_trees_to_join.at(ii)->v_pattern;*/


		pvGenePatternsTrees->push_back(v_trees_to_join.at(ii));
	}//for (int ii = 0; ii < v_trees_to_join.size(); ii++)
}//void  C3LO::v_build_trees_from_dsm(double  **pdDSM, vector<C3LOPattern  *>  *pvGenePatternsTrees)



void  C3LOSingle::v_build_gene_relations_from_dsm(double  **pdDSM, int  **piGeneRelations)
{
	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int iy = 0; iy < i_templ_length; iy++)
			piGeneRelations[ix][iy] = -1;
	}//for (int ix = 0; ix < i_templ_length; ix++)


	int  *pi_already_used_buffer;
	pi_already_used_buffer = new int[i_templ_length];

	double  d_largest_relation;
	int  i_largest_relation_gene_offset;
	int  i_counter;

	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_already_used_buffer[ii] = 0;
		d_largest_relation = 1;
		i_counter = 0;

		while (d_largest_relation > 0)
		{
			d_largest_relation = 0;
			for (int iy = 0; iy < i_templ_length; iy++)
			{
				if (
					(d_largest_relation < pdDSM[ix][iy])&&
					(pi_already_used_buffer[iy] == 0)
					)
				{
					d_largest_relation = pdDSM[ix][iy];
					i_largest_relation_gene_offset = iy;
				}//if (d_largest_relation < pdDSM[ix][iy])
			}//for (int iy = 0; iy < i_templ_length; iy++)

			if (d_largest_relation > 0)
			{
				pi_already_used_buffer[i_largest_relation_gene_offset] = 1;
				piGeneRelations[ix][i_counter] = i_largest_relation_gene_offset;
				i_counter++;
			}//if (d_largest_relation > 0)
		}//while (d_largest_relation > 0)
		
	}//for (int ix = 0; ix < i_templ_length; ix++)


	CString  s_buf;
	s_buf.Format("zz_DSM_%.2d_relations.txt", i_dsm_reporting_enumeration_tool - 1);
	//v_report_dsm(piGeneRelations, s_buf);


}//void  C3LO::v_build_gene_relations_from_dsm(double  **pdDSM, int  **piGeneRalations)



void  C3LOSingle::v_build_dsm(double  **pdDSM, vector<C3LOPattern  *>  *pvGenePatternsParts)
{
	CString  s_buf;


	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int iy = 0; iy < i_templ_length; iy++)
			pdDSM[ix][iy] = 0;
	}//for (int ix = 0; ix < i_templ_length; ix++)


	for (int i_pattern_scrap = 0; i_pattern_scrap < pvGenePatternsParts->size(); i_pattern_scrap++)
	{
		for (int i_gene_first = 0; i_gene_first < pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.size() - 1; i_gene_first++)
		{
			for (int i_gene_second = i_gene_first + 1; i_gene_second < pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.size(); i_gene_second++)
			{
				pdDSM[pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_first).iGenePos()][pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_second).iGenePos()] += 1;
				pdDSM[pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_second).iGenePos()][pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_first).iGenePos()] += 1;
			}//for (int i_gene_second = i_gene_first + 1; i_gene_second < pvGenePatternsParts->size(); i_gene_second++)			
		}//for (int i_first_gene = 0; i_first_gene < pvGenePatternsParts->size(); i_first_gene++)
	}//for (int i_pattern_scrap = 0; i_pattern_scrap < pvGenePatternsParts->size(); i_pattern_scrap++)

	//double  **pd_dsm_isis;
	//v_load_isis_dsm(&pd_dsm_isis, "IsingSpinGlass_pm_784_0.txt");
	//v_report_dsm(pd_dsm_isis, "zz_DSM_isis.txt");


	//print dsm...
	s_buf.Format("zz_DSM_%.2d.txt", i_dsm_reporting_enumeration_tool);
	//v_report_dsm(pd_dsm_matrix, s_buf);
	s_buf.Format("zz_DSM_%.2d_comp.txt", i_dsm_reporting_enumeration_tool);
	//v_report_dsm_comparison(pd_dsm_matrix, pd_dsm_isis, s_buf);
	//i_dsm_get_miss_number(pd_dsm_matrix, pd_dsm_isis);
	i_dsm_reporting_enumeration_tool++;


	//first flush current trees
	//for (int ii = 0; ii < pvGenePatternsTrees->size(); ii++)
	//delete  pvGenePatternsTrees->at(ii);
	//pvGenePatternsTrees->clear();

	//for (int ii = 0; ii < v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).size(); ii++)
	//delete  v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).at(ii);
	//v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).clear();



	s_buf.Format("zzz_pats_%.2d.txt", i_dsm_reporting_enumeration_tool - 1);
	//v_save_trees_and_patterns(s_buf);

	s_buf.Format("LastM: %d  ", i_last_m_from_dsm);
	pc_log->vPrintLine(s_buf, true);


	//for (int ii = 0; ii < i_templ_length; ii++)
		//delete  pd_dsm_isis[ii];
	//delete  pd_dsm_isis;


	return;
}//void  C3LO::v_build_dsm(double  ***pdDSM, vector<C3LOPattern  *>  *pvGenePatternsParts)




void  C3LOSingle::v_build_pattern_tree_dsm_like(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel)
{
	CString  s_buf;

	
	for (int ix = 0; ix < i_templ_length; ix++)
	{
		for (int iy = 0; iy < i_templ_length; iy++)
			pd_dsm_matrix[ix][iy] = 0;
	}//for (int ix = 0; ix < i_templ_length; ix++)


	for (int i_pattern_scrap = 0; i_pattern_scrap < pvGenePatternsParts->size(); i_pattern_scrap++)
	{
		for (int i_gene_first = 0; i_gene_first < pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.size() - 1; i_gene_first++)
		{
			for (int i_gene_second = i_gene_first + 1; i_gene_second < pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.size(); i_gene_second++)
			{
				pd_dsm_matrix[pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_first).iGenePos()][pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_second).iGenePos()] += 1;
				pd_dsm_matrix[pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_second).iGenePos()][pvGenePatternsParts->at(i_pattern_scrap)->v_pattern.at(i_gene_first).iGenePos()] += 1;
			}//for (int i_gene_second = i_gene_first + 1; i_gene_second < pvGenePatternsParts->size(); i_gene_second++)			
		}//for (int i_first_gene = 0; i_first_gene < pvGenePatternsParts->size(); i_first_gene++)
	}//for (int i_pattern_scrap = 0; i_pattern_scrap < pvGenePatternsParts->size(); i_pattern_scrap++)

	//double  **pd_dsm_isis;
	//v_load_isis_dsm(&pd_dsm_isis, "IsingSpinGlass_pm_784_0.txt");
	//v_report_dsm(pd_dsm_isis, "zz_DSM_isis.txt");

	
	//print dsm...
	s_buf.Format("zz_DSM_%d_%.2d.txt", iLevel, i_dsm_reporting_enumeration_tool);
	//v_report_dsm(pd_dsm_matrix, s_buf);
	s_buf.Format("zz_DSM_%.2d_comp.txt", i_dsm_reporting_enumeration_tool);
	//v_report_dsm_comparison(pd_dsm_matrix, pd_dsm_isis, s_buf);
	//i_dsm_get_miss_number(pd_dsm_matrix, pd_dsm_isis);
	i_dsm_reporting_enumeration_tool++;

	
	//ALREADY DONE IN bRunIteration
	//first flush current trees
	//for (int ii = 0; ii < pvGenePatternsTrees->size(); ii++)
		//delete  pvGenePatternsTrees->at(ii);
	//pvGenePatternsTrees->clear();

	/*for (int ii = 0; ii < v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).size(); ii++)
		delete  v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).at(ii);*/
	//v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).clear();
	
	
	v_build_trees_from_dsm(pd_dsm_matrix, pvGenePatternsTrees, iLevel);


	s_buf.Format("zzz_pats_%.2d.txt", iLevel, i_dsm_reporting_enumeration_tool - 1);
	//v_save_trees_and_patterns(s_buf);

	s_buf.Format("LastM: %d  tree patterns number: %d", i_last_m_from_dsm, v_gene_patterns_parts_interceptions_and_differencies.at(iLevel).size());
	pc_log->vPrintLine(s_buf, true);


	//for (int ii = 0; ii < i_templ_length; ii++)
		//delete  pd_dsm_isis[ii];
	//delete  pd_dsm_isis;
	

	return;
}//void  C3LO::v_build_pattern_tree_dsm_like(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel)



void  C3LOSingle::v_build_pattern_tree(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees, int  iLevel)
{
	vector<C3LOPattern  *>  v_parts_buffer;
	v_parts_buffer = *pvGenePatternsParts;

	//first flush current trees
	for (int ii = 0; ii < pvGenePatternsTrees->size(); ii++)
		delete  pvGenePatternsTrees->at(ii);
	pvGenePatternsTrees->clear();

	C3LOPattern  *pc_new_tree;
	bool  b_enhanced;
	int  i_common_pos, i_uncommon_pos, i_unmarked_pos;
	while (v_parts_buffer.size() > 0)
	{
		pc_new_tree = new C3LOPattern(i_templ_length);

		pvGenePatternsTrees->push_back(pc_new_tree);
		pc_new_tree->v_pattern = v_parts_buffer.at(0)->v_pattern;
		v_parts_buffer.erase(v_parts_buffer.begin());

		b_enhanced = true;
		while (b_enhanced == true)
		{
			b_enhanced = false;
			for (int i_current_pattern_part = 0; i_current_pattern_part < v_parts_buffer.size(); i_current_pattern_part++)
			{
				pc_new_tree->vGetCommonAndUNcommonPositions(v_parts_buffer.at(i_current_pattern_part), &i_common_pos, &i_uncommon_pos, &i_unmarked_pos);

				//if two pattern have ANY common position and ANY different position - we glue them
				if (
					(i_common_pos > 0) //&&
					//((i_uncommon_pos > 0) || (i_unmarked_pos > 0))
					)
				{
					pc_new_tree->vJoinSinglePatternAdd(v_parts_buffer.at(i_current_pattern_part));

					v_parts_buffer.erase(v_parts_buffer.begin() + i_current_pattern_part);
					i_current_pattern_part--;
					b_enhanced = true;
				}//if (
			
			
			}//for (int i_current_pattern_part = 0; i_current_pattern_part < v_parts_buffer.size(); i_current_pattern_part++)
		}//while (b_enhanced == true)
	
	
	}//while (v_parts_buffer.size() > 0)
	

	return;//old version below...


	//vector<C3LOPattern  *>  v_parts_buffer;

	v_parts_buffer = *pvGenePatternsParts;

	//first build pattern pairs similarity
	vector<C3LOPatternPair>  v_pattern_pairs;
	int  i_common_positions;
	C3LOPatternPair  c_buffer;

	for (int i_first = 0; i_first < v_parts_buffer.size(); i_first++)
	{
		for (int i_second = i_first + 1; i_second < v_parts_buffer.size(); i_second++)
		{
			i_common_positions =
				v_parts_buffer.at(i_first)->iCommonPositions
				(
					v_parts_buffer.at(i_second), false, false
				);

			if (i_common_positions > 0)
			{
				c_buffer.d_similarity = i_common_positions;
				c_buffer.pc_first = v_parts_buffer.at(i_first);
				c_buffer.pc_second = v_parts_buffer.at(i_second);

				v_pattern_pairs.push_back(c_buffer);
			}//if  (i_common_positions > 0)

		}//for  (int  i_second = i_first + 1; i_second < v_parts_buffer.size(); i_second++)
	}//for  (int  i_first = 0; i_first < v_parts_buffer.size(); i_first++)


	double  d_most_similar;
	int  i_most_similar_offset = -1;
	d_most_similar = d_get_most_similar_pattern_pair_offset(&i_most_similar_offset, &v_pattern_pairs);

	//int  i_common_pos, i_uncommon_pos, i_unmarked_pos;
	C3LOPattern  *pc_joined_pattern;
	while (i_most_similar_offset >= 0)
	{
		v_pattern_pairs.at(i_most_similar_offset).pc_first->vGetCommonAndUNcommonPositions(v_pattern_pairs.at(i_most_similar_offset).pc_second, &i_common_pos, &i_uncommon_pos, &i_unmarked_pos);


		if ((i_unmarked_pos == 0) || (i_uncommon_pos == 0))
		{
			if (i_unmarked_pos == 0)
			{
				pc_joined_pattern = v_pattern_pairs.at(i_most_similar_offset).pc_first;
				pc_joined_pattern->vJoinSinglePatternAdd(v_pattern_pairs.at(i_most_similar_offset).pc_second);


				//v_remove_pattern_from_buffer(&v_parts_buffer, v_pattern_pairs.at(i_most_similar_offset).pc_first);
				//v_remove_pattern_from_pairs_buffer(&v_pattern_pairs, v_pattern_pairs.at(i_most_similar_offset).pc_first);
			}//if  (i_unmarked_pos  ==  0)
			else
			{
				pc_joined_pattern = v_pattern_pairs.at(i_most_similar_offset).pc_second;
				pc_joined_pattern->vJoinSinglePatternAdd(v_pattern_pairs.at(i_most_similar_offset).pc_first);


				//v_remove_pattern_from_buffer(&v_parts_buffer, v_pattern_pairs.at(i_most_similar_offset).pc_second);
				//v_remove_pattern_from_pairs_buffer(&v_pattern_pairs, v_pattern_pairs.at(i_most_similar_offset).pc_second);
			}//else  if  (i_unmarked_pos  ==  0)

			v_remove_pattern_from_buffer(&v_parts_buffer, pc_joined_pattern);
			v_remove_pattern_from_pairs_buffer(&v_pattern_pairs, pc_joined_pattern);
		}//if  ((i_unmarked_pos  ==  0)||(i_uncommon_pos  ==  0))
		else
		{
			pc_joined_pattern = new C3LOPattern(i_templ_length);
			pc_joined_pattern->vJoinTwoPatterns(v_pattern_pairs.at(i_most_similar_offset).pc_first, v_pattern_pairs.at(i_most_similar_offset).pc_second);
		}//else  if  ((i_unmarked_pos  ==  0)||(i_uncommon_pos  ==  0))

		 //PRW FINISHED HERE


		 //now remove joined patterns from buffer (but DO NOT DELETE, the parts buffer IS NOT the owner)


		 /*FILE  *pf_test;
		 pf_test = fopen("tesdt.txt", "w+");

		 //fprintf(pf_test,"similarity: %.2lf\n\n\n", (v_pattern_pairs.at(i_most_similar_offset).d_similarity));
		 //v_pattern_pairs.at(i_most_similar_offset).pc_first->vSavePattern(pf_test, pc_fitness);
		 //v_pattern_pairs.at(i_most_similar_offset).pc_second->vSavePattern(pf_test, pc_fitness);

		 //fprintf(pf_test,"AFTER JOIN: \n\n\n");

		 pc_joined_pattern->vSavePattern(pf_test, pc_fitness);

		 fclose(pf_test);

		 ::Tools::vShow(i_most_similar_offset);//*/


		for (int i_child = 0; i_child < pc_joined_pattern->pvGetChildren()->size(); i_child++)
		{
			v_remove_pattern_from_buffer(&v_parts_buffer, pc_joined_pattern->pvGetChildren()->at(i_child));
			v_remove_pattern_from_pairs_buffer(&v_pattern_pairs, pc_joined_pattern->pvGetChildren()->at(i_child));
		}//for  (int  i_child = 0; i_child < pc_joined_pattern->pvGetChildren()->size(); i_child++)


		 //now add pairs for a new part
		for (int i_other = 0; i_other < v_parts_buffer.size(); i_other++)
		{
			i_common_positions =
				v_parts_buffer.at(i_other)->iCommonPositions
				(
					pc_joined_pattern, false, false
				);

			if (i_common_positions > 0)
			{
				c_buffer.d_similarity = i_common_positions;
				c_buffer.pc_first = pc_joined_pattern;
				c_buffer.pc_second = v_parts_buffer.at(i_other);

				v_pattern_pairs.push_back(c_buffer);
			}//if  (i_common_positions > 0)
		}//for  (int  i_other = 0; i_other < v_parts_buffer.size(); i_other++)

		v_parts_buffer.push_back(pc_joined_pattern);

		d_most_similar = d_get_most_similar_pattern_pair_offset(&i_most_similar_offset, &v_pattern_pairs);
	}//while  (i_most_similar_offset >= 0)




	 /*FILE  *pf_test;
	 pf_test = fopen("tesdt.txt", "w+");


	 for  (int  ii = 0; ii < v_parts_buffer.size(); ii++)
	 {
	 v_parts_buffer.at(ii)->vSavePattern(pf_test, pc_fitness);
	 }//for  (int  ii = 0; ii < v_parts_buffer.size(); ii++)

	 fclose(pf_test);

	 ::Tools::vShow(v_parts_buffer.size());//*/

	*pvGenePatternsTrees = v_parts_buffer;

	v_add_linkage_covering_at_level(iLevel);
	for (int i_tree = 0; i_tree < pvGenePatternsTrees->size(); i_tree++)
	{
		for (int i_gene_offset = 0; i_gene_offset < pvGenePatternsTrees->at(i_tree)->pvGetPattern()->size(); i_gene_offset++)
		{
			v_genes_marked_by_linkage.at(iLevel)[pvGenePatternsTrees->at(i_tree)->pvGetPattern()->at(i_gene_offset).iGenePos()] = 1;
			//pi_genes_marked_by_linkage[pvGenePatternsTrees->at(i_tree)->pvGetPattern()->at(i_gene_offset).iGenePos()] = 1;
		}//for  (int i_gene_offset = 0; i_gene_offset < pvGenePatternsTrees->at(i_tree)->pvGetPattern()->size(); i_gene_offset++)
	}//for  (int i_tree = 0; i_tree < pvGenePatternsTrees->size(); i_tree++)
}//void  C3LO::v_build_pattern_tree(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsTrees)





void  C3LOSingle::v_update_high_lvl_genes_and_trees()
{
	vector<C3LOPattern *> v_gene_candidates;
	C3LOPattern  *pc_pattern;


	for (int i_lvl = 1; i_lvl < v_gene_patterns_trees.size(); i_lvl++)
	{
		if (i_lvl >= v_higher_level_trees_genes.size())  b_add_higher_level_gene(NULL, i_lvl);

		if (i_lvl == 2)
		{
			int ig = 0;
			ig++;
		}//if  (i_lvl  == 2)


		v_gene_candidates.clear();

		for (int i_tree = 0; i_tree < v_gene_patterns_trees.at(i_lvl - 1).size(); i_tree++)
			v_gene_patterns_trees.at(i_lvl - 1).at(i_tree)->vGetBestCovering(&v_gene_candidates, NULL, i_lvl - 1);

		for (int i_candidate = 0; i_candidate < v_gene_candidates.size(); i_candidate++)
		{
			pc_pattern = new C3LOPattern(i_templ_length);
			pc_pattern->i_pattern_level = v_gene_candidates.at(i_candidate)->iGetPatternLevel() + 1;
			pc_pattern->v_pattern = v_gene_candidates.at(i_candidate)->v_pattern;

			b_add_higher_level_gene(pc_pattern, i_lvl);
		}//for  (int i_candidate = 0; i_candidate < v_gene_candidates.size(); i_candidate++)
	}//for  (int  i_lvl = 1; i_lvl < v_gene_patterns_trees.size(); i_lvl++)


	v_count_higher_genes_freq();
	/*
	vector<C3LOPattern *> v_diffr_covering_patterns;
	C3LOPattern *pc_best_covering_pattern;


	for  (int i_tree = 0; i_tree < v_gene_patterns_trees.at(pcCrossingTree->iGetPatternLevel()).size(); i_tree++)
	v_gene_patterns_trees.at(pcCrossingTree->iGetPatternLevel()).at(i_tree)->vGetBestCovering(&v_diffr_covering_patterns, pi_differences, pcCrossingTree->iGetPatternLevel());
	*/

}//void  C3LO::v_update_high_lvl_genes_and_trees()




void  C3LOSingle::vShufflePatternsToCrossSingleLevel(int iLevel)
{
	//vector<vector<C3LOPattern *>>  v_gene_patterns_parts_interceptions_and_differencies;

	vector<C3LOPattern *>  v_patterns_original;
	v_patterns_original = v_gene_patterns_parts_interceptions_and_differencies.at(iLevel);

	vector<C3LOPattern *>  *pv_dest;
	pv_dest = &(v_gene_patterns_parts_interceptions_and_differencies.at(iLevel));
	
	pv_dest->clear();


	int  i_current_pattern_length;
	vector<C3LOPattern *>  v_shuffling_buffer;

	int i_chosen_pattern;
	i_current_pattern_length = -1;
	bool  b_move_patterns;
	while (v_patterns_original.size() > 0)
	{
		i_current_pattern_length = v_patterns_original.at(0)->v_pattern.size();
		b_move_patterns = true;
		v_shuffling_buffer.clear();

		while (b_move_patterns == true)
		{
			if (v_patterns_original.size() > 0)
			{
				if (v_patterns_original.at(0)->v_pattern.size() == i_current_pattern_length)
				{
					v_shuffling_buffer.push_back(v_patterns_original.at(0));
					v_patterns_original.erase(v_patterns_original.begin());
				}//if (v_patterns_original.at(0)->v_pattern.size() == i_current_pattern_length)
				else
					b_move_patterns = false;

			}//if (v_patterns_original.size() > 0)
		}//while (b_move_patterns == true)


		while (v_shuffling_buffer.size() > 0)
		{
			i_chosen_pattern = RandUtils::iRandNumber(0, v_shuffling_buffer.size() - 1);
			pv_dest->push_back(v_shuffling_buffer.at(i_chosen_pattern));
			v_shuffling_buffer.erase(v_shuffling_buffer.begin() + i_chosen_pattern);
		}//while (v_shuffling_buffer.size() > 0)

	}//while (v_patterns_original.size() > 0)


}//void  C3LOSingle::vShufflePatternsToCrossSingleLevel(int iLevel)



void  C3LOSingle::vShufflePatternsToCross()
{
	for (int i_level = 0; i_level < v_gene_patterns_parts_interceptions_and_differencies.size(); i_level++)
	{
		vShufflePatternsToCrossSingleLevel(i_level);
	}//for (int i_level = 0; i_level < v_gene_patterns_parts_interceptions_and_differencies.size(); i_level++)
}//void  C3LOSingle::vShufflePatternsToCross()



void  C3LOSingle::vShufflePatternsAndFreqs()
{
	vector<C3LOPattern *>  v_trees_backup;
	int  i_random_offset;

	//first shuffle the trees
	for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)
	{
		v_trees_backup = v_gene_patterns_trees.at(i_pattern_level);

		v_gene_patterns_trees.at(i_pattern_level).clear();

		while (v_trees_backup.size() > 0)
		{
			i_random_offset = RandUtils::iRandNumber(0, v_trees_backup.size() - 1);

			v_trees_backup.at(i_random_offset)->vShuffleFreqs();

			v_gene_patterns_trees.at(i_pattern_level).push_back(v_trees_backup.at(i_random_offset));
			v_trees_backup.erase(v_trees_backup.begin() + i_random_offset);
		}//while (v_trees_backup.size() > 0)
	}//for (int i_pattern_level = 0; i_pattern_level < v_gene_patterns_trees.size(); i_pattern_level++)

	
}//void  C3LOSingle::vShufflePatternsAndFreqs()



void  C3LOSingle::v_create_trees_from_ind(C3LOIndividual  *pcIndiv, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel)
{
	if (pc_parent->b_use_dsm == true)
	{
		return;
	}//if (pc_parent->b_use_dsm == true)

	if (pcIndiv->bGetTreesGenerated(iPatternLevel) == false)
	{
		b_linkage_generated = true;

		double  d_ffe_before, d_time_before;
		double  d_ffe_after, d_time_after;
		d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		c_time_counter.bGetTimePassed(&d_time_before);

		//vShufflePatternsAndFreqs();
		pcIndiv->vSetOptimized(false);

		i_linkage_generations++;

		v_insert_gene_pattern_part_level(iPatternLevel);
		
		CString  s_buf;
		s_buf.Format("Linkage generation start...");
		pc_log->vPrintLine(s_buf, true);

		//pcIndiv->vGeneratePatternSingle(&(v_gene_patterns_parts.at(iPatternLevel)), &(v_gene_patterns_parts_original.at(iPatternLevel)), iPatternLevel);
		pcIndiv->vGeneratePatternSingle(&(v_gene_patterns_parts.at(iPatternLevel)), NULL, iPatternLevel);

		


		v_add_linkage_covering_at_level(iPatternLevel);
		/*
		for (int i_pattern = 0; i_pattern < v_gene_patterns_parts_original.at(iPatternLevel).size(); i_pattern++)
		{
			for (int i_gene_off = 0; i_gene_off < v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern)->pvGetPattern()->size(); i_gene_off++)
				v_genes_marked_by_linkage.at(iPatternLevel)[v_gene_patterns_parts_original.at(iPatternLevel).at(i_pattern)->pvGetPattern()->at(i_gene_off).iGenePos()] = 1;
			//pi_genes_marked_by_linkage[pvGenePatternsTrees->at(i_tree)->pvGetPattern()->at(i_gene_offset).iGenePos()] = 1;
		}//for  (int i_gene_offset = 0; i_gene_offset < pvGenePatternsTrees->at(i_tree)->pvGetPattern()->size(); i_gene_offset++)
		*/

		for (int i_pattern = 0; i_pattern < pcIndiv->v_gene_patterns_parts.at(iPatternLevel).size(); i_pattern++)
		{
			for (int i_gene_off = 0; i_gene_off < pcIndiv->v_gene_patterns_parts.at(iPatternLevel).at(i_pattern)->pvGetPattern()->size(); i_gene_off++)
				v_genes_marked_by_linkage.at(iPatternLevel)[pcIndiv->v_gene_patterns_parts.at(iPatternLevel).at(i_pattern)->pvGetPattern()->at(i_gene_off).iGenePos()] = 1;
			//pi_genes_marked_by_linkage[pvGenePatternsTrees->at(i_tree)->pvGetPattern()->at(i_gene_offset).iGenePos()] = 1;
		}//for  (int i_gene_offset = 0; i_gene_offset < pvGenePatternsTrees->at(i_tree)->pvGetPattern()->size(); i_gene_offset++)

		
		
		d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		c_time_counter.bGetTimePassed(&d_time_after);

		s_buf.Format("...linkage generation end COST(%.2lf sec  %.0lf ffe)\n", d_time_after - d_time_before, d_ffe_after - d_ffe_before);
		pc_log->vPrintLine(s_buf, true);

		
		
		v_build_dsm(pd_dsm_matrix, &(v_gene_patterns_parts.at(iPatternLevel)));
		/*v_build_gene_relations_from_dsm(pd_dsm_matrix, pi_gene_realtions);

		int i_relation_length_cur;

		i_relation_length_max = 0;
		for (int ix = 0; ix < i_templ_length; ix++)
		{
			i_relation_length_cur = 0;
			for (int iy = 0; iy < i_templ_length; iy++)
			{
				if (pi_gene_realtions[ix][iy] >= 0)  i_relation_length_cur++;
			}//for (int iy = 0; iy < i_templ_length; iy++)

			if (i_relation_length_max < i_relation_length_cur)  i_relation_length_max = i_relation_length_cur;
		}//for (int ix = 0; ix < i_templ_length; ix++)

		s_buf.Format("max relation length: %d", i_relation_length_max);
		pc_log->vPrintLine(s_buf, true);//*/
		

		

		//v_build_pattern_tree_dsm_like(&(pcIndiv->v_gene_patterns_parts.at(iPatternLevel)), &(pcIndiv->v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);
		v_build_pattern_tree_dsm_like(&(v_gene_patterns_parts.at(iPatternLevel)), &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);//*/

		/*if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)
		{
			//if (v_population_levels.size() > iPatternLevel)
			{
				for (int i_hl_gene = 0; i_hl_gene < v_gene_patterns_trees.at(iPatternLevel).size(); i_hl_gene++)
				{
					v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vZeroHighLvlGeneFreq();

					for (int i_3lo_pop = 0; i_3lo_pop < pc_parent->v_3lo_pops.size(); i_3lo_pop++)
					{
						v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
						(
							&(pc_parent->v_3lo_pops.at(i_3lo_pop)->v_population_levels.at(0)), iPatternLevel
						);
					}//for (int ii = 0; pc_parent->v_3lo_pops.size(); ii++)


					v_gene_patterns_trees.at(iPatternLevel).at(i_hl_gene)->vCountHighLvlGeneFreq
					(
						&(pc_parent->pc_3lo_single_tool->v_population_levels.at(0)), iPatternLevel + 1
					);
				}//for  (int  i_hl_gene = 0; i_hl_gene < v_higher_level_trees_genes.at(i_lvl).size(); i_hl_gene++)
			}//if  (v_population_levels.size() > iPatternLevel)
		}//if (v_gene_patterns_trees.at(iPatternLevel).size() > 1)*/


					

		//v_build_pattern_tree(&(v_gene_patterns_parts.at(iPatternLevel)), &(v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);
		//v_build_create_intercept_patterns(&(v_gene_patterns_parts.at(iPatternLevel)), &(v_gene_patterns_parts_interceptions_and_differencies.at(iPatternLevel)));
		

		//v_build_pattern_tree(&(pcIndiv->v_gene_patterns_parts.at(iPatternLevel)), &(pcIndiv->v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);		

		d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		c_time_counter.bGetTimePassed(&d_time_after);


		d_linkage_ffe += d_ffe_after - d_ffe_before;
		d_linkage_time += d_time_after - d_time_before;


		s_buf.Format("...linkage COST: %.2lfs %.0lf ", d_time_after - d_time_before, d_ffe_after - d_ffe_before);
		pc_log->vPrintLine(s_buf, true);


		return;//to remove
	}//if  (pv_population->at(ii)->bGetTreesGenerated() ==  false)

}//void  C3LO::v_create_trees_from_ind_high_lvl(C3LOIndividual  *pcIndiv, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel)


void C3LOSingle::v_create_trees_from_dsm(vector<C3LOPattern *>  *pvTreesDest)
{
	b_linkage_generated = true;

	int pi_ones_or_zeros_occurrences[4];
	int i_total;

	double d_joint, d_sum_of_separated;

	for (int ii = 0; ii < i_templ_length - 1; ii++)
	{
		pd_dsm_matrix[ii][ii] = 0;

		for (int jj = ii + 1; jj < i_templ_length; jj++)
		{
			pd_dsm_matrix[ii][jj] = 0;

			pi_ones_or_zeros_occurrences[0] = pppi_genes_values_occurrences[ii][jj][0] + pppi_genes_values_occurrences[ii][jj][1];
			pi_ones_or_zeros_occurrences[1] = pppi_genes_values_occurrences[ii][jj][2] + pppi_genes_values_occurrences[ii][jj][3];
			pi_ones_or_zeros_occurrences[2] = pppi_genes_values_occurrences[ii][jj][0] + pppi_genes_values_occurrences[ii][jj][2];
			pi_ones_or_zeros_occurrences[3] = pppi_genes_values_occurrences[ii][jj][1] + pppi_genes_values_occurrences[ii][jj][3];

			i_total = pi_ones_or_zeros_occurrences[0] + pi_ones_or_zeros_occurrences[1];

			d_joint = MathUtils::dComputeEntropy((uint32_t*)pppi_genes_values_occurrences[ii][jj], (uint32_t)4, (uint32_t)i_total);

			if (d_joint > 0)
			{
				d_sum_of_separated = MathUtils::dComputeEntropy((uint32_t*)pi_ones_or_zeros_occurrences, (uint32_t)4, (uint32_t)i_total);
				pd_dsm_matrix[ii][jj] = 2 - d_sum_of_separated / d_joint;
			}//if (d_joint > 0)

			pd_dsm_matrix[ii][jj] = 1 - pd_dsm_matrix[ii][jj] + 1;
			pd_dsm_matrix[jj][ii] = pd_dsm_matrix[ii][jj];
		}//for (int jj = ii + 1; jj < i_templ_length; jj++)
	}//for (int ii = 0; ii < i_templ_length - 1; ii++)

	v_build_trees_from_dsm(pd_dsm_matrix, pvTreesDest, 0);
}//void C3LOSingle::v_create_trees_from_dsm(vector<C3LOPattern *>  *pvTreesDest)


void  C3LOSingle::v_create_trees_from_ind_2(C3LOIndividual  *pcIndSource, C3LOIndividual  *pcIndDest, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel)
{
	v_insert_gene_pattern_part_level(iPatternLevel);

	v_gene_patterns_parts.at(iPatternLevel).clear();
	v_gene_patterns_trees.at(iPatternLevel).clear();


	CString  s_buf;
	s_buf.Format("Linkage generation start...");
	pc_log->vPrintLine(s_buf, true);

	//pcIndSource->vGeneratePatternSingle(pcIndDest, &(v_gene_patterns_parts.at(iPatternLevel)), &(v_gene_patterns_parts_original.at(iPatternLevel)), iPatternLevel);
	pcIndSource->vGeneratePatternSingle(pcIndDest, &(v_gene_patterns_parts.at(iPatternLevel)), NULL, iPatternLevel);
	v_build_pattern_tree(&(pcIndSource->v_gene_patterns_parts.at(iPatternLevel)), &(pcIndSource->v_gene_patterns_trees.at(iPatternLevel)), iPatternLevel);

	for (int i_new_tree = 0; i_new_tree < pcIndSource->v_gene_patterns_trees.at(iPatternLevel).size(); i_new_tree++)
		pvTreesDest->push_back(pcIndSource->v_gene_patterns_trees.at(iPatternLevel).at(i_new_tree));

	s_buf.Format("...linkage generation end");
	pc_log->vPrintLine(s_buf, true);

}//void  C3LO::v_create_trees_from_ind_2(C3LOIndividual  *pcIndSource, C3LOIndividual  *pcIndDest, vector<C3LOPattern *>  *pvTreesDest, int iPatternLevel)


void C3LOSingle::v_update_genes_values_occurrences(C3LOIndividual *pcInd)
{
	int *pi_genotype = pcInd->piGetGenotype();

	for (int ii = 0; ii < i_templ_length - 1; ii++)
	{
		for (int jj = ii + 1; jj < i_templ_length; jj++)
		{
			if (pi_genotype[ii] == 0 && pi_genotype[jj] == 0)
			{
				pppi_genes_values_occurrences[ii][jj][0]++;
				pppi_genes_values_occurrences[jj][ii][0]++;
			}//if (pi_genotype[ii] == 0 && pi_genotype[jj] == 0)
			else if (pi_genotype[ii] == 0 && pi_genotype[jj] == 1)
			{
				pppi_genes_values_occurrences[ii][jj][1]++;
				pppi_genes_values_occurrences[jj][ii][1]++;
			}//else if (pi_genotype[ii] == 0 && pi_genotype[jj] == 1)
			else if (pi_genotype[ii] == 1 && pi_genotype[jj] == 0)
			{
				pppi_genes_values_occurrences[ii][jj][2]++;
				pppi_genes_values_occurrences[jj][ii][2]++;
			}//else if (pi_genotype[ii] == 1 && pi_genotype[jj] == 0)
			else if (pi_genotype[ii] == 1 && pi_genotype[jj] == 1)
			{
				pppi_genes_values_occurrences[ii][jj][3]++;
				pppi_genes_values_occurrences[jj][ii][3]++;
			}//else if (pi_genotype[ii] == 1 && pi_genotype[jj] == 1)
		}//for (int jj = ii + 1; jj < i_templ_length; jj++)
	}//for (int ii = 0; ii < i_templ_length - 1; ii++)
}//void C3LOSingle::v_update_genes_values_occurrences(C3LOIndividual *pcInd)


void  C3LOIndividual::vGeneratePatternSingle(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel)
{
	//dComputeFitness(iPatternLevel);//force phenotype creation

	CString  s_buf;
	s_buf.Format("00 vGeneratePatternSingle_2 level=%d", iPatternLevel);
	pc_parent->pc_log->vPrintLine(s_buf, true);


	if (iPatternLevel  <  v_gene_patterns_parts.size())
	{
		if (v_gene_patterns_parts.at(iPatternLevel).size() > 0)  return;
	}//if  (iPatternLevel  <  v_gene_patterns_parts.size())
	else
	{
		while (iPatternLevel >= v_gene_patterns_parts.size())  v_gene_patterns_parts.push_back(vector<C3LOPattern *>());
		while (iPatternLevel >= v_gene_patterns_trees.size())  v_gene_patterns_trees.push_back(vector<C3LOPattern *>());
	}//else  if  (iPatternLevel  <  v_gene_patterns_parts.size())

	s_buf.Format("01 vGeneratePatternSingle_2 level=%d", iPatternLevel);
	pc_parent->pc_log->vPrintLine(s_buf, true);

	vGetPatternsParts(pcIndDest, &v_gene_patterns_parts.at(iPatternLevel), iPatternLevel);

	/*C3LOPattern  *pc_pattern;
	for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)
	{
	pc_pattern = new C3LOPattern(i_templ_length);
	pc_pattern->v_pattern = v_gene_patterns_parts.at(iPatternLevel).at(ii)->v_pattern;
	pvGenePatternsPartsOriginal->push_back(pc_pattern);
	}//for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)*/



	for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)
	{
		if (v_gene_patterns_parts.at(iPatternLevel).at(ii)->bAmIThere(pvGenePatternsParts) == false)
			pvGenePatternsParts->push_back(v_gene_patterns_parts.at(iPatternLevel).at(ii));
	}//for  (int  ii = 0; ii < v_gene_patterns_parts.size(); ii++)



	while (iPatternLevel >= v_pattern_generated.size())  v_pattern_generated.push_back(false);
	v_pattern_generated.at(iPatternLevel) = true;
}//void  C3LOIndividual::vGeneratePatternSingle(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel)






void  C3LOIndividual::vGeneratePatternSingle(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel)
{
	dComputeFitness(iPatternLevel);//force phenotype creation

	CString  s_buf;
	s_buf.Format("00 vGeneratePatternSingle level=%d", iPatternLevel);
	pc_parent->pc_log->vPrintLine(s_buf, true);


	if (iPatternLevel  <  v_gene_patterns_parts.size())
	{
		if (v_gene_patterns_parts.at(iPatternLevel).size() > 0)  return;
	}//if  (iPatternLevel  <  v_gene_patterns_parts.size())
	else
	{
		while (iPatternLevel >= v_gene_patterns_parts.size())  v_gene_patterns_parts.push_back(vector<C3LOPattern *>());
		while (iPatternLevel >= v_gene_patterns_trees.size())  v_gene_patterns_trees.push_back(vector<C3LOPattern *>());
	}//else  if  (iPatternLevel  <  v_gene_patterns_parts.size())

	s_buf.Format("01 vGeneratePatternSingle level=%d", iPatternLevel);
	pc_parent->pc_log->vPrintLine(s_buf, true);

	vGetPatternsParts(&v_gene_patterns_parts.at(iPatternLevel), iPatternLevel);
	
	/*C3LOPattern  *pc_pattern;
	for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)
	{
		pc_pattern = new C3LOPattern(i_templ_length);
		pc_pattern->v_pattern = v_gene_patterns_parts.at(iPatternLevel).at(ii)->v_pattern;
		pvGenePatternsPartsOriginal->push_back(pc_pattern);
	}//for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)*/
	


	for (int ii = 0; ii < v_gene_patterns_parts.at(iPatternLevel).size(); ii++)
	{
		//if (v_gene_patterns_parts.at(iPatternLevel).at(ii)->bAmIThere(pvGenePatternsParts) == false)
			pvGenePatternsParts->push_back(v_gene_patterns_parts.at(iPatternLevel).at(ii));
	}//for  (int  ii = 0; ii < v_gene_patterns_parts.size(); ii++)*/



	while (iPatternLevel >= v_pattern_generated.size())  v_pattern_generated.push_back(false);
	v_pattern_generated.at(iPatternLevel) = true;
}//void  C3LOIndividual::vGeneratePatternSingle(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<C3LOPattern  *>  *pvGenePatternsPartsOriginal, int iPatternLevel)



void  C3LOIndividual::vGetPatternsParts(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)
{
	if (iPatternLevel == 0)
		v_get_patterns_parts_zero_level(pcIndDest, pvGenePatternsParts);
	else
		v_get_patterns_parts_high_level(pcIndDest, pvGenePatternsParts, iPatternLevel);
}//void  C3LOIndividual::vGetPatternsParts(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)




void  C3LOIndividual::vGetPatternsParts(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)
{
	if (iPatternLevel == 0)
		v_get_patterns_parts_zero_level(pvGenePatternsParts);
	else
		v_get_patterns_parts_high_level(pvGenePatternsParts, iPatternLevel);
}//void  C3LOIndividual::vGetMaskedPatternsParts(vector<C3LOPattern  *>  *pvGenePatternsParts, vector<CMessyGene>  *pvMask)






bool  C3LOIndividual::b_phenotypes_the_same(int iPatternLevel, int iOrder, int  iPosition, C3LOIndividual *pcOther)
{
	if (pv_optimization_orders == NULL)  return(false);
	if (iOrder >= pv_optimization_orders->size())  return(false);
	if (iOrder >= pcOther->pv_optimization_orders->size())  return(false);
	if (iPatternLevel >= vpi_optimized_genotypes.size())  return(false);
	if (iPatternLevel >= pcOther->vpi_optimized_genotypes.size())  return(false);

	if (vpi_optimized_genotypes.at(iPatternLevel)[iOrder][iPosition] != pcOther->vpi_optimized_genotypes.at(iPatternLevel)[iOrder][iPosition])  return(false);

	return(true);
}//bool  C3LOIndividual::b_phenotypes_the_same(int iPatternLevel, int iOrder, int  iPosition, C3LOIndividual *pcOther)


/*
bool  C3LOIndividual::b_phenotypes_the_same(int  iPosition, C3LOIndividual *pcOther)
{
	CString  s_buf;

	//first make all phenotypes in both individuals
	int  i_phenotype_levels_max;
	i_phenotype_levels_max = vpi_optimized_genotypes.size();
	if (i_phenotype_levels_max < pcOther->vpi_optimized_genotypes.size())  i_phenotype_levels_max = pcOther->vpi_optimized_genotypes.size();

	//create all phenotypes if they are not made yet
	for (int i_pheno = 0; i_pheno < i_phenotype_levels_max; i_pheno++)
	{
		dComputeFitness(i_pheno);
		pcOther->dComputeFitness(i_pheno);
	}//for (int ii = 0; ii < i_phenotype_levels_max; ii++)


	for (int i_pheno = 0; i_pheno < i_phenotype_levels_max; i_pheno++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			//s_buf.Format("ORDERS: %d order %d   pos: %d", pv_optimization_orders->size(), i_order, iPosition);
			//::Tools::vReportInFile("test.txt", s_buf);

			if (vpi_optimized_genotypes.at(i_pheno)[i_order][iPosition] != pcOther->vpi_optimized_genotypes.at(i_pheno)[i_order][iPosition])  return(false);
		}//for  (int  i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
	}//for (int i_pheno = 0; i_pheno < i_phenotype_levels_max; i_pheno++)

	return(true);
}//bool  C3LOIndividual::b_phenotypes_the_same(int  iPosition)
*/




void  C3LOIndividual::v_get_patterns_parts_zero_level(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts)
{
	//v_optimized.at(0)[0] = 0;
	//dComputeFitness(0);//force phenotype creation

	vector<CMessyGene>  v_new_pattern;
	C3LOPattern  *pc_pattern;
	CString  s_buf;


	double  d_fitnes_buf;
	int i_gene_val;
	int  *pi_genotype_source, *pi_phenotype_source, *pi_phenotype_source_buffer, *pi_phenotype_destination;
	pi_genotype_source = new int[i_templ_length];
	pi_phenotype_source = new int[i_templ_length];
	pi_phenotype_source_buffer = new int[i_templ_length];
	double  d_fitness_source, d_fitness_source_buffer;

		
	pi_phenotype_destination = pcIndDest->vpi_optimized_genotypes_for_linkage.at(0)[0];

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_source[ii] = pi_genotype[ii];
	d_fitness_source = d_compute_fitness_at_level_zero(0, pi_genotype_source, pi_phenotype_source);
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_phenotype_source_buffer[ii] = pi_phenotype_source[ii];



	for (int i_gene_offset = 0; i_gene_offset < pv_optimization_orders->at(0).size(); i_gene_offset++)
	{
		if (pi_genotype_source[pv_optimization_orders->at(0).at(i_gene_offset)] != pcIndDest->pi_genotype[pv_optimization_orders->at(0).at(i_gene_offset)])
		{
			pi_genotype_source[pv_optimization_orders->at(0).at(i_gene_offset)] = pcIndDest->pi_genotype[pv_optimization_orders->at(0).at(i_gene_offset)];
			d_fitness_source_buffer = d_compute_fitness_at_level_zero(0, pi_genotype_source, pi_phenotype_source);
		
			v_new_pattern.clear();
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				//if (b_phenotypes_the_same(i_gene_offset, &c_disturbed) == false)
				if (pi_phenotype_source_buffer[ii] != pi_phenotype_source[ii])
				{
					if (d_fitness_source > d_fitness_source_buffer)
						v_new_pattern.push_back(CMessyGene(pi_genotype_source[ii], ii));
					else
						v_new_pattern.push_back(CMessyGene(pi_genotype_source[ii] * (-1) + 1, ii));
				}//if  (b_phenotypes_the_same(i_gene_pos_for_check, &c_disturbed)  ==  false)
			}//for  (int  i_gene_offset = 0; i_gene_offset < pvMask->size(); i_gene_offset++)


			if (v_new_pattern.size() > 1)
			{
				pc_pattern = new C3LOPattern(i_templ_length);
				pc_pattern->v_pattern = v_new_pattern;
				pc_pattern->i_pattern_level = 0;
				pvGenePatternsParts->push_back(pc_pattern);
				
				/*FILE  *pf_test;
				pf_test = fopen("test.txt", "w+");

				pc_pattern->vSavePattern(pf_test, pc_fitness);
				this->vSave(pf_test);
				c_disturbed.vSave(pf_test);
				fclose(pf_test);

				::vShow(0);//*/
			}//if  (v_new_pattern.size() > 1)
					
		}//if (pi_genotype_source[pv_optimization_orders->at(0).at(i_gene_offset)] != pcIndDest->pi_genotype[pv_optimization_orders->at(0).at(i_gene_offset)])

		if (i_gene_offset % 500 == 0)
		{
			double  d_buf;
			d_buf = i_gene_offset;
			d_buf = d_buf / i_templ_length;
			s_buf.Format("patterns parts (%d/%d) %.4lf", i_gene_offset, i_templ_length, d_buf);
			pc_parent->pc_log->vPrintLine(s_buf, true);
		}//if (i_gene_pos % 50)
	}//for (int i_gene_offset = 0; i_gene_offset < pv_optimization_orders->size(); i_gene_offset++)


	delete pi_genotype_source;
	delete pi_phenotype_source;
	delete pi_phenotype_source_buffer;
}//void  C3LOIndividual::v_get_patterns_parts_zero_level(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)




void  C3LOIndividual::v_get_patterns_parts_zero_level(vector<C3LOPattern  *>  *pvGenePatternsParts)
{
	v_add_phenotype_level(1);
	v_optimized.at(1)[0] = 0;
	dComputeFitness(1);//force phenotype creation

	vector<CMessyGene>  v_new_pattern;
	C3LOPattern  *pc_pattern;
	CString  s_buf;


	double  d_fitnes_buf;
	int i_gene_val;
	int  *pi_genotype_buffer, *pi_phenotype_buffer;
	pi_genotype_buffer = new int[i_templ_length];
	pi_phenotype_buffer = new int[i_templ_length];
	//C3LOIndividual  c_disturbed(i_templ_length, pc_problem, pv_optimization_orders, pc_parent);


	for (int i_gene_pos = 0; i_gene_pos < i_templ_length; i_gene_pos++)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			pi_genotype_buffer[ii] = pi_genotype[ii];

		i_gene_val = pi_genotype[i_gene_pos];

		if (i_gene_val == 0)
			i_gene_val = 1;
		else
			i_gene_val = 0;


		double  d_ffe_before, d_ffe_after;
		double  d_time_before, d_time_after;

		//d_ffe_before = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		//pc_parent->c_time_counter.bGetTimePassed(&d_time_before);

		
		pi_genotype_buffer[i_gene_pos] = i_gene_val;
		d_fitnes_buf = d_compute_fitness_at_level_zero(0, pi_genotype_buffer, pi_phenotype_buffer);


		//d_ffe_after = (double)pc_problem->pcGetEvaluation()->iGetFFE();
		//pc_parent->c_time_counter.bGetTimePassed(&d_time_after);

		//s_buf.Format("SINGLE ZERO LEVEL COST(%.2lf sec  %.0lf ffe)\n", d_time_after - d_time_before, d_ffe_after - d_ffe_before);
		//pc_parent->pc_log->vPrintLine(s_buf, true);
		//Tools::vShow(0);


		v_new_pattern.clear();
		for (int i_gene_offset = 0; i_gene_offset < i_templ_length; i_gene_offset++)
		{
			//if (b_phenotypes_the_same(i_gene_offset, &c_disturbed) == false)
			if  (pi_phenotype_buffer[i_gene_offset] != vpi_optimized_genotypes.at(1)[0][i_gene_offset])
			{
				if (d_fitnes_buf > dComputeFitness(1))
					v_new_pattern.push_back(CMessyGene(pi_genotype_buffer[i_gene_offset], i_gene_offset));
				else
					v_new_pattern.push_back(CMessyGene(pi_genotype[i_gene_offset], i_gene_offset));
			}//if  (b_phenotypes_the_same(i_gene_pos_for_check, &c_disturbed)  ==  false)
		}//for  (int  i_gene_offset = 0; i_gene_offset < pvMask->size(); i_gene_offset++)

		 /*if  (ii  ==  15)
		 {
		 FILE  *pf_test;
		 pf_test = fopen("test.txt", "w+");
		 c_disturbed.vSave(pf_test);
		 fclose(pf_test);
		 }//if  (ii  ==  15)*/

		if (v_new_pattern.size() > 1)
		{
			pc_pattern = new C3LOPattern(i_templ_length);
			pc_pattern->v_pattern = v_new_pattern;
			pc_pattern->i_pattern_level = 0;
			pvGenePatternsParts->push_back(pc_pattern);


			/*FILE  *pf_test;
			pf_test = fopen("test.txt", "w+");

			pc_pattern->vSavePattern(pf_test, pc_fitness);
			this->vSave(pf_test);
			c_disturbed.vSave(pf_test);
			fclose(pf_test);

			::vShow(0);//*/
		}//if  (v_new_pattern.size() > 1)


		if (i_gene_pos % 500 == 0)
		{
			double  d_buf;
			d_buf = i_gene_pos;
			d_buf = d_buf / i_templ_length;
			s_buf.Format("patterns parts (%d/%d) %.4lf", i_gene_pos, i_templ_length, d_buf);
			pc_parent->pc_log->vPrintLine(s_buf, true);
		}//if (i_gene_pos % 50)

	}//for  (int  ii = 0; ii < pvMask->size(); ii++)


	delete  pi_genotype_buffer;
	delete  pi_phenotype_buffer;
}//void  C3LOIndividual::v_get_patterns_parts_zero_level(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)



void  C3LOIndividual::v_get_patterns_parts_high_level(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)
{

}//void  C3LOIndividual::v_get_patterns_parts_high_level(C3LOIndividual  *pcIndDest, vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)



void  C3LOIndividual::v_print_genotype(int  *piGenotypeToPrint, FILE  *pfDest)
{
	CString  s_genotype, s_buf;;

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		s_buf.Format("%d", piGenotypeToPrint[ii]);
		s_genotype += s_buf;
	}//for (int ii = 0; ii < i_templ_length; ii++)

	fprintf(pfDest, s_genotype);
	fprintf(pfDest, "\n");
}//void  C3LOIndividual::v_print_genotype(int  *piGenotypeToPrint, FILE  *pfDest)


void  C3LOIndividual::vGetBestGenes
	(
		int  iBestGenesNumber, vector<C3LOHighLvlGeneFreq *> *pvResult, vector<C3LOHighLvlGeneFreq  *>  *pvGenesPool,
		int  *piOriginalGenotype, int  *piGenotypeBuffer, 
		int  iPatternLevel
	)
{
	vector<double>  v_high_lvl_genes_fit;
	double  d_fit_buf;
	int  i_position;
	bool  b_pos_found;

	for (int i_gene_pool_offset = 0; i_gene_pool_offset < pvGenesPool->size(); i_gene_pool_offset++)
	{
		for (int ii = 0; ii < i_templ_length; ii++)
			piGenotypeBuffer[ii] = piOriginalGenotype[ii];//we consider the effect of dComputeFitness(iPatternLevel + 1);

		pvGenesPool->at(i_gene_pool_offset)->vInfectPhenotype(piGenotypeBuffer);
		d_fit_buf = pc_parent->dComputeFitness(piGenotypeBuffer);

		i_position = 0;
		b_pos_found = false;
		for (int i_result_gene = 0; (i_result_gene < v_high_lvl_genes_fit.size())&&(b_pos_found == false); i_result_gene++)
		{
			i_position = i_result_gene;
			if (v_high_lvl_genes_fit.at(i_result_gene) < d_fit_buf)  b_pos_found = true;
		}//for (int i_result_gene = 0; (i_result_gene < v_high_lvl_genes_fit.size())&&(b_pos_found == false); i_result_gene++)

		if (b_pos_found == false)
		{
			if (pvResult->size() < iBestGenesNumber)
			{
				pvResult->push_back(pvGenesPool->at(i_gene_pool_offset));
				v_high_lvl_genes_fit.push_back(d_fit_buf);
			}//if (pvResult->size() < iBestGenesNumber)
		}//if (b_pos_found == false)
		else
		{
			pvResult->insert(pvResult->begin() + i_position, pvGenesPool->at(i_gene_pool_offset));
			v_high_lvl_genes_fit.insert(v_high_lvl_genes_fit.begin() + i_position, d_fit_buf);
		}//else  if (b_pos_found == false)

		//control if not too long
		while (v_high_lvl_genes_fit.size() > iBestGenesNumber)  v_high_lvl_genes_fit.erase(v_high_lvl_genes_fit.begin() + v_high_lvl_genes_fit.size() - 1);
		while (pvResult->size() > iBestGenesNumber)  pvResult->erase(pvResult->begin() + pvResult->size() - 1);

	}//for (int i_gene_pool_offset = 0; i_gene_pool_offset < pvGenesPool->size(); i_gene_pool_offset++)

}//void  C3LOIndividual::v_get_best_genes




void  C3LOIndividual::v_get_patterns_parts_high_level(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)
{
	if (iPatternLevel - 1 >= pc_parent->v_gene_patterns_trees.size())  return;
	if (pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).size() < 2)  return;

	dComputeFitness(iPatternLevel + 1);//force phenotype creation - for zero level we use PatternLevel = 1(bitwise optimization), for level 1 we use the first (zero) level trees and so on

	C3LOPattern  *pc_pattern;
	int  *pi_pattern;
	bool  b_pat_with_at_least_one_val;

	int i_gene_pos_for_check;
	int i_gene_val;
	C3LOPattern  *pc_current_high_lvl_gene;
	vector<C3LOHighLvlGeneFreq  *>  *pv_gene_freq;
	vector<C3LOHighLvlGeneFreq *>  v_chosen_genes;


	double  d_fitnes_buf;
	int  *pi_genotype_buffer, *pi_phenotype_buffer;
	pi_genotype_buffer = new int[i_templ_length];
	pi_phenotype_buffer = new int[i_templ_length];

	pi_pattern = new int[i_templ_length];

	/*FILE  *pf_test;
	pf_test = fopen("zzz_test.txt", "w+");
	
	//::MessageBox(NULL, "", "", MB_OK);//*/

	for (int i_high_gene_offset = 0; i_high_gene_offset < pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).size(); i_high_gene_offset++)
	{
		pc_current_high_lvl_gene = pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).at(i_high_gene_offset);
		pv_gene_freq = pc_current_high_lvl_gene->pvGetHighLevelGeneFreq();

		v_chosen_genes.clear();
		/*if (pc_current_high_lvl_gene->dCheckHighLevelGeneValues() > 0)
		{
			Tools::vShow(pc_current_high_lvl_gene->dCheckHighLevelGeneValues());
		}//if (pc_current_high_lvl_gene->dCheckHighLevelGeneValues() > 0)*/

		vGetBestGenes
			(
			2, &v_chosen_genes, pv_gene_freq, 
			vpi_optimized_genotypes.at(iPatternLevel + 1)[0], pi_genotype_buffer, 
			iPatternLevel + 1
			);

	
		/*fprintf(pf_test, "\n\n\n [%d|%.1lf] GENE VALUES (%d):\n", i_high_gene_offset, pc_current_high_lvl_gene->dCheckHighLevelGeneValues(), v_chosen_genes.size());
		for (int ii = 0; ii < v_chosen_genes.size(); ii++)
		{
			v_chosen_genes.at(ii)->vSave(pf_test);
			fprintf(pf_test, "\n");
		}//for  (int ii = 0; ii < v_chosen_genes.size(); ii++)*/

	

		//for (int i_gene_freq = 0; (i_gene_freq < pv_gene_freq->size()) && (i_gene_freq < C3LO_INDIVIDUAL_MAX_HIGH_LEVEL_GENES_VALUES); i_gene_freq++)
		for (int i_gene_freq = 0; i_gene_freq < v_chosen_genes.size(); i_gene_freq++)
		{
			//vCopyTo(&c_disturbed);
			for (int ii = 0; ii < i_templ_length; ii++)
				pi_genotype_buffer[ii] = vpi_optimized_genotypes.at(iPatternLevel + 1)[0][ii];//we consider the effect of dComputeFitness(iPatternLevel + 1);

			b_pat_with_at_least_one_val = false;
			for (int ii = 0; ii < i_templ_length; ii++)
				pi_pattern[ii] = 0;

			
			//pv_gene_freq->at(i_gene_freq)->vInfectPhenotype(pi_genotype_buffer);
			v_chosen_genes.at(i_gene_freq)->vInfectPhenotype(pi_genotype_buffer);
			d_fitnes_buf = d_compute_fitness_at_level(iPatternLevel + 1, 0, pi_genotype_buffer, pi_phenotype_buffer);

			/*fprintf(pf_test, "USING:\n");
			v_chosen_genes.at(i_gene_freq)->vSave(pf_test);
			fprintf(pf_test, "\n");
			v_print_genotype(vpi_optimized_genotypes.at(iPatternLevel + 1)[0], pf_test);
			v_print_genotype(pi_genotype_buffer, pf_test);
			v_print_genotype(pi_phenotype_buffer, pf_test);

			fprintf(pf_test, "\n\n");//*/


			pc_current_high_lvl_gene->vGeneratePatternTable();//generates pi_pattern_table
			

			for (int i_gene_offset = 0; i_gene_offset < i_templ_length; i_gene_offset++)
			{
				if (pc_current_high_lvl_gene->pi_pattern_table[i_gene_offset] == 1)
				{
					pi_pattern[i_gene_offset] = 1;
				}//if (pc_current_high_lvl_gene->pi_pattern_table[i_gene_offset] == 1)
				else
				{
					//if (b_phenotypes_the_same(i_gene_offset, &c_disturbed) == false)
					if (vpi_optimized_genotypes.at(iPatternLevel + 1)[0][i_gene_offset] != pi_phenotype_buffer[i_gene_offset])
					{
						pi_pattern[i_gene_offset] = 1;
						b_pat_with_at_least_one_val = true;
					}//if  (b_phenotypes_the_same(i_gene_pos_for_check, &c_disturbed)  ==  false)
				}//else  if (pc_current_high_lvl_gene->pi_pattern_table[i_gene_offset] == 1)
			}//for  (int  i_gene_offset = 0; i_gene_offset < pvMask->size(); i_gene_offset++)

			/*if  (b_pat_with_at_least_one_val == true) 
				fprintf(pf_test, "CHANGE\n");
			else
				fprintf(pf_test, "no change\n");
			v_print_genotype(pi_pattern, pf_test);

			fprintf(pf_test, "\n\n");//*/


			if (b_pat_with_at_least_one_val == true)
			{

				for (int i_high_gene_offset_2 = 0; i_high_gene_offset_2 < pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).size(); i_high_gene_offset_2++)
				{
					if (pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).at(i_high_gene_offset_2)->bDoITouchGeneGroup(pi_pattern) == true)
					{
						for (int ii = 0; ii < pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).at(i_high_gene_offset_2)->v_pattern.size(); ii++)
							pi_pattern[pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).at(i_high_gene_offset_2)->v_pattern.at(ii).iGenePos()] = 1;
					}//if (pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).at(ii)->bDoITouchGeneGroup(pi_pattern) == true)
				}//for (int ii = 0; ii < pc_parent->v_gene_patterns_trees.at(iPatternLevel - 1).size(); ii++)


				pc_pattern = new C3LOPattern(i_templ_length);
				pc_pattern->i_pattern_level = iPatternLevel;


				
				for (int ii = 0; ii < i_templ_length; ii++)
				{
					if (pi_pattern[ii] == 1)  pc_pattern->v_pattern.push_back(CMessyGene(1, ii));
				}//for  (int  ii = 0; ii < i_templ_length; ii++)	


				pc_pattern->vGeneratePatternTable();
				pvGenePatternsParts->push_back(pc_pattern);

				//pc_pattern->vSavePattern(pf_test, NULL);
			}//if  (b_pat_with_at_least_one_val == true)
		}//for  (int  i_gene_freq = 0; i_gene_freq < pv_gene_freq->size(); i_gene_freq++)		

	}//for  (int  i_high_gene_offset = 0; i_high_gene_offset < pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).size(); i_high_gene_offset++)

	 //fclose(pf_test);
	 //::MessageBox(NULL, "", "", MB_OK);


	delete pi_genotype_buffer;
	delete  pi_phenotype_buffer;
	delete  pi_pattern;
}//void  C3LOIndividual::v_get_patterns_parts_zero_level(vector<C3LOPattern  *>  *pvGenePatternsParts, int iPatternLevel)





bool bGeneFreqGreater(C3LOHighLvlGeneFreq *elem1, C3LOHighLvlGeneFreq *elem2)
{
	return elem1->iGetGeneNumber() > elem2->iGetGeneNumber();
}//bool bGeneFreqGreater (C3LOHighLvlGeneFreq *elem1, C3LOHighLvlGeneFreq *elem2)



void  C3LOPattern::vZeroHighLvlGeneFreq()
{
	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
		delete  v_high_level_gene_freq.at(ii);
	v_high_level_gene_freq.clear();
}//void  C3LOPattern::vZeroHighLvlGeneFreq()



void  C3LOPattern::v_brutal_values_generate_high_lvl(vector<C3LOPattern *>  *pvHighLvlGenes, int  iGeneOffset, int  *piGenotypeBuf)
{
	if (iGeneOffset < pvHighLvlGenes->size())
	{
		for (int i_high_lvl_value = 0; i_high_lvl_value < pvHighLvlGenes->at(iGeneOffset)->v_high_level_gene_freq.size(); i_high_lvl_value++)
		{
			pvHighLvlGenes->at(iGeneOffset)->v_high_level_gene_freq.at(i_high_lvl_value)->vInfectPhenotype(piGenotypeBuf);
			v_brutal_values_generate_high_lvl(pvHighLvlGenes, iGeneOffset + 1, piGenotypeBuf);
		}//for (int i_high_lvl_value = 0; i_high_lvl_value < pvHighLvlGenes->at(iGeneOffset)->v_high_level_gene_freq.size(); i_high_lvl_value++)
		
		return;
	}//if (iGeneOffset < pvGenesChosen->size())


	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	{
		if (v_high_level_gene_freq.at(ii)->bCheckHighLevelGene(piGenotypeBuf) == true)  return;
	}//for  (int  ii = 0; ii < v_high_level_gene_freq.size(); ii++)

	 //create high level gene
	C3LOHighLvlGeneFreq  *pc_new_gene_val_freq;

	pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
	pc_new_gene_val_freq->vCreate(this, piGenotypeBuf);

	v_high_level_gene_freq.push_back(pc_new_gene_val_freq);

}//void  C3LOPattern::v_brutal_values_generate_high_lvl(vector<C3LOPattern *>  *&v_pattern_genes_high_lvl, int  iGeneOffset, int  *piGenotypeBuf)



void  C3LOPattern::v_brutal_values_generate_zero_lvl(vector<CMessyGene>  *pvGenesChosen, int  iGeneOffset, int  *piGenotypeBuf)
{
	if (iGeneOffset < pvGenesChosen->size())
	{
		piGenotypeBuf[pvGenesChosen->at(iGeneOffset).iGenePos()] = 0;
		v_brutal_values_generate_zero_lvl(pvGenesChosen, iGeneOffset + 1, piGenotypeBuf);

		piGenotypeBuf[pvGenesChosen->at(iGeneOffset).iGenePos()] = 1;
		v_brutal_values_generate_zero_lvl(pvGenesChosen, iGeneOffset + 1, piGenotypeBuf);

		return;
	}//if (iGeneOffset < pvGenesChosen->size())


	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	{
		if (v_high_level_gene_freq.at(ii)->bCheckHighLevelGene(piGenotypeBuf) == true)  return;
	}//for  (int  ii = 0; ii < v_high_level_gene_freq.size(); ii++)

	 //create high level gene
	C3LOHighLvlGeneFreq  *pc_new_gene_val_freq;

	pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
	pc_new_gene_val_freq->vCreate(this, piGenotypeBuf);

	v_high_level_gene_freq.push_back(pc_new_gene_val_freq);


}//void  C3LOPattern::v_brutal_values_generate_zero_lvl(vector<CMessyGene>  *pvGenesChosen, int  iGeneOffset, int  *piGenotypeBuf)




void  C3LOPattern::vBrutalValuesGenerationHighLvl(C3LOIndividual *pcIndForContext, int iPatternLevel, vector<C3LOPattern *> *pvHighLvlGenesFroLowerLevel)
{
	int  i_old_level;
	i_old_level = pcIndForContext->i_level;
	pcIndForContext->i_level = iPatternLevel + 1;
	pcIndForContext->dComputeFitness(iPatternLevel + 1);
	pcIndForContext->i_level = i_old_level;

	if (pcIndForContext->piGetPhenotype(iPatternLevel + 1, 0) == NULL)  return;

	int  *pi_genotype_buf;
	pi_genotype_buf = new int[i_templ_length];

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_buf[ii] = pcIndForContext->piGetPhenotype(iPatternLevel + 1, 0)[ii];


	vector<C3LOPattern *>  v_pattern_genes_high_lvl;

	vGeneratePatternTable();
	for (int ii = 0; ii < pvHighLvlGenesFroLowerLevel->size(); ii++)
	{
		if (pvHighLvlGenesFroLowerLevel->at(ii)->bAmIInsideGeneGroup(pi_pattern_table) == true)  
			v_pattern_genes_high_lvl.push_back(pvHighLvlGenesFroLowerLevel->at(ii));
	}//for (int ii = 0; ii < pvHighLvlGenesFroLowerLevel->size(); ii++)
	
	//throw out some random genes if necessary
	int  i_rand_offset;
	while (v_pattern_genes_high_lvl.size() > _3LO_HIGH_LVL_GENE_VALUE_GENERATION_SUBGENE_MAX)
	{
		i_rand_offset = RandUtils::iRandNumber(0, v_pattern_genes_high_lvl.size() - 1);
		v_pattern_genes_high_lvl.erase(v_pattern_genes_high_lvl.begin() + i_rand_offset);
	}//if (v_pattern.size() > _3LO_HIGH_LVL_GENE_VALUE_GENERATION_SUBGENE_MAX)

	v_brutal_values_generate_high_lvl(&v_pattern_genes_high_lvl, 0, pi_genotype_buf);


	delete  pi_genotype_buf;
}//void  C3LOPattern::vBrutalValuesGenerationHighLvl(C3LOIndividual *pcIndForContext, int iPatternLevel, vector<C3LOPattern *> *pvHighLvlGenesFroLowerLevel)




void  C3LOPattern::vBrutalValuesGenerationZeroLvl(C3LOIndividual *pcIndForContext, int iPatternLevel)
{
	int  i_old_level;
	i_old_level = pcIndForContext->i_level;
	pcIndForContext->i_level = iPatternLevel + 1;
	pcIndForContext->dComputeFitness(iPatternLevel + 1);
	pcIndForContext->i_level = i_old_level;

	if (pcIndForContext->piGetPhenotype(iPatternLevel + 1, 0) == NULL)  return;

	int  *pi_genotype_buf;
	pi_genotype_buf = new int[i_templ_length];

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype_buf[ii] = pcIndForContext->piGetPhenotype(iPatternLevel + 1, 0)[ii];


	vector<CMessyGene>  v_genes_chosen;
	v_genes_chosen = v_pattern;

	//throw out some random genes if necessary
	int  i_rand_offset;
	while (v_genes_chosen.size() > _3LO_HIGH_LVL_GENE_VALUE_GENERATION_SUBGENE_MAX)
	{
		i_rand_offset = RandUtils::iRandNumber(0, v_genes_chosen.size() - 1);
		v_genes_chosen.erase(v_genes_chosen.begin() + i_rand_offset);
	}//if (v_pattern.size() > _3LO_HIGH_LVL_GENE_VALUE_GENERATION_SUBGENE_MAX)
		
	v_brutal_values_generate_zero_lvl(&v_genes_chosen, 0, pi_genotype_buf);


	delete  pi_genotype_buf;
}//void  C3LOPattern::vBrutalValuesGeneration()


void  C3LOPattern::vCountHighLvlGeneFreq(vector<C3LOIndividual *>  *pvPopoulation, int iPatternLevel)
{
	
	for (int ii = 0; ii < pvPopoulation->size(); ii++)
	{
		v_check_high_level_gene_for_ind(pvPopoulation->at(ii), iPatternLevel);
	}//for  (int  ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	
	//std::sort(v_high_level_gene_freq.begin(), v_high_level_gene_freq.end(), bGeneFreqGreater);
}//void  C3LOPattern::vCountHighLvlGeneFreq(vector<C3LOIndividual *>  *pvPopoulation)



void  C3LOPattern::vShuffleFreqs()
{
	vector<C3LOHighLvlGeneFreq *>  v_high_level_gene_freq_backup;

	v_high_level_gene_freq_backup = v_high_level_gene_freq;
	v_high_level_gene_freq.clear();

	int  i_random_offset;
	while (v_high_level_gene_freq_backup.size() > 0)
	{
		i_random_offset = RandUtils::iRandNumber(0, v_high_level_gene_freq_backup.size() - 1);
		v_high_level_gene_freq.push_back(v_high_level_gene_freq_backup.at(i_random_offset));
		v_high_level_gene_freq_backup.erase(v_high_level_gene_freq_backup.begin() + i_random_offset);
	}//while (v_high_level_gene_freq_backup.size() > 0)

}//void  C3LOPattern::vShuffleFreqs()



double  C3LOPattern::dCheckHighLevelGeneValues()
{
	bool  b_all_zeroes = false;
	bool  b_all_ones = false;

	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	{
		if (v_high_level_gene_freq.at(ii)->bCheckAll(0) == true)  b_all_zeroes = true;
		if (v_high_level_gene_freq.at(ii)->bCheckAll(1) == true)  b_all_ones = true;

		if ((b_all_zeroes == true) && (b_all_ones == true)) return(1);
	}//for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)	

	if ( (b_all_zeroes == true) || (b_all_ones == true)) return(0.5);

	return(0);
}//bool  C3LOPattern::bCheckHighLevelGeneValues()




bool bScrapHitNumber(C3LOPattern *elem1, C3LOPattern *elem2)
{
	return elem1->iGetHitNumber() > elem2->iGetHitNumber();
}//bool bScrapHitNumber(C3LOPattern *elem1, C3LOPattern *elem2)


void  C3LOSingle::v_save_trees_and_patterns(CString  sPatternFile)
{
	FILE  *pf_test;
	
	pf_test = fopen(sPatternFile, "w+");

	fprintf(pf_test, "\n\nINTERCEPTING SCRAPS:\n");

	
	for (int i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_parts_interceptions_and_differencies.size(); i_pat_lvl++)
	{
		fprintf(pf_test, "\n\nSCRAP LEVEL: %d\n", i_pat_lvl);

		std::sort(v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).begin(), v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).end(), bScrapHitNumber);

		for (int ii = 0; ii < v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).size(); ii++)
			v_gene_patterns_parts_interceptions_and_differencies.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", true);
	}//for (int i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_parts_interceptions_and_differencies.size(); i_pat_lvl++)

	

	fprintf(pf_test, "\n\nPATTERNS:\n");
	for (int i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_trees.size(); i_pat_lvl++)
	{
		fprintf(pf_test, "\n\nPATTERN LEVEL: %d\n", i_pat_lvl);

		//fprintf(pf_test, "\n\nPATTERN SCRAPS ORIGINAL: %d\n", i_pat_lvl);

		//for (int ii = 0; ii < v_gene_patterns_parts_original.at(i_pat_lvl).size(); ii++)
			//v_gene_patterns_parts_original.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", false);

		/*fprintf(pf_test, "\n\nPATTERN SCRAPS: %d\n", i_pat_lvl);

		for (int ii = 0; ii < v_gene_patterns_parts.at(i_pat_lvl).size(); ii++)
			v_gene_patterns_parts.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", false);*/



		fprintf(pf_test, "\n\nPATTERN TREES: %d\n", i_pat_lvl);
		for (int ii = 0; ii < v_gene_patterns_trees.at(i_pat_lvl).size(); ii++)
			v_gene_patterns_trees.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", false);



		if (i_pat_lvl < v_higher_level_trees_genes.size())
		{
			fprintf(pf_test, "\n\nPATTERN LEVEL GENES: %d\n", i_pat_lvl);

			for (int ii = 0; ii < v_higher_level_trees_genes.at(i_pat_lvl).size(); ii++)
			{
				v_higher_level_trees_genes.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", true);
			}//for  (int  ii = 0; ii < v_higher_level_trees_genes.at(i_pat_lvl).size(); ii++)
		}//if  (i_pat_lvl < v_higher_level_trees_genes.size())
	}//for  (int  i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_trees.size(); i_pat_lvl++)


	fprintf(pf_test, "\n\nHIGH LEVEL GENES FREQ:\n");
	for (int i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_trees.size(); i_pat_lvl++)
	{
		for (int ii = 0; ii < v_gene_patterns_trees.at(i_pat_lvl).size(); ii++)
		{
			v_gene_patterns_trees.at(i_pat_lvl).at(ii)->vSaveHighLevelGene(pf_test, NULL);
			fprintf(pf_test, "\n");
		}//for  (int  ii = 0; ii < v_higher_level_trees_genes.at(i_pat_lvl).size(); ii++)
	}//for  (int  i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_trees.size(); i_pat_lvl++)



	/*fprintf(pf_test, "\n\nPOPULATION: %d\n", v_population_levels.at(0).size());
	for (int i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)
	{
		v_population_levels.at(0).at(i_ind)->vSave(pf_test, NULL, true);
	}//for  (int  i_ind = 0; i_ind < v_population_levels.at(0).size(); i_ind++)*/

	fclose(pf_test);//*/
}//void  C3LO::v_save_trees_and_patterns(CString  sPatternFile)





double C3LOSingle::dComputeFitness(int32_t *piBits)
{
	double d_fitness_value = d_compute_fitness_value(piBits);
	return(d_fitness_value);
}//double C3LO::dComputeFitness(int32_t *piBits)




 //---------------------------------------------C3LOIndividual-------------------------------------------------------
C3LOIndividual::C3LOIndividual(int iLevel, int  iTemplLength, CProblem<CBinaryCoding, CBinaryCoding > *pcProblem, vector<vector<int>> *pvOptimizationOrders, C3LOSingle  *pcParent)
{
	i_templ_length = iTemplLength;
	pc_parent = pcParent;
	i_level = iLevel;
	i_level_inner = 0;

	pi_genotype = new int[i_templ_length];


	pc_problem = pcProblem;

	v_optimization_order.push_back(vector<int>());
	vGenerateOrder();
	
	pv_optimization_orders = &v_optimization_order;
	//pv_optimization_orders = pvOptimizationOrders;

	v_add_phenotype_level(0);
	

	b_trees_generated = false;

}//C3LOIndividual::C3LOIndividual(int  iTemplLength, gcroot<Random*> pcRandomGen, CConcatDecFunc  *pcFitness, vector<vector<int>> *pvOptimizationOrders)



C3LOIndividual::~C3LOIndividual()
{
	for (int i_pheno = 0; i_pheno < vpi_optimized_genotypes.size(); i_pheno++)
	{
		for (int ii = 0; ii < pv_optimization_orders->size(); ii++)
		{
			delete[] vpi_optimized_genotypes_for_linkage.at(i_pheno)[ii];
			delete[] vpi_optimized_genotypes.at(i_pheno)[ii];
		}//for  (int  ii = 0; ii < pvOptimizationOrders->size(); ii++)

		delete[] vpi_optimized_genotypes_for_linkage.at(i_pheno);
		delete[] vpi_optimized_genotypes.at(i_pheno);

		delete[] vpd_fit_values.at(i_pheno);
			
	}//for (int i_pheno = 0; i_pheno < vpi_optimized_genotypes.size(); i_pheno++)

}//C3LOIndividual::~C3LOIndividual()



void  C3LOIndividual::vGenerateOrder()
{
	pc_parent->v_add_order_random(&(v_optimization_order.at(0)), NULL, false);
}//void  C3LOIndividual::vGenerateOrder()



void  C3LOIndividual::v_add_phenotype_level(int iPhenotypeLevel)
{
	int  **pi_optimized_genotypes;
	double  *pd_fit_values;
	int  *pi_optimized;


	//NOT: while (iPhenotypeLevel + 1> vpi_optimized_genotypes.size() - 1)
	while (iPhenotypeLevel + 1> vpi_optimized_genotypes.size())
	{
		pi_optimized_genotypes = new int*[pv_optimization_orders->size()];
		for (int ii = 0; ii < pv_optimization_orders->size(); ii++)
			pi_optimized_genotypes[ii] = new int[i_templ_length];
		vpi_optimized_genotypes.push_back(pi_optimized_genotypes);

		pi_optimized_genotypes = new int*[pv_optimization_orders->size()];
		for (int ii = 0; ii < pv_optimization_orders->size(); ii++)
			pi_optimized_genotypes[ii] = new int[i_templ_length];
		vpi_optimized_genotypes_for_linkage.push_back(pi_optimized_genotypes);


		
		pd_fit_values = new double[pv_optimization_orders->size()];

		vpd_fit_values.push_back(pd_fit_values);

		pi_optimized = new int[pv_optimization_orders->size()];

		for (int ii = 0; ii < pv_optimization_orders->size(); ii++)
			pi_optimized[ii] = 0;

		v_optimized.push_back(pi_optimized);
	}//while (iPhenotypeLevel + 1> vpi_optimized_genotypes.size())

}//void  C3LOIndividual::v_add_phenotype_level(int iPhenotypeLevel)




void  C3LOIndividual::v_rem_phenotype_level(int iPhenotypeLevel)
{
	int  **pi_optimized_genotypes;
	double  *pd_fit_values;
	int  *pi_optimized;

	while (iPhenotypeLevel < vpi_optimized_genotypes.size() - 1)
	{
		for (int ii = 0; ii < pv_optimization_orders->size(); ii++)
		{
			delete vpi_optimized_genotypes.at(vpi_optimized_genotypes.size() - 1)[ii];
		}//for  (int  ii = 0; ii < pvOptimizationOrders->size(); ii++)

		delete vpi_optimized_genotypes.at(vpi_optimized_genotypes.size() - 1);		
		vpi_optimized_genotypes.erase(vpi_optimized_genotypes.begin() + vpi_optimized_genotypes.size() -1);
		
		delete vpd_fit_values.at(vpd_fit_values.size() - 1);
		vpd_fit_values.erase(vpd_fit_values.begin() + vpd_fit_values.size() - 1);

		delete v_optimized.at(v_optimized.size() - 1);
		v_optimized.erase(v_optimized.begin() + v_optimized.size() - 1);
	}//while (iPhenotypeLevel > vpi_optimized_genotypes.size() - 1)

}//void  C3LOIndividual::v_rem_phenotype_level(int iPhenotypeLevel)




void C3LOIndividual::vRandomInit()
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		pi_genotype[ii] = RandUtils::iRandNumber(0, 1);
	}//for  (int ii = 0; ii < i_templ_length; ii++)
	
	vSetOptimized(false);
};//C3LOIndividual::vRandomInit()




/*double  C3LOIndividual::dComputeFitnessOfCurrentPhenotype(int  iPhenotype)
{
	if (b_optimized == true)  return(d_fit_value);


	d_fit_value = pc_parent->dComputeFitness(pi_optimized_genotypes[iPhenotype]);
	//pc_fitness->eGetFuncValue(pi_optimized_genotypes[iPhenotype], i_templ_length, &d_fit_value);

	pd_fit_values[iPhenotype] = d_fit_value;
	b_optimized = true;

	return(d_fit_value);
}//double  C3LOIndividual::dComputeFitnessOfCurrentPhenotype(int  iPhenotype)*/




double  C3LOIndividual::dGetBestComputedFitness(int  **piBestPhenotype)
{
	dComputeFitness(-1);

	double  d_fit_best, d_fit_buf;

	*piBestPhenotype = NULL;
	d_fit_best = 0;
	for (int i_pattern_level = 0; i_pattern_level < v_optimized.size(); i_pattern_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			if (v_optimized.at(i_pattern_level)[i_order] == 1)
			{
				d_fit_buf = vpd_fit_values.at(i_pattern_level)[i_order];

				if (d_fit_best < d_fit_buf)
				{
					d_fit_best = d_fit_buf;
					*piBestPhenotype = vpi_optimized_genotypes.at(i_pattern_level)[i_order];
				}//if (d_fit_buf > d_fit_best)
			}//if (v_optimized.at(i_pattern_level)[i_order] == 1)
		}//for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		
	}//for (int i_pattern_level = 0; i_pattern_level < v_optimized.size(); i_pattern_level++)

	return(d_fit_best);
}//double  C3LOIndividual::dGetBestComputedFitness(int  **piBestPhenotype)





double  C3LOIndividual::dComputeFitness(int  iPatternLevel)
{
	double  d_fitnes_buf;

	//if (iPatternLevel < 0) iPatternLevel = vpi_optimized_genotypes.size() - 1;
	if (iPatternLevel < 0) iPatternLevel = i_level;
	if (iPatternLevel > i_level) iPatternLevel = i_level;

	v_add_phenotype_level(iPatternLevel);

	for (int i_pattern_level = 0; i_pattern_level <= iPatternLevel; i_pattern_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			if (v_optimized.at(i_pattern_level)[i_order] == 0)
			{
				if (i_pattern_level == 0)
				{
					for (int ii = 0; ii < i_templ_length; ii++)
						vpi_optimized_genotypes.at(0)[i_order][ii] = pi_genotype[ii];

					d_fitnes_buf = pc_parent->dComputeFitness(vpi_optimized_genotypes.at(0)[i_order]);
				}//if (i_pattern_level == 0)

				if (i_pattern_level == 1)
					d_fitnes_buf = d_compute_fitness_at_level_zero(i_order, pi_genotype, vpi_optimized_genotypes.at(1)[i_order]);
				
				if (i_pattern_level > 1)
					d_fitnes_buf = d_compute_fitness_at_level(i_pattern_level, i_order, vpi_optimized_genotypes.at(i_pattern_level - 1)[i_order], vpi_optimized_genotypes.at(i_pattern_level)[i_order]);

				vpd_fit_values.at(i_pattern_level)[i_order] = d_fitnes_buf;
				v_optimized.at(i_pattern_level)[i_order] = 1;
			}//if (v_optimized.at(i_pattern_level)[i_order] == 0)
		}//for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)	
	}//for (int i_pattern_level = 0; i_pattern_level <= iPatternLevel; i_pattern_level++)

	//copying for linkage generation
	for (int i_pattern_level = 0; i_pattern_level <= iPatternLevel; i_pattern_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			for (int i_gene = 0; i_gene < i_templ_length; i_gene++)
				vpi_optimized_genotypes_for_linkage.at(i_pattern_level)[i_order][i_gene] = vpi_optimized_genotypes.at(i_pattern_level)[i_order][i_gene];
		}//for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)	
	}//for (int i_pattern_level = 0; i_pattern_level <= iPatternLevel; i_pattern_level++)

	


	double  d_fit_result;

	d_fit_result = vpd_fit_values.at(iPatternLevel)[0];

	for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
	{
		if (d_fit_result < vpd_fit_values.at(iPatternLevel)[i_order])  d_fit_result = vpd_fit_values.at(iPatternLevel)[i_order];		
	}//for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)	


	return(d_fit_result);
}//double  C3LOIndividual::dComputeFitness()




double  C3LOIndividual::d_compute_fitness_at_level_zero(int  iOrder, int  *piGenotypeToOptimize, int *piPhenotype)
{
	//if (v_optimized.at(0)[iOrder] == 1)  return(vpd_fit_values.at(0)[iOrder]);




	double  d_best_fitness, d_fit_buf;
	int  i_best_fitnes_gene_value;
	int  i_genes_number_optimized;

	//int  *pi_optimized_genotype;
	//pi_optimized_genotype = piPhenotype;//vpi_optimized_genotypes.at(0)[iOrder];

	vector  <int>  v_possible_gene_values;

	for (int ii = 0; ii < i_templ_length; ii++)
		piPhenotype[ii] = piGenotypeToOptimize[ii];

	d_best_fitness = pc_parent->dComputeFitness(piPhenotype);


	i_genes_number_optimized = 0;
	int  i_gene_pos;
	bool  b_at_least_one_gene_optimized;

	b_at_least_one_gene_optimized = true;
	while (b_at_least_one_gene_optimized == true)
	{
		b_at_least_one_gene_optimized = false;

		for (int i_ordered_gene_pos = 0; i_ordered_gene_pos < pv_optimization_orders->at(iOrder).size(); i_ordered_gene_pos++)
		{
			i_gene_pos = pv_optimization_orders->at(iOrder).at(i_ordered_gene_pos);

			//i_original_gene_value = pi_genotype[i_gene_pos];
			i_best_fitnes_gene_value = piPhenotype[i_gene_pos];
			//first get other possible gene values from template and genotype
			v_possible_gene_values.clear();

			//this is prepared for more than binary coding...
			if (piPhenotype[i_gene_pos] == 0)
				v_possible_gene_values.push_back(1);
			else
				v_possible_gene_values.push_back(0);



			for (int ii = 0; ii < v_possible_gene_values.size(); ii++)
			{
				piPhenotype[i_gene_pos] = v_possible_gene_values.at(ii);
				d_fit_buf = pc_parent->dComputeFitness(piPhenotype);
				//pc_fitness->eGetFuncValue(pi_optimized_genotypes[iOrder], i_templ_length, &d_fit_buf);

				if (d_fit_buf > d_best_fitness)
				{
					b_at_least_one_gene_optimized = true;

					i_genes_number_optimized++;

					d_best_fitness = d_fit_buf;
					i_best_fitnes_gene_value = v_possible_gene_values.at(ii);
				}//if  (d_fit_buf  >  d_best_fitness)
			}//for  (int  ii = 0; ii < v_possible_gene_values.size(); ii++)

			piPhenotype[i_gene_pos] = i_best_fitnes_gene_value;
		}//for  (int  i_genotype_iterator = v_genotype.size() - 1; i_genotype_iterator >= 0; i_genotype_iterator--)
	}//while  (b_at_least_one_gene_optimized == true)


	//v_optimized.at(0)[iOrder] = 1;
	//vpd_fit_values.at(0)[iOrder] = d_best_fitness;

	return(d_best_fitness);
}//double  C3LOIndividual::d_compute_fitness_at_level_zero(bool  bOptimize)





double  C3LOIndividual::d_compute_fitness_at_level(int  iPatternLevel, int  iOrder, int  *piGenotypeToOptimize, int *piPhenotype)
{
	double  d_result;
	C3LOPattern  *pc_high_level_gene;

	//if (v_optimized.at(iPatternLevel)[iOrder] == 1)  return(vpd_fit_values.at(iPatternLevel)[iOrder]);


	//int  *pi_optimized_genotype;
	//pi_optimized_genotype = vpi_optimized_genotypes.at(iPatternLevel)[iOrder];

	for (int ii = 0; ii < i_templ_length; ii++)
		piPhenotype[ii] = piGenotypeToOptimize[ii];//vpi_optimized_genotypes.at(iPatternLevel - 1)[iOrder][ii];

	double  d_best_fitness, d_fit_buf;
	d_best_fitness = pc_parent->dComputeFitness(piPhenotype);

	int  *pi_optimized_genotype_buffer;
	pi_optimized_genotype_buffer = new int[i_templ_length];

	bool  b_at_least_one_gene_optimized;

	b_at_least_one_gene_optimized = true;
	
	if (iPatternLevel - 2  <  pc_parent->v_gene_patterns_trees.size())//0- no optimization, 1-genewise optimization, 2 and higher - tree-genes optimization
	{
		/*FILE  *pf_test;
		int i_got_in = 0;
		pf_test = fopen("zz_d_compute_fitness_at_level.txt", "w+");

		fprintf(pf_test, "START: \n");
		vSave(pf_test, true);*/

		while (b_at_least_one_gene_optimized == true)
		{
			b_at_least_one_gene_optimized = false;


			for (int i_high_gene_offset = 0; i_high_gene_offset < pc_parent->v_gene_patterns_trees.at(iPatternLevel - 2).size(); i_high_gene_offset++)
			{
				pc_high_level_gene = pc_parent->v_gene_patterns_trees.at(iPatternLevel - 2).at(i_high_gene_offset);

				//for (int i_high_gene_value = 0; (i_high_gene_value < pc_high_level_gene->pvGetHighLevelGeneFreq()->size()) && (i_high_gene_value < C3LO_INDIVIDUAL_MAX_HIGH_LEVEL_GENES_VALUES); i_high_gene_value++)
				for (int i_high_gene_value = 0; i_high_gene_value < pc_high_level_gene->pvGetHighLevelGeneFreq()->size(); i_high_gene_value++)
				{
					/*i_got_in++;
					fprintf(pf_test, "\nGENE %d VALUE:%d: \n", i_high_gene_offset, i_high_gene_value);
					pc_high_level_gene->vSavePattern(pf_test, pc_fitness, " ", true)
					pc_high_level_gene->pvGetHighLevelGeneFreq()->at(i_high_gene_value)->vSave(pf_test);*/

					//make a copy for a takeback
					for (int ii = 0; ii < i_templ_length; ii++)
						pi_optimized_genotype_buffer[ii] = piPhenotype[ii];
					pc_high_level_gene->pvGetHighLevelGeneFreq()->at(i_high_gene_value)->vInfectPhenotype(piPhenotype);

					d_fit_buf = pc_parent->dComputeFitness(piPhenotype);
					if (d_best_fitness < d_fit_buf)
					{
						d_best_fitness = d_fit_buf;
						b_at_least_one_gene_optimized = true;
					}//if (d_best_fitness < d_fit_buf)
					else
					{
						for (int ii = 0; ii < i_templ_length; ii++)
							piPhenotype[ii] = pi_optimized_genotype_buffer[ii];
					}//else  if (d_best_fitness < d_fit_buf)

				}//for  (i_high_gene_value = 0; i_high_gene_value < pc_high_level_gene->pvGetHighLevelGeneFreq()->size(); i_high_gene_value++)

			}//for  (int  i_high_gene_offset = 0; i_high_gene_offset < pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).size();  i_high_gene_offset++)	
		}//while (b_at_least_one_gene_optimized == true)

		
		/*fprintf(pf_test, "\n\nFINAL:\n");
		vSave(pf_test, true);

		fclose(pf_test);
		if  (i_got_in > 4)  ::vShow("high lvl opt");*/
	}//if  (iPatternLevel  <  pc_parent->pvGetHigherLevelTreesGenes()->size())
	
	delete  pi_optimized_genotype_buffer;
		
	return(d_best_fitness);
}//double  C3LOIndividual::d_compute_fitness_at_level(int  iPatternLevel)








bool  C3LOIndividual::bGetTreesGenerated(int  iPatternLevel)
{
	if (iPatternLevel >= v_pattern_generated.size())  return(false);

	return(v_pattern_generated.at(iPatternLevel));
}//bool  C3LOIndividual::bGetTreesGenerated(int  iPatternLevel)




/*void  C3LOIndividual::vCopyPhenotypeTo(C3LOIndividual  *pcOther)
{
	for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
	{
		pcOther->pd_fit_values[i_order] = pd_fit_values[i_order];

		for (int ii = 0; ii < i_templ_length; ii++)
		{
			pcOther->pi_optimized_genotypes[i_order][ii] = pi_optimized_genotypes[i_order][ii];
		}//for  (int  ii = 0; ii < i_templ_length; ii++)
	}//for  (int  i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
}//void  C3LOIndividual::vCopyPhenotypeTo(C3LOIndividual  *pcOther)*/



void  C3LOIndividual::vCopyTo(C3LOIndividual  *pcOther)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		pcOther->pi_genotype[ii] = pi_genotype[ii];
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	//vCopyPhenotypeTo(pcOther);
	if (vpi_optimized_genotypes.size() < pcOther->vpi_optimized_genotypes.size())
	{
		pcOther->v_rem_phenotype_level(vpi_optimized_genotypes.size());
	}//while (vpi_optimized_genotypes.size() < pcOther->vpi_optimized_genotypes.size())

	if (vpi_optimized_genotypes.size() > pcOther->vpi_optimized_genotypes.size())
	{
		pcOther->v_add_phenotype_level(vpi_optimized_genotypes.size() - 1);
	}//while (vpi_optimized_genotypes.size() < pcOther->vpi_optimized_genotypes.size())



	for (int i_pheno = 0; i_pheno < vpi_optimized_genotypes.size(); i_pheno++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			for  (int  ii = 0; ii < i_templ_length; ii++)
				pcOther->vpi_optimized_genotypes.at(i_pheno)[i_order][ii] = vpi_optimized_genotypes.at(i_pheno)[i_order][ii];

			pcOther->vpd_fit_values.at(i_pheno)[i_order] = vpd_fit_values.at(i_pheno)[i_order];
			pcOther->v_optimized.at(i_pheno)[i_order] = v_optimized.at(i_pheno)[i_order];
		}//for  (int  ii = 0; ii < pvOptimizationOrders->size(); ii++)
	}//for (int i_pheno = 0; i_pheno < vpi_optimized_genotypes.size(); i_pheno++)

	pcOther->pv_optimization_orders = pv_optimization_orders;
}//void  C3LOIndividual::vCopyTo(C3LOIndividual  *pcOther)


void  C3LOIndividual::vSetGenotype(int  *piNewGenotype)
{
	for  (int ii = 0; ii < i_templ_length; ii++)
		pi_genotype[ii] = piNewGenotype[ii];

	vSetOptimized(false);
}//void  C3LOIndividual::vSetGenotype(int  *piNewGenotype)


void  C3LOIndividual::vSetOptimized(bool  bOptimized)
{
	int  i_value = 0;
	if (bOptimized == true)  i_value = 1;

	for (int i_level = 0; i_level < v_optimized.size(); i_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
			v_optimized.at(i_level)[i_order] = i_value;
	}//for (int i_level = 0; i_level < v_optimized.size(); i_level++)
}//void  C3LOIndividual::vSetOptimized(bool  bOptimized)



bool  C3LOIndividual::bComparePhenotype(int  iPatternLevel, int iOrder, int  *piGeneString)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGeneString[ii] != this->vpi_optimized_genotypes.at(iPatternLevel)[iOrder][ii])  return(false);
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOIndividual::bComparePhenotype(C3LOIndividual  *pcOther)



bool  C3LOIndividual::bComparePhenotype(int  iPatternLevel, int iOrder, C3LOIndividual  *pcOther)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pcOther->vpi_optimized_genotypes.at(iPatternLevel)[iOrder][ii] != vpi_optimized_genotypes.at(iPatternLevel)[iOrder][ii])  return(false);
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOIndividual::bComparePhenotype(C3LOIndividual  *pcOther)



bool  C3LOIndividual::bCompareGenotype(C3LOIndividual  *pcOther)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pcOther->pi_genotype[ii] != pi_genotype[ii])  return(false);
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOIndividual::bCompareGenotype(C3LOIndividual  *pcOther)



bool  C3LOIndividual::bCompareGenotype(int  *piGeneString)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (piGeneString[ii] != pi_genotype[ii])  return(false);
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOIndividual::bCompareGenotype(C3LOIndividual  *pcOther)


int  C3LOIndividual::iGetOptDistance(int  iPatternLevel, int iOrder)
{
	if (iPatternLevel < 0)  return(-1);
	if (iPatternLevel >= vpi_optimized_genotypes.size())  return(-1);
	if (iOrder >= pv_optimization_orders->size()) return(-1);

	dComputeFitness(iPatternLevel);

	int  *pi_considered_phenotype;
	pi_considered_phenotype = vpi_optimized_genotypes.at(iPatternLevel)[iOrder];

	//1011011011101001010011011111111010111010000001101110000011010010110011011001011001001010010011010101000100000100111110010011111001101001111001010111100100111000100011001011101111110101100100111110010001101001100101000001110100010111011111000101111111100000110000010101011110010001000100101111001010111110101000100110101001001101001100101001001010011011101010101111111001000001010100011010111001000010111110101001101010101100110111011010101010101010011100100110000101111111010110100011100011010110100101110100001011110000000000010101110101010011010110011010101101010100000101100100101011010000100100011111111100010001000111000010100011010111010100010100110001010110010100101011001111101110110000100110001101000110001101011001001100000100010010110011011100001111101000010010100110101010
	CString  s_opt_as_string;
	s_opt_as_string = "1011011011101001010011011111111010111010000001101110000011010010110011011001011001001010010011010101000100000100111110010011111001101001111001010111100100111000100011001011101111110101100100111110010001101001100101000001110100010111011111000101111111100000110000010101011110010001000100101111001010111110101000100110101001001101001100101001001010011011101010101111111001000001010100011010111001000010111110101001101010101100110111011010101010101010011100100110000101111111010110100011100011010110100101110100001011110000000000010101110101010011010110011010101101010100000101100100101011010000100100011111111100010001000111000010100011010111010100010100110001010110010100101011001111101110110000100110001101000110001101011001001100000100010010110011011100001111101000010010100110101010";
	int  i_differences = 0;
	int  i_bit_from_opt;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (s_opt_as_string.GetAt(ii) == '0')
		{
			i_bit_from_opt = 0;
			//::MessageBox(NULL, "0", "0", MB_OK);
		}//if (s_opt_as_string.GetAt(ii) == '0')
		else
		{
			i_bit_from_opt = 1;
			//::MessageBox(NULL, "1", "1", MB_OK);
		}//else if (s_opt_as_string.GetAt(ii) == '0')

		

		if (pi_considered_phenotype[ii] != i_bit_from_opt)  i_differences++;
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(i_differences);
}//double  C3LOIndividual::dGetOptDistance(int  iPatternLevel, int iOrder)



int*  C3LOIndividual::piGetPhenotype(int  iPatternLevel, int iOrder)
{
	if (iPatternLevel < 0)  return(NULL);
	if (iPatternLevel >= vpi_optimized_genotypes.size())  return(NULL);
	if (iOrder >= pv_optimization_orders->size()) return(NULL);

	dComputeFitness(iPatternLevel);
	//d_optimize_genotype(iPhenotype);

	return(vpi_optimized_genotypes.at(iPatternLevel)[iOrder]);
}//int*  C3LOIndividual::piGetPhenotype(int  iPhenotype)



int*  C3LOIndividual::piGetBestPhenotype()
{
	dComputeFitness(vpi_optimized_genotypes.size() - 1);

	int  *pi_result;
	double  d_best_fit = 0;
	pi_result = NULL;

	for (int i_pattern_level = 0; i_pattern_level < vpi_optimized_genotypes.size(); i_pattern_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			if (d_best_fit < vpd_fit_values.at(i_pattern_level)[i_order])
			{
				pi_result = vpi_optimized_genotypes.at(i_pattern_level)[i_order];
				d_best_fit = vpd_fit_values.at(i_pattern_level)[i_order];
			}//if (vpd_fit_values.at(i_pattern_level)[i_order] > d_best_fit)
		}//for (int i_order = 0; i_order <  pv_optimization_orders->size(); i_order++)
	}//for (int i_pattern_level = 0; i_pattern_level < vpi_optimized_genotypes.size(); i_pattern_level++)

	return(pi_result);
}//int*  C3LOIndividual::piGetBestPhenotype()



double  C3LOIndividual::dUnitation(int  iPatternLevel, int iOrder)
{
	if (iPatternLevel < -1)  return(-1);
	if (iPatternLevel >= vpi_optimized_genotypes.size())  return(-1);
	if (iOrder >= pv_optimization_orders->size()) return(-1);

	dComputeFitness(iPatternLevel);

	double  d_unitation;
	d_unitation = 0;

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (iPatternLevel >= 0)
		{
			if (vpi_optimized_genotypes.at(iPatternLevel)[iOrder][ii] == 1)  d_unitation++;
		}//if  (bOptimizedGenotype  ==  true)
		else
		{
			if (pi_genotype[ii] == 1)  d_unitation++;
		}//else  if  (bOptimizedGenotype  ==  true)

	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	
	d_unitation = d_unitation / i_templ_length;

	return(d_unitation);
}//double  C3LOIndividual::dUnitation()




void  C3LOIndividual::vSaveGenePatternParts(CString  sDebugName)
{
	FILE  *pf_test;


	pf_test = fopen(sDebugName, "w+");
	fprintf(pf_test, "\n\nPATTERNS:\n");
	for (int i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_parts.size(); i_pat_lvl++)
	{
		fprintf(pf_test, "\n\nPATTERN LEVEL: %d\n", i_pat_lvl);

		for (int ii = 0; ii < v_gene_patterns_parts.at(i_pat_lvl).size(); ii++)
			v_gene_patterns_parts.at(i_pat_lvl).at(ii)->vSavePattern(pf_test, NULL, "", true);
	}//for  (int  i_pat_lvl = 0; i_pat_lvl < v_gene_patterns_trees.size(); i_pat_lvl++)


	fclose(pf_test);//*/
}//void  C3LOIndividual::vSaveGenePatternParts(CString  sDebugName)




void  C3LOIndividual::vSave(FILE  *pfReport, int  *piProblemBlockStructure, bool bDoNotForcePhenotype)
{
	double  d_fitness;
	d_fitness = dComputeFitness(-1);

	fprintf(pfReport, "%.8lf \t lvl: %d %.5lf \t %.5lf\n", d_fitness, i_level, dUnitation(0), dUnitation());

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		fprintf(pfReport, "%d", pi_genotype[ii]);
		//if (pc_fitness->bLastBitOfSubfunction(ii) == true)  fprintf(pfReport, " ");
		if (piProblemBlockStructure != NULL)
		{
			if (piProblemBlockStructure[ii] == true)  fprintf(pfReport, " ");
		}//if (piProblemBlockStructure != NULL)
	}//for  (int ii = 0; ii < i_templ_length; ii++)
	fprintf(pfReport, "\n");


	for (int i_pat_level = 0; i_pat_level < vpi_optimized_genotypes.size(); i_pat_level++)
	{
		for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
		{
			for (int ii = 0; ii < i_templ_length; ii++)
			{
				fprintf(pfReport, "%d", vpi_optimized_genotypes.at(i_pat_level)[i_order][ii]);
				//if (pc_fitness->bLastBitOfSubfunction(ii) == true)  fprintf(pfReport, " ");
				if (piProblemBlockStructure != NULL)
				{
					if (piProblemBlockStructure[ii] == true)  fprintf(pfReport, " ");
				}//if (piProblemBlockStructure != NULL)
			}//for  (int ii = 0; ii < i_templ_length; ii++)
			fprintf(pfReport, "\n");
		}//for (int i_order = 0; i_order < pv_optimization_orders->size(); i_order++)
	}//for  (int  i_pat_level = 0; i_pat_level < vpi_optimized_genotypes.size(); i_pat_level++)


	 //for  (int i_pat = 0; i_pat < v_gene_patterns_parts.size(); i_pat++)
	 //	v_gene_patterns_parts.at(i_pat)->vSavePattern(pfReport, pc_fitness);


}//void  C3LOIndividual::vSave(FILE  *pfReport)




void  C3LOIndividual::vIncrementalBrutalUpdate(vector<C3LOPattern *> *pvGenePatterns, vector<C3LOPattern *> *pvGenePatternPartsOriginal, int iPatternLevel)
{
	vector<C3LOPattern *> v_patterns_buf, v_patterns_shuffled;

	v_patterns_buf = *pvGenePatterns;
	int  i_pattern_offset;
	while (v_patterns_buf.size() > 0)
	{
		
		//i_pattern_offset = pc_random_gen->Next(0, v_patterns_buf.size());
		i_pattern_offset = RandUtils::iRandNumber(0, v_patterns_buf.size() - 1);
		v_patterns_shuffled.push_back(v_patterns_buf.at(i_pattern_offset));
		v_patterns_buf.erase(v_patterns_buf.begin() + i_pattern_offset);
	}//while  (v_patterns_buf.size() > 0)


	C3LOIndividual  *pc_best_update;
	C3LOIndividual  *pc_best_different_ind;
	C3LOPattern  *pc_pattern;

	/*FILE  *pf_test;
	pf_test = fopen("zz_init.txt", "w+");
	fprintf(pf_test, "\nSTART:\n");
	vSave(pf_test);//*/

	while (v_patterns_shuffled.size()  >  0)
	{
		pc_pattern = v_patterns_shuffled.at(0);
		v_patterns_shuffled.erase(v_patterns_shuffled.begin());

		pc_best_update = NULL;
		pc_best_different_ind = NULL;
		if (bLinkedGenesBrutalSearch(pc_pattern, iPatternLevel, pvGenePatternPartsOriginal, &pc_best_update, &pc_best_different_ind) == true)  pc_best_update->vCopyTo(this);

		if (pc_best_update != NULL)  delete  pc_best_update;
		if (pc_best_different_ind != NULL)  delete  pc_best_different_ind;

	}//while  (v_patterns_shuffled.size()  >  0)


	 /*fprintf(pf_test, "\nAFTER:\n");
	 vSave(pf_test);
	 fclose(pf_test);//*/
}//void  C3LOIndividual::vIncrementalBrutalUpdate(vector<C3LOPattern> *pvGenePatterns, int iPatternLevel, C3LOIndividual  **pcBestUpdate)



bool  C3LOIndividual::bLinkedGenesBrutalSearch(C3LOPattern *pcLinkedGenes, int iPatternLevel, vector<C3LOPattern *> *pvGenePatternPartsOriginal, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent)
{
	//return(false);
	/*FILE  *pf_test;
	pf_test = fopen("z_bLinkedGenesBrutalSearch.txt", "w+");
	fprintf(pf_test, "%d\n", i_global_counter++);
	fclose(pf_test);*/

	dComputeFitness(iPatternLevel);


	*pcBestButDifferent = NULL;
	*pcBestUpdate = NULL;
	bool  b_result;

	if (iPatternLevel == 0)
	{
		b_result = bLinkedGenesBrutalSearchZero(pcLinkedGenes, pvGenePatternPartsOriginal, pcBestUpdate, pcBestButDifferent);

		/*pf_test = fopen("z_bLinkedGenesBrutalSearch.txt", "a");
		fprintf(pf_test, "ENDED\n", i_global_counter++);
		fclose(pf_test);*/

		pc_parent->i_brutal_random_search++;
		if (b_result == true)  pc_parent->i_brutal_random_search_effective++;


		return(b_result);
	}//if  (iPatternLevel == 0)  

	b_result = bLinkedGenesBrutalSearchHigher(pcLinkedGenes, iPatternLevel, pcBestUpdate, pcBestButDifferent);

	/*pf_test = fopen("z_bLinkedGenesBrutalSearch.txt", "a");
	fprintf(pf_test, "ENDED\n", i_global_counter++);
	fclose(pf_test);*/

	pc_parent->i_brutal_random_search++;
	if (b_result == true)  pc_parent->i_brutal_random_search_effective++;
	return(b_result);

	return(false);
}//void  C3LOIndividual::vLinkedGenesBrutalSearch(C3LOPattern *pcLinkedGenes, int iPatternLevel)



C3LOPattern*  C3LOIndividual::pc_get_pattern_part(int  iGenePosition, int *piFilledGenes, vector<C3LOPattern *> *pvGenePatternPartsOriginal)
{
	vector<C3LOPattern*>  v_possible_patterns;
	bool  b_pattern_possible;

	for (int i_pattern = 0; i_pattern < pvGenePatternPartsOriginal->size(); i_pattern++)
	{
		if (pvGenePatternPartsOriginal->at(i_pattern)->bGeneThere(iGenePosition) == true)
		{
			b_pattern_possible = false;
			for (int ii = 0; (ii < pvGenePatternPartsOriginal->at(i_pattern)->v_pattern.size())&&(b_pattern_possible == false); ii++)
			{
				if  (piFilledGenes[pvGenePatternPartsOriginal->at(i_pattern)->v_pattern.at(ii).iGenePos()] < 1)  b_pattern_possible = true;
			}//for (int ii = 0; (ii < pvGenePatternPartsOriginal->at(i_pattern)->v_pattern.size())&&(b_pattern_possible == false); ii++)

			if (b_pattern_possible == true)  v_possible_patterns.push_back(pvGenePatternPartsOriginal->at(i_pattern));		
		}//if (pvGenePatternPartsOriginal->at(i_pattern)->bGeneThere(iGenePosition) == true)
	}//for (int i_pattern = 0; i_pattern < pvGenePatternPartsOriginal->size(); i_pattern++)

	if (v_possible_patterns.size() == 0)  return(NULL);

	int  i_pat_offset;
	i_pat_offset = RandUtils::iRandNumber(0, v_possible_patterns.size() - 1);

	return(v_possible_patterns.at(i_pat_offset));
}//C3LOPattern*  C3LOIndividual::pc_get_pattern_part(int  iGenePosition, int *piFilledGenes, vector<C3LOPattern *> *pvGenePatternPartsOriginal)



bool  C3LOIndividual::bLinkedGenesBrutalSearchZero(C3LOPattern *pcLinkedGenes, vector<C3LOPattern *> *pvGenePatternPartsOriginal, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent)
{
	if (pcLinkedGenes->pvGetPattern()->size() == 0)  return(false);



	int  *pi_genes_included;
	pi_genes_included = new  int[i_templ_length];

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genes_included[ii] = 0;

	for (int ii = 0; ii < pcLinkedGenes->pvGetPattern()->size(); ii++)
	{
		pi_genes_included[pcLinkedGenes->pvGetPattern()->at(ii).iGenePos()] = 1;
	}//for  (int  ii = 0; ii < pcLinkedGenes->v_pattern.size(); ii++)


	vector<int>  v_genes_included;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_genes_included[ii] > 0)  v_genes_included.push_back(ii);
	}//for  (int  ii = 0; ii < i_templ_length; ii++)


	//use it again for a different purpose - gene marking in brutal search pattern set
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_genes_included[ii] = 0;

	C3LOPattern  *pc_gene_pattern_part;
	vector<int>  v_genes_shuffled;
	vector<CMessyGene>  v_genes_directly_connected;
	int  i_gene_offset, i_gene_position, i_gene_pos_new;
	int  i_current_gene_offset = -1;
	while ((v_genes_included.size() > 0) && (v_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))
	{
		//i_gene_offset = pc_random_gen->Next(0, v_genes_included.size());
		i_gene_offset = RandUtils::iRandNumber(0, v_genes_included.size() - 1);
		i_gene_position = v_genes_included.at(i_gene_offset);
		v_genes_included.erase(v_genes_included.begin() + i_gene_offset);


		if (pi_genes_included[i_gene_position] == 0)
		{
			v_genes_shuffled.push_back(i_gene_position);
			i_current_gene_offset = v_genes_shuffled.size() - 1;
			pi_genes_included[i_gene_position] = 1;

			//now get genes directly dependent on the selected gene
			pc_gene_pattern_part = pc_get_pattern_part(i_gene_position, pi_genes_included, pvGenePatternPartsOriginal);

			while ((pc_gene_pattern_part != NULL)&&(v_genes_directly_connected.size() > 0) && (v_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))
			{
				v_genes_directly_connected = pc_gene_pattern_part->v_pattern;

				while ((v_genes_directly_connected.size() > 0) && (v_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))
				{
					i_gene_offset = RandUtils::iRandNumber(0, v_genes_directly_connected.size() - 1);
					i_gene_pos_new = v_genes_directly_connected.at(i_gene_offset).iGenePos();
					v_genes_directly_connected.erase(v_genes_directly_connected.begin() + i_gene_offset);

					if (pi_genes_included[i_gene_pos_new] == 0)
					{
						v_genes_shuffled.push_back(i_gene_pos_new);
						pi_genes_included[i_gene_pos_new] = 1;
					}//if (pi_genes_included[i_gene_pos_new] == 0)
				}//while ((v_genes_directly_connected.size() > 0) && (v_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))

				pc_gene_pattern_part = pc_get_pattern_part(i_gene_position, pi_genes_included, pvGenePatternPartsOriginal);

				while ((pc_gene_pattern_part == NULL) && (i_current_gene_offset >= 0))
				{
					i_current_gene_offset++;
					if (i_current_gene_offset < v_genes_shuffled.size())
					{
						i_gene_position = v_genes_shuffled.at(i_current_gene_offset);
						pc_gene_pattern_part = pc_get_pattern_part(i_gene_position, pi_genes_included, pvGenePatternPartsOriginal);
					}//if (i_current_gene_offset < v_genes_shuffled.size())
					else
						i_current_gene_offset = -1;
				}//while ((pc_gene_pattern_part == NULL) && (i_current_gene_offset >= 0))
			}//if  (pc_gene_pattern_part != NULL)
		}//if  (pi_genes_included[i_gene_position] == 0)
	}//while (v_genes_included.size() > 0)	


	



	delete  pi_genes_included;//DONT REMOVE during updates!
	
	/*vector<int>  v_genes_shuffled;
	int  i_gene_offset;
	while ((v_genes_included.size() > 0) && (v_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))
	{
		//i_gene_offset = pc_random_gen->Next(0, v_genes_included.size());
		i_gene_offset = RandUtils::iRandNumber(0, v_genes_included.size()-1);
		v_genes_shuffled.push_back(v_genes_included.at(i_gene_offset));
		v_genes_included.erase(v_genes_included.begin() + i_gene_offset);
	}//while (v_genes_included.size() > 0)*/


	::MessageBox(NULL, "0 level at ind creation 0", "", MB_OK);
	C3LOIndividual  c_coptimized(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
	C3LOIndividual  c_best_and_different(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
	vCopyTo(&c_coptimized);
	vCopyTo(&c_best_and_different);
	int  i_best_but_different_found;

	i_best_but_different_found = 0;
	if (b_brutal_check(&v_genes_shuffled, &c_coptimized, &c_best_and_different, &i_best_but_different_found, 0) == true)
	{
		::MessageBox(NULL, "0 level at ind creation 1", "", MB_OK);
		*pcBestUpdate = new  C3LOIndividual(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
		c_coptimized.vCopyTo(*pcBestUpdate);
		(*pcBestUpdate)->vSetOptimized(false);

		double  d_opt_fitness, d_opt_fitness_bald, d_curr_fit;
		d_opt_fitness = c_coptimized.dComputeFitness(-1);
		d_opt_fitness_bald = (*pcBestUpdate)->dComputeFitness(-1);
		d_curr_fit = dComputeFitness(-1);

		if (d_opt_fitness > d_opt_fitness_bald)
		{
			for (int ii = 0; ii < i_templ_length; ii++)
				((*pcBestUpdate)->pi_genotype)[ii] = c_coptimized.vpi_optimized_genotypes.at(0)[0][ii];

			if (pc_parent->bCanIndBeAddedAtAnyLevel((*pcBestUpdate)) == true)
			{
				(*pcBestUpdate)->vSetOptimized(false);
				return(true);
			}//if (pc_parent->bCanIndBeAddedAtAnyLevel((*pcBestUpdate)) == true)


			return(false);
			//CString  s_buf;

			//s_buf.Format("opt: %.4lf > %.4lf  (curr: %.4lf)", d_opt_fitness, d_opt_fitness_bald, d_curr_fit);
			//::MessageBox(NULL, s_buf, s_buf, MB_OK);
		}//if (d_opt_fitness < d_opt_fitness_bald)

		return(true);
	}//if  (b_brutal_check(&v_genes_shuffled, &c_coptimized, 0) == true)
	else
	{
		if (i_best_but_different_found != 0)
		{
			if (bComparePhenotype(0, 0, &c_best_and_different) == false)
			{
				::MessageBox(NULL, "0 level at ind creation 2", "", MB_OK);
				*pcBestButDifferent = new  C3LOIndividual(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
				c_best_and_different.vCopyTo(*pcBestButDifferent);
				(*pcBestButDifferent)->vSetOptimized(false);
			}//if  (bComparePhenotype(&c_coptimized)  ==  false)
		}//if  (i_best_but_different_found  !=  0)
	}//else  if  (b_brutal_check(&v_genes_shuffled, &c_coptimized, 0) == true)



	return(false);
}//void  C3LOIndividual::vLinkedGenesBrutalSearch(C3LOPattern *pcHighLvlGene)



bool  C3LOIndividual::b_brutal_check(vector<int>  *pvGenesShuffled, C3LOIndividual  *pcOptimizationBuf, C3LOIndividual  *pcBestThatFar, int  *piBestThatFarInitialized, int  iPos)
{
	if (iPos >= pvGenesShuffled->size())
	{
		pcOptimizationBuf->vSetOptimized(false);
		pcOptimizationBuf->dComputeFitness(-1);
		//pcOptimizationBuf->dComputeFitnessOfCurrentPhenotype(0);

		if (pc_parent->bCanIndBeAddedAtAnyLevel(pcOptimizationBuf) == true)
		{
			if (*piBestThatFarInitialized == 0)
			{
				*piBestThatFarInitialized = 1;
				pcOptimizationBuf->vCopyTo(pcBestThatFar);
			}//if  (*piBestThatFarInitialized  ==  0)
			else
			{
				if (pcBestThatFar->dComputeFitness(-1)  <  pcOptimizationBuf->dComputeFitness(-1))  pcOptimizationBuf->vCopyTo(pcBestThatFar);
			}//else  if  (*piBestThatFarInitialized  ==  0)

			if (dComputeFitness(-1)  <  pcOptimizationBuf->dComputeFitness(-1))  return(true);
		}//if  (b_can_ind_be_added_at_any_level(pcOptimizationBuf) == true)

		return(false);
	}//if  (iPos  >= pvGenesShuffled->size())

	int  i_val;
	//i_val = pc_random_gen->Next(0, 2);
	i_val = RandUtils::iRandNumber(0,1);

	pcOptimizationBuf->pi_genotype[pvGenesShuffled->at(iPos)] = i_val;
	pcOptimizationBuf->vpi_optimized_genotypes.at(0)[0][pvGenesShuffled->at(iPos)] = i_val;
	if (b_brutal_check(pvGenesShuffled, pcOptimizationBuf, pcBestThatFar, piBestThatFarInitialized, iPos + 1) == true)  return(true);

	if (i_val == 1)
		i_val = 0;
	else
		i_val = 1;

	pcOptimizationBuf->pi_genotype[pvGenesShuffled->at(iPos)] = i_val;
	pcOptimizationBuf->vpi_optimized_genotypes.at(0)[0][pvGenesShuffled->at(iPos)] = i_val;
	if (b_brutal_check(pvGenesShuffled, pcOptimizationBuf, pcBestThatFar, piBestThatFarInitialized, iPos + 1) == true)  return(true);


	return(false);
}//bool  C3LOIndividual::b_brutal_check(vector<int>  *pvGenesShuffled, C3LOIndividual  *pcOptimizationBuf, int  iPos)




bool  C3LOIndividual::bLinkedGenesBrutalSearchHigher(C3LOPattern *pcLinkedGenes, int iPatternLevel, C3LOIndividual  **pcBestUpdate, C3LOIndividual  **pcBestButDifferent)
{
	if (pc_parent->pvGetHigherLevelTreesGenes()->size() <= iPatternLevel)  return(false);
	if (pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).size() == 0)  return(false);



	vector<C3LOPattern  *>  v_high_level_genes_included;
	int  i_common, i_uncommon, i_unmarked;
	for (int ii = 0; ii < pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).size(); ii++)
	{
		pcLinkedGenes->vGetCommonAndUNcommonPositions
		(
			pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).at(ii),
			&i_common, &i_uncommon, &i_unmarked
		);

		if (
			(i_common  >  0) && (i_unmarked == 0)
			)
			v_high_level_genes_included.push_back(pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).at(ii));
	}//for  (int  ii = 0; ii < pc_parent->pvGetHigherLevelTreesGenes()->at(iPatternLevel).size(); ii++)


	vector<C3LOPattern  *>  v_high_level_genes_shuffled;
	int  i_gene_offset;
	while ((v_high_level_genes_included.size() > 0) && (v_high_level_genes_shuffled.size() < C3LO_INDIVIDUAL_MAX_GENES_FOR_BS))
	{
		//i_gene_offset = pc_random_gen->Next(0, v_high_level_genes_included.size());
		i_gene_offset = RandUtils::iRandNumber(0, v_high_level_genes_included.size() - 1);
		v_high_level_genes_shuffled.push_back(v_high_level_genes_included.at(i_gene_offset));
		v_high_level_genes_included.erase(v_high_level_genes_included.begin() + i_gene_offset);
	}//while (v_genes_included.size() > 0)


	 /*FILE  *pf_test;
	 pf_test = fopen("z_brutal_test.txt", "a");

	 fprintf(pf_test, "\n\nBEFORE:\n");
	 vSave(pf_test);

	 fprintf(pf_test, "\n\nPATTERN:\n");
	 pcLinkedGenes->vSavePattern(pf_test, pc_fitness, "", true);//*/


	::MessageBox(NULL, "0 level at ind creation 3", "", MB_OK);
	C3LOIndividual  c_coptimized(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
	C3LOIndividual  c_best_and_different(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
	vCopyTo(&c_coptimized);
	vCopyTo(&c_best_and_different);
	int  i_best_but_different_found;



	i_best_but_different_found = 0;
	if (b_brutal_check_high_level(&v_high_level_genes_shuffled, &c_coptimized, &c_best_and_different, &i_best_but_different_found, iPatternLevel, 0) == true)
	{
		::MessageBox(NULL, "0 level at ind creation 4", "", MB_OK);
		*pcBestUpdate = new  C3LOIndividual(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);

		c_coptimized.vCopyTo(*pcBestUpdate);

		/*fprintf(pf_test, "\n\nOPTIMIZED AFTER:\n");
		vSave(pf_test);
		fclose(pf_test);*/

		return(true);
	}//if  (b_brutal_check(&v_genes_shuffled, &c_coptimized, 0) == true)
	else
	{
		if (i_best_but_different_found != 0)
		{
			if (bComparePhenotype(0,0, &c_best_and_different) == false)
			{
				::MessageBox(NULL, "0 level at ind creation 5", "", MB_OK);
				*pcBestButDifferent = new  C3LOIndividual(0, i_templ_length, pc_problem, pv_optimization_orders, pc_parent);
				c_best_and_different.vCopyTo(*pcBestButDifferent);
			}//if  (bComparePhenotype(&c_coptimized)  ==  false)
		}//if  (i_best_but_different_found  !=  0)
	}//else  if  (b_brutal_check(&v_genes_shuffled, &c_coptimized, 0) == true)

	return(false);
}//void  C3LOIndividual::vLinkedGenesBrutalSearchHigher(C3LOPattern *pcLinkedGenes, int iPatternLevel)




bool  C3LOIndividual::b_brutal_check_high_level(vector<C3LOPattern  *>  *pvHighLevelGenesShuffled, C3LOIndividual  *pcOptimizationBuf, C3LOIndividual  *pcBestThatFar, int  *piBestThatFarInitialized, int iPatternLevel, int  iPos)
{
	if (iPos >= pvHighLevelGenesShuffled->size())
	{
		pcOptimizationBuf->vSetOptimized(false);
		pcOptimizationBuf->dComputeFitness(iPatternLevel);

		if (pc_parent->bCanIndBeAddedAtAnyLevel(pcOptimizationBuf) == true)
		{
			if (*piBestThatFarInitialized == 0)
			{
				*piBestThatFarInitialized = 1;
				pcOptimizationBuf->vCopyTo(pcBestThatFar);
			}//if  (*piBestThatFarInitialized  ==  0)
			else
			{
				if (pcBestThatFar->dComputeFitness(iPatternLevel)  <  pcOptimizationBuf->dComputeFitness(iPatternLevel))  pcOptimizationBuf->vCopyTo(pcBestThatFar);
			}//else  if  (*piBestThatFarInitialized  ==  0)

			if (dComputeFitness(iPatternLevel)  <  pcOptimizationBuf->dComputeFitness(iPatternLevel))  return(true);
		}//if  (b_can_ind_be_added_at_any_level(pcOptimizationBuf) == true)

		return(false);
	}//if  (iPos  >= pvGenesShuffled->size())





	if (pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->size() == 0)  return(false);


	if (pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->size() == 1)
	{
		pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->at(0)->vInfect(pcOptimizationBuf);
		if (b_brutal_check_high_level(pvHighLevelGenesShuffled, pcOptimizationBuf, pcBestThatFar, piBestThatFarInitialized, iPatternLevel, iPos + 1) == true)  return(true);
	}//if  (pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->size() == 1)


	if (pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->size() > 1)
	{
		int  i_val;
		//i_val = pc_random_gen->Next(0, 2);
		i_val = RandUtils::iRandNumber(0,1);

		pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->at(i_val)->vInfect(pcOptimizationBuf);
		if (b_brutal_check_high_level(pvHighLevelGenesShuffled, pcOptimizationBuf, pcBestThatFar, piBestThatFarInitialized, iPatternLevel, iPos + 1) == true)  return(true);


		if (i_val == 1)
			i_val = 0;
		else
			i_val = 1;


		pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->at(i_val)->vInfect(pcOptimizationBuf);
		if (b_brutal_check_high_level(pvHighLevelGenesShuffled, pcOptimizationBuf, pcBestThatFar, piBestThatFarInitialized, iPatternLevel, iPos + 1) == true)  return(true);

	}//if  (pvHighLevelGenesShuffled->at(iPos)->pvGetHighLevelGeneFreq()->size() == 1)



	return(false);
}//bool  C3LOIndividual::b_brutal_check(vector<int>  *pvGenesShuffled, C3LOIndividual  *pcOptimizationBuf, int  iPos)








 //---------------------------------------------C3LOPattern-------------------------------------------------------

C3LOPattern::C3LOPattern(int iTemplLength)
{
	i_templ_length = iTemplLength;
	i_pattern_level = 0;

	pi_pattern_table = NULL;

	i_number_of_hits = 0;
	d_dsm_value = 0;
}//C3LOPattern::C3LOPattern()


C3LOPattern::~C3LOPattern()
{
	/*FILE  *pf_test;
	pf_test = fopen("zz_pat_del", "a+");
	vSavePattern(pf_test, NULL, "", true);
	fclose(pf_test);*/

	for  (int ii = 0; ii < v_children.size(); ii++)
	  delete  v_children.at(ii);


	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
		delete  v_high_level_gene_freq.at(ii);

	if (pi_pattern_table != NULL)  delete  pi_pattern_table;
}//C3LOPattern::~C3LOPattern()



 //common: we have both the same; uncommon: this has it, pcOther doues not; Unmarked: this does not have it, pcOther does have it
void  C3LOPattern::vGetCommonAndUNcommonPositions(C3LOPattern  *pcOther, int *piCommonPos, int *piUncommonPos, int *piUnmarkedPos)
{
	int  *pi_differences;
	pi_differences = new int[i_templ_length];

	for (int ii = 0; ii < i_templ_length; ii++)  pi_differences[ii] = 0;

	for (int ii = 0; ii < pcOther->v_pattern.size(); ii++)
	{
		pi_differences[pcOther->v_pattern.at(ii).iGenePos()] = 1;
	}//for  (int  ii = 0; ii < v_pattern.size(); ii++)

	vGetCommonAndUNcommonPositions(pi_differences, piCommonPos, piUncommonPos, piUnmarkedPos);

	delete  pi_differences;
}//void  C3LOPattern::vGetCommonAndUNcommonPositions(C3LOPattern  *pcOther, int *piCommonPos, int *piUncommonPos, int *piUnmarkedPos)


void  C3LOPattern::vGetCommonAndUNcommonPositions(int  *piDifferences, int *piCommonPos, int *piUncommonPos, int *piUnmarkedPos)
{
	int  *pi_tool;
	int  i_uncommon_pos, i_common_pos;

	pi_tool = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_tool[ii] = piDifferences[ii];



	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (pi_tool[v_pattern.at(ii).iGenePos()] == 0)
		{
			pi_tool[v_pattern.at(ii).iGenePos()] = -1;
		}//if  (pi_tool[v_pattern.at(ii).iGenePos()] == 0)


		if (pi_tool[v_pattern.at(ii).iGenePos()] == 1)
		{
			pi_tool[v_pattern.at(ii).iGenePos()] = 2;
		}//if  (pi_tool[v_pattern.at(ii).iGenePos()] == 1)

	}//for  (int  i_my_offset = 0; i_my_offset < v_pattern.size(); i_my_offset++)


	*piCommonPos = 0;
	*piUncommonPos = 0;
	*piUnmarkedPos = 0;

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_tool[ii] == -1) (*piUncommonPos)++;
		if (pi_tool[ii] == 1) (*piUnmarkedPos)++;
		if (pi_tool[ii] == 2) (*piCommonPos)++;
	}//for  (int ii = 0; ii < i_templ_length; ii++)


	delete  pi_tool;

}//int  C3LOPattern::iGetUNcommonPositions(int  *piDifferences)



int  C3LOPattern::iCommonPositions(C3LOPattern  *pcOther, bool  bCheckValues /*= false*/, bool  bAllMustBeTheSame /*= true*/)
{
	if (bAllMustBeTheSame == true)
	{
		if (v_pattern.size() != pcOther->v_pattern.size())  return(0);
	}//if  (bBreakAtTheFirstDifference == true)

	int  i_result;
	bool  b_found;

	i_result = 0;
	for (int i_my_offset = 0; i_my_offset < v_pattern.size(); i_my_offset++)
	{
		b_found = false;
		for (int i_other_offset = 0; (i_other_offset < pcOther->v_pattern.size()) && (b_found == false); i_other_offset++)
		{
			if (bCheckValues == true)
			{
				if (
					(v_pattern.at(i_my_offset).iGenePos() == pcOther->v_pattern.at(i_other_offset).iGenePos()) &&
					(v_pattern.at(i_my_offset).iGeneVal() == pcOther->v_pattern.at(i_other_offset).iGeneVal())
					)
					b_found = true;
			}//if  (bCheckValues == true)
			else
			{
				if (
					(v_pattern.at(i_my_offset).iGenePos() == pcOther->v_pattern.at(i_other_offset).iGenePos())
					)
					b_found = true;
			}//else  if  (bCheckValues == true)
		}//for  (int  i_other_offset = 0; (i_other_offset < pcOther->v_pattern.size())&&(b_found = false); i_other_offset++)

		if (b_found == false)
		{
			if (bAllMustBeTheSame == true)
			{
				return(0);
			}//if  (bBreakAtTheFirstDifference == true)
		}//if  (b_found == false)
		else
		{
			i_result++;
		}//else  if  (b_found == false)
	}//for  (int  i_my_offset = 0; i_my_offset < v_pattern.size(); i_my_offset++)


	return(i_result);
}//bool  C3LOPattern::bIsTheSame(C3LOPattern  *pcOther)




int  C3LOPattern::iDifferentPositions(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int iPatternLevel)
{
	int  i_result;


	pcInd0->dComputeFitness(iPatternLevel);//phenotype creation forcing...
	pcInd1->dComputeFitness(iPatternLevel);//phenotype creation forcing...

	int  *pi_phenotype_0, *pi_phenotype_1;
	pi_phenotype_0 = pcInd0->piGetPhenotype(0, 0);
	pi_phenotype_1 = pcInd1->piGetPhenotype(0, 0);

	i_result = 0;
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (pi_phenotype_0[v_pattern.at(ii).iGenePos()] != pi_phenotype_1[v_pattern.at(ii).iGenePos()])  i_result++;
	}//for  (int  ii = 0; ii < v_pattern.size(); ii++)

	return(i_result);
}//int  C3LOPattern::iDifferentPositions(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1)






void  C3LOPattern::v_check_high_level_gene_for_ind(C3LOIndividual *pcInd, int iPatternLevel)
{
	/*if (iPatternLevel == 1)
	{
		if (v_high_level_gene_freq.size() > 0)  return;

		C3LOHighLvlGeneFreq  *pc_new_gene_val_freq;

		int  i_first_val = RandUtils::iRandNumber(0,1);

		pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
		pc_new_gene_val_freq->vCreateTest(this, i_first_val, iPatternLevel);
		v_high_level_gene_freq.push_back(pc_new_gene_val_freq);

		for (int ii = 0; ii < 20; ii++)
		{
			pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
			pc_new_gene_val_freq->vCreateTest(this, -1, iPatternLevel);
			v_high_level_gene_freq.push_back(pc_new_gene_val_freq);
		}//for (int ii = 0; ii < 20; ii++)

		pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
		pc_new_gene_val_freq->vCreateTest(this, i_first_val*(-1) + 1, iPatternLevel);
		v_high_level_gene_freq.push_back(pc_new_gene_val_freq);

		return;
	}//if (iPatternLevel == 1)*/


	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	{
		if (v_high_level_gene_freq.at(ii)->bCheckHighLevelGene(pcInd, iPatternLevel) == true)  return;
	}//for  (int  ii = 0; ii < v_high_level_gene_freq.size(); ii++)

	 //create high level gene
	C3LOHighLvlGeneFreq  *pc_new_gene_val_freq;

	pc_new_gene_val_freq = new C3LOHighLvlGeneFreq(i_templ_length);
	pc_new_gene_val_freq->vCreate(this, pcInd, iPatternLevel);

	v_high_level_gene_freq.push_back(pc_new_gene_val_freq);
}//void  C3LOPattern::v_check_high_level_gene_for_ind(C3LOIndividual *pcInd)



void  C3LOPattern::vGetBestCovering(vector<C3LOPattern*> *pvDiffrCoveringPatterns, int  *piDifferences, int iAcceptableLevel)
{
	if (i_pattern_level <= iAcceptableLevel)
	{
		int  i_common_pos, i_uncommon_pos, i_unmarked_pos;

		if (piDifferences != NULL)
		{
			vGetCommonAndUNcommonPositions(piDifferences, &i_common_pos, &i_uncommon_pos, &i_unmarked_pos);

			//if  (i_uncommon_pos == 0)
			if (i_common_pos > 0)
			{
				pvDiffrCoveringPatterns->push_back(this);
				return;
			}//if  (i_uncommon_pos == 0)
		}//	if  (piDifferences  !=  NULL)
		else
		{
			pvDiffrCoveringPatterns->push_back(this);
			return;
		}//else  if  (piDifferences  !=  NULL)
	}//if  (i_pattern_level <= iAcceptableLevel)


	 //for  (int i_child = 0; i_child < v_children.size(); i_child++)
	 //	v_children.at(i_child)->vGetBestCovering(pvDiffrCoveringPatterns, piDifferences, iAcceptableLevel);	

}//void  C3LOPattern::vGetBestCovering(vector<C3LOPattern*> *pvDiffrCoveringPatterns, int  *piDifferences, int iAcceptableLevel)




C3LOPattern*  C3LOPattern::pcGetLongestReasonableBranch(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int *piDifferenceSearchTool, int iPatternLevel)
{
	if (bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == false)  return(NULL);
	if (bUnMarkedPhenotypeDifferences(pcInd0->piGetPhenotype(iPatternLevel, 0), pcInd1->piGetPhenotype(iPatternLevel, 0), piDifferenceSearchTool) == true)  return(this);


	vector<C3LOPattern *>  v_candidates;
	C3LOPattern *pc_candidate;


	for (int i_child = 0; i_child < v_children.size(); i_child++)
	{
		pc_candidate = v_children.at(i_child)->pcGetLongestReasonableBranch(pcInd0, pcInd1, piDifferenceSearchTool, iPatternLevel);
		if (pc_candidate != NULL)  v_candidates.push_back(pc_candidate);
	}//for  (int i_child = 0; i_child < v_children.size(); i_child++)


	/*vector<C3LOPattern *>  v_candidates;

	for (int i_child = 0; i_child < v_children.size(); i_child++)
	{
		if (v_children.at(i_child)->bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == true)
		{
			if (v_children.at(i_child)->bUnMarkedPhenotypeDifferences(pcInd0, pcInd1, piDifferenceSearchTool, iPatternLevel) == true)  v_candidates.push_back(v_children.at(i_child));
		}//if (v_children.at(i_child)->bMarkedPhenotypeDifferences(pcInd0, pcInd1, iPatternLevel) == true)		
	}//for  (int i_child = 0; i_child < v_children.size(); i_child++)*/

	if (v_candidates.size() > 0)
	{
		int  i_chosen_candidate;
		i_chosen_candidate = RandUtils::iRandNumber(0, v_candidates.size() - 1);

		return(v_candidates.at(i_chosen_candidate));
	}//if (v_candidates.size() > 0)


	return(NULL);
}//C3LOPattern  C3LOPattern::pcGetLongestReasonableBranch(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int iPatternLevel)



void  C3LOPattern::vGetRandomTreeBranch(int  *piPhenotype0, int *piPhenotype1, vector<C3LOPattern  *>  *pvTreeBranch, int iPatternLevel)
{
	if (pvTreeBranch == NULL)  return;

	if (pvTreeBranch->size() == 0)
	{
		if (bMarkedPhenotypeDifferences(piPhenotype0, piPhenotype1) == false)  return;
		pvTreeBranch->push_back(this);
	}//if  (pvTreeBranch->size() ==  0)

	vector<C3LOPattern *>  v_allowed_children;
	for (int i_child = 0; i_child < v_children.size(); i_child++)
	{
		if (v_children.at(i_child)->bMarkedPhenotypeDifferences(piPhenotype0, piPhenotype1) == true)  v_allowed_children.push_back(v_children.at(i_child));
	}//for  (int i_child = 0; i_child < v_children.size(); i_child++)

	if (v_allowed_children.size() == 0)  return;

	int  i_chosen_node;
	//i_chosen_node = pc_random_gen->Next(0, v_allowed_children.size());
	i_chosen_node = RandUtils::iRandNumber(0, v_allowed_children.size() - 1);
	if (i_chosen_node >= v_allowed_children.size())  i_chosen_node = v_allowed_children.size() - 1;

	pvTreeBranch->push_back(v_allowed_children.at(i_chosen_node));
	v_allowed_children.at(i_chosen_node)->vGetRandomTreeBranch(piPhenotype0, piPhenotype1, pvTreeBranch, iPatternLevel);
}//void  C3LOPattern::vGetRandomTreeBranch(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, vector<C3LOPattern  *>  *pvTreeBranch)




C3LOPattern  *C3LOPattern::pcGetBestCrossingNode(int iPatternLevel, C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int *piDiffrentPositions, double  *pdDifrPosPerc)
{
	int  i_diffr_pos_my, i_diffr_pos_other, i_diffr_pos_best;
	double  d_diffr_pos_perc_my, d_diffr_pos_perc_other, d_diffr_pos_perc_best;

	i_diffr_pos_my = iDifferentPositions(pcInd0, pcInd1, iPatternLevel);
	d_diffr_pos_perc_my = i_diffr_pos_my;
	d_diffr_pos_perc_my = d_diffr_pos_perc_my / v_pattern.size();


	C3LOPattern  *pc_node_best, *pc_node_other;

	pc_node_best = this;
	i_diffr_pos_best = i_diffr_pos_my;
	d_diffr_pos_perc_best = d_diffr_pos_perc_my;


	for (int i_child = 0; i_child < v_children.size(); i_child++)
	{
		pc_node_other = v_children.at(i_child)->pcGetBestCrossingNode(iPatternLevel, pcInd0, pcInd1, &i_diffr_pos_other, &d_diffr_pos_perc_other);

		if (d_diffr_pos_perc_other > d_diffr_pos_perc_best)
		{
			pc_node_best = pc_node_other;
			i_diffr_pos_best = i_diffr_pos_other;
			d_diffr_pos_perc_best = d_diffr_pos_perc_other;
		}//if  (d_diffr_pos_perc_other > d_diffr_pos_perc_best)		
	}//for  (int i_child = 0; i_child < v_children.size(); i_child++)


	if (piDiffrentPositions != NULL)  *piDiffrentPositions = i_diffr_pos_best;
	if (pdDifrPosPerc != NULL)  *pdDifrPosPerc = d_diffr_pos_perc_best;


	return(pc_node_best);
}//C3LOPattern  *C3LOPattern::pcGetBestCrossingNode(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1)




bool  C3LOPattern::bUnMarkedPhenotypeDifferences(int  *piPhenotype0, int *piPhenotype1, int *piDifferenceSearchTool)
{

	for (int ii = 0; ii < i_templ_length; ii++)
		piDifferenceSearchTool[ii] = -1;

	for (int ii = 0; ii < v_pattern.size(); ii++)
		piDifferenceSearchTool[v_pattern.at(ii).iGenePos()] = 1;
	
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		//we consider genes unmarked by tree
		if (piDifferenceSearchTool[ii] < 0)
		{
			if (piPhenotype0[ii] != piPhenotype1[ii])  return(true);
		}//if (piDifferenceSearchTool[ii] > 0)
	}//for (int ii = 0; ii < i_templ_length; ii++)

	return(false);
}//bool  C3LOPattern::bUnMarkedPhenotypeDifferences(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int *piDifferenceSearchTool, int iPatternLevel)



bool  C3LOPattern::bMarkedPhenotypeDifferences(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1, int iPatternLevel)
{
	return
		(
			bMarkedPhenotypeDifferences(pcInd0->piGetPhenotype(iPatternLevel, 0), pcInd1->piGetPhenotype(iPatternLevel, 0))
		);
}//bool  C3LOPattern::bMarkedPhenotypeDifferences(int  *piGenotype, C3LOIndividual  *pcInd1, int iPatternLevel);


bool  C3LOPattern::bMarkedPhenotypeDifferences(int  *piPhenotype0, int *piPhenotype1)
{
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piPhenotype0[v_pattern.at(ii).iGenePos()] != piPhenotype1[v_pattern.at(ii).iGenePos()])  return(true);
	}//for  (int  ii = 0; ii < v_pattern.size(); ii++)

	return(false);
}//bool  C3LOPattern::bMarkedPhenotypeDifferences(C3LOIndividual *pcInd0, C3LOIndividual  *pcInd1)





void  C3LOPattern::vJoinSinglePatternAdd(C3LOPattern  *pcPatOther)
{
	int  i_my_size_on_start;
	i_my_size_on_start = v_pattern.size();

	bool  b_found;
	for (int i_pat_other_offset = 0; i_pat_other_offset < pcPatOther->v_pattern.size(); i_pat_other_offset++)
	{
		b_found = false;

		for (int i_my_offset = 0; (i_my_offset < i_my_size_on_start) && (b_found == false); i_my_offset++)
		{
			if (v_pattern.at(i_my_offset).iGenePos() == pcPatOther->v_pattern.at(i_pat_other_offset).iGenePos())  b_found = true;
		}//for  (int  i_pat1_offset = 0; (i_pat1_offset < pcPat1->v_pattern.size())&&(b_found == false); i_pat1_offset++)

		if (b_found == false)  v_pattern.push_back(pcPatOther->v_pattern.at(i_pat_other_offset));

	}//for  (int  i_pat0_offset = 0; i_pat0_offset < pcPat0->v_pattern.size(); i_pat0_offset++)

	
	if (pi_pattern_table != NULL)
	{
		delete  pi_pattern_table;
		pi_pattern_table = NULL;
	}//if (pi_pattern_table != NULL)


	v_children.push_back(pcPatOther);
	pcPatOther->v_parents.push_back(this);
}//void  C3LOPattern::vJoinPattern(C3LOPattern  *pcPatOther)



void  C3LOPattern::vJoinTwoPatterns(C3LOPattern  *pcPat0, C3LOPattern  *pcPat1)
{
	if (v_pattern.size() == 0)
	{
		for (int i_pat0_offset = 0; i_pat0_offset < pcPat0->v_pattern.size(); i_pat0_offset++)
		{
			v_pattern.push_back(pcPat0->v_pattern.at(i_pat0_offset));
		}//for  (int  i_pat0_offset = 0; i_pat0_offset < pcPat0->v_pattern.size(); i_pat0_offset++)

		v_children.push_back(pcPat0);
		pcPat0->v_parents.push_back(this);
	}//if (v_pattern.size() == 0)
	else
		vJoinSinglePatternAdd(pcPat0);


	vJoinSinglePatternAdd(pcPat1);

	i_pattern_level = pcPat0->i_pattern_level;
	if (i_pattern_level < pcPat1->i_pattern_level)  i_pattern_level = pcPat1->i_pattern_level;

}//void  C3LOPattern::vJoinTwoPatterns(C3LOPattern  *pcPat0, C3LOPattern  *pcPat1)



bool  C3LOPattern::bDoITouchGeneGroup(int  *piGeneGroup)
{
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piGeneGroup[v_pattern.at(ii).iGenePos()] == 1)  return(true);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bDoITouchGeneGroup(int  *piGeneGroup)




bool  C3LOPattern::bAmIInsideGeneGroup(int  *piGeneGroup)
{
	bool  b_inside;

	b_inside = false;

	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piGeneGroup[v_pattern.at(ii).iGenePos()] == 1)
			b_inside = true;
		else
			return(false);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	if (b_inside == true)  return(true);
	return(false);
}//bool  C3LOPattern::bDoITouchGeneGroup(int  *piGeneGroup)




bool  C3LOPattern::bDoIExtendGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask)
{
	bool  b_any_pos_different, b_any_pos_the_same;
	b_any_pos_different = false;
	b_any_pos_the_same = false;

	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piDifferenceMask[v_pattern.at(ii).iGenePos()] == 1)
		{
			if (piGeneGroup[v_pattern.at(ii).iGenePos()] == 1)
				b_any_pos_the_same = true;
			else
				b_any_pos_different = true;

			if ((b_any_pos_the_same == true) && (b_any_pos_different == true))  return(true);
		}//if (piDifferenceMask[v_pattern.at(ii).iGenePos()] == 1)
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bDoITouchGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask)*/


bool  C3LOPattern::bDoITouchGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask)
{
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piDifferenceMask[v_pattern.at(ii).iGenePos()] == 1)
		{
			if (piGeneGroup[v_pattern.at(ii).iGenePos()] == 1)  return(true);
		}//if (piDifferenceMask[v_pattern.at(ii).iGenePos()] == 1)
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bDoITouchGeneGroupWithMask(int  *piGeneGroup, int *piDifferenceMask)




bool  C3LOPattern::bDoIExtendGeneGroup(int  *piGeneGroup)
{
	bool  b_any_pos_different, b_any_pos_the_same;
	b_any_pos_different = false;
	b_any_pos_the_same = false;

	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piGeneGroup[v_pattern.at(ii).iGenePos()] == 1)  
			b_any_pos_the_same = true;
		else
			b_any_pos_different = true;

		if ((b_any_pos_the_same == true) && (b_any_pos_different == true))  return(true);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bDoIExtendGeneGroup(int  *piGeneGroup)




bool  C3LOPattern::bAmIThere(vector<C3LOPattern  *>  *pvGenePatternsParts, bool  bCheckValues /*= false*/)
{
	for (int ii = 0; ii < pvGenePatternsParts->size(); ii++)
	{
		if (pvGenePatternsParts->at(ii)->iCommonPositions(this) > 0)  return(true);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bAmIThere(vector<C3LOPattern  *>  *pvGenePatternsParts)



bool  C3LOPattern::bGeneThere(int  iGenePosition)
{
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (v_pattern.at(ii).iGenePos() == iGenePosition)  return(true);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(false);
}//bool  C3LOPattern::bGeneThere(int  iGenePosition)




bool  C3LOPattern::bDoesIndividualContainMarkedBlock(int  *piGenotype0, int  *piControlledGenotype)
{
	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (piGenotype0[v_pattern.at(ii).iGenePos()] != piControlledGenotype[v_pattern.at(ii).iGenePos()])  return(false);
	}//for  (int ii = 0; ii < pvGenePatternsParts->size(); ii++)

	return(true);
}//bool  C3LOPattern::bDoesIndividualContainMarkedBlock(int  *piGenotype0, int  *piControlledGenotype)


bool  C3LOPattern::bDoesIndividualContainMarkedBlock(int  *piGenotype0, C3LOIndividual  *pcInd)
{
	return(bDoesIndividualContainMarkedBlock(piGenotype0, pcInd->piGetGenotype()));
}//bool  C3LOPattern::bDoesIndividualContainMarkedBlock(int  *piGenotype0, C3LOIndividual  *pcInd)





void  C3LOPattern::vGeneratePatternTable()
{
	if (pi_pattern_table != NULL)  return;

	pi_pattern_table = new int[i_templ_length];

	for (int ii = 0; ii < i_templ_length; ii++)
		pi_pattern_table[ii] = 0;

	for (int ii = 0; ii < v_pattern.size(); ii++)
		pi_pattern_table[v_pattern.at(ii).iGenePos()] = 1;

}//void  C3LOPattern::vGeneratePatternTable()



bool  C3LOPattern::bEqualToDefinition(int  *piPatternDefinition)
{
	vGeneratePatternTable();

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_pattern_table[ii] != piPatternDefinition[ii])  return(false);
	}//for (int ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOPattern::bEqualToDefinition(int  *piPatternDefinition)



int  C3LOPattern::iGetSubtraction(C3LOPattern *pcOther, int  *piResult)
{
	vGeneratePatternTable();
	pcOther->vGeneratePatternTable();

	int  i_marked_genes_counter;

	i_marked_genes_counter = 0;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		piResult[ii] = 0;


		if ( (pi_pattern_table[ii] == 1)&&(pcOther->pi_pattern_table[ii] == 0) )
		{
			piResult[ii] = 1;
			i_marked_genes_counter++;
		}//if ( (pi_pattern_table[ii] == 1)&&(pcOther->pi_pattern_table[ii] == 0) )
		else
			piResult[ii] = 0;		
	}//for (int ii = 0; ii < i_templ_length; ii++)

	return(i_marked_genes_counter);
}//int  C3LOPattern::iGetSubtraction(C3LOPattern *pcOther, int  *piResult)



int  C3LOPattern::iGetAND(C3LOPattern *pcOther, int  *piResult)
{
	vGeneratePatternTable();
	pcOther->vGeneratePatternTable();

	int  i_marked_genes_counter;

	i_marked_genes_counter = 0;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		piResult[ii] = 0;


		if ((pi_pattern_table[ii] == 1) && (pcOther->pi_pattern_table[ii] == 1))
		{
			piResult[ii] = 1;
			i_marked_genes_counter++;
		}//if ( (pi_pattern_table[ii] == 1)&&(pcOther->pi_pattern_table[ii] == 0) )
		else
			piResult[ii] = 0;
	}//for (int ii = 0; ii < i_templ_length; ii++)

	return(i_marked_genes_counter);

}//int  C3LOPattern::iGetAND(C3LOPattern *pcOther, int  *piResult)



double  C3LOPattern::dGetDSM_Similarity(C3LOPattern  *pcOther, double  **pdDSM)
{
	for (int ii = 0; ii < v_pattern_distances.size(); ii++)
	{
		if (v_pattern_distances.at(ii).pc_first == pcOther)  return(v_pattern_distances.at(ii).d_similarity);
		if (v_pattern_distances.at(ii).pc_second == pcOther)  return(v_pattern_distances.at(ii).d_similarity);
	}//for (int ii = 0; v_pattern_distances.size(); ii++)

	double  d_similarity;

	d_similarity = 0;

	for (int ix = 0; ix < v_pattern.size(); ix++)
	{
		for (int iy = 0; iy < pcOther->v_pattern.size(); iy++)
		{
			if (pdDSM[v_pattern.at(ix).iGenePos()][pcOther->v_pattern.at(iy).iGenePos()] > d_similarity)
				d_similarity = pdDSM[v_pattern.at(ix).iGenePos()][pcOther->v_pattern.at(iy).iGenePos()];
		}//if ((pi_pattern_table[ii] == 1) && (pcOther->pi_pattern_table[ii] == 1))
	}//for (int ii = 0; ii < i_templ_length; ii++)

	C3LOPatternPair  c_new_pair;
	c_new_pair.pc_first = this;
	c_new_pair.pc_second = pcOther;
	c_new_pair.d_similarity = d_similarity;

	v_pattern_distances.push_back(c_new_pair);
	pcOther->v_pattern_distances.push_back(c_new_pair);


	return(d_similarity);
}//double  C3LOPattern::dGetDSM_Distance(C3LOPattern  *pcOther, double  **pdDSM)





double  C3LOPatternPair::dGetDSM_Similarity()
{
	double  d_result;

	//compute pattern length after glueing
	d_result = pc_first->v_pattern.size() + pc_second->v_pattern.size();
	d_result = 1.0 / d_result;

	d_result += d_similarity;

	return(d_result);
}//double  C3LOPatternPair::dGetDSM_Similarity()



double  C3LOPattern::dGetMaxDSM()
{
	double  d_dsm_max, d_dsm_cur;
	d_dsm_max = 0;

	for (int ii = 0; ii < v_pattern_distances.size(); ii++)
	{
		//v_pattern_distances.at(ii).pc_first->vGeneratePatternTable();
		//v_pattern_distances.at(ii).pc_second->vGeneratePatternTable();

		d_dsm_cur = v_pattern_distances.at(ii).dGetDSM_Similarity();
		if (d_dsm_max < d_dsm_cur)  d_dsm_max = d_dsm_cur;
	}//for (int ii = 0; v_pattern_distances.size(); ii++)

	return(d_dsm_max);
}//double  C3LOPattern::dGetMaxDSM()




void  C3LOPattern::vRemoveFromDSM_SimilarityList(C3LOPattern  *pcOther)
{
	for (int ii = 0; ii < v_pattern_distances.size(); ii++)
	{
		if (ii >= 0)
		{
			if (v_pattern_distances.at(ii).pc_first == pcOther)
			{
				v_pattern_distances.erase(v_pattern_distances.begin() + ii);
				ii--;
			}//if (v_pattern_distances.at(ii).pc_first == pcOther)
		}//if  (ii >= 0)

		if (ii >= 0)
		{
			if (v_pattern_distances.at(ii).pc_second == pcOther)
			{
				v_pattern_distances.erase(v_pattern_distances.begin() + ii);
				ii--;
			}//if (v_pattern_distances.at(ii).pc_first == pcOther)
		}//if  (ii >= 0)
	}//for (int ii = 0; v_pattern_distances.size(); ii++)

}//void  C3LOPattern::vRemoveFromDSM_SimilarityList(C3LOPattern  *pcOther)




void  C3LOPattern::vGetAll_DSM_PairsWithSimilarity(double  dDSM_Similarity, vector<C3LOPatternPair> *pvDest)
{
	for (int ii = 0; ii < v_pattern_distances.size(); ii++)
	{
		if (v_pattern_distances.at(ii).dGetDSM_Similarity() == dDSM_Similarity)
		{
			pvDest->push_back(v_pattern_distances.at(ii));
		}//if (v_pattern_distances.at(ii).dGetDSM_Similarity() == dDSM_Similarity)
	}//for (int ii = 0; v_pattern_distances.size(); ii++)

}//void  C3LOPattern::vGetAll_DSM_PairsWithSimilarity(double  dDSM_Max, vector<C3LOPatternPair> *pvDest)



void  C3LOPattern::vCopyPatternAndFreq(C3LOPattern  *pcOther)
{
	v_pattern = pcOther->v_pattern;

	C3LOHighLvlGeneFreq  *pc_new_freq;
	for (int ii = 0; ii < pcOther->v_high_level_gene_freq.size(); ii++)
	{
		pc_new_freq = new C3LOHighLvlGeneFreq(i_templ_length);
		pc_new_freq->vCopy(pcOther->v_high_level_gene_freq.at(ii));
		v_high_level_gene_freq.push_back(pc_new_freq);
	}//for (int ii = 0; ii < pcOther->v_high_level_gene_freq.size(); ii++)

}//void  C3LOPattern::vCopyPatternAndFreq(C3LOPattern  *pcOther)


void  C3LOPattern::vSavePattern(FILE  *pfDest, int  *piProblemBlockStructure, CString  sTab /*=""*/, bool bOnlyRoot /*= false*/)
{
	bool  b_gene_pos_found;
	int  i_gene_val;

	fprintf(pfDest, "%s", sTab);

	fprintf(pfDest, "[%d](%.2lf)", i_pattern_level, d_dsm_value);

	int  i_offset_min, i_offset_max;

	i_offset_min = -1;
	i_offset_max = -1;

	for (int ii = 0; ii < v_pattern.size(); ii++)
	{
		if (i_offset_min == -1)  i_offset_min = v_pattern.at(ii).iGenePos();
		if (i_offset_max == -1)  i_offset_max = v_pattern.at(ii).iGenePos();

		if (i_offset_min  >  v_pattern.at(ii).iGenePos())  i_offset_min = v_pattern.at(ii).iGenePos();
		if (i_offset_max  <  v_pattern.at(ii).iGenePos())  i_offset_max = v_pattern.at(ii).iGenePos();

	}//for  (int  ii = 0; ii < v_pattern.size(); ii++)

	int  *pi_gene_count_tool;
	pi_gene_count_tool = new int[i_templ_length];
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_count_tool[ii] = 0;


	for (int ii = 0; ii < v_pattern.size(); ii++)
		pi_gene_count_tool[v_pattern.at(ii).iGenePos()] = 1;		

	int  i_gene_count = 0;
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_count_tool[ii] == 1)  i_gene_count++;
	}//for  (int  ii = 0; ii < v_pattern.size(); ii++)

	delete pi_gene_count_tool;

	if  (i_gene_count == v_pattern.size())
		fprintf(pfDest, " (%d->%d [%d])", i_offset_min, i_offset_max, v_pattern.size());
	else
		fprintf(pfDest, " DIFFERENT (%d->%d)[%d!=%d]", i_offset_min, i_offset_max, i_gene_count, v_pattern.size());


	fprintf(pfDest, "[HitNo:%d]", i_number_of_hits);




	for (int i_gene_pos = 0; i_gene_pos < i_templ_length; i_gene_pos++)
	{
		b_gene_pos_found = false;
		for (int ii = 0; (ii < v_pattern.size()) && (b_gene_pos_found == false); ii++)
		{
			if (i_gene_pos == v_pattern.at(ii).iGenePos())
			{
				i_gene_val = v_pattern.at(ii).iGeneVal();
				b_gene_pos_found = true;
			}//if  (i_gene_pos == v_disturbance_for_this.at(ii).iGenePos())  
		}//for  (int  ii = 0; (ii < v_distrub_candidates.size())&&(b_gene_pos_found == false); ii++)

		if (b_gene_pos_found == true)
			fprintf(pfDest, "%d", i_gene_val);
		else
			fprintf(pfDest, "*");

		/*if (pcFitness != NULL)
		{
			if (pcFitness->bLastBitOfSubfunction(i_gene_pos) == true)  fprintf(pfDest, " ");
		}//if  (pcFitness != NULL)*/


		if (piProblemBlockStructure != NULL)
		{
			if (piProblemBlockStructure[i_gene_pos] == 1)  fprintf(pfDest, " ");
		}//if (piProblemBlockStructure != NULL)

		
	}//for  (int ii = 0; ii < i_templ_length; ii++)

	fprintf(pfDest, "\n");

	if (bOnlyRoot == true)  return;


	sTab += C3LO_PATTERN_TAB_INCREASE;
	for (int ii = 0; ii < v_children.size(); ii++)
	{
		v_children.at(ii)->vSavePattern(pfDest, piProblemBlockStructure, sTab);
	}//for  (int  ii = 0; ii < v_children.size();  ii++)

}//void  C3LO::v_save_pattern(FILE  *pfDest, vector<CMessyGene> *pvPattern)



void  C3LOPattern::vSaveHighLevelGene(FILE  *pfDest, int  *piProblemBlockStructure)
{
	vSavePattern(pfDest, piProblemBlockStructure);

	for (int ii = 0; ii < v_high_level_gene_freq.size(); ii++)
	{
		v_high_level_gene_freq.at(ii)->vSave(pfDest, " ");
	}//for  (int  ii = 0; ii < v_high_level_gene_freq.size(); ii++)
}//void  C3LOPattern::vSaveHighLevelGene(FILE  *pfDest, CConcatDecFunc *pcFitness)





 //---------------------------------------------C3LOHighLvlGeneFreq-------------------------------------------------------
C3LOHighLvlGeneFreq::C3LOHighLvlGeneFreq(int  iTemplLength)
{
	pi_gene_def = NULL;
	i_templ_length = -1;
	i_gene_number = 0;

	vSetTemplLen(iTemplLength);
}//C3LOHighLvlGeneFreq::C3LOHighLvlGeneFreq()


C3LOHighLvlGeneFreq::~C3LOHighLvlGeneFreq()
{
	if (pi_gene_def != NULL)  delete  pi_gene_def;
}//C3LOHighLvlGeneFreq::~C3LOHighLvlGeneFreq()



void  C3LOHighLvlGeneFreq::vSetTemplLen(int  iTemplLength)
{
	if (iTemplLength <= 0)  return;
	i_templ_length = iTemplLength;

	if (pi_gene_def != NULL)  delete  pi_gene_def;
	pi_gene_def = new int[i_templ_length];
}//void  C3LOHighLvlGeneFreq::vSetTemplLen(int  iTemplLength)



void  C3LOHighLvlGeneFreq::vCopy(C3LOHighLvlGeneFreq  *pcOther)
{
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_def[ii] = pcOther->pi_gene_def[ii];
}//void  C3LOHighLvlGeneFreq::vCopy(C3LOHighLvlGeneFreq  *pcOther)



void  C3LOHighLvlGeneFreq::vSave(FILE  *pfDest, CString  sTab /*= " "*/)
{
	CString  s_line, s_genotype, s_buf;

	s_genotype = "";

	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii]  <  0)
			s_buf = "*";
		else
			s_buf.Format("%d", pi_gene_def[ii]);

		s_genotype += s_buf;
	}//for  (int  ii = 0; ii < i_templ_length; ii++)


	//s_line.Format("%d\t %s\n", i_gene_number, s_genotype);
	s_line.Format("%s\n", s_genotype);
	fprintf(pfDest, s_line);

}//void  C3LOHighLvlGeneFreq::vSave(FILE  *pfDest, CString  sTab /*= " "*/)



void  C3LOHighLvlGeneFreq::vInfect(C3LOIndividual  *pcInd)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii] >= 0)  pcInd->pi_genotype[ii] = pi_gene_def[ii];
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

}//void  C3LOHighLvlGeneFreq::vInfect(C3LOIndividual  *pcInd)



void  C3LOHighLvlGeneFreq::vInfectPhenotype(C3LOIndividual  *pcInd, int iPatternLevel, int iOrder)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii] >= 0)  pcInd->vpi_optimized_genotypes.at(iPatternLevel)[iOrder][ii] = pi_gene_def[ii];
	}//for  (int  ii = 0; ii < i_templ_length; ii++)
}//void  C3LOHighLvlGeneFreq::vInfectPhenotype(C3LOIndividual  *pcInd, int iPhenotype)


void  C3LOHighLvlGeneFreq::vInfectPhenotype(int  *piGenotypeToInfcet)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii] >= 0)  piGenotypeToInfcet[ii] = pi_gene_def[ii];
	}//for  (int  ii = 0; ii < i_templ_length; ii++)
}//void  C3LOHighLvlGeneFreq::vInfectPhenotype(int  *piGenotypeToInfcet)



bool  C3LOHighLvlGeneFreq::bCheckAll(int iValue)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii] >= 0)
		{
			if (pi_gene_def[ii] != iValue)  return(false);
		}//if  (pi_gene_def[ii]  >=  0)
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	return(true);
}//bool  C3LOHighLvlGeneFreq::bCheckAll(int iValue)


bool  C3LOHighLvlGeneFreq::bCheckHighLevelGene(C3LOIndividual *pcInd, int iPatternLevel)
{
	int  i_old_level;
	i_old_level = pcInd->i_level;
	pcInd->i_level = iPatternLevel;
	pcInd->dComputeFitness(iPatternLevel);//force phenotype creation
	pcInd->i_level = i_old_level;

	if (pcInd->piGetPhenotype(iPatternLevel, 0) != NULL)
	{
		bool  b_result;
		b_result = bCheckHighLevelGene(pcInd->piGetPhenotype(iPatternLevel, 0));

		return(b_result);
	}//if (pcInd->piGetPhenotype(iPatternLevel, 0) != NULL)

	return(false);
}//bool  C3LOHighLvlGeneFreq::bCheckHighLevelGene(C3LOIndividual *pcInd)



bool  C3LOHighLvlGeneFreq::bCheckHighLevelGene(int *piGenotypeToCheck)
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		if (pi_gene_def[ii] >= 0)
		{
			if (piGenotypeToCheck[ii] != pi_gene_def[ii])  return(false);
		}//if  (pi_gene_def[ii]  >=  0)
	}//for  (int  ii = 0; ii < i_templ_length; ii++)

	i_gene_number++;
	return(true);
}//bool  C3LOHighLvlGeneFreq::bCheckHighLevelGene(int *piGenotypeToCheck)




void  C3LOHighLvlGeneFreq::vCreateTest(C3LOPattern  *pcPattern, int  iVal, int iPatternLevel)
{
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_def[ii] = -1;

	if ((iVal == 0) || (iVal == 1))
	{
		for (int ii = 0; ii < pcPattern->pvGetPattern()->size(); ii++)
			pi_gene_def[pcPattern->pvGetPattern()->at(ii).iGenePos()] = iVal;
	}//if ((iVal == 0) || (iVal == 1))
	else
	{
		for (int ii = 0; ii < pcPattern->pvGetPattern()->size(); ii++)
			pi_gene_def[pcPattern->pvGetPattern()->at(ii).iGenePos()] = 0;

		int  i_ones, i_rand_off;
		
		i_ones = 0;
		while (i_ones < pcPattern->pvGetPattern()->size() / 2)
		{
			i_rand_off = RandUtils::iRandNumber(0, pcPattern->pvGetPattern()->size() - 1);
			if (pi_gene_def[pcPattern->pvGetPattern()->at(i_rand_off).iGenePos()] == 0)
			{
				pi_gene_def[pcPattern->pvGetPattern()->at(i_rand_off).iGenePos()] = 1;
				i_ones++;
			}//if (pi_gene_def[pcPattern->pvGetPattern()->at(i_rand_off).iGenePos()] == 0)
		}//while (i_ones < pcPattern->pvGetPattern()->size() / 2)
	}//else  if ((iVal == 0) || (iVal == 1))

	i_gene_number = 1;
}//void  C3LOHighLvlGeneFreq::vCreateTest(C3LOPattern  *pcPattern, int  iVal, int iPatternLevel)



void  C3LOHighLvlGeneFreq::vCreate(C3LOPattern  *pcPattern, C3LOIndividual *pcInd, int iPatternLevel)
{
	int  i_old_level;
	i_old_level = pcInd->i_level;
	pcInd->i_level = iPatternLevel;
	pcInd->dComputeFitness(iPatternLevel);//force phenotype creation
	pcInd->i_level = i_old_level;

	vCreate(pcPattern, pcInd->piGetPhenotype(iPatternLevel, 0));
	
	i_gene_number = 1;
}//void  C3LOHighLvlGeneFreq::vCreate(C3LOPattern  *pcPattern, C3LOIndividual)



void  C3LOHighLvlGeneFreq::vCreate(C3LOPattern  *pcPattern, int  *piGenotype)
{
	for (int ii = 0; ii < i_templ_length; ii++)
		pi_gene_def[ii] = -1;
	i_gene_number = 1;

	if (piGenotype == NULL)  return;

	for (int ii = 0; ii < pcPattern->pvGetPattern()->size(); ii++)
		pi_gene_def[pcPattern->pvGetPattern()->at(ii).iGenePos()] = piGenotype[pcPattern->pvGetPattern()->at(ii).iGenePos()];

}//void  C3LOHighLvlGeneFreq::vCreate(C3LOPattern  *pcPattern, int  *piGenotype)









 //---------------------------class  CMessyGene--------------------------------------

CMessyGene::CMessyGene()
{
	i_gene_pos = 0;
	i_gene_val = 0;
}//CMessyGene::CMessyGene()


CMessyGene::CMessyGene(int  iGeneVal, int  iGenePos)
{
	i_gene_pos = iGenePos;
	i_gene_val = iGeneVal;
}//CMessyGene::CMessyGene(int  iGeneVal,  int  iGenePos)


CMessyGene::CMessyGene(CMessyGene  *pcOther)
{
	i_gene_pos = pcOther->i_gene_pos;
	i_gene_val = pcOther->i_gene_val;
}//CMessyGene::CMessyGene(CMessyGene  *pcOther)


CMessyGene::~CMessyGene()
{

};//CMessyGene::~CMessyGene()



int i_func_call(int iArg0, int iArg1)
{
	int *pi_result;

	pi_result = new int();
	*pi_result = iArg0 + iArg1;

	return(*pi_result);
};//int i_func_call(int iArg0, int iArg1)






  //---------------------------------------------CGenotypeInfoNode-------------------------------------------------------
CGenotypeInfoNode::CGenotypeInfoNode()
{
	int  i_value = -1;
	i_succ_opt = 0;
};//CGenotypeInfoNode::CGenotypeInfoNode()

CGenotypeInfoNode::~CGenotypeInfoNode()
{
	for (int ii = 0; ii < v_children.size(); ii++)
		delete  v_children.at(ii);
};//CGenotypeInfoNode::CGenotypeInfoNode()



C3LOIndividual *CGenotypeInfoNode::pcCheckIndividual(C3LOIndividual *pcIndToCheck)
{
	CGenotypeInfoNode  *pc_node;
	pc_node = pc_get_child_node(pcIndToCheck->pi_genotype[0]);

	C3LOIndividual *pc_result;
	pc_result = pc_node->pc_check_individual(pcIndToCheck, 1);

	if (pc_result != NULL)  i_succ_opt++;

	return(pc_result);
};//C3LOIndividual *CGenotypeInfoNode::pcCheckIndividual(C3LOIndividual *pcIndToCheck)



C3LOIndividual *CGenotypeInfoNode::pc_check_individual(C3LOIndividual *pcIndToCheck, int  iPos)
{
	//if (iPos >= pcIndToCheck->pc_fitness->iGetProblemBitLength())
	if (iPos >= pcIndToCheck->pc_problem->pcGetEvaluation()->iGetNumberOfElements())
	{
		if (v_individuals.size()  >  0)
			return(v_individuals.at(0));
		else
		{
			v_individuals.push_back(pcIndToCheck);
			return(NULL);
		}//else  if  (v_individuals.size()  >  0)  
	}//if  (iPos >= pc_fitness->iGetProblemBitLength())

	CGenotypeInfoNode  *pc_node;
	pc_node = pc_get_child_node(pcIndToCheck->pi_genotype[iPos]);

	return(pc_node->pc_check_individual(pcIndToCheck, ++iPos));
};//C3LOIndividual *CGenotypeInfoNode::pc_check_individual(C3LOIndividual *pcIndToCheck, int  iPos)



CGenotypeInfoNode  *CGenotypeInfoNode::pc_get_child_node(int  iVal)
{
	for (int ii = 0; ii < v_children.size(); ii++)
	{
		if (v_children.at(ii)->i_value == iVal)  return(v_children.at(ii));
	}//for  (int  ii = 0; ii < v_children.size(); ii++)

	CGenotypeInfoNode  *pc_new_node;
	pc_new_node = new CGenotypeInfoNode();
	pc_new_node->i_value = iVal;

	v_children.push_back(pc_new_node);

	return(pc_new_node);
};//CGenotypeInfoNode  *CGenotypeInfoNode::pc_get_child_node(int  iVal)

