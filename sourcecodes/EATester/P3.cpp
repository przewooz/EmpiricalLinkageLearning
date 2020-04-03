#include "P3.h"

#include "CommandParam.h"
#include "FloatCommandParam.h"
#include "StringUtils.h"

#include <algorithm>
#include <functional>
#include <vector>

CP3::CP3(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: CP3BasedOptimizer("Pyramid", pcProblem, pcLog, iRandomSeed)
{

}//CP3::CP3(CP3Params *pcParams, CLog *pcLog, uint32_t iRandomSeed)

CP3::CP3(CP3 *pcOther)
	: CP3BasedOptimizer(pcOther)
{
	b_with_feedback = pcOther->b_with_feedback;
	f_feedback_probability = pcOther->f_feedback_probability;
}//CP3::CP3(CP3 *pcOther)

CError CP3::eConfigure(istream *psSettings)
{
	CError c_error = CP3BasedOptimizer::eConfigure(psSettings);

	if (!c_error)
	{
		CBoolCommandParam p_with_feedback("with_feedback");
		b_with_feedback = p_with_feedback.bGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error && b_with_feedback)
	{
		CFloatCommandParam p_feedback_probability("feedback_probability", 0.0f, 1.0f);
		f_feedback_probability = p_feedback_probability.fGetValue(psSettings, &c_error);
	}//if (!c_error && b_with_feedback)

	return c_error;
}//CError CP3::eConfigure(istream *psSettings)

void CP3::vInitialize(time_t tStartTime)
{
	CP3BasedOptimizer::vInitialize(tStartTime);
	pc_p3 = static_pointer_cast<Pyramid>(pc_optimizer);

	i_feedback_counter = 0;
	i_effective_feedback_counter = 0;
}//void CP3::vInitialize(time_t tStartTime)

#include "RandUtils.h"

bool CP3::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_updated = CP3BasedOptimizer::bRunIteration(iIterationNumber, tStartTime);

	//if (iIterationNumber % 1 == 0)
	{
		//v_log_linkage();
	}//if (iIterationNumber % 1 == 0)

	if (b_with_feedback && RandUtils::bSatisfyProbability(f_feedback_probability))
	{
		CBinaryCoding *pc_genotype = pc_best_individual->pcGetGenotype();//new CBinaryCoding(pc_best_individual->pcGetGenotype());

		//for (uint16_t i = 0; i < pc_genotype->iGetNumberOfBits(); i++)
		//{
		//	if (RandUtils::bSatisfyProbability(0.1))
		//	{
		//		*(pc_genotype->piGetBits() + i) ^= 1;
		//	}//if (RandUtils::bSatisfyProbability(0.5))
		//}//for (uint16_t i = 0; i < pc_genotype->iGetNumberOfBits(); i++)

		pc_log->vPrintLine("FEEDBACK");

		b_updated = bClimb(pc_genotype, iIterationNumber, tStartTime) || b_updated;

		pc_log->vPrintLine("END FEEDBACK");

		i_feedback_counter++;
	}//if (b_with_feedback)

	return b_updated;
}//bool CP3::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)

void CP3::vRun()
{
	CP3BasedOptimizer::vRun();

	CString s_feedback_report;
	s_feedback_report.Format("all feedback: %d; effective feedback: %d (%f)", i_feedback_counter, 
		i_effective_feedback_counter, (float)i_effective_feedback_counter / (float)i_feedback_counter);

	pc_log->vPrintLine(s_feedback_report, true);
}//void CP3::vRun()

bool CP3::bClimb(CBinaryCoding *pcGenotype, uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_updated = false;
	
	vector<bool> v_p3_individual;

	v_p3_individual.reserve(pcGenotype->iGetNumberOfBits());

	for (uint16_t i = 0; i < pcGenotype->iGetNumberOfBits(); i++)
	{
		v_p3_individual.push_back(*(pcGenotype->piGetBits() + i) == 1);
	}//for (uint16_t i = 0; i < pcGenotype->iGetNumberOfBits(); i++)

	float f_p3_individual_fitness_value = pc_recorder->evaluate(v_p3_individual);

	vector<size_t> v_added_to_levels;

	if (pc_p3->climb(v_p3_individual, f_p3_individual_fitness_value, v_added_to_levels))
	{
		i_effective_feedback_counter++;

		CIndividual<CBinaryCoding, CBinaryCoding> *pc_previous_best_individual = new CIndividual<CBinaryCoding, CBinaryCoding>(pc_best_individual);

		b_updated = b_update_best_individual(iIterationNumber, tStartTime);

		//if (b_updated)
		{
			pc_log->vPrintEmptyLine();
			pc_log->vPrintEmptyLine();
			pc_log->vPrintLine("EFFECTIVE FEEDBACK");
			pc_log->vPrintLine("PREVIOUS BEST");
			pc_log->vPrintLine(pc_previous_best_individual->pcGetGenotype()->sToString(10));
			pc_log->vPrintLine("NEW BEST");
			pc_log->vPrintLine(pc_best_individual->pcGetGenotype()->sToString(10));
			
			CString s_number_of_levels;
			s_number_of_levels.Format("NUMBER OF LEVELS: %d", pc_p3->pops.size());
			pc_log->vPrintLine(s_number_of_levels);

			CString s_added_to_levels;
			s_added_to_levels.Append("ADDED TO LEVELS:");

			for (size_t i = 0; i < v_added_to_levels.size(); i++)
			{
				s_added_to_levels.AppendFormat(" %d", v_added_to_levels.at(i));
			}//for (size_t i = 0; i < v_added_to_levels.size(); i++)

			pc_log->vPrintLine(s_added_to_levels);

			pc_log->vPrintEmptyLine();
			pc_log->vPrintEmptyLine();
		}//if (b_updated)
	}//if (pc_p3->climb(v_p3_individual, f_p3_individual_fitness_value, v_added_to_levels))

	//CString s_linkage_before_climg = s_create_linkage_report();

	//if (pc_p3->climb(v_p3_individual, f_p3_individual_fitness_value))
	//{
	//	pc_log->vPrintEmptyLine();
	//	pc_log->vPrintEmptyLine();
	//	pc_log->vPrintLine("EFFECTIVE FEEDBACK");
	//	pc_log->vPrintLine("CURRENT BEST");
	//	pc_log->vPrintLine(pc_best_individual->pcGetGenotype()->sToString(10));
	//	pc_log->vPrintLine("LINKAGE BEFORE");
	//	pc_log->vPrintLine(s_linkage_before_climg);
	//	pc_log->vPrintLine("LINKAGE AFTER");
	//	pc_log->vPrintLine(s_create_linkage_report());
	//	pc_log->vPrintEmptyLine();
	//	pc_log->vPrintEmptyLine();

	//	CIndividual<CBinaryCoding, CBinaryCoding> *pc_previous_best_individual = new CIndividual<CBinaryCoding, CBinaryCoding>(pc_best_individual);

	//	b_updated = b_update_best_individual(iIterationNumber, tStartTime);

	//	if (b_updated)
	//	{
	//		pc_log->vPrintEmptyLine();
	//		pc_log->vPrintEmptyLine();
	//		pc_log->vPrintLine("PREVIOUS BEST");
	//		pc_log->vPrintLine(pc_previous_best_individual->pcGetGenotype()->sToString(10));
	//		pc_log->vPrintLine("NEW BEST");
	//		pc_log->vPrintLine(pc_best_individual->pcGetGenotype()->sToString(10));
	//		pc_log->vPrintEmptyLine();
	//		pc_log->vPrintEmptyLine();
	//	}//if (b_updated)

	//	i_effective_feedback_counter++;
	//}//if (pc_p3->climb(v_p3_individual, f_p3_individual_fitness_value))

	return b_updated;
}//bool CP3::bClimb(CBinaryCoding *pcGenotype, time_t tStartTime, uint32_t iIterationNumber)

CString CP3::s_create_linkage_report()
{
	CString s_linkage_report;

	Population *pc_population;
	vector<int> *pv_linkage;

	uint16_t i_linkage_block_size = 10;
	uint16_t i_number_of_linkage_blocks = pc_problem->pcGetEvaluation()->iGetNumberOfElements() / i_linkage_block_size;
	uint16_t i_good_linkage_counter;

	size_t i_linkage_length_sum = 0;

	uint16_t i_linkage_value;
	uint16_t i_linkage_block_start, i_linkage_block_end;

	function<uint16_t(uint16_t)> fi_get_block_start = [=](uint16_t iBlockIndex) { return iBlockIndex * i_linkage_block_size; };
	function<uint16_t(uint16_t)> fi_get_block_end = [=](uint16_t iBlockIndex) { return (iBlockIndex + 1) * i_linkage_block_size - 1; };

	size_t i_number_of_linkages = 0;
	size_t i_number_of_linkages_with_at_least_one_whole_block = 0;

	size_t *pi_only_block_counters = new size_t[i_number_of_linkage_blocks];
	size_t *pi_block_with_other_blocks_counters = new size_t[i_number_of_linkage_blocks];
	size_t *pi_block_with_other_genes_counters = new size_t[i_number_of_linkage_blocks];

	vector<uint16_t> v_whole_blocks_indexes;
	bool b_at_least_one_not_whole_block;

	vector<CString> v_pure_linkages;
	vector<CString> v_dirty_linkages;

	function<void(size_t*)> fv_increase_counter = [&](size_t *piCounter)
	{
		for (size_t i = 0; i < v_whole_blocks_indexes.size(); i++)
		{
			(*(piCounter + v_whole_blocks_indexes.at(i)))++;
		}//for (size_t i = 0; i < v_whole_blocks_indexes.size(); i++)
	};//function<void(size_t*)> fv_increase_counter = [&](size_t *piCounter)

	for (uint16_t i = 0; i < i_number_of_linkage_blocks; i++)
	{
		*(pi_only_block_counters + i) = 0;
		*(pi_block_with_other_blocks_counters + i) = 0;
		*(pi_block_with_other_genes_counters + i) = 0;
	}//for (uint16_t i = 0; i < i_number_of_linkage_blocks; i++)

	for (size_t i = 0; i < pc_p3->pops.size(); i++)
	{
		pc_population = &pc_p3->pops.at(i);

		for (size_t j = 0; j < pc_population->cluster_ordering.size(); j++)
		{
			pv_linkage = &pc_population->clusters.at(pc_population->cluster_ordering.at(j));

			sort(pv_linkage->begin(), pv_linkage->end());

			v_whole_blocks_indexes.clear();
			b_at_least_one_not_whole_block = false;

			for (uint16_t k = 0; k < i_number_of_linkage_blocks; k++)
			{
				i_good_linkage_counter = 0;

				i_linkage_block_start = fi_get_block_start(k);
				i_linkage_block_end = fi_get_block_end(k);

				for (size_t m = 0; m < pv_linkage->size(); m++)
				{
					i_linkage_value = (uint16_t)pv_linkage->at(m);

					if (i_linkage_value >= i_linkage_block_start && i_linkage_value <= i_linkage_block_end)
					{
						i_good_linkage_counter++;
					}//if (i_linkage_value >= i_linkage_block_start && i_linkage_value <= i_linkage_block_end)
				}//for (size_t m = 0; k < pv_linkage->size(); m++)

				if (i_good_linkage_counter == i_linkage_block_size)
				{
					v_whole_blocks_indexes.push_back(k);
				}//if (i_good_linkage_counter == i_linkage_block_size)
				else
				{
					if (i_good_linkage_counter > 0)
					{
						b_at_least_one_not_whole_block = true;
					}//if (i_good_linkage_counter > 0)
				}//else if (i_good_linkage_counter == i_linkage_block_size)
			}//for (uint16_t k = 0; k < i_number_of_linkage_blocks; k++)

			if (b_at_least_one_not_whole_block)
			{
				if (!v_whole_blocks_indexes.empty())		
				{
					fv_increase_counter(pi_block_with_other_genes_counters);

					CString s_linkage;

					for (size_t k = 0; k < v_whole_blocks_indexes.size(); k++)
					{
						s_linkage.AppendFormat("[%d-%d] ", fi_get_block_start(v_whole_blocks_indexes.at(k)), fi_get_block_end(v_whole_blocks_indexes.at(k)));
					}//for (size_t k = 0; k < v_whole_blocks_indexes.size(); k++)

					for (size_t k = 0; k < pv_linkage->size(); k++)
					{
						i_linkage_value = pv_linkage->at(k);

						bool b_within_any_block = false;

						for (size_t m = 0; m < v_whole_blocks_indexes.size() && !b_within_any_block; m++)
						{
							i_linkage_block_start = fi_get_block_start(v_whole_blocks_indexes.at(m));
							i_linkage_block_end = fi_get_block_end(v_whole_blocks_indexes.at(m));
							b_within_any_block = i_linkage_value >= i_linkage_block_start && i_linkage_value <= i_linkage_block_end;
						}//for (size_t m = 0; m < v_whole_blocks_indexes.size() && !b_within_any_block; m++)

						if (!b_within_any_block)
						{
							s_linkage.AppendFormat("%d ", i_linkage_value);
						}//if (!b_within_any_block)
					}//for (size_t k = 0; k < pv_linkage->size(); k++)

					s_linkage.AppendFormat("(population: %d; order: %d; all: %d)", i, j, pc_population->cluster_ordering.size());
					v_dirty_linkages.push_back(s_linkage);
				}//if (!v_whole_blocks_indexes.empty())
			}//if (b_at_least_one_not_whole_block)
			else
			{
				if (v_whole_blocks_indexes.size() == 1)
				{
					fv_increase_counter(pi_only_block_counters);
				}//if (v_whole_blocks_indexes.size() == 1)
				else
				{
					fv_increase_counter(pi_block_with_other_blocks_counters);
				}//else if (v_whole_blocks_indexes.size() == 1)

				CString s_linkage_blocks;

				for (size_t k = 0; k < v_whole_blocks_indexes.size(); k++)
				{
					s_linkage_blocks.AppendFormat("[%d-%d] ", fi_get_block_start(v_whole_blocks_indexes.at(k)), fi_get_block_end(v_whole_blocks_indexes.at(k)));
				}//for (size_t k = 0; k < v_whole_blocks_indexes.size(); k++)

				s_linkage_blocks.AppendFormat("(population: %d; order: %d; all: %d)", i, j, pc_population->cluster_ordering.size());
				v_pure_linkages.push_back(s_linkage_blocks);
			}//else if (b_at_least_one_not_whole_block)

			if (!v_whole_blocks_indexes.empty())
			{
				i_number_of_linkages_with_at_least_one_whole_block++;
			}//if (!v_whole_blocks_indexes.empty())

			i_linkage_length_sum += pv_linkage->size();

			i_number_of_linkages++;
		}//for (size_t j = 0; j < pc_population->cluster_ordering.size(); j++)
	}//for (size_t i = 0; i < pc_p3->pops.size(); i++)

	s_linkage_report.AppendFormat("number of linkages: %d", i_number_of_linkages);
	StringUtils::vAppendEndLine(&s_linkage_report);

	s_linkage_report.AppendFormat("average linkage length: %f", (float)i_linkage_length_sum / (float)i_number_of_linkages);
	StringUtils::vAppendEndLine(&s_linkage_report);

	s_linkage_report.AppendFormat("number of linkages with at least one whole block: %d (%f)", i_number_of_linkages_with_at_least_one_whole_block,
		(float)i_number_of_linkages_with_at_least_one_whole_block / (float)i_number_of_linkages);
	StringUtils::vAppendEndLine(&s_linkage_report);

	s_linkage_report.AppendFormat("number of linkages without blocks: %d (%f)", i_number_of_linkages - i_number_of_linkages_with_at_least_one_whole_block,
		(float)(i_number_of_linkages - i_number_of_linkages_with_at_least_one_whole_block) / (float)i_number_of_linkages);
	StringUtils::vAppendEndLine(&s_linkage_report);

	s_linkage_report.AppendFormat("number of populations: %d", pc_p3->pops.size());
	StringUtils::vAppendEndLine(&s_linkage_report);

	s_linkage_report.AppendFormat("linkage per population: %f", (float)i_number_of_linkages / (float)pc_p3->pops.size());
	StringUtils::vAppendEndLine(&s_linkage_report);

	for (uint16_t i = 0; i < i_number_of_linkage_blocks; i++)
	{
		i_linkage_block_start = fi_get_block_start(i);
		i_linkage_block_end = fi_get_block_end(i);

		s_linkage_report.AppendFormat("[%d-%d] - only: %d (%f); with other whole blocks: %d (%f); with other genes: %d (%f)", i_linkage_block_start,
			i_linkage_block_end, *(pi_only_block_counters + i), (float)*(pi_only_block_counters + i) / (float)i_number_of_linkages, 
			*(pi_block_with_other_blocks_counters + i), (float)*(pi_block_with_other_blocks_counters + i) / (float)i_number_of_linkages,
			*(pi_block_with_other_genes_counters + i), (float)*(pi_block_with_other_genes_counters + i) / (float)i_number_of_linkages);
		StringUtils::vAppendEndLine(&s_linkage_report);
	}//for (uint16_t i = 0; i < i_number_of_linkage_blocks; i++)

	s_linkage_report.Append("pure linkage:");
	StringUtils::vAppendEndLine(&s_linkage_report);

	for (size_t i = 0; i < v_pure_linkages.size(); i++)
	{
		s_linkage_report.Append(v_pure_linkages.at(i));
		StringUtils::vAppendEndLine(&s_linkage_report);
	}//for (size_t i = 0; i < v_pure_linkages.size(); i++)

	s_linkage_report.Append("dirty linkage:");
	StringUtils::vAppendEndLine(&s_linkage_report);

	for (size_t i = 0; i < v_dirty_linkages.size(); i++)
	{
		s_linkage_report.Append(v_dirty_linkages.at(i));
		StringUtils::vAppendEndLine(&s_linkage_report);
	}//for (size_t i = 0; i < v_dirty_linkages.size(); i++)

	delete pi_only_block_counters;
	delete pi_block_with_other_blocks_counters;
	delete pi_block_with_other_genes_counters;

	return s_linkage_report;
}//CString CP3::s_create_linkage_report()

void CP3::v_log_linkage(bool bEcho)
{
	pc_log->vPrintLine(s_create_linkage_report(), bEcho);

	/*Population *pc_population;
	vector<int> *pv_linkage;

	uint16_t i_linkage_block_size = 10;
	uint16_t i_number_of_linkage_blocks = pc_problem->pcGetEvaluation()->iGetNumberOfElements() / i_linkage_block_size;
	uint16_t i_good_linkage_counter;

	uint16_t i_linkage_value;
	uint16_t i_linkage_block_start, i_linkage_block_end;*/

	//vector<CString> v_not_whole_blocks;
	//v_not_whole_blocks.resize(i_number_of_linkage_blocks);

	//for (size_t i = 0; i < pc_p3->pops.size(); i++)
	//{
	//	pc_population = &pc_p3->pops.at(i);

	//	CString s_population_level;
	//	s_population_level.Format("population level: %d", i);

	//	pc_log->vPrintLine(s_population_level, bEcho);

	//	for (size_t j = 0; j < pc_population->cluster_ordering.size(); j++)
	//	{
	//		pv_linkage = &pc_population->clusters.at(pc_population->cluster_ordering.at(j));

	//		v_not_whole_blocks.clear();

	//		CString s_linkage;

	//		for (uint16_t k = 0; k < i_number_of_linkage_blocks; k++)
	//		{
	//			i_good_linkage_counter = 0;

	//			i_linkage_block_start = k * i_linkage_block_size;
	//			i_linkage_block_end = (k + 1) * i_linkage_block_size - 1;

	//			CString s_block_linkage;

	//			for (size_t m = 0; m < pv_linkage->size(); m++)
	//			{
	//				i_linkage_value = (uint16_t)pv_linkage->at(m);

	//				if (i_linkage_value >= i_linkage_block_start && i_linkage_value <= i_linkage_block_end)
	//				{
	//					s_block_linkage.AppendFormat("%d ", i_linkage_value);
	//					i_good_linkage_counter++;
	//				}//if (i_linkage_value >= i_linkage_block_start && i_linkage_value <= i_linkage_block_end)
	//			}//for (size_t m = 0; k < pv_linkage->size(); m++)

	//			if (i_good_linkage_counter == i_linkage_block_size)
	//			{
	//				s_linkage.AppendFormat("%d-%d ", i_linkage_block_start, i_linkage_block_end);
	//			}//if (i_good_linkage_counter == i_linkage_block_size)
	//			else
	//			{
	//				v_not_whole_blocks.push_back(s_block_linkage);
	//			}//else if (i_good_linkage_counter == i_linkage_block_size)
	//		}//for (uint16_t k = 0; k < i_number_of_linkage_blocks; k++)

	//		if (!s_linkage.IsEmpty())
	//		{
	//			for (uint16_t k = 0; k < (uint16_t)v_not_whole_blocks.size(); k++)
	//			{
	//				s_linkage.Append(v_not_whole_blocks.at(k));
	//			}//for (uint16_t k = 0; k < (uint16_t)v_not_whole_blocks.size(); k++)

	//			pc_log->vPrintLine(s_linkage, bEcho);
	//		}//if (!s_linkage.IsEmpty())
	//	}//for (size_t j = 0; j < pc_population->cluster_ordering.size(); j++)
	//}//for (size_t i = 0; i < pc_p3->pops.size(); i++)
}//void CP3::v_log_linkage(bool bEcho)