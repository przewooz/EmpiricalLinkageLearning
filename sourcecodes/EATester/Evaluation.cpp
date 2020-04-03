#include "Evaluation.h"

#include "BinaryCoding.h"
#include "EvaluationUtils.h"
#include "UIntCommandParam.h"

#include <algorithm>
#include <atlstr.h>
#include <cfloat>
#include <utility>

template <class TFenotype>
uint32_t CEvaluation<TFenotype>::iERROR_PARENT_CEVALUATION = CError::iADD_ERROR_PARENT("CEvaluation");

template <class TFenotype>
CEvaluation<TFenotype>::CEvaluation()
{
	v_init(0, DBL_MAX);
}//CEvaluation<TFenotype>::CEvaluation()

template <class TFenotype>
CEvaluation<TFenotype>::CEvaluation(uint16_t iNumberOfElements, double dMaxValue)
{
	v_init(iNumberOfElements, dMaxValue);
}//CEvaluation<TFenotype>::CEvaluation(uint16_t iNumberOfElements, double dMaxValue)

template <class TFenotype>
CEvaluation<TFenotype>::~CEvaluation()
{

}//CEvaluation<TFenotype>::~CEvaluation()

template <class TFenotype>
double CEvaluation<TFenotype>::dEvaluate(TFenotype *pcFenotype)
{
	i_ffe++;

	return d_evaluate(pcFenotype, 0);
}//double CEvaluation<TFenotype>::dEvaluate(TFenotype *pcFenotype)

template <class TFenotype>
void CEvaluation<TFenotype>::v_init(uint16_t iNumberOfElements, double dMaxValue)
{
	i_ffe = 0;
	i_number_of_elements = iNumberOfElements;
	d_max_value = dMaxValue;
}//void CEvaluation<TFenotype>::v_init(uint16_t iNumberOfElements, double dMaxValue)


template <class TFenotype>
CInversedEvaluation<TFenotype>::CInversedEvaluation(CEvaluation<TFenotype> *pcEvaluation)
	: CEvaluation<TFenotype>(pcEvaluation->iGetNumberOfElements(), -pcEvaluation->dGetMaxValue())
{
	pc_evaluation = pcEvaluation;
}//CInversedEvaluation<TFenotype>::CInversedEvaluation(CEvaluation<TFenotype> *pcEvaluation)

template <class TFenotype>
CInversedEvaluation<TFenotype>::~CInversedEvaluation()
{
	delete pc_evaluation;
}//CInversedEvaluation<TFenotype>::~CInversedEvaluation()

template <class TFenotype>
double CInversedEvaluation<TFenotype>::d_evaluate(TFenotype * pcFenotype, uint16_t iShift)
{
	return -pc_evaluation->d_evaluate(pcFenotype, iShift);
}//double CInversedEvaluation<TFenotype>::d_evaluate(TFenotype * pcFenotype, uint16_t iShift)


template <class TFenotype>
CConcatenationEvaluation<TFenotype>::CConcatenationEvaluation()
	: CEvaluation<TFenotype>(0, 0)
{
	pc_sample_fenotype = nullptr;
}//CConcatenationEvaluation<TFenotype>::CConcatenationEvaluation()

template <class TFenotype>
CConcatenationEvaluation<TFenotype>::~CConcatenationEvaluation()
{
	v_clear();
}//CConcatenationEvaluation<TFenotype>::~CConcatenationEvaluation()

template <class TFenotype>
CError CConcatenationEvaluation<TFenotype>::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CEvaluation<TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_components(EVALUATION_ARGUMENT_NUMBER_OF_COMPONENTS, 1, UINT16_MAX);
		uint16_t i_number_of_components = p_number_of_components.iGetValue(psSettings, &c_error);

		CEvaluation<TFenotype> *pc_component = nullptr;

		for (uint16_t i = 0; i < i_number_of_components && !c_error; i++)
		{
			CUIntCommandParam p_component_multiplier(EVALUATION_ARGUMENT_COMPONENT_MULTIPLIER, 1, UINT16_MAX);
			uint16_t i_component_multiplier = p_component_multiplier.iGetValue(psSettings, &c_error);

			if (!c_error)
			{
				pc_component = EvaluationUtils::pcGetEvaluation<TFenotype>(psSettings, &c_error);
			}//if (!c_error)

			if (!c_error)
			{
				v_add_component(pc_component, i_component_multiplier);
			}//if (!c_error)
		}//for (uint16_t i = 0; i < i_number_of_components && !c_error; i++)
	}//if (!c_error)

	return c_error;
}//CError CConcatenationEvaluation<TFenotype>::eConfigure(istream *psSettings)

template <class TFenotype>
TFenotype * CConcatenationEvaluation<TFenotype>::pcCreateSampleFenotype()
{
	if (!pc_sample_fenotype)
	{
		if (!v_components.empty())
		{
			pc_sample_fenotype = v_components.front()->pcCreateSampleFenotype();

			TFenotype *pc_other_sample_fenotype;

			for (uint16_t i = 1; i < (uint16_t)v_components.size(); i++)
			{
				pc_other_sample_fenotype = v_components.at(i)->pcCreateSampleFenotype();
				pc_sample_fenotype->vMerge(pc_other_sample_fenotype);
				delete pc_other_sample_fenotype;
			}//for (uint16_t i = 1; i < (uint16_t)v_components.size(); i++)
		}//if (!v_components.empty())
	}//if (!pc_sample_fenotype)

	return pc_sample_fenotype ? new TFenotype(pc_sample_fenotype) : nullptr;
}//TFenotype * CConcatenationEvaluation<TFenotype>::pcCreateSampleFenotype()

template <class TFenotype>
double CConcatenationEvaluation<TFenotype>::d_evaluate(TFenotype *pcFenotype, uint16_t iShift)
{
	double d_fitness_value = 0;

	CEvaluation<TFenotype> *pc_component;

	for (uint16_t i = 0; i < (uint16_t)v_components.size(); i++)
	{
		pc_component = v_components.at(i);

		d_fitness_value += pc_component->d_evaluate(pcFenotype, iShift);
		iShift += pc_component->iGetNumberOfElements();
	}//for (uint16_t i = 0; i < (uint16_t)v_components.size(); i++)

	return d_fitness_value;
}//double CConcatenationEvaluation<TFenotype>::d_evaluate(TFenotype *pcFenotype, uint16_t iShift)

template <class TFenotype>
void CConcatenationEvaluation<TFenotype>::v_clear_sample_fenotype()
{
	delete pc_sample_fenotype;
	pc_sample_fenotype = nullptr;
}//void CConcatenationEvaluation<TFenotype>::v_clear_sample_fenotype()

template <class TFenotype>
void CConcatenationEvaluation<TFenotype>::v_clear()
{
	if (!v_components.empty())
	{
		CEvaluation<TFenotype> *pc_previous_component = v_components.at(0);

		for (uint16_t i = 1; i < (uint16_t)v_components.size(); i++)
		{
			if (pc_previous_component != v_components.at(i))
			{
				delete pc_previous_component;
				pc_previous_component = v_components.at(i);
			}//if (pc_previous_component != v_components.at(i))
		}//for (uint16_t i = 1; i < (uint16_t)v_components.size(); i++)

		delete pc_previous_component;
	}//if (!v_components.empty())

	v_components.clear();

	v_clear_sample_fenotype();
}//void CConcatenationEvaluation<TFenotype>::v_clear()

template <class TFenotype>
void CConcatenationEvaluation<TFenotype>::v_add_component(CEvaluation<TFenotype> *pcComponent, uint16_t iMultiplier)
{
	i_number_of_elements += iMultiplier * pcComponent->iGetNumberOfElements();
	d_max_value += iMultiplier * pcComponent->dGetMaxValue();

	for (uint16_t i = 0; i < iMultiplier; i++)
	{
		v_components.push_back(pcComponent);
	}//for (uint16_t i = 0; i < iMultiplier; i++)

	v_clear_sample_fenotype();
}//void CConcatenationEvaluation<TFenotype>::v_add_component(CEvaluation<TFenotype> *pcComponent, uint16_t iMultiplier)


template <class TFenotype>
CCyclicFenotypeEvaluation<TFenotype>::~CCyclicFenotypeEvaluation()
{
	v_clear();
}//CCyclicFenotypeEvaluation<TFenotype>::~CCyclicFenotypeEvaluation()

template <class TFenotype>
CError CCyclicFenotypeEvaluation<TFenotype>::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CEvaluation<TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		pc_actual_evaluation = EvaluationUtils::pcGetEvaluation<TFenotype>(psSettings, &c_error);

		if (!c_error)
		{
			i_number_of_elements = pc_actual_evaluation->i_number_of_elements;
			d_max_value = pc_actual_evaluation->d_max_value;

			pc_evaluation_fenotype = pc_actual_evaluation->pcCreateSampleFenotype();
		}//if (!c_error)
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_shift(EVALUATION_ARGUMENT_SHIFT, 1, i_number_of_elements - 1);
		i_cyclic_shift = p_shift.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CCyclicFenotypeEvaluation<TFenotype>::eConfigure(istream *psSettings)

template <class TFenotype>
double CCyclicFenotypeEvaluation<TFenotype>::d_evaluate(TFenotype *pcFenotype, uint16_t iShift)
{
	v_prepare_evaluation_fenotype(pcFenotype, iShift);
	return pc_actual_evaluation->d_evaluate(pc_evaluation_fenotype, iShift);
}//double CCyclicFenotypeEvaluation<TFenotype>::d_evaluate(TFenotype *pcFenotype, uint16_t iShift)

template <class TFenotype>
void CCyclicFenotypeEvaluation<TFenotype>::v_clear()
{
	delete pc_actual_evaluation;
	pc_actual_evaluation = nullptr;

	delete pc_evaluation_fenotype;
	pc_evaluation_fenotype = nullptr;
}//void CCyclicFenotypeEvaluation<TFenotype>::v_clear()


template class CEvaluation<CBinaryCoding>;
template class CInversedEvaluation<CBinaryCoding>;
template class CConcatenationEvaluation<CBinaryCoding>;
template class CCyclicFenotypeEvaluation<CBinaryCoding>;