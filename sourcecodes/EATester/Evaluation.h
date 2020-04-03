#ifndef EVALUATION_H
#define EVALUATION_H

#define EVALUATION_ARGUMENT_TYPE "evaluation_type"
#define EVALUATION_ARGUMENT_TYPE_CONCATENATION "concatenation"
#define EVALUATION_ARGUMENT_TYPE_BINARY_CYCLIC_FENOTYPE "binary_cyclic_fenotype"
#define EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE "binary_deceptive"
#define EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE_CONCATENATION "binary_deceptive_concatenation"
#define EVALUATION_ARGUMENT_TYPE_BINARY_KNAPSACK "binary_knapsack"
#define EVALUATION_ARGUMENT_TYPE_BINARY_MAX_SAT "binary_max_sat"
#define EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE_STEP_TRAP "binary_deceptive_step_trap"
#define EVALUATION_ARGUMENT_TYPE_BINARY_NEAREST_NEIGHBOR_NK "binary_nearest_neighbor_nk"
#define EVALUATION_ARGUMENT_TYPE_BINARY_ISING_SPIN_GLASS "binary_ising_spin_glass"
#define EVALUATION_ARGUMENT_TYPE_BINARY_MAX_3_SAT "binary_max_3_sat"
#define EVALUATION_ARGUMENT_TYPE_BINARY_DISCRETIZED_RASTRIGIN "binary_discretized_rastrigin"
#define EVALUATION_ARGUMENT_TYPE_BINARY_DISCRETIZED_ROSENBROCK "binary_discretized_rosenbrock"
#define EVALUATION_ARGUMENT_TYPE_BINARY_HIFF "binary_hiff"
#define EVALUATION_ARGUMENT_TYPE_BINARY_CYCLIC_TRAP "binary_cyclic_trap"

#define EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS "evaluation_number_of_elements"

#define EVALUATION_ARGUMENT_SHIFT "evaluation_shift"

#define EVALUATION_ARGUMENT_NUMBER_OF_COMPONENTS "evaluation_number_of_components"
#define EVALUATION_ARGUMENT_COMPONENT_MULTIPLIER "evaluation_component_multiplier"

#define EVALUATION_ARGUMENT_CONFIG_FILE_PATH "evaluation_config_file_path"

#include "Error.h"

#include <cstdint>
#include <istream>
#include <vector>

using namespace std;

template <class TFenotype> class CInversedEvaluation;
template <class TFenotype> class CConcatenationEvaluation;
template <class TFenotype> class CCyclicFenotypeEvaluation;


enum EEvaluationType
{
	EVALUATION_CONCATENATION = 0,
	EVALUATION_BINARY_CYCLIC_FENOTYPE = 1,
	EVALUATION_BINARY_DECEPTIVE = 2,
	EVALUATION_BINARY_DECEPTIVE_CONCATENATION = 3,
	EVALUATION_BINARY_KNAPSACK = 4,
	EVALUATION_BINARY_MAX_SAT = 5,
	EVALUATION_BINARY_DECEPTIVE_STEP_TRAP = 6,
	EVALUATION_BINARY_NEAREST_NEIGHBOR_NK = 7,
	EVALUATION_BINARY_ISING_SPIN_GLASS = 8,
	EVALUATION_BINARY_MAX_3_SAT = 9,
	EVALUATION_BINARY_DISCRETIZED_RASTRIGIN = 10,
	EVALUATION_BINARY_DISCRETIZED_ROSENBROCK = 11,
	EVALUATION_BINARY_HIFF = 12,
	EVALUATION_BINARY_CYCLIC_TRAP = 13
};//enum EEvaluationType


template <class TFenotype>
class CEvaluation
{
friend class CInversedEvaluation<TFenotype>;
friend class CConcatenationEvaluation<TFenotype>;
friend class CCyclicFenotypeEvaluation<TFenotype>;

public:
	CEvaluation();
	CEvaluation(uint16_t iNumberOfElements, double dMaxValue);

	virtual ~CEvaluation();

	virtual CError eConfigure(istream *psSettings) { return CError(iERROR_PARENT_CEVALUATION); };

	double dEvaluate(TFenotype *pcFenotype);

	virtual TFenotype *pcCreateSampleFenotype() = 0;

	virtual bool bIsMaxValue(double dValue) { return dValue == d_max_value; };

	uint16_t iGetNumberOfElements() { return i_number_of_elements; };

	double dGetMaxValue() { return d_max_value; };

	uint64_t iGetFFE() { return i_ffe; };

	void vSetFFE(uint64_t iFFE) { i_ffe = iFFE; };
	void vSetZeroFFE() { vSetFFE(0); };

	void vIncreaseFFE() { i_ffe++; };
	void vDecreaseFFE() { i_ffe--; };

protected:
	static uint32_t iERROR_PARENT_CEVALUATION;

	virtual double d_evaluate(TFenotype *pcFenotype, uint16_t iShift) = 0;

	void v_init(uint16_t iNumberOfElements, double dMaxValue);

	uint16_t i_number_of_elements;

	double d_max_value;

	uint64_t i_ffe;
};//class CEvaluation


template <class TFenotype>
class CInversedEvaluation : public CEvaluation<TFenotype>
{
public:
	CInversedEvaluation(CEvaluation<TFenotype> *pcEvaluation);

	virtual ~CInversedEvaluation();

protected:
	virtual double d_evaluate(TFenotype *pcFenotype, uint16_t iShift);

private:
	CEvaluation<TFenotype> *pc_evaluation;
};//class CInversedEvaluation : public CEvaluation<TFenotype>


template <class TFenotype>
class CConcatenationEvaluation : public CEvaluation<TFenotype>
{
public:
	CConcatenationEvaluation();

	virtual ~CConcatenationEvaluation();

	virtual CError eConfigure(istream *psSettings);

	virtual TFenotype *pcCreateSampleFenotype();

protected:
	virtual double d_evaluate(TFenotype *pcFenotype, uint16_t iShift);

private:
	void v_clear_sample_fenotype();
	void v_clear();

	void v_add_component(CEvaluation<TFenotype> *pcComponent, uint16_t iMultiplier);

	vector<CEvaluation<TFenotype>*> v_components;

	TFenotype *pc_sample_fenotype;
};//class CConcatenationEvaluation : public CEvaluation<TFenotype>


template <class TFenotype>
class CCyclicFenotypeEvaluation : public CEvaluation<TFenotype>
{
public:
	virtual ~CCyclicFenotypeEvaluation();

	virtual CError eConfigure(istream *psSettings);

	virtual TFenotype *pcCreateSampleFenotype() { return pc_actual_evaluation->pcCreateSampleFenotype(); };

protected:
	virtual double d_evaluate(TFenotype *pcFenotype, uint16_t iShift);

	virtual void v_prepare_evaluation_fenotype(TFenotype *pcFenotype, uint16_t iShift) = 0;

	uint16_t i_cyclic_shift;

	TFenotype *pc_evaluation_fenotype;

private:
	void v_clear();

	CEvaluation<TFenotype> *pc_actual_evaluation;
};//class CCyclicFenotypeEvaluation : public CEvaluation<TFenotype>

#endif//EVALUATION_H