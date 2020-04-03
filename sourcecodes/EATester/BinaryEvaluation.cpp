#include "BinaryEvaluation.h"

#include "FloatCommandParam.h"
#include "StringCommandParam.h"
#include "UIntCommandParam.h"

#include <atlstr.h>
#include <cfloat>
#include <functional>


void CBinaryCyclicFenotypeEvaluation::v_prepare_evaluation_fenotype(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	uint16_t i_index_without_shift;

	for (uint16_t i = 0; i < i_number_of_elements; i++)
	{
		i_index_without_shift = (i + i_cyclic_shift) % i_number_of_elements;
		*(pc_evaluation_fenotype->piGetBits() + i + iShift) = *(pcFenotype->piGetBits() + i_index_without_shift + iShift);
	}//for (uint16_t i = 0; i < i_number_of_elements; i++)
}//void CBinaryCyclicFenotypeEvaluation::v_prepare_evaluation_fenotype(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryEvaluation::CBinaryEvaluation()
	: CEvaluation<CBinaryCoding>()
{

}//CBinaryEvaluation::CBinaryEvaluation()

CBinaryEvaluation::CBinaryEvaluation(uint16_t iNumberOfElements, double dMaxValue)
	: CEvaluation<CBinaryCoding>(iNumberOfElements, dMaxValue)
{

}//CBinaryEvaluation::CBinaryEvaluation(uint16_t iNumberOfElements, double dMaxValue)

CBinaryCoding * CBinaryEvaluation::pcCreateSampleFenotype()
{
	return new CBinaryCoding(i_number_of_elements);
}//CBinaryCoding * CBinaryEvaluation::pcCreateSampleFenotype()


CError CBinaryFileConfigEvaluation::eConfigure(istream *psSettings)
{
	CError c_error = CEvaluation<CBinaryCoding>::eConfigure(psSettings);

	if (!c_error)
	{
		CStringCommandParam p_config_file_path(EVALUATION_ARGUMENT_CONFIG_FILE_PATH);
		CString s_config_file_path = p_config_file_path.sGetValue(psSettings, &c_error);

		FILE *pf_config = fopen(s_config_file_path, "r");

		if (pf_config)
		{
			c_error = e_init(pf_config);
		}//if (pf_config)
		else
		{
			c_error.vSetError(CError::iERROR_CODE_SYSTEM_FILE_NOT_FOUND, s_config_file_path);
		}//else if (pf_config)
	}//if (!c_error)

	return c_error;
}//CError CBinaryFileConfigEvaluation::eConfigure(istream *psSettings)


CBinaryDeceptiveEvaluation::CBinaryDeceptiveEvaluation()
{
	pc_deceptive_function = nullptr;
}//CBinaryDeceptiveEvaluation::CBinaryDeceptiveEvaluation()

CBinaryDeceptiveEvaluation::CBinaryDeceptiveEvaluation(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryDeceptiveEvaluation::CBinaryDeceptiveEvaluation(FILE *pfConfig, CError *pcError)

CBinaryDeceptiveEvaluation::~CBinaryDeceptiveEvaluation()
{
	delete pc_deceptive_function;
}//CBinaryDeceptiveEvaluation::~CBinaryDeceptiveEvaluation()

CError CBinaryDeceptiveEvaluation::e_init(FILE *pfConfig)
{
	delete pc_deceptive_function;
	pc_deceptive_function = nullptr;
	
	CError c_error = eLoadSingleCompProblem(pfConfig, &pc_deceptive_function);

	if (pc_deceptive_function)
	{
		i_number_of_elements = (uint16_t)pc_deceptive_function->iGetBitLength();
		d_max_value = pc_deceptive_function->dGetMaxFuncVal();
	}//if (pc_deceptive_function)

	fclose(pfConfig);

	return c_error;
}//CError CBinaryDeceptiveEvaluation::e_init(FILE *pfConfig)

double CBinaryDeceptiveEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	return pc_deceptive_function->dGetFuncValue(iShift, pcFenotype->piGetBits(), pcFenotype->iGetNumberOfBits());
}//double CBinaryDeceptiveEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryDeceptiveConcatenationEvaluation::CBinaryDeceptiveConcatenationEvaluation()
{
	pc_deceptive_concatenation_function = nullptr;
}//CBinaryDeceptiveConcatenationEvaluation::CBinaryDeceptiveConcatenationEvaluation()

CBinaryDeceptiveConcatenationEvaluation::CBinaryDeceptiveConcatenationEvaluation(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryDeceptiveConcatenationEvaluation::CBinaryDeceptiveConcatenationEvaluation(FILE *pfConfig, CError *pcError)

CBinaryDeceptiveConcatenationEvaluation::~CBinaryDeceptiveConcatenationEvaluation()
{
	delete pc_deceptive_concatenation_function;
}//CBinaryDeceptiveConcatenationEvaluation::~CBinaryDeceptiveConcatenationEvaluation()

CError CBinaryDeceptiveConcatenationEvaluation::e_init(FILE *pfConfig)
{
	delete pc_deceptive_concatenation_function;
	pc_deceptive_concatenation_function = new CConcatCompProblem();
	
	CError c_error = pc_deceptive_concatenation_function->eLoad(pfConfig);

	if (!c_error)
	{
		i_number_of_elements = (uint16_t)pc_deceptive_concatenation_function->iGetProblemBitLength();
		d_max_value = pc_deceptive_concatenation_function->dGetMaxFuncVal();
	}//if (!c_error)

	fclose(pfConfig);

	return c_error;
}//CError CBinaryDeceptiveConcatenationEvaluation::e_init(FILE *pfConfig)

double CBinaryDeceptiveConcatenationEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	double d_fitness_value;
	pc_deceptive_concatenation_function->eGetFuncValue(pcFenotype->piGetBits(), pcFenotype->iGetNumberOfBits(), &d_fitness_value);

	return d_fitness_value;
}//void CBinaryDeceptiveConcatenationEvaluation::vEvaluate(CBinaryCoding *pcFenotype, double *pdFitnessValue)

CError CBinaryDeceptiveConcatenationEvaluation::eReport(FILE *pfReport)
{
	return pc_deceptive_concatenation_function->eCreateReport(pfReport);
}//CError CBinaryDeceptiveConcatenationEvaluation::eReport(FILE *pfReport)


CBinaryKnapsackEvaluation::CBinaryKnapsackEvaluation()
{
	
}//CBinaryKnapsackEvaluation::CBinaryKnapsackEvaluation()

CBinaryKnapsackEvaluation::CBinaryKnapsackEvaluation(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryKnapsackEvaluation::CBinaryKnapsackEvaluation(uint16_t iNumberOfElements, double dCapacity, double *pdValues, double *pdWeights)

CError CBinaryKnapsackEvaluation::e_init(FILE *pfConfig)
{
	CError c_error = c_knapsack.eLoadSettings(pfConfig);

	if (!c_error)
	{
		i_number_of_elements = (uint16_t)c_knapsack.iGetBitLength();
		d_max_value = c_knapsack.dGetMaxFuncVal();
	}//if (!c_error)

	fclose(pfConfig);

	return c_error;
}//CError CBinaryKnapsackEvaluation::e_init(FILE *pfConfig)

double CBinaryKnapsackEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	//c_knapsack.bRepairGreedy(iShift, pcFenotype->piGetBits(), pcFenotype->iGetNumberOfBits());
	return c_knapsack.dGetFuncValue(iShift, pcFenotype->piGetBits(), pcFenotype->iGetNumberOfBits());
}//double CBinaryKnapsackEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryMaxSatEvaluation::CBinaryMaxSatEvaluation()
{

}//CBinaryMaxSatEvaluation::CBinaryMaxSatEvaluation()

CBinaryMaxSatEvaluation::CBinaryMaxSatEvaluation(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryMaxSatEvaluation::CBinaryMaxSatEvaluation(FILE *pfConfig, CError *pcError)

CError CBinaryMaxSatEvaluation::e_init(FILE *pfConfig)
{
	CError c_error = c_max_sat.eLoadSettings(pfConfig);

	if (!c_error)
	{
		i_number_of_elements = (uint16_t)c_max_sat.iGetBitLength();
		d_max_value = c_max_sat.dGetMaxFuncVal();
	}//if (!c_error)

	fclose(pfConfig);

	return c_error;
}//CError CBinaryMaxSatEvaluation::e_init(FILE *pfConfig)

double CBinaryMaxSatEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	return c_max_sat.dGetFuncValue(iShift, pcFenotype->piGetBits(), pcFenotype->iGetNumberOfBits());
}//double CBinaryMaxSatEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryDeceptiveStepTrapEvaluation::CBinaryDeceptiveStepTrapEvaluation()
{
	d_max_value = 1.0;
}//CBinaryDeceptiveStepTrapEvaluation::CBinaryDeceptiveStepTrapEvaluation()

CError CBinaryDeceptiveStepTrapEvaluation::eConfigure(istream *psSettings)
{
	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_trap_size(EVALUATION_BINARY_DECEPTIVE_STEP_TRAP_ARGUMENT_TRAP_SIZE, 1, UINT16_MAX);
		i_trap_size = p_trap_size.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_step_size(EVALUATION_BINARY_DECEPTIVE_STEP_TRAP_ARGUMENT_STEP_SIZE, 1, i_trap_size);
		i_step_size = p_step_size.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		i_offset = (i_trap_size - i_step_size) % i_step_size;
	}//if (!c_error)

	return c_error;
}//CError CBinaryDeceptiveStepTrapEvaluation::eConfigure(istream *psSettings)

double CBinaryDeceptiveStepTrapEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	uint16_t i_total = 0;

	uint16_t i_trap_maximum = (i_offset + i_trap_size) / i_step_size;

	uint16_t i_partial;

	for (uint16_t i = 0; i < i_number_of_elements; i += i_trap_size)
	{
		i_partial = 0;

		for (uint16_t j = i; j < i + i_trap_size; j++)
		{
			i_partial += (uint16_t)(*(pcFenotype->piGetBits() + j + iShift));
		}//for (uint16_t j = 0; j < i + i_trap_size; j++)

		if (i_partial < i_trap_size)
		{
			i_partial = i_trap_size - i_partial - 1;
		}//if (i_partial < i_trap_size)

		i_total += (i_offset + i_partial) / i_step_size;
	}//for (uint16_t i = 0; i < i_number_of_elements; i += i_trap_size)

	return (double)(i_total * i_trap_size) / (double)(i_number_of_elements * i_trap_maximum);
}//double CBinaryDeceptiveStepTrapEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryNearestNeighborNKEvaluation::CBinaryNearestNeighborNKEvaluation()
{
	d_max_value = 1.0;
	pc_native_nearest_neighbor_nk = nullptr;
}//CBinaryNearestNeighborNKEvaluation::CBinaryNearestNeighborNKEvaluation()

CBinaryNearestNeighborNKEvaluation::~CBinaryNearestNeighborNKEvaluation()
{
	v_clear();
}//CBinaryNearestNeighborNKEvaluation::~CBinaryNearestNeighborNKEvaluation()

CError CBinaryNearestNeighborNKEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = (uint16_t)p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_problem_seed(EVALUATION_BINARY_P3_ARGUMENT_PROBLEM_SEED);
		i_problem_seed = p_problem_seed.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_precision(EVALUATION_BINARY_P3_ARGUMENT_PRECISION);
		i_precision = p_precision.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		CUIntCommandParam p_k(EVALUATION_BINARY_NEAREST_NEIGHBOR_NK_ARGUMENT_K, 1, UINT16_MAX);
		i_k = (uint16_t)p_k.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		Configuration c_config;

		c_config.set("length", (int)i_number_of_elements);
		
		c_config.set("problem_seed", (int)i_problem_seed);
		
		c_config.set("precision", (int)i_precision);
		c_config.set("k", (int)i_k);

		c_config.set("problem_folder", string(""));

		pc_native_nearest_neighbor_nk = new NearestNeighborNK(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryNearestNeighborNKEvaluation::eConfigure(istream *psSettings)

double CBinaryNearestNeighborNKEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_nearest_neighbor_nk->evaluate(v_native_solution);
}//double CBinaryNearestNeighborNKEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)

void CBinaryNearestNeighborNKEvaluation::v_clear()
{
	delete pc_native_nearest_neighbor_nk;
	pc_native_nearest_neighbor_nk = nullptr;
}//void CBinaryNearestNeighborNKEvaluation::v_clear()


CBinaryIsingSpinGlassEvaluation::CBinaryIsingSpinGlassEvaluation()
{
	d_max_value = 1.0;
	pc_native_ising_spin_glass = nullptr;
}//CBinaryIsingSpinGlassEvaluation::CBinaryIsingSpinGlassEvaluation()

CBinaryIsingSpinGlassEvaluation::~CBinaryIsingSpinGlassEvaluation()
{
	v_clear();
}//CBinaryIsingSpinGlassEvaluation::~CBinaryIsingSpinGlassEvaluation()

CError CBinaryIsingSpinGlassEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = (uint16_t)p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	uint32_t i_problem_seed = 0;

	if (!c_error)
	{
		CUIntCommandParam p_problem_seed(EVALUATION_BINARY_P3_ARGUMENT_PROBLEM_SEED);
		i_problem_seed = p_problem_seed.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		Configuration c_config;

		c_config.set("length", (int)i_number_of_elements);
		c_config.set("problem_seed", (int)i_problem_seed);
		c_config.set("precision", 65536);
		c_config.set("ising_type", string("pm"));
		c_config.set("problem_folder", string(""));

		pc_native_ising_spin_glass = new IsingSpinGlass(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryIsingSpinGlassEvaluation::eConfigure(istream *psSettings)

double CBinaryIsingSpinGlassEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_ising_spin_glass->evaluate(v_native_solution);
}//double CBinaryIsingSpinGlassEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)

void CBinaryIsingSpinGlassEvaluation::v_clear()
{
	delete pc_native_ising_spin_glass;
	pc_native_ising_spin_glass = nullptr;
}//void CBinaryIsingSpinGlassEvaluation::v_clear()


CBinaryMax3SatEvaluation::CBinaryMax3SatEvaluation()
{
	d_max_value = 1.0;
	pc_native_max_sat = nullptr;
}//CBinaryMax3SatEvaluation::CBinaryMax3SatEvaluation()

CBinaryMax3SatEvaluation::~CBinaryMax3SatEvaluation()
{
	v_clear();
}//CBinaryMax3SatEvaluation::~CBinaryMax3SatEvaluation()

CError CBinaryMax3SatEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = (uint16_t)p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	uint32_t i_problem_seed = 0;

	if (!c_error)
	{
		CUIntCommandParam p_problem_seed(EVALUATION_BINARY_P3_ARGUMENT_PROBLEM_SEED);
		i_problem_seed = p_problem_seed.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	float f_clause_ratio = 0;

	if (!c_error)
	{
		CFloatCommandParam p_clause_ratio(EVALUATION_BINARY_MAX_3_SAT_ARGUMENT_CLAUSE_RATIO, 0, FLT_MAX);
		f_clause_ratio = p_clause_ratio.fGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		Configuration c_config;

		c_config.set("length", (int)i_number_of_elements);
		c_config.set("problem_seed", (int)i_problem_seed);
		c_config.set("precision", 65536);
		c_config.set("clause_ratio", f_clause_ratio);

		pc_native_max_sat = new MAXSAT(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryMax3SatEvaluation::eConfigure(istream *psSettings)

void CBinaryMax3SatEvaluation::v_clear()
{
	delete pc_native_max_sat;
	pc_native_max_sat = nullptr;
}//void CBinaryMax3SatEvaluation::v_clear()

double CBinaryMax3SatEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_max_sat->evaluate(v_native_solution);
}//double CBinaryMax3SatEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)


CBinaryDiscretizedRastriginEvaluation::CBinaryDiscretizedRastriginEvaluation()
{
	d_max_value = 1.0;
	pc_native_rastrigin = nullptr;
}//CBinaryDiscretizedRastriginEvaluation::CBinaryDiscretizedRastriginEvaluation()

CBinaryDiscretizedRastriginEvaluation::~CBinaryDiscretizedRastriginEvaluation()
{
	v_clear();
}//CBinaryDiscretizedRastriginEvaluation::~CBinaryDiscretizedRastriginEvaluation()

CError CBinaryDiscretizedRastriginEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		Configuration c_config;

		c_config.set("length", (int)i_number_of_elements);
		c_config.set("precision", 65536);
		c_config.set("bits_per_float", 10);

		pc_native_rastrigin = new Rastrigin(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryDiscretizedRastriginEvaluation::eConfigure(istream *psSettings)

double CBinaryDiscretizedRastriginEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_rastrigin->evaluate(v_native_solution);
}//double CBinaryDiscretizedRastriginEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)

void CBinaryDiscretizedRastriginEvaluation::v_clear()
{
	delete pc_native_rastrigin;
	pc_native_rastrigin = nullptr;
}//void CBinaryDiscretizedRastriginEvaluation::v_clear()


CBinaryDiscretizedRosenbrockEvaluation::CBinaryDiscretizedRosenbrockEvaluation()
{
	d_max_value = 1.0;
	pc_native_rosenbrock = nullptr;
}//CBinaryDiscretizedRosenbrockEvaluation::CBinaryDiscretizedRosenbrockEvaluation()

CBinaryDiscretizedRosenbrockEvaluation::~CBinaryDiscretizedRosenbrockEvaluation()
{
	v_clear();
}//CBinaryDiscretizedRosenbrockEvaluation::~CBinaryDiscretizedRosenbrockEvaluation()

CError CBinaryDiscretizedRosenbrockEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = (uint16_t)p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	uint8_t i_bits_per_float = 0;

	if (!c_error)
	{
		CUIntCommandParam p_bits_per_float(EVALUATION_BINARY_P3_ARGUMENT_BITS_PER_FLOAT, 1, UINT8_MAX);
		i_bits_per_float = (uint8_t)p_bits_per_float.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		Configuration c_config;

		c_config.set("length", (int)i_number_of_elements);
		c_config.set("precision", 65536);
		c_config.set("bits_per_float", (int)i_bits_per_float);

		pc_native_rosenbrock = new Rosenbrock(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryDiscretizedRosenbrockEvaluation::eConfigure(istream *psSettings)

double CBinaryDiscretizedRosenbrockEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_rosenbrock->evaluate(v_native_solution);
}//double CBinaryDiscretizedRosenbrockEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)

void CBinaryDiscretizedRosenbrockEvaluation::v_clear()
{
	delete pc_native_rosenbrock;
	pc_native_rosenbrock = nullptr;
}//void CBinaryDiscretizedRosenbrockEvaluation::v_clear()


CBinaryHIFFEvaluation::CBinaryHIFFEvaluation()
{
	d_max_value = 1.0;
	pc_native_hiff = nullptr;
}//CBinaryHIFFEvaluation::CBinaryHIFFEvaluation()

CBinaryHIFFEvaluation::~CBinaryHIFFEvaluation()
{
	v_clear();
}//CBinaryHIFFEvaluation::~CBinaryHIFFEvaluation()

CError CBinaryHIFFEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		Configuration c_config;

		i_number_of_elements = 2048;

		c_config.set("length", (int)i_number_of_elements);
		c_config.set("precision", 65536);

		pc_native_hiff = new HIFF(c_config, 0);

		v_native_solution.resize(i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryHIFFEvaluation::eConfigure(istream *psSettings)

double CBinaryHIFFEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)
	{
		v_native_solution.at(i) = *(pcFenotype->piGetBits() + i + iShift) == 1;
	}//for (uint16_t i = 0; i < pcFenotype->iGetNumberOfBits(); i++)

	return (double)pc_native_hiff->evaluate(v_native_solution);
}//double CBinaryHIFFEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)

void CBinaryHIFFEvaluation::v_clear()
{
	delete pc_native_hiff;
	pc_native_hiff = nullptr;
}//void CBinaryHIFFEvaluation::v_clear()


CBinaryDSMGA2BasedEvaluation::CBinaryDSMGA2BasedEvaluation()
{
	pc_native_chromosome = nullptr;
}//CBinaryDSMGA2BasedEvaluation::CBinaryDSMGA2BasedEvaluation()

CBinaryDSMGA2BasedEvaluation::~CBinaryDSMGA2BasedEvaluation()
{
	v_clear();
}//CBinaryDSMGA2BasedEvaluation::~CBinaryDSMGA2BasedEvaluation()

CError CBinaryDSMGA2BasedEvaluation::eConfigure(istream *psSettings)
{
	v_clear();

	CError c_error = CBinaryEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_number_of_elements(EVALUATION_ARGUMENT_NUMBER_OF_ELEMENTS, 1, UINT16_MAX);
		i_number_of_elements = p_number_of_elements.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	if (!c_error)
	{
		pc_native_chromosome = new Chromosome((int)i_number_of_elements);
	}//if (!c_error)

	return c_error;
}//CError CBinaryDSMGA2BasedEvaluation::eConfigure(istream *psSettings)

void CBinaryDSMGA2BasedEvaluation::v_clear()
{
	delete pc_native_chromosome;
	pc_native_chromosome = nullptr;
}//void CBinaryDSMGA2BasedEvaluation::v_clear()


CError CBinaryCyclicTrapEvaluation::eConfigure(istream *psSettings)
{
	CError c_error = CBinaryDSMGA2BasedEvaluation::eConfigure(psSettings);

	if (!c_error)
	{
		Chromosome::function = Chromosome::Function::CYCTRAP;
		d_max_value = pc_native_chromosome->getMaxFitness() + 1e-8;
	}//if (!c_error)

	return c_error;
}//CError CBinaryCyclicTrapEvaluation::eConfigure(istream *psSettings)

double CBinaryCyclicTrapEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	for (uint16_t i = 0; i < i_number_of_elements; i++)
	{
		pc_native_chromosome->setVal((int)i, (int)*(pcFenotype->piGetBits() + i + iShift));
	}//for (uint16_t i = 0; i < i_number_of_elements; i++)

	return pc_native_chromosome->cycTrap(1, 0.8);
}//double CBinaryCyclicTrapEvaluation::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)