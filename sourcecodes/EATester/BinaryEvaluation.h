#ifndef BINARY_EVALUATION_H
#define BINARY_EVALUATION_H

#define EVALUATION_BINARY_P3_ARGUMENT_PROBLEM_SEED "evaluation_problem_seed"
#define EVALUATION_BINARY_P3_ARGUMENT_PRECISION "evaluation_precision"
#define EVALUATION_BINARY_P3_ARGUMENT_BITS_PER_FLOAT "evaluation_bits_per_float"

#define EVALUATION_BINARY_DECEPTIVE_STEP_TRAP_ARGUMENT_TRAP_SIZE "evaluation_trap_size"
#define EVALUATION_BINARY_DECEPTIVE_STEP_TRAP_ARGUMENT_STEP_SIZE "evaluation_step_size"

#define EVALUATION_BINARY_NEAREST_NEIGHBOR_NK_ARGUMENT_K "evaluation_k"

#define EVALUATION_BINARY_MAX_3_SAT_ARGUMENT_CLAUSE_RATIO "evaluation_clause_ratio"

#include "BinaryCoding.h"
#include "CompProblem.h"
#include "Error.h"
#include "Evaluation.h"
#include "FiberNet.h"

#include "../DSMGA2/chromosome.h"
#include "../P3/Evaluation.h"

#include <cstdint>
#include <iostream>
#include <vector>

using namespace FiberNets;
using namespace ProblemTools;
using namespace std;


class CBinaryCyclicFenotypeEvaluation : public CCyclicFenotypeEvaluation<CBinaryCoding>
{
protected:
	virtual void v_prepare_evaluation_fenotype(CBinaryCoding *pcFenotype, uint16_t iShift);
};//class CBinaryCyclicFenotypeEvaluation : public CCyclicFenotypeEvaluation<CBinaryCoding>


class CBinaryEvaluation : public CEvaluation<CBinaryCoding>
{
public:
	CBinaryEvaluation();
	CBinaryEvaluation(uint16_t iNumberOfElements, double dMaxValue);

	virtual CBinaryCoding *pcCreateSampleFenotype();
};//class CBinaryEvaluation : public CEvaluation<CBinaryCoding>


class CBinaryFileConfigEvaluation : public CBinaryEvaluation
{
public:
	virtual CError eConfigure(istream *psSettings);

protected:
	virtual CError e_init(FILE *pfConfig) = 0;
};//class CBinaryFileConfigEvaluation : public CBinaryEvaluation


class CBinaryDeceptiveEvaluation : public CBinaryFileConfigEvaluation
{
public:
	CBinaryDeceptiveEvaluation();
	CBinaryDeceptiveEvaluation(FILE *pfConfig, CError *pcError);
	
	virtual ~CBinaryDeceptiveEvaluation();

protected:
	virtual CError e_init(FILE *pfConfig);

	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	CCompProblem *pc_deceptive_function;
};//class CBinaryDeceptiveEvaluation


class CBinaryDeceptiveConcatenationEvaluation : public CBinaryFileConfigEvaluation
{
public:
	CBinaryDeceptiveConcatenationEvaluation();
	CBinaryDeceptiveConcatenationEvaluation(FILE *pfConfig, CError *pcError);

	virtual ~CBinaryDeceptiveConcatenationEvaluation();

	CError eReport(FILE *pfReport);

protected:
	virtual CError e_init(FILE *pfConfig);

	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	CConcatCompProblem *pc_deceptive_concatenation_function;
};//class CBinaryDeceptiveConcatenationEvaluation


class CBinaryKnapsackEvaluation : public CBinaryFileConfigEvaluation
{
public:
	CBinaryKnapsackEvaluation();
	CBinaryKnapsackEvaluation(FILE *pfConfig, CError *pcError);

protected:
	virtual CError e_init(FILE *pfConfig);

	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	CKnapsackProblem c_knapsack;
};//class CBinaryKnapsackEvaluation


class CBinaryMaxSatEvaluation : public CBinaryFileConfigEvaluation
{
public:
	CBinaryMaxSatEvaluation();
	CBinaryMaxSatEvaluation(FILE *pfConfig, CError *pcError);


protected:
	virtual CError e_init(FILE *pfConfig);

	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	CMaxSatProblem c_max_sat;
};//class CBinaryMaxSatEvaluation


class CBinaryDeceptiveStepTrapEvaluation : public CBinaryEvaluation
{
public:
	CBinaryDeceptiveStepTrapEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	uint16_t i_trap_size;
	uint16_t i_step_size;
	uint16_t i_offset;
};//class CBinaryDeceptiveStepTrapEvaluation : public CBinaryEvaluation


class CBinaryNearestNeighborNKEvaluation : public CBinaryEvaluation
{
public:
	CBinaryNearestNeighborNKEvaluation();

	virtual ~CBinaryNearestNeighborNKEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	uint16_t i_k;

	uint32_t i_precision;
	uint32_t i_problem_seed;

	vector<bool> v_native_solution;

	NearestNeighborNK *pc_native_nearest_neighbor_nk;
};//class CBinaryNearestNeighborNKEvaluation : public CBinaryEvaluation


class CBinaryIsingSpinGlassEvaluation : public CBinaryEvaluation
{
public:
	CBinaryIsingSpinGlassEvaluation();

	virtual ~CBinaryIsingSpinGlassEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	vector<bool> v_native_solution;

	IsingSpinGlass *pc_native_ising_spin_glass;
};//class CBinaryIsingSpinGlassEvaluation : public CBinaryEvaluation


class CBinaryMax3SatEvaluation : public CBinaryEvaluation
{
public:
	CBinaryMax3SatEvaluation();

	virtual ~CBinaryMax3SatEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	vector<bool> v_native_solution;

	MAXSAT *pc_native_max_sat;
};//class CBinaryMax3SatEvaluation : public CBinaryEvaluation


class CBinaryDiscretizedRastriginEvaluation : public CBinaryEvaluation
{
public:
	CBinaryDiscretizedRastriginEvaluation();

	virtual ~CBinaryDiscretizedRastriginEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	vector<bool> v_native_solution;

	Rastrigin *pc_native_rastrigin;
};//class CBinaryDiscretizedRastriginEvaluation : public CBinaryEvaluation


class CBinaryDiscretizedRosenbrockEvaluation : public CBinaryEvaluation
{
public:
	CBinaryDiscretizedRosenbrockEvaluation();

	virtual ~CBinaryDiscretizedRosenbrockEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	vector<bool> v_native_solution;

	Rosenbrock *pc_native_rosenbrock;
};//class CBinaryDiscretizedRosenbrockEvaluation : public CBinaryEvaluation


class CBinaryHIFFEvaluation : public CBinaryEvaluation
{
public:
	CBinaryHIFFEvaluation();

	virtual ~CBinaryHIFFEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);

private:
	void v_clear();

	vector<bool> v_native_solution;

	HIFF *pc_native_hiff;
};//class CBinaryHIFFEvaluation : public CBinaryEvaluation


class CBinaryDSMGA2BasedEvaluation : public CBinaryEvaluation
{
public:
	CBinaryDSMGA2BasedEvaluation();

	virtual ~CBinaryDSMGA2BasedEvaluation();

	virtual CError eConfigure(istream *psSettings);

protected:
	Chromosome *pc_native_chromosome;

private:
	void v_clear();
};//class CBinaryDSMGA2BasedEvaluation : public CBinaryEvaluation


class CBinaryCyclicTrapEvaluation : public CBinaryDSMGA2BasedEvaluation
{
public:
	virtual CError eConfigure(istream *psSettings);

protected:
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);
};//class CBinaryCyclicTrapEvaluation : public CBinaryDSMGA2BasedEvaluation

#endif//BINARY_EVALUATION_H