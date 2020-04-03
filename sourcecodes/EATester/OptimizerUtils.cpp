#include "OptimizerUtils.h"

#include "BinaryCoding.h"
#include "DSMGA2.h"
#include "3LO.h"
#include "EnumCommandParam.h"
#include "LTGA.h"
#include "LTGAOriginal.h"
#include "LTGAOriginalPopulationSizing.h"
#include "P3.h"
#include "PopulationSizingOptimizer.h"
#include "StringUtils.h"

#include <atlstr.h>
#include <unordered_map>
#include <utility>

using namespace ThreeLO;


template <class TGenotype, class TFenotype>
COptimizer<TGenotype, TFenotype> * OptimizerUtils::pcGetOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed, istream *psSettings, CError *pcError, bool bIsObligatory)
{
	COptimizer<TGenotype, TFenotype> *pc_optimizer = nullptr;

	size_t i_genotype_type_hash_code = typeid(TGenotype).hash_code();
	size_t i_fenotype_type_hash_code = typeid(TFenotype).hash_code();

	unordered_map<CString, EOptimizerType> m_optimizer_types;

	m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_POPULATION_SIZING, OPTIMIZER_POPULATION_SIZING));

	if (i_genotype_type_hash_code == typeid(CBinaryCoding).hash_code() && i_fenotype_type_hash_code == typeid(CBinaryCoding).hash_code())
	{
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_P3, OPTIMIZER_P3));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_LTGA, OPTIMIZER_LTGA));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_LTGA_ORIGINAL, OPTIMIZER_LTGA_ORIGINAL));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_LTGA_ORIGINAL_POPULATION_SIZING, OPTIMIZER_LTGA_ORIGINAL_POPULATION_SIZING));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_DSMGA2, OPTIMIZER_DSMGA2));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_3LO, OPTIMIZER_3LO));
		m_optimizer_types.insert(pair<const CString, EOptimizerType>(OPTIMIZER_ARGUMENT_TYPE_3LO_SINGLE, OPTIMIZER_3LO_SINGLE));
	}//if (i_genotype_type_hash_code == typeid(CBinaryCoding).hash_code() && i_fenotype_type_hash_code == typeid(CBinaryCoding).hash_code())

	CEnumCommandParam<EOptimizerType> p_type(OPTIMIZER_ARGUMENT_TYPE, &m_optimizer_types, bIsObligatory);
	EOptimizerType e_type = p_type.eGetValue(psSettings, pcError);

	if (!*pcError && p_type.bHasValue())
	{
		switch (e_type)
		{
			case OPTIMIZER_P3:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new CP3((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_P3
			case OPTIMIZER_LTGA:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new CLTGA((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_LTGA
			case OPTIMIZER_DSMGA2:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new CDSMGA2((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_DSMGA2
			case OPTIMIZER_LTGA_ORIGINAL:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new CLTGAOriginal((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_LTGA_ORIGINAL
			case OPTIMIZER_LTGA_ORIGINAL_POPULATION_SIZING:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new CLTGAOriginalPopulationSizing((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_LTGA_ORIGINAL_POPULATION_SIZING
			case OPTIMIZER_3LO:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new C3LO((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_3LO
			case OPTIMIZER_3LO_SINGLE:
			{
				pc_optimizer = (COptimizer<TGenotype, TFenotype>*)new C3LOSingle((CProblem<CBinaryCoding, CBinaryCoding>*)pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_3LO
			case OPTIMIZER_POPULATION_SIZING:
			{
				pc_optimizer = new CPopulationSizingOptimizer<TGenotype, TFenotype>(pcProblem, pcLog, iRandomSeed);
				break;
			}//case OPTIMIZER_POPULATION_SIZING
			default:
			{
				pcError->vSetError(CError::iERROR_CODE_OPERATOR_NOT_FOUND, "optimizer");
				break;
			}//default
		}//switch (e_type)
	}//if (!*pcError && p_type.bHasValue())

	if (!*pcError)
	{
		*pcError = pc_optimizer->eConfigure(psSettings);
	}//if (!*pcError)

	return pc_optimizer;
}//COptimizer<TGenotype, TFenotype> * OptimizerUtils::pcGetOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed, istream *psSettings, COptimizerParams<TGenotype, TFenotype> **ppcParams, CError *pcError, bool bIsObligatory)

template<class TGenotype, class TFenotype>
CPopulationSizingSingleOptimizer<TGenotype, TFenotype>* OptimizerUtils::pcGetPopulationSizingSingleOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed, istream *psSettings, CError *pcError, bool bIsObligatory)
{
	COptimizer<TGenotype, TFenotype> *pc_optimizer = pcGetOptimizer(pcProblem, pcLog, iRandomSeed, psSettings, pcError, bIsObligatory);

	if (!*pcError && pc_optimizer)
	{
		if (!dynamic_cast<CPopulationSizingSingleOptimizer<TGenotype, TFenotype>*>(pc_optimizer))
		{
			pcError->vSetError(CError::iERROR_CODE_OPERATOR_NOT_FOUND, "optimizer instead of population sizing single optimizer");
		}//if (!dynamic_cast<CPopulationSizingSingleOptimizer<TGenotype, TFenotype>*>(pc_optimizer))
	}//if (!*pcError && pc_optimizer)

	return (CPopulationSizingSingleOptimizer<TGenotype, TFenotype>*)pc_optimizer;
}//CPopulationSizingSingleOptimizer<TGenotype, TFenotype>* OptimizerUtils::pcGetPopulationSizingSingleOptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed, istream *psSettings, CError *pcError, bool bIsObligatory)


template COptimizer<CBinaryCoding, CBinaryCoding> * OptimizerUtils::pcGetOptimizer(CProblem<CBinaryCoding, CBinaryCoding>*, CLog*, uint32_t, istream*, CError*, bool);
template CPopulationSizingSingleOptimizer<CBinaryCoding, CBinaryCoding> * OptimizerUtils::pcGetPopulationSizingSingleOptimizer(CProblem<CBinaryCoding, CBinaryCoding>*, CLog*, uint32_t, istream*, CError*, bool);