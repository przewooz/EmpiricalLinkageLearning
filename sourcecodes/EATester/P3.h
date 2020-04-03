#ifndef P3_H
#define P3_H

#include "BinaryCoding.h"
#include "Error.h"
#include "Log.h"
#include "P3BasedOptimizer.h"
#include "Optimizer.h"
#include "Problem.h"

#include "../P3/Pyramid.h"

#include <atlstr.h>
#include <cstdint>
#include <ctime>
#include <istream>

using namespace std;

class CIncrementalP3;
class CIslandModelP3;

namespace DynamicOptimizer
{
	class CDynamicP3;
}//namespace DynamicOptimizer

class CP3 : public CP3BasedOptimizer
{
friend class CIncrementalP3;
friend class CIslandModelP3;
friend class DynamicOptimizer::CDynamicP3;

public:
	CP3(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CP3(CP3 *pcOther);

	virtual CError eConfigure(istream *psSettings);

	virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CP3(this); };

	virtual void vInitialize(time_t tStartTime);
	virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

	virtual void vRun();

	bool bClimb(CBinaryCoding *pcGenotype, uint32_t iIterationNumber, time_t tStartTime);

private:
	CString s_create_linkage_report();
	void v_log_linkage(bool bEcho = false);

	bool b_with_feedback;
	float f_feedback_probability;

	uint32_t i_feedback_counter;
	uint32_t i_effective_feedback_counter;

	shared_ptr<Pyramid> pc_p3;
};//class CP3 : public COptimizer<CBinaryCoding, CBinaryCoding>

#endif//P3_H