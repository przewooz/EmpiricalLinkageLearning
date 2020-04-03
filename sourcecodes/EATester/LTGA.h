#ifndef LTGA_H
#define LTGA_H

#include "BinaryCoding.h"
#include "Log.h"
#include "P3BasedOptimizer.h"
#include "Optimizer.h"
#include "Problem.h"

#include "../P3/LTGA.h"

#include <cstdint>
#include <ctime>

using namespace std;

namespace DynamicOptimizer
{
	class CDynamicLTGA;
}//namespace DynamicOptimizer

class CLTGA : public CP3BasedOptimizer
{
friend class DynamicOptimizer::CDynamicLTGA;

public:
	CLTGA(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CLTGA(CLTGA *pcOther);

	virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CLTGA(this); };

	virtual void vInitialize(time_t tStartTime);
	virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

private:
	shared_ptr<LTGA> pc_ltga;
};//class CLTGA : public CP3BasedOptimizer

#endif//LTGA_H