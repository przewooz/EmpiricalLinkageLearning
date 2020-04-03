#ifndef P3_BASED_OPTIMIZER
#define P3_BASED_OPTIMIZER

#include "BinaryCoding.h"
#include "Error.h"
#include "Individual.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"

#include "../P3/Configuration.h"
#include "../P3/MiddleLayer.h"
#include "../P3/OptimizationCollection.h"
#include "../P3/Util.h"

#include <ctime>
#include <cstdint>
#include <cstring>
#include <istream>

using namespace std;

class CP3BasedOptimizer : public COptimizer<CBinaryCoding, CBinaryCoding>
{
public:
	CP3BasedOptimizer(string sOptimizerName, CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CP3BasedOptimizer(CP3BasedOptimizer *pcOther);

	~CP3BasedOptimizer();

	virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CP3BasedOptimizer(this); };

	virtual CError eConfigure(istream *psSettings);

	virtual void vInitialize(time_t tStartTime);
	virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

protected:
	bool b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime);

	shared_ptr<Optimizer> pc_optimizer;
	Middle_Layer *pc_recorder;

private:
	string s_optimizer_name;

	Configuration *pc_config;
	Random *pc_rand;
};//class CP3BasedOptimizer : public COptimizer<CBinaryCoding, CBinaryCoding>

#endif//P3_BASED_OPTIMIZER