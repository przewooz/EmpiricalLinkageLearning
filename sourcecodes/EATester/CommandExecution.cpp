#include "CommandExecution.h"

#include "BinaryCoding.h"
#include "CompProblem.h"
#include "EnumCommandParam.h"
#include "Evaluation.h"
#include "EvaluationUtils.h"
#include "FloatCommandParam.h"
#include "MathUtils.h"
#include "Problem.h"
#include "RandUtils.h"
#include "StringCommandParam.h"
#include "StringUtils.h"
#include "System.h"
#include "Transformation.h"
#include "UIntCommandParam.h"

#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>


CCommandExecution::CCommandExecution(CString sName, istream *psSettings, CLog *pcLog, ostream *psSummary)
{
	s_name = sName;
	ps_settings = psSettings;
	pc_log = pcLog;
	ps_summary = psSummary;
}//CCommandExecution::CCommandExecution(CString sName, CLog *pcLog, ostream *psSummary)

CCommandExecution::~CCommandExecution()
{

}//CCommandExecution::~CCommandExecution()