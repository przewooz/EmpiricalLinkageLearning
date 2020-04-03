#include "LTGA.h"

#include <atlstr.h>
#include <vector>

CLTGA::CLTGA(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
	: CP3BasedOptimizer("LTGA", pcProblem, pcLog, iRandomSeed)
{

}//CLTGA::CLTGA(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

CLTGA::CLTGA(CLTGA *pcOther)
	: CP3BasedOptimizer(pcOther)
{

}//CLTGA::CLTGA(CLTGA *pcOther)

void CLTGA::vInitialize(time_t tStartTime)
{
	CP3BasedOptimizer::vInitialize(tStartTime);
	pc_ltga = static_pointer_cast<LTGA>(pc_optimizer);
}//void CLTGA::vInitialize(time_t tStartTime)

bool CLTGA::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	bool b_updated = CP3BasedOptimizer::bRunIteration(iIterationNumber, tStartTime);

	//if (iIterationNumber < 5)
	//{
	//	vector<vector<bool>> v_population(pc_ltga->pop.solutions);

	//	CString s_log_message;

	//	pc_log->vPrintEmptyLine();

	//	for (size_t i = 0; i < v_population.size(); i++)
	//	{
	//		s_log_message.Empty();

	//		for (size_t j = 0; j < v_population.at(i).size(); j++)
	//		{
	//			s_log_message.AppendFormat("%d", v_population.at(i).at(j) ? 1 : 0);
	//		}//for (size_t j = 0; j < v_population.at(i).size(); j++)

	//		pc_log->vPrintLine(s_log_message);
	//	}//for (size_t i = 0; i < v_population.size(); i++)

	//	pc_log->vPrintEmptyLine();
	//}//if (iIterationNumber < 5)

	return b_updated;
}//bool CLTGA::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)