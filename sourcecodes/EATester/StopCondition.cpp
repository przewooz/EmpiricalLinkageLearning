#include "StopCondition.h"

#include "BinaryCoding.h"
#include "UIntCommandParam.h"

template <class TGenotype, class TFenotype>
CStopCondition<TGenotype, TFenotype>::CStopCondition(CEvaluation<TFenotype> *pcEvaluation)
{
	pc_evaluation = pcEvaluation;
}//CStopCondition<TGenotype, TFenotype>::CStopCondition(CEvaluation<TFenotype> *pcEvaluation)

template <class TGenotype, class TFenotype>
CStopCondition<TGenotype, TFenotype>::~CStopCondition()
{

}//CStopCondition<TGenotype, TFenotype>::~CStopCondition()

template <class TGenotype, class TFenotype>
bool CStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)
{
	bool b_stop = false;

	if (pcBestIndividual)
	{
		b_stop = pc_evaluation->bIsMaxValue(pcBestIndividual->dGetFitnessValue());
	}//if (pcBestIndividual)

	return b_stop;
}//bool CStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)


template <class TGenotype, class TFenotype>
CRunningTimeStopCondition<TGenotype, TFenotype>::CRunningTimeStopCondition(CEvaluation<TFenotype> *pcEvaluation)
	: CStopCondition<TGenotype, TFenotype>(pcEvaluation)
{

}//CRunningTimeStopCondition<TGenotype, TFenotype>::CRunningTimeStopCondition(CEvaluation<TFenotype> *pcEvaluation)

template <class TGenotype, class TFenotype>
CRunningTimeStopCondition<TGenotype, TFenotype>::CRunningTimeStopCondition(CEvaluation<TFenotype> *pcEvaluation, time_t tRunningTime)
	: CStopCondition<TGenotype, TFenotype>(pcEvaluation)
{
	t_running_time = tRunningTime;
}//CRunningTimeStopCondition<TGenotype, TFenotype>::CRunningTimeStopCondition(CEvaluation<TFenotype> *pcEvaluation, time_t tRunningTime)

template <class TGenotype, class TFenotype>
CError CRunningTimeStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	CError c_error = CStopCondition<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_running_time(STOP_CONDITION_ARGUMENT_TIME);
		t_running_time = (time_t)p_running_time.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CRunningTimeStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)

template <class TGenotype, class TFenotype>
bool CRunningTimeStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)
{
	bool b_stop = CStopCondition::bStop(tStart, iIterationNumber, iFFE, pcBestIndividual);

	if (!b_stop)
	{
		time_t t_now = time(nullptr);
		b_stop = tStart + t_running_time <= t_now;
	}//if (!b_stop)

	return b_stop;
}//bool CRunningTimeStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)


template <class TGenotype, class TFenotype>
CIterationsStopCondition<TGenotype, TFenotype>::CIterationsStopCondition(CEvaluation<TFenotype> *pcEvaluation)
	: CStopCondition<TGenotype, TFenotype>(pcEvaluation)
{

}//CIterationsStopCondition<TGenotype, TFenotype>::CIterationsStopCondition(CEvaluation<TFenotype> *pcEvaluation)

template <class TGenotype, class TFenotype>
CError CIterationsStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	CError c_error = CStopCondition<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_max_number_of_iterations(STOP_CONDITION_ARGUMENT_MAX_NUMBER_OF_ITERATIONS);
		i_max_number_of_iterations = p_max_number_of_iterations.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CIterationsStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)

template <class TGenotype, class TFenotype>
bool CIterationsStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)
{
	return CStopCondition<TGenotype, TFenotype>::bStop(tStart, iIterationNumber, iFFE, pcBestIndividual) || iIterationNumber >= i_max_number_of_iterations;
}//bool CIterationsStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)


template <class TGenotype, class TFenotype>
CProgressStopCondition<TGenotype, TFenotype>::CProgressStopCondition(CEvaluation<TFenotype> *pcEvaluation)
	: CStopCondition<TGenotype, TFenotype>(pcEvaluation)
{
	pc_last_best = nullptr;
	i_last_best_iteration = 0;
}//CProgressStopCondition<TGenotype, TFenotype>::CProgressStopCondition(CEvaluation<TFenotype> *pcEvaluation)


template <class TGenotype, class TFenotype>
CError CProgressStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	CError c_error = CStopCondition<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_max_number_of_iterations_without_progress(STOP_CONDITION_ARGUMENT_MAX_NUMBER_OF_ITERATIONS_WITHOUT_PROGRESS);
		i_max_number_of_iterations_without_progress = p_max_number_of_iterations_without_progress.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CProgressStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)

template <class TGenotype, class TFenotype>
bool CProgressStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)
{
	if (i_last_best_iteration > iIterationNumber || pc_last_best != pcBestIndividual)
	{
		pc_last_best = pcBestIndividual;
		i_last_best_iteration = iIterationNumber;
	}//if (i_last_best_iteration > iIterationNumber || pc_last_best != pcBestIndividual)

	return CStopCondition<TGenotype, TFenotype>::bStop(tStart, iIterationNumber, iFFE, pcBestIndividual) 
		|| iIterationNumber - i_last_best_iteration > i_max_number_of_iterations_without_progress;
}//bool CProgressStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)


template <class TGenotype, class TFenotype>
CFFEStopCondition<TGenotype, TFenotype>::CFFEStopCondition(CEvaluation<TFenotype> *pcEvaluation)
	: CStopCondition<TGenotype, TFenotype>(pcEvaluation)
{

}//CFFEStopCondition<TGenotype, TFenotype>::CFFEStopCondition(CEvaluation<TFenotype> *pcEvaluation)

template <class TGenotype, class TFenotype>
CError CFFEStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	CError c_error = CStopCondition<TGenotype, TFenotype>::eConfigure(psSettings);

	if (!c_error)
	{
		CUIntCommandParam p_max_number_of_ffe(STOP_CONDITION_ARGUMENT_MAX_NUMBER_OF_FFE);
		i_max_number_of_ffe = (uint64_t)p_max_number_of_ffe.iGetValue(psSettings, &c_error);
	}//if (!c_error)

	return c_error;
}//CError CFFEStopCondition<TGenotype, TFenotype>::eConfigure(istream *psSettings)


template <class TGenotype, class TFenotype>
bool CFFEStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)
{
	return CStopCondition<TGenotype, TFenotype>::bStop(tStart, iIterationNumber, iFFE, pcBestIndividual) || iFFE >= i_max_number_of_ffe;
}//bool CFFEStopCondition<TGenotype, TFenotype>::bStop(time_t tStart, uint32_t iIterationNumber, uint64_t iFFE, CIndividual<TGenotype, TFenotype> *pcBestIndividual)


template class CStopCondition<CBinaryCoding, CBinaryCoding>;
template class CRunningTimeStopCondition<CBinaryCoding, CBinaryCoding>;
template class CIterationsStopCondition<CBinaryCoding, CBinaryCoding>;
template class CProgressStopCondition<CBinaryCoding, CBinaryCoding>;
template class CFFEStopCondition<CBinaryCoding, CBinaryCoding>;