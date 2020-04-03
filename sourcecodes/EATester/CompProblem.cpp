//#include "stdafx.h"

#include "CompProblem.h"



using namespace ProblemTools;
using namespace std;
using namespace Tools;



CError  ProblemTools::eLoadSingleUnitationProblem(FILE *pfSource, CUnitationBased **ppcUnitationProblem)
{
	CError  c_err;
	CString  s_problem_type, s_setting_name;
	CString  s_buf;

	*ppcUnitationProblem = NULL;


	s_problem_type  =  Tools::sReadLine(pfSource, &s_setting_name);
	//::vShow(s_problem_type);
	if  (s_setting_name  !=  COMP_PROBLEM_TYPE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_TYPE, s_setting_name);
		c_err.vSetError(s_buf);
		fclose(pfSource);
		return(c_err);	
	}//if  (s_setting_name  !=  COMP_PROBLEM_TYPE)


	if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)
	{
		*ppcUnitationProblem = new CDecFunc1Modal(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)

	if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)
	{
		*ppcUnitationProblem = new CDecFunc2Modal(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)

	if  (s_problem_type == DEC_FUNC_EASY_NAME)
	{
		*ppcUnitationProblem = new CDecFuncEasy(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_EASY_NAME)

	if  (s_problem_type == DEC_FUNC_ANY_NAME)
	{
		*ppcUnitationProblem = new CDecFuncAny(1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_ANY_NAME)

	if  (*ppcUnitationProblem  ==  NULL)
	{
		s_buf.Format("Unknown function type: (%s)", s_problem_type);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (*pcFunc  ==  NULL)

	c_err = (*ppcUnitationProblem)->eLoadSettings(pfSource);

	return(c_err);
}//CError  ProblemTools::eLoadSingleUnitationProblem(FILE *pfSource, CUnitationBased **ppcUnitationProblem)



CError  ProblemTools::eLoadSingleCompProblem(FILE *pfSource, CCompProblem **ppcCompProblem)
{
	CError  c_err;
	CString  s_problem_type, s_setting_name;
	CString  s_buf;

	*ppcCompProblem = NULL;


	s_problem_type  =  Tools::sReadLine(pfSource, &s_setting_name);
	//::vShow(s_problem_type);
	if  (s_setting_name  !=  COMP_PROBLEM_TYPE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_TYPE, s_setting_name);
		c_err.vSetError(s_buf);
		fclose(pfSource);
		return(c_err);	
	}//if  (s_setting_name  !=  COMP_PROBLEM_TYPE)


	if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)
	{
		*ppcCompProblem = new CDecFunc1Modal(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)

	if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)
	{
		*ppcCompProblem = new CDecFunc2Modal(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)

	if  (s_problem_type == DEC_FUNC_EASY_NAME)
	{
		*ppcCompProblem = new CDecFuncEasy(1, 1, 0, 1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_EASY_NAME)

	if  (s_problem_type == DEC_FUNC_ANY_NAME)
	{
		*ppcCompProblem = new CDecFuncAny(1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_ANY_NAME)

	if (s_problem_type == DEC_FUNC_ANY_SCALED_NAME)
	{
		*ppcCompProblem = new CDecFuncAnyScaled(1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_ANY_NAME)
	
	if (s_problem_type == DEC_FUNC_ANY_SCALED_NOISED_NAME)
	{
		*ppcCompProblem = new CDecFuncAnyScaledNoised(1, &c_err);
	}//if  (s_problem_type == DEC_FUNC_ANY_NAME)
				
	if  (s_problem_type == HIFF_PROBLEM_NAME)
	{
		*ppcCompProblem = new CHiffProblem();
	}//if  (s_problem_type == HIFF_PROBLEM_NAME)

	if  (s_problem_type == KNAPSACK_PROBLEM_NAME)
	{
		*ppcCompProblem = new CKnapsackProblem();
	}//if  (s_problem_type == KNAPSACK_PROBLEM_NAME)

	if (s_problem_type == MAX_SAT_PROBLEM_NAME)
	{
		*ppcCompProblem = new CMaxSatProblem();
	}//if (s_problem_type == MAX_SAT_PROBLEM_NAME)

	if  (*ppcCompProblem  ==  NULL)
	{
		s_buf.Format("Unknown function type: (%s)", s_problem_type);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (*pcFunc  ==  NULL)

	c_err = (*ppcCompProblem)->eLoadSettings(pfSource);

	return(c_err);
}//CError  ProblemTools::eLoadSingleCompProblem(FILE *pfSource, CCompProblem **ppcCompProblem)


//---------------------------class  CConcatCompProblem--------------------------------------


CConcatCompProblem::CConcatCompProblem()
{

}//CConcatCompProblem::CConcatCompProblem()


CConcatCompProblem::~CConcatCompProblem()
{
	for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.size();  ii++)
		delete  v_COMP_PROBLEMs.at(ii);
}//CConcatCompProblem::~CConcatCompProblem()



void  CConcatCompProblem::operator=(CConcatCompProblem  &cOther)
{
	CError  c_err;

	for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.size();  ii++)
		delete  v_COMP_PROBLEMs.at(ii);

	for  (int  ii = 0;  ii < (int) cOther.v_COMP_PROBLEMs.size();  ii++)
	{
		c_err  =  eAddCompProblem(cOther.v_COMP_PROBLEMs.at(ii));
		if  (c_err)  return;
	}//for  (int  ii = 0;  ii < (int) cOther.v_COMP_PROBLEMs.size();  ii++)

}//void  CConcatCompProblem::operator=(CConcatCompProblem  &cOther)



int  CConcatCompProblem::iGetProblemBitLength()
{
	int  i_problem_bit_len;

	i_problem_bit_len = 0;

	for  (int  ii = 0; ii < v_COMP_PROBLEMs.size(); ii++)
		i_problem_bit_len += v_COMP_PROBLEMs.at(ii)->iGetBitLength();

	return(i_problem_bit_len);
}//int  CConcatCompProblem::iGetProblemBitLength()




CError  CConcatCompProblem::eLoad(FILE *pfSource)
{
	CError  c_err;
	CString  s_problem_type, s_setting_name;
	CString  s_buf;
	int  i_problem_multiplier;


	while  (feof(pfSource) != true)
	{
		i_problem_multiplier  =  Tools::iReadLine(pfSource, &s_setting_name);
		//::vShow(i_problem_multiplier);
		if  (s_setting_name == COMP_PROBLEM_DEF_END)  return(c_err);

		
		if  (s_setting_name  !=  COMP_PROBLEM_MULTIPLIER)
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_MULTIPLIER, s_setting_name);
			c_err.vSetError(s_buf);
			fclose(pfSource);
			return(c_err);	
		}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)

		
		/*s_problem_type  =  Tools::sReadLine(pfSource, &s_setting_name);
		//::vShow(s_problem_type);
		if  (s_setting_name  !=  DEC_FUNC)
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", DEC_FUNC, s_setting_name);
			c_err.vSetError(s_buf);
			fclose(pfSource);
			return(c_err);	
		}//if  (s_setting_name  !=  DEC_FUNC_1_MODAL_NAME)*/

		CCompProblem  *pc_fitness = NULL;
		c_err = eLoadSingleCompProblem(pfSource, &pc_fitness);

		//FILE  *pf_single;
		//pf_single = fopen("zz_sin_test.txt", "w+");
		//pc_fitness->eCreateReportSummary(pf_single);
		//fclose(pf_single);

		if  (c_err)
		{
			fclose(pfSource);				
			return(c_err);	
		}//if  (c_err)

		if  (pc_fitness == NULL)
		{
			c_err.vSetError("No fitness function read");
			fclose(pfSource);				
			return(c_err);			
		}//if  (pc_fitness == NULL)

		for  (int  ij = 0; ij < i_problem_multiplier; ij++)
		{
			c_err  =  eAddCompProblem(pc_fitness);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)
		}//for  (int  ii = 0; ii < 4; ii++)

		delete  pc_fitness;


		/*if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)
		{
			CDecFunc1Modal  c_fitness(1, 1, 0, 1, &c_err);
			c_err = c_fitness.eLoadSettings(pfSource);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)

			for  (int  ij = 0; ij < i_problem_multiplier; ij++)
			{
				c_err  =  eAddDecFunc(&c_fitness);
				if  (c_err)
				{
					fclose(pfSource);				
					return(c_err);	
				}//if  (c_err)
			}//for  (int  ii = 0; ii < 4; ii++)
		}//if  (s_problem_type == DEC_FUNC_1_MODAL_NAME)



		if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)
		{
			CDecFunc2Modal  c_fitness(1, 1, 0, 1, &c_err);
			c_err = c_fitness.eLoadSettings(pfSource);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)

			for  (int  ij = 0; ij < i_problem_multiplier; ij++)
			{
				c_err  =  eAddDecFunc(&c_fitness);
				if  (c_err)
				{
					fclose(pfSource);				
					return(c_err);	
				}//if  (c_err)
			}//for  (int  ii = 0; ii < 4; ii++)
		}//if  (s_problem_type == DEC_FUNC_2_MODAL_NAME)



		if  (s_problem_type == DEC_FUNC_EASY_NAME)
		{
			CDecFuncEasy  c_easy_fitness(1, 1, 0, 1, &c_err);
			c_err = c_easy_fitness.eLoadSettings(pfSource);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)

			for  (int  ij = 0; ij < i_problem_multiplier; ij++)
			{
				c_err  =  eAddDecFunc(&c_easy_fitness);
				if  (c_err)
				{
					fclose(pfSource);				
					return(c_err);	
				}//if  (c_err)
			}//for  (int  ii = 0; ii < 4; ii++)
		}//if  (s_problem_type == DEC_FUNC_EASY_NAME)


		if  (s_problem_type == DEC_FUNC_ANY_NAME)
		{
			CDecFuncAny  c_any_fitness(1, &c_err);
			c_err = c_any_fitness.eLoadSettings(pfSource);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)

			for  (int  ij = 0; ij < i_problem_multiplier; ij++)
			{
				c_err  =  eAddDecFunc(&c_any_fitness);
				if  (c_err)
				{
					fclose(pfSource);				
					return(c_err);	
				}//if  (c_err)
			}//for  (int  ii = 0; ii < 4; ii++)
		}//if  (s_problem_type == DEC_FUNC_EASY_NAME)


		if  (s_problem_type == HIFF_PROBLEM_NAME)
		{
			CDecFuncAny  c_any_fitness(1, &c_err);
			c_err = c_any_fitness.eLoadSettings(pfSource);
			if  (c_err)
			{
				fclose(pfSource);				
				return(c_err);	
			}//if  (c_err)

			for  (int  ij = 0; ij < i_problem_multiplier; ij++)
			{
				c_err  =  eAddDecFunc(&c_any_fitness);
				if  (c_err)
				{
					fclose(pfSource);				
					return(c_err);	
				}//if  (c_err)
			}//for  (int  ii = 0; ii < 4; ii++)
		}//if  (s_problem_type == DEC_FUNC_EASY_NAME)*/


	}//while  (feof(pfSource) != true)



	/*CString  s_setting_name;
	CString  s_buf;
	CDecFunc1Modal  c_fitness(1, 1, 0, 1, &c_err);
	CDecFuncEasy  c_easy_fitness(1, 1, 0, 1, &c_err);
	int  i_problem_multiplier;
	int  i_problem_bit_length = 0;
	for  (int  ii = 0; (s_setting_name !=  DEC_KIND_CONCAT_END)&&(feof(pfSource) != true);ii++)
	{
		i_problem_multiplier  =  Tools::iReadLine(pfSource, &s_setting_name);
		if  ( (s_setting_name  !=  COMP_PROBLEM_MULTIPLIER)&&(s_setting_name  !=  DEC_KIND_CONCAT_END)&&(s_setting_name  !=  DEC_FUNC) )
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_MULTIPLIER, s_setting_name);
			c_err.vSetError(s_buf);
			fclose(pfSource);
			return(c_err);	
		}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


		if  (s_setting_name  !=  DEC_KIND_CONCAT_END)
		{
			if  (s_setting_name  !=  DEC_FUNC)
			{
				c_err  =  c_fitness.eLoadSettings(pfSource);

				if  (c_err)
				{
					fclose(pfSource);
					return(c_err);		
				}//if  (c_err)

				for  (int  ij = 0; ij < i_problem_multiplier; ij++)
				{
					i_problem_bit_length  +=  c_fitness.iGetBitLength();
					c_err  =  eAddDecFunc(&c_fitness);
					if  (c_err)
					{
						fclose(pfSource);				
						return(c_err);	
					}//if  (c_err)
				}//for  (int  ii = 0; ii < 4; ii++)
			}//if  (s_setting_name  !=  SYSTEM_PARAM_DEC_FUNC)
			else
			{
				i_problem_multiplier  =  Tools::iReadLine(pfSource, &s_setting_name);
				if  ( (s_setting_name  !=  COMP_PROBLEM_MULTIPLIER)&&(s_setting_name  !=  DEC_KIND_CONCAT_END)&&(s_setting_name  !=  DEC_FUNC) )
				{
					s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_MULTIPLIER, s_setting_name);
					c_err.vSetError(s_buf);
					fclose(pfSource);
					return(c_err);	
				}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)

				c_err  =  c_easy_fitness.eLoadSettings(pfSource);

				if  (c_err)
				{
					fclose(pfSource);
					return(c_err);		
				}//if  (c_err)

				for  (int  ij = 0; ij < i_problem_multiplier; ij++)
				{
					i_problem_bit_length  +=  c_easy_fitness.iGetBitLength();
					c_err  =  eAddDecFunc(&c_easy_fitness);
					if  (c_err)
					{
						fclose(pfSource);				
						return(c_err);	
					}//if  (c_err)
				}//for  (int  ii = 0; ii < 4; ii++)
				
			}//else  if  (s_setting_name  !=  SYSTEM_PARAM_DEC_FUNC)

		
		}//if  (s_setting_name  !=  SYSTEM_PARAM_DEC_END)
	}//for  (int  ii = 0; (s_setting_name !=  SYSTEM_PARAM_DEC_END)&&(feof(pf_buf) != true);ii++)*/


	return(c_err);
}//CError  CConcatCompProblem::eLoad(FILE *pfSource)


CError  CConcatCompProblem::eAddCompProblem(CCompProblem  *pcCompProblem)  
{
	CError  c_err;
	CCompProblem  *pc_new_comp_problem;

	pc_new_comp_problem  =  pcCompProblem->pcClone();
	if  (pc_new_comp_problem  ==  NULL)
	{
		c_err.vSetError("Could not add computation problem - probably unknown type");
		return(c_err);
	}//if  (pc_new_comp_problem  ==  NULL)

	v_COMP_PROBLEMs.push_back(pc_new_comp_problem);
	return(c_err);
};//void  CConcatCompProblem::eAddCompProblem(CCompProblem  *pcCompProblem)  



CError  CConcatCompProblem::eCreateReportSeparateFuncSummary(FILE  *pfReport)
{
	CError  c_err;

	fprintf(pfReport, "FITNESS FUNCTION VALUES\n\n");

	for  (int  ii=0; ii < (int) v_COMP_PROBLEMs.size();  ii++)
	{
		fprintf(pfReport, "Subfunction (%d/%d)\n", ii, v_COMP_PROBLEMs.size());
		v_COMP_PROBLEMs.at(ii)->eCreateReportSummary(pfReport);
	}//for  (int  ii=0; ii < (int) v_COMP_PROBLEMs.size();  ii++)

	return(c_err);
}//CError  CConcatCompProblem::eCreateReportSeparateFunc(FILE  *pfReport)




CError  CConcatCompProblem::eCreateReportSeparateFunc(FILE  *pfReport)
{
	CError  c_err;

	fprintf(pfReport, "FITNESS FUNCTION VALUES\n\n");

	for  (int  ii=0; ii < (int) v_COMP_PROBLEMs.size();  ii++)
	{
		fprintf(pfReport, "Subfunction (%d/%d)\n", ii, v_COMP_PROBLEMs.size());
		v_COMP_PROBLEMs.at(ii)->eCreateReport(pfReport);
	}//for  (int  ii=0; ii < (int) v_COMP_PROBLEMs.size();  ii++)

	return(c_err);
}//CError  CConcatCompProblem::eCreateReportSeparateFunc(FILE  *pfReport)



CError  CConcatCompProblem::eCreateReport(FILE  *pfReport)
{
	CError  c_err;
	int  *pi_argument;
	int  i_bit_length  =  0;

	for  (int  ii=0; ii < (int) v_COMP_PROBLEMs.size();  ii++)
		i_bit_length  +=  v_COMP_PROBLEMs.at(ii)->iGetBitLength();


	pi_argument  =  new  int[i_bit_length];
	vector  <CReportGrid *>  v_result_grids;

	for  (int  ii=0; ii < i_bit_length;  ii++)  pi_argument[ii]  =  0;


	c_err  =  e_create_report(&v_result_grids,  pi_argument,  i_bit_length,  0);


	sort(v_result_grids.begin(), v_result_grids.end(), SReportGridSort());
	for  (int  ii = 0;  ii < (int) v_result_grids.size();  ii++)
	{
		v_result_grids.at(ii)->vPrint(pfReport);
		fprintf(pfReport,  "\n");	
	}//for  (int  ii = 0;  ii < v_result_grids.size();  ii++)

    for  (int  ii = 0; ii < (int) v_result_grids.size();  ii++)
		delete  v_result_grids.at(ii);

	if  (c_err)  return(c_err);

	delete  pi_argument;
	return(c_err);
};//void  CConcatCompProblem::vCreateReport(FILE  *pfReport)




CError  CConcatCompProblem::e_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  iBitLength,  int  iCurrentBit)
{
	CError  c_err;

	if  (iCurrentBit  <  iBitLength)
	{
		piArgument[iCurrentBit]  =  0;
		c_err  =  e_create_report(pvReportGrids,  piArgument,  iBitLength,  iCurrentBit+1);
		if  (c_err)  return(c_err);


		piArgument[iCurrentBit]  =  1;
		c_err  =  e_create_report(pvReportGrids,  piArgument,  iBitLength,  iCurrentBit+1);
		if  (c_err)  return(c_err);
	
	}//if  (iCurrentBit  <  i_bit_length)
	else
	{
		double  d_func_val;
		CReportGrid  *pc_rep_grid;

		c_err  =  eGetFuncValue(piArgument,  iBitLength,  &d_func_val);
		if  (c_err)  return(c_err);
		pc_rep_grid  =  new  CReportGrid(iBitLength, piArgument, d_func_val);
		pvReportGrids->push_back(pc_rep_grid);
	}//else  if  (iCurrentBit  <  i_bit_length)


	return(c_err);
}//void  CConcatCompProblem::v_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  ICurrentBit)




double  CConcatCompProblem::dGetMaxFuncVal()
{
	double  d_result;

	d_result  =  0;
	for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.size();  ii++)
	{
		d_result  +=  v_COMP_PROBLEMs.at(ii)->dGetMaxFuncVal();
	}//for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.size();  ii++)

	//for overlapping
	//d_result += 100;


	return(d_result);
}//double  CConcatCompProblem::dGetMaxFuncVal()



CError  CConcatCompProblem::eGetFuncValue(int  *piArgument,  int  iArgLen,  double  *pdFuncVal)
{
	CError  c_err;

	double d_single_problem_value;
	double  d_sum  =  0;
	int  i_unitation;
	int  iCurrentBit  =  0;

	for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.size();  ii++)
	{
		d_single_problem_value  =  v_COMP_PROBLEMs.at(ii)->dGetFuncValue(iCurrentBit, piArgument, iArgLen);

		if (d_single_problem_value == -1)
		{
			*pdFuncVal = -1;
			return(c_err);
		}//if (d_single_problem_value == -1)

		d_sum  +=  d_single_problem_value;


		//CString  s_buf;
		//s_buf.Format("Shift: %d  Length: %d  Bit len: %d", iCurrentBit, iArgLen, v_COMP_PROBLEMs.at(ii)->iGetBitLength());
		//Tools::vReportInFile("test.txt", s_buf);
		//Tools::vShow(0);


		iCurrentBit += v_COMP_PROBLEMs.at(ii)->iGetBitLength();

		/*if  (v_COMP_PROBLEMs.at(ii)->iGetFuncType() != DEC_FUNC_HIFF)
		{
			i_unitation  =  0;
			for  (int  ik = 0;  ik < (int) v_COMP_PROBLEMs.at(ii)->iGetBitLength();  ik++)
			{
				if  (iCurrentBit  >=  iArgLen)
				{
					c_err.vSetError("Argument legnth below function argument legnth");
					return(c_err);
				}//if  (iCurrentBit  >=  iArgLen)

				i_unitation  +=  piArgument[iCurrentBit];
				iCurrentBit++;		
			}//for  (int  ii = 0;  ii < (int) v_COMP_PROBLEMs.at(ii)->iGetBitLength();  ii++)

			d_sum  +=  v_COMP_PROBLEMs.at(ii)->dGetFuncValue(i_unitation);	
		}//if  (v_COMP_PROBLEMs.at(ii)->iGetFuncType() != DEC_FUNC_HIFF)
		else
		{
		
		}//else  if  (v_COMP_PROBLEMs.at(ii)->iGetFuncType() != DEC_FUNC_HIFF)*/

		
	}//for  (int  ii = 0;  ii < v_COMP_PROBLEMs.size();  ii++)


	if  (iCurrentBit  <  iArgLen)
	{
		c_err.vSetError("Argument legnth exceeds function argument legnth");
		return(c_err);
	}//if  (iCurrentBit  >=  iArgLen)


	/*//for overlapping problems :)
	int  i_ones_num = 0;
	for (int i_bit = 2; i_bit < iArgLen; i_bit+=5)
	{
		if (piArgument[i_bit] == 1)  i_ones_num++;
	}//for (int i_bit = 1; i_bit < iArgLen; i_bit++)
	if (i_ones_num % 4 == 0)  d_sum += 50;

	for (int i_bit = 3; i_bit < iArgLen; i_bit += 5)
	{
		if (piArgument[i_bit] == 1)  i_ones_num++;
	}//for (int i_bit = 1; i_bit < iArgLen; i_bit++)
	if (i_ones_num % 3 == 0)  d_sum += 50;*/



	*pdFuncVal  =  d_sum;
	d_func_evaluations++;
	/*CString  s_buf;
	s_buf.Format("%.0lf", d_func_evaluations);
	::MessageBox(NULL, s_buf, "", MB_OK);*/

	return(c_err);
}//CError  CConcatCompProblem::eGetFuncValue(int  *piArgument,  int  iArgLen,  double  *pdFuncVal)

double CConcatCompProblem::d_func_evaluations = 0;



//---------------------------class  CCompProblem--------------------------------------


CCompProblem::CCompProblem()
{
	i_bit_length  =  0;
	d_max_val  =  0;
	d_min_val  =  0;
	d_min_max_dist_grid  =  0;
}//CCompProblem::CCompProblem()


CCompProblem::CCompProblem(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale,  CError  *pcErr)
{
	i_bit_length  =  0;
	d_max_val  =  0;
	d_min_val  =  0;
	d_min_max_dist_grid  =  0;
	d_lower_than_max_scaling  =  dLowerThanMaxScale;
	*pcErr  =  eConfigure(iBitLength,  dMaxVal,  dMinVal, dLowerThanMaxScale);
}//CCompProblem::CCompProblem(int  iBitLength,  double  dMaxVal,  double  dMinVal)


CCompProblem::~CCompProblem()
{

}//CCompProblem::~CCompProblem()


CError  CCompProblem::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale)
{
	CError  c_err;

	if  (iBitLength  <  0)
	{
		c_err.vSetError("Bit length below 0");
		return(c_err);
	}//if  (iBitLength  <  0)

	if  (d_max_val  <  d_min_val)
	{
		c_err.vSetError("Min function value greater than max function value");
		return(c_err);
	}//if  (iBitLength  <  0)

	i_bit_length  =  iBitLength;
	d_max_val  =  dMaxVal;
	d_min_val  =  dMinVal;
	d_lower_than_max_scaling  =  dLowerThanMaxScale;
	//d_min_max_dist_grid  =  (d_max_val  -  d_min_val)/(i_bit_length/2);

	return(c_err);
}//CError  CCompProblem::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal)



CError  CCompProblem::eCreateReportUnitation(FILE  *pfReport)
{
	CError  c_err;
	int  *pi_argument;

	
	pi_argument = new  int[i_bit_length];

	fprintf(pfReport, "Unitation\tFitness\n\n");

	for (int ii = 0; ii < i_bit_length; ii++)  pi_argument[ii] = 0;

	double  d_val;
	for (int ii = -1; ii < i_bit_length; ii++)
	{
		if  (ii >= 0)  pi_argument[ii] = 1;

		d_val = dGetFuncValue(0, pi_argument, i_bit_length);

		fprintf(pfReport, "%d\t%.2lf\n", ii+1, d_val);
	}//for (int ii = 0; ii < i_bit_length; ii++)
	
	

	delete  pi_argument;
	return(c_err);
}//CError  CCompProblem::eCreateReportUnitation(FILE  *pfReport)



CError  CCompProblem::eCreateReportSummary(FILE  *pfReport)
{
	CError  c_err;
	int  *pi_argument;

	//fprintf(pfReport, "d_lower_than_max_scaling = %.2lf\n", d_lower_than_max_scaling);
	//fprintf(pfReport, "d_max_val = %.2lf\n", d_max_val);
	//fprintf(pfReport, "d_min_max_dist_grid = %.2lf\n", d_min_max_dist_grid);
	//fprintf(pfReport, "d_min_val = %.2lf\n", d_min_val);
	//fprintf(pfReport, "i_bit_length = %d\n", i_bit_length);

	pi_argument  =  new  int[i_bit_length];
	vector  <CReportGrid *>  v_result_grids;

	for  (int  ii=0; ii < i_bit_length;  ii++)  pi_argument[ii]  =  0;


	c_err  =  e_create_report(&v_result_grids,  pi_argument,  0, true);


	sort(v_result_grids.begin(), v_result_grids.end(), SReportGridSort());
	for  (int  ii = 0;  ii < (int) v_result_grids.size();  ii++)
	{
		v_result_grids.at(ii)->vPrint(pfReport);
		fprintf(pfReport,  "\n");	
	}//for  (int  ii = 0;  ii < v_result_grids.size();  ii++)

    for  (int  ii = 0; ii < (int) v_result_grids.size();  ii++)
		delete  v_result_grids.at(ii);

	if  (c_err)  return(c_err);

	delete  pi_argument;
	return(c_err);
}//CError  CCompProblem::eCreateReportSummary(FILE  *pfReport)




CError  CCompProblem::eCreateReport(FILE  *pfReport)
{
	CError  c_err;
	int  *pi_argument;

	pi_argument  =  new  int[i_bit_length];
	vector  <CReportGrid *>  v_result_grids;

	for  (int  ii=0; ii < i_bit_length;  ii++)  pi_argument[ii]  =  0;


	c_err  =  e_create_report(&v_result_grids,  pi_argument,  0);


	sort(v_result_grids.begin(), v_result_grids.end(), SReportGridSort());
	for  (int  ii = 0;  ii < (int) v_result_grids.size();  ii++)
	{
		v_result_grids.at(ii)->vPrint(pfReport);
		fprintf(pfReport,  "\n");	
	}//for  (int  ii = 0;  ii < v_result_grids.size();  ii++)

    for  (int  ii = 0; ii < (int) v_result_grids.size();  ii++)
		delete  v_result_grids.at(ii);

	if  (c_err)  return(c_err);

	delete  pi_argument;
	return(c_err);
};//void  CCompProblem::vCreateReport(FILE  *pfReport)



CError  CCompProblem::e_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  iCurrentBit, bool  bLinear, int iCurrentValForLinear)
{
	CError  c_err;

	if  (iCurrentBit  <  i_bit_length)
	{
		if  ( (bLinear  ==  false)||(iCurrentValForLinear  ==  0) )
		{
			piArgument[iCurrentBit]  =  0;
			c_err  =  e_create_report(pvReportGrids,  piArgument,  iCurrentBit+1, bLinear, 0);
			if  (c_err)  return(c_err);
		}//if  ( (bLinear  ==  false)||(iCurrentValForLinear  ==  0) )


		piArgument[iCurrentBit]  =  1;
		c_err  =  e_create_report(pvReportGrids,  piArgument,  iCurrentBit+1, bLinear, 1);
		if  (c_err)  return(c_err);
	
	}//if  (iCurrentBit  <  i_bit_length)
	else
	{
		double  d_func_val;
		CReportGrid  *pc_rep_grid;

		//Tools::vShow(i_bit_length);
		d_func_val = dGetFuncValue(0, piArgument,  i_bit_length);
		if  (c_err)  return(c_err);
		pc_rep_grid  =  new  CReportGrid(i_bit_length, piArgument, d_func_val);
		pvReportGrids->push_back(pc_rep_grid);
	}//else  if  (iCurrentBit  <  i_bit_length)


	return(c_err);
}//void  CCompProblem::v_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  ICurrentBit)



//---------------------------class  CUnitationBased--------------------------------------


CUnitationBased::CUnitationBased() : CCompProblem()
{

}//CUnitationBased::CUnitationBased()


CUnitationBased::CUnitationBased(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale,  CError  *pcErr) : CCompProblem(iBitLength, dMaxVal, dMinVal, dLowerThanMaxScale, pcErr)
{

}//CUnitationBased::CUnitationBased(int  iBitLength,  double  dMaxVal,  double  dMinVal)



/*
CError  CUnitationBased::eGetFuncValue(int  *piArgument,  int  iArgLen,  double  *pdFuncVal)
{
	CError  c_err;

	if  (iArgLen  !=  i_bit_length)
	{
		c_err.vSetError("Argument length not equals the function argument length");
		return(c_err);	
	}//if  (iArgLen  !=  i_bit_length)

	int  i_unitation  =  0;
    for  (int  ii = 0;  ii < i_bit_length;  ii++)
	{
		if  (piArgument[ii]  >  0)  i_unitation++;	
	}//for  (int  ii = 0;  ii < i_bit_length;  ii++)

	*pdFuncVal  =  dGetFuncValue(i_unitation);
	if  (c_err)  return(c_err);

	return(c_err);
};//CError  CUnitationBased::eGetFuncValue(int  *piArgument,  int  iArgLen,  double  *pdFuncVal)
*/


double  CUnitationBased::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	if  (iShift < 0)  return(-1);
	if  (iShift + iGetBitLength() >= iArgLen)  return(-1);

	int  i_unitation;
	int i_cur_pos;

	i_unitation = 0;
	for  (int i_pos = 0; i_pos < i_bit_length; i_pos++)
	{
		i_cur_pos = i_pos + iShift;
		if  (piArgument[i_cur_pos] == 1)  i_unitation++;
	}//for  (int i_pos = 0; i_pos < i_bit_length; i_pos++)

	return(dGetFuncValue(i_unitation));
}//double  CUnitationBased::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)


int  CUnitationBased::i_get_unitation(int  iShift, int  *piArgument,  int  iArgLen)
{
	//CString  s_buf;
	//s_buf.Format("Shift: %d  Length: %d  Bit len: %d", iShift, iArgLen, iGetBitLength());
	//Tools::vReportInFile("test.txt", s_buf);

	if  (iShift < 0)  return(-1);
	if  (iShift + iGetBitLength() > iArgLen)  return(-1);

	int  i_unitation;
	int i_cur_pos;

	i_unitation = 0;
	for  (int i_pos = 0; i_pos < i_bit_length; i_pos++)
	{
		i_cur_pos = i_pos + iShift;
		if  (piArgument[i_cur_pos] == 1)  i_unitation++;
	}//for  (int i_pos = 0; i_pos < i_bit_length; i_pos++)

	return(i_unitation);
}//double  CUnitationBased::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)



//---------------------------class  CConstraintBased--------------------------------------


CConstraintBased::CConstraintBased() : CCompProblem()
{

}//CConstraintBased::CConstraintBased()


CConstraintBased::CConstraintBased(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale,  CError  *pcErr) : CCompProblem(iBitLength, dMaxVal, dMinVal, dLowerThanMaxScale, pcErr)
{

}//CConstraintBased::CConstraintBased(int  iBitLength,  double  dMaxVal,  double  dMinVal)


bool CConstraintBased::bSatisfyConstraints(int  iShift, int  *piArgument,  int  iArgLen)
{
	int i_constraint_count = iGetConstraintCount();
	double *pd_constraints = pdGetConstraints();
	double *pd_constraints_values = new double[i_constraint_count];

	if (!bGetConstraintsValues(iShift, piArgument, iArgLen, pd_constraints_values)) return(false);

	bool b_satisfy = true;

	for (int i = 0; i < i_constraint_count && b_satisfy; i++)
		b_satisfy = *(pd_constraints_values + i) <= *(pd_constraints + i);

	delete pd_constraints_values;

	return(b_satisfy);
}//bool CConstraintBased::bSatisfyConstraints(int  iShift, int  *piArgument,  int  iArgLen)



//---------------------------class  CDecFuncEasy--------------------------------------

CDecFuncEasy::CDecFuncEasy()  :  CUnitationBased()
{

}//CDecFuncEasy::CDecFuncEasy()


CDecFuncEasy::CDecFuncEasy(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr)  :  CUnitationBased(iBitLength,  dMaxVal,  dMinVal,  dLowerThanMaxScale, pcErr)
{
	*pcErr  =  eConfigure(iBitLength,  dMaxVal,  dMinVal, dLowerThanMaxScale);
}//CDecFuncEasy::CDecFuncEasy(int  iBitLength,  double  dMaxVal,  double  dMinVal,  CError  *pcErr)


CDecFuncEasy::~CDecFuncEasy()
{

}//CDecFuncEasy::~CDecFuncEasy()



CError  CDecFuncEasy::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale)
{
	CError  c_err;

	c_err  =  CUnitationBased::eConfigure(iBitLength,  dMaxVal,  dMinVal, dLowerThanMaxScale);

	d_min_max_dist_grid  =  (d_max_val  -  d_min_val)/(i_bit_length);
	//d_min_max_dist_grid  =  d_min_max_dist_grid/dLowerThanMaxScale;

	return(c_err);
}//CError  CDecFunc2Modal::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal)




CError  CDecFuncEasy::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	CString  s_setting_name,  s_buf;


	
	i_bit_length  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	
	d_max_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_min_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MIN_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MIN_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_lower_than_max_scaling  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)



	c_err  =  eConfigure(i_bit_length,  d_max_val,  d_min_val, d_lower_than_max_scaling);
		

	return(c_err);
}//CError  CDecFuncEasy::eLoadSettings(FILE  *pfSettings)




double  CDecFuncEasy::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument,  iArgLen);
	if  (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFuncEasy::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)





double  CDecFuncEasy::dGetFuncValue(int  iUnitation)
{
	int  i_mirror_uni;


	if  (iUnitation  ==  i_bit_length)  return(d_max_val);

	
	
	double  d_func_val;
	d_func_val  =  d_min_val  +  ((double) iUnitation)  *  d_min_max_dist_grid;
	
	return(d_func_val);
}//double  CDecFuncEasy::dGetFuncValue(int  iUnitation)



CCompProblem* CDecFuncEasy::pcClone()
{
	CError  c_err;
	CDecFuncEasy  *pc_result;

	pc_result = new CDecFuncEasy(1, 1, 0, 1, &c_err);

	*pc_result = *this;

	return(pc_result);
}//CCompProblem* CDecFuncEasy::pcClone()




//---------------------------class  CDecFuncAny--------------------------------------

CDecFuncAny::CDecFuncAny()  :  CUnitationBased()
{

}//CDecFuncAny::CDecFuncAny()



CDecFuncAny::CDecFuncAny(int  iBitLength,  CError  *pcErr)
{
	*pcErr  =  eConfigure(iBitLength,  1,  -1, 0);
	d_max_val = 0;
	d_min_val = 0;
}//CDecFuncAny::CDecFuncAny(int  iBitLength,  CError  *pcErr)


CDecFuncAny::~CDecFuncAny()
{

}//CDecFuncAny::~CDecFuncAny()


CError  CDecFuncAny::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	CString  s_buf;
	CString  s_values_for_unitation, s_setting_name;

	int  i_index;
	double  d_value;


	i_bit_length  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)

	c_err  =  eConfigure(i_bit_length,  1,  -1, 0);
	

	i_index = 0;
	s_values_for_unitation  =  Tools::sReadLine(pfSettings, &s_setting_name);
	//::vShow(s_values_for_unitation);
	int  i_unitation = 0;
	while  (i_index  <  s_values_for_unitation.GetLength())
	{
		d_value = dExtractFromString(s_values_for_unitation, i_index, &i_index);

		c_err = eAddUnitationLevel(i_unitation, d_value);
		if  (c_err)  return(c_err);
		i_unitation++;
		//s_buf.Format("%lf  ind:%d   max:%d", d_value, i_index, s_values_for_unitation.GetLength());
		//::vShow(s_buf);
	}//while  (i_index  <  s_buf.GetLength())

	if  (i_unitation  <  i_bit_length)
	{
		s_buf.Format("Not enough unitation levels defined. Levels: %d  Defined: %d", i_bit_length, i_unitation);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if  (i_unitation  <  i_bit_length)


	for (int ii = 0; ii < v_unitation_values.size(); ii++)
	{
		if (ii == 0)
		{
			d_min_val = v_unitation_values.at(ii);
			d_max_val = v_unitation_values.at(ii);
		}//if (ii == 0)

		if (d_min_val > v_unitation_values.at(ii))  d_min_val = v_unitation_values.at(ii);
		if (d_max_val < v_unitation_values.at(ii))  d_max_val = v_unitation_values.at(ii);
	}//for (int ii = 0; ii < v_unitation_values.size(); ii++)

	return(c_err);
}//CError  CDecFuncAny::eLoadSettings(FILE  *pfSettings)



CError  CDecFuncAny::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale)
{
	CError  c_err;

	i_bit_length = iBitLength;
	v_unitation_values.clear();

	for  (int  ii = 0; ii < i_bit_length + 1; ii++)
		v_unitation_values.push_back(0);	

	d_max_val = 0;
	d_min_val = 0;

	return(c_err);
}//CError  CDecFuncAny::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale)



CError  CDecFuncAny::eAddUnitationLevel(int  iUnitation, double dUnitValue)
{
	CError  c_err;

	if  (iUnitation < 0)
	{
		CString  s_buf;

		s_buf.Format("Unitation below zero: %d", iUnitation);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if  (iUnitation < 0)


	if  (iUnitation >= v_unitation_values.size())
	{
		CString  s_buf;

		s_buf.Format("Unitation too high: %d >= %d", iUnitation, v_unitation_values.size());
		c_err.vSetError(s_buf);

		return(c_err);	
	}//if  (iUnitation >= v_unitation_values.size())


	v_unitation_values.at(iUnitation) = dUnitValue;

	for  (int  ii = 0; ii < v_unitation_values.size(); ii++)
	{
		if  (d_min_val > v_unitation_values.at(ii))  d_min_val = v_unitation_values.at(ii);
		if  (d_max_val < v_unitation_values.at(ii))  d_max_val = v_unitation_values.at(ii);	
	}//for  (int  ii = 0; ii < v_unitation_values.size(); ii++)

	return(c_err);
}//CError  CDecFuncAny::eAddUnitationLevel(double dUnitValue)




double  CDecFuncAny::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument,  iArgLen);
	if  (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFuncAny::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)




double  CDecFuncAny::dGetFuncValue(int  iUnitation)
{
	if  (iUnitation < 0)  return(0);
	if  (iUnitation >= v_unitation_values.size())  return(0);

	return(v_unitation_values.at(iUnitation));
}//double  CDecFuncAny::dGetFuncValue(int  iUnitation)


CCompProblem* CDecFuncAny::pcClone()
{
	CError  c_err;
	CDecFuncAny  *pc_result;

	pc_result = new CDecFuncAny(1, &c_err);
	(*pc_result) = (*this);

	return(pc_result);
}//CCompProblem* CDecFuncAny::pcClone()




 //---------------------------class  CHiffProblem--------------------------------------

CDecFuncAnyScaled::CDecFuncAnyScaled() : CDecFuncAny()
{
	i_bit_number_multiplier = 1;
}//CDecFuncAnyScaled::CDecFuncAnyScaled()



CDecFuncAnyScaled::CDecFuncAnyScaled(int  iBitLength, CError  *pcErr) : CDecFuncAny(iBitLength, pcErr)
{
	i_bit_number_multiplier = 1;
}//CDecFuncAnyScaled::CDecFuncAnyScaled(int  iBitLength,  CError  *pcErr)


CDecFuncAnyScaled::~CDecFuncAnyScaled()
{

}//CDecFuncAnyScaled::~CDecFuncAnyScaled()


CError  CDecFuncAnyScaled::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	
	c_err = CDecFuncAny::eLoadSettings(pfSettings);
	if (c_err)  return(c_err);
	
	
	CString  s_buf;
	CString  s_setting_name;
		

	i_bit_number_multiplier = Tools::iReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_BIT_NUMBER_MULTIPLIER)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_NUMBER_MULTIPLIER, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (s_setting_name  !=  COMP_PROBLEM_PARAM_SCALED_BIT_NUMBER)


	i_bit_length = i_bit_number_multiplier * i_bit_length;


	/*FILE  *pf_test;
	pf_test = fopen("test.txt", "w+");
	eCreateReportUnitation(pf_test);
	//eCreateReport(pf_test);
	fclose(pf_test);//*/

	return(c_err);
}//CError  CDecFuncAnyScaled::eLoadSettings(FILE  *pfSettings)



CCompProblem* CDecFuncAnyScaled::pcClone()
{
	CError  c_err;
	CDecFuncAnyScaled  *pc_result;

	pc_result = new CDecFuncAnyScaled(1, &c_err);
	(*pc_result) = (*this);

	return(pc_result);
}//CCompProblem* CDecFuncAnyScaled::pcClone()


CError  CDecFuncAnyScaled::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)
{
	CError  c_err;
	c_err = CDecFuncAny::eConfigure(iBitLength, dMaxVal, dMinVal, dLowerThanMaxScale);
	return(c_err);
}//CError  CDecFuncAnyScaled::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)


double  CDecFuncAnyScaled::dGetFuncValue(int  iShift, int  *piArgument, int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument, iArgLen);
	if (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFuncAnyScaled::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)




double  CDecFuncAnyScaled::dGetFuncValue(int  iUnitation)
{
	if (iUnitation < 0)  return(0);
	if (iUnitation >= v_unitation_values.size() * i_bit_number_multiplier)  return(0);

	for (int i_level = 0; i_level < v_unitation_values.size(); i_level++)
	{
		if (i_level * i_bit_number_multiplier == iUnitation)  return(v_unitation_values.at(i_level));
		if (i_level * i_bit_number_multiplier > iUnitation)
		{
			double  d_value, d_buf;
			d_buf = i_bit_number_multiplier - (i_level * i_bit_number_multiplier - iUnitation);
			d_buf = d_buf / i_bit_number_multiplier;

			d_value = d_buf * (v_unitation_values.at(i_level) - v_unitation_values.at(i_level - 1));
			d_value = v_unitation_values.at(i_level - 1) + d_value;

			return(d_value);
		}//if (i_level * i_bit_number_multiplier < iUnitation)
	}//for (int i_level = 0; i_level < v_unitation_values.size(); i_level++)

	return(0);
}//double  CDecFuncAnyScaled::dGetFuncValue(int  iUnitation)



//---------------------------class  CDecFuncAnyScaledNoised--------------------------------------

CDecFuncAnyScaledNoised::CDecFuncAnyScaledNoised() :CDecFuncAnyScaled()
{
	//d_unitation_noise_size = 0;
	i_unitation_noise_step_in_bits = 1;
}//CDecFuncAnyScaled::CDecFuncAnyScaled()



CDecFuncAnyScaledNoised::CDecFuncAnyScaledNoised(int  iBitLength, CError  *pcErr) : CDecFuncAnyScaled(iBitLength, pcErr)
{
	//d_unitation_noise_size = 0;
	i_unitation_noise_step_in_bits = 1;
}//CDecFuncAnyScaled::CDecFuncAnyScaled(int  iBitLength,  CError  *pcErr)


CDecFuncAnyScaledNoised::~CDecFuncAnyScaledNoised()
{

}//CDecFuncAnyScaledNoised::~CDecFuncAnyScaledNoised()


CError  CDecFuncAnyScaledNoised::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;

	c_err = CDecFuncAnyScaled::eLoadSettings(pfSettings);
	if (c_err)  return(c_err);


	CString  s_buf;
	CString  s_setting_name;



	i_unitation_noise_step_in_bits = Tools::iReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_UNITATION_NOISE_STEP)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_UNITATION_NOISE_STEP, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (s_setting_name  !=  COMP_PROBLEM_PARAM_UNITATION_NOISE_STEP)


	/*d_unitation_noise_size = Tools::dReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_UNITATION_NOISE_SIZE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_UNITATION_NOISE_SIZE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (s_setting_name  !=  COMP_PROBLEM_PARAM_UNITATION_NOISE_STEP)*/


	double  d_val_before, d_val_before_before, d_my_val;
	double  d_modification;
	for (int i_unitation = 0; i_unitation <= i_bit_length; i_unitation++)
	{
		CDecFuncAnyScaled::dGetFuncValue(i_unitation);

		d_modification = 0;
		if ((1 < i_unitation) && (i_unitation < i_bit_length))
		{
			if (i_unitation % i_unitation_noise_step_in_bits == 0)
			{
				d_val_before_before = CDecFuncAnyScaled::dGetFuncValue(i_unitation - 2);
				d_val_before = CDecFuncAnyScaled::dGetFuncValue(i_unitation - 1);
				d_my_val = CDecFuncAnyScaled::dGetFuncValue(i_unitation);

				d_modification = (d_val_before_before + d_val_before) / 2 - d_my_val;
			}//if (iUnitation % i_unitation_noise_step_in_bits == 0)			
		}//if ((ii > 1) && (ii < i_bit_length - 1))

		v_modifications.push_back(d_modification);
	}//for (int ii = 0; ii <= i_bit_length; ii++)
	

	 /*FILE  *pf_test;
	 pf_test = fopen("test.txt", "w+");
	 eCreateReportUnitation(pf_test);
	 //eCreateReport(pf_test);
	 fclose(pf_test);//*/


	return(c_err);
}//CError  CDecFuncAnyScaledNoised::eLoadSettings(FILE  *pfSettings)




CError  CDecFuncAnyScaledNoised::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)
{
	CError  c_err;
	c_err = CDecFuncAnyScaled::eConfigure(iBitLength, dMaxVal, dMinVal, dLowerThanMaxScale);
	return(c_err);
}//CError  CDecFuncAnyScaledNoised::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)


CCompProblem* CDecFuncAnyScaledNoised::pcClone()
{
	CError  c_err;
	CDecFuncAnyScaledNoised  *pc_result;

	pc_result = new CDecFuncAnyScaledNoised(1, &c_err);
	(*pc_result) = (*this);

	return(pc_result);
}//CCompProblem* CDecFuncAnyScaledNoised::pcClone()



double  CDecFuncAnyScaledNoised::dGetFuncValue(int  iShift, int  *piArgument, int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument, iArgLen);
	if (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFuncAnyScaledNoised::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)




double  CDecFuncAnyScaledNoised::dGetFuncValue(int  iUnitation)
{
	double  d_result;

	d_result = CDecFuncAnyScaled::dGetFuncValue(iUnitation);
	d_result +=	v_modifications.at(iUnitation);

	return(d_result);
}//double  CDecFuncAnyScaledNoised::dGetFuncValue(int  iUnitation)




//---------------------------class  CHiffProblem--------------------------------------

CHiffProblem::CHiffProblem()  :  CCompProblem()
{
	pc_root = NULL;
}//CHiffProblem::CHiffProblem()



CHiffProblem::~CHiffProblem()
{

}//CHiffProblem::~CHiffProblem()



CError  CHiffProblem::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	CString  s_buf;
	CString  s_setting_name;

	
	int  i_index;
	double  d_value;


	i_demanded_tree_level  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  HIFF_PROBLEM_PARAM_TREE_DEPTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", HIFF_PROBLEM_PARAM_TREE_DEPTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	if  (i_demanded_tree_level  <  0)
	{
		s_buf.Format("Demanded tree level below 0:%d", i_demanded_tree_level);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (i_demanded_tree_level  <  0)

	//HIFF_PROBLEM_PARAM_MULTIPLIER


	CUnitationBased  *pc_fit_func;
	double  d_multiplier_buf;
	for  (int  i_level = 0; i_level < i_demanded_tree_level; i_level++)
	{
		d_multiplier_buf  =  Tools::dReadLine(pfSettings, &s_setting_name);
		if  (s_setting_name  !=  HIFF_PROBLEM_PARAM_MULTIPLIER)
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", HIFF_PROBLEM_PARAM_MULTIPLIER, s_setting_name);
			c_err.vSetError(s_buf);
			return(c_err);	
		}//if  (s_setting_name  !=  HIFF_PROBLEM_PARAM_MULTIPLIER)


		c_err = eLoadSingleUnitationProblem(pfSettings, &pc_fit_func);
		if  (c_err)  return(c_err);

		if  (i_level == 0)
		{
			v_all_hiff_nodes.push_back(CHiffNode(this));
			v_all_hiff_nodes.at(0).pc_parent = this;
			v_all_hiff_nodes.at(0).d_multiplier = d_multiplier_buf;
			v_all_hiff_nodes.at(0).pc_unitation_problem = (CUnitationBased*)pc_fit_func->pcClone();

			pc_root = &(v_all_hiff_nodes.at(0));
		}//if  (i_level == 0)
		else
		{
			int  i_first_node_to_fill_offset = -1;
			for  (int  ii = 0; (ii < v_all_hiff_nodes.size())&&(i_first_node_to_fill_offset  ==  -1); ii++)
			{
				if  (v_all_hiff_nodes.at(ii).v_child_funcs_offsets.size() == 0)  i_first_node_to_fill_offset = ii;
			}//for  (int  ii = 0; (ii < v_all_hiff_nodes.size())&&(i_first_node_to_fill_offset  ==  -1); ii++)

			if  (i_first_node_to_fill_offset  <  0)
			{
				c_err.vSetError("fist node to fill offset not found");
				return(c_err);
			}//if  (i_first_node_to_fill_offset  <  0)

			int  i_last_node_to_fill_offset = v_all_hiff_nodes.size()-1;

			//now we add the next level;
			for  (int  i_node_to_fill_with_children = i_first_node_to_fill_offset; i_node_to_fill_with_children <= i_last_node_to_fill_offset; i_node_to_fill_with_children++)
			{
				for  (int i_gene_num = 0; i_gene_num < v_all_hiff_nodes.at(i_node_to_fill_with_children).iGetBitLength(); i_gene_num++)
				{
					v_all_hiff_nodes.push_back(CHiffNode(this));
					v_all_hiff_nodes.at(v_all_hiff_nodes.size() - 1).d_multiplier = d_multiplier_buf;
					v_all_hiff_nodes.at(v_all_hiff_nodes.size() - 1).pc_unitation_problem = (CUnitationBased*)pc_fit_func->pcClone();

					v_all_hiff_nodes.at(i_node_to_fill_with_children).v_child_funcs_offsets.push_back(v_all_hiff_nodes.size() - 1);
				}//for  (int i_gene_num = 0; i_gene_num < v_all_hiff_nodes.at(i_node_to_fill_with_children).iGetBitLength(); i_gene_num++)			
			}//for  (int  i_node_to_fill_with_children = i_first_node_to_fill_offset; i_node_to_fill_with_children < i_last_node_to_fill_offset; i_node_to_fill_with_children++)
			
		
		}//else  if  (i_level == 0)*/

		delete  pc_fit_func;
		
	}//for  (int  i_level = 0; i_level < i_demanded_tree_level; i_level++)

	
	//now fill the last level with gene offsets
	int  i_first_node_to_fill_offset = -1;
	for  (int  ii = 0; (ii < v_all_hiff_nodes.size())&&(i_first_node_to_fill_offset  ==  -1); ii++)
	{
		if  (v_all_hiff_nodes.at(ii).v_child_funcs_offsets.size() == 0)  i_first_node_to_fill_offset = ii;
	}//for  (int  ii = 0; (ii < v_all_hiff_nodes.size())&&(i_first_node_to_fill_offset  ==  -1); ii++)

	if  (i_first_node_to_fill_offset  <  0)
	{
		c_err.vSetError("fist node to fill offset not found (genes)");
		return(c_err);
	}//if  (i_first_node_to_fill_offset  <  0)

	int  i_last_node_to_fill_offset = v_all_hiff_nodes.size()-1;

	//now we add the next level;
	int  i_gene_offset = 0;
	for  (int  i_node_to_fill_with_children = i_first_node_to_fill_offset; i_node_to_fill_with_children <= i_last_node_to_fill_offset; i_node_to_fill_with_children++)
	{
		for  (int i_gene_num = 0; i_gene_num < v_all_hiff_nodes.at(i_node_to_fill_with_children).iGetBitLength(); i_gene_num++)
		{
			v_all_hiff_nodes.at(i_node_to_fill_with_children).v_genes.push_back(i_gene_offset++);
		}//for  (int i_gene_num = 0; i_gene_num < v_all_hiff_nodes.at(i_node_to_fill_with_children).iGetBitLength(); i_gene_num++)			
	}//for  (int  i_node_to_fill_with_children = i_first_node_to_fill_offset; i_node_to_fill_with_children < i_last_node_to_fill_offset; i_node_to_fill_with_children++)
	


	i_bit_length = i_gene_offset;

	s_buf  =  Tools::sReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  HIFF_PROBLEM_PARAM_TREE_END)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", HIFF_PROBLEM_PARAM_TREE_END, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  HIFF_PROBLEM_PARAM_TREE_END)


	int  *pi_arg;
	pi_arg = new int[i_bit_length];
	for (int ii = 0; ii < i_bit_length; ii++)
		pi_arg[ii] = 1;
	d_max_val = dGetFuncValue(0, pi_arg, i_bit_length);
	delete pi_arg;

	//s_buf.Format("%.2lf", d_max_val);
	//::MessageBox(NULL, s_buf, s_buf, MB_OK);
		
	
	return(c_err);
}//CError  CHiffProblem::eLoadSettings(FILE  *pfSettings)



CError  CHiffProblem::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)
{
	CError  c_err;

	if (iBitLength  <  0)
	{
		c_err.vSetError("Bit length below 0");
		return(c_err);
	}//if  (iBitLength  <  0)

	if (d_max_val  <  d_min_val)
	{
		c_err.vSetError("Min function value greater than max function value");
		return(c_err);
	}//if  (iBitLength  <  0)

	i_bit_length = iBitLength;
	d_max_val = dMaxVal;
	d_min_val = dMinVal;
	d_lower_than_max_scaling = dLowerThanMaxScale;
	

	int  *pi_arg;
	pi_arg = new int[i_bit_length];
	for (int ii = 0; ii < i_bit_length; ii++)
		pi_arg[ii] = 1;
	d_max_val = dGetFuncValue(0, pi_arg, i_bit_length);
	delete pi_arg;

	CString  s_buf;
	s_buf.Format("%.2lf", d_max_val);
	::MessageBox(NULL, s_buf, s_buf, MB_OK);

	return(c_err);
}//CError  CHiffProblem::eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale)



double  CHiffProblem::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	double  d_result, d_hiff_node_res;

	d_result = 0;

	for  (int  ii = 0; ii < v_all_hiff_nodes.size(); ii++)
	{
		d_hiff_node_res = v_all_hiff_nodes.at(ii).dGetFuncValue(iShift, piArgument,  iArgLen);
		if  (d_hiff_node_res < 0)  
		{
			return(-1);
		}//if  (d_hiff_node_res < 0)  
		d_result += d_hiff_node_res;
	}//for  (int  ii = 0; ii < v_all_hiff_nodes.size(); ii++)

	return(d_result);
}//double  CHiffProblem::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)


CCompProblem* CHiffProblem::pcClone()
{
	CError  c_err;
	CHiffProblem  *pc_result;


	pc_result = new CHiffProblem();
	
	for  (int ii = 0; ii < v_all_hiff_nodes.size(); ii++)
	{
		pc_result->v_all_hiff_nodes.push_back(CHiffNode(v_all_hiff_nodes.at(ii)));
		pc_result->v_all_hiff_nodes.at(ii).pc_parent = pc_result;	
	}//for  (int ii = 0; ii < v_all_hiff_nodes.size(); ii++)

	if  (pc_result->v_all_hiff_nodes.size()  >  0)
		pc_result->pc_root = &(pc_result->v_all_hiff_nodes.at(0));
	else
		pc_result->pc_root = NULL;

	pc_result->i_demanded_tree_level = i_demanded_tree_level;
	pc_result->i_bit_length = i_bit_length;
	pc_result->d_min_val = d_min_val;
	pc_result->d_max_val = d_max_val;

	return(pc_result);
}//CCompProblem* CHiffProblem::pcClone()




//---------------------------class  CHiffNode--------------------------------------
CHiffNode::CHiffNode()
{
	pc_unitation_problem = NULL;
	pc_parent = NULL;
}//CHiffNode::CHiffNode()


CHiffNode::CHiffNode(CHiffProblem  *pcParent)
{
	pc_unitation_problem = NULL;
	pc_parent = pcParent;
}//CHiffNode::CHiffNode(CHiffProblem  *pcParent)


CHiffNode::~CHiffNode()
{
	if  (pc_unitation_problem != NULL)  delete  pc_unitation_problem;
}//CHiffNode::CHiffNode()


CHiffNode::CHiffNode(const CHiffNode  &pcOther)
{
	if  (pcOther.pc_unitation_problem  !=  NULL)
		pc_unitation_problem = (CUnitationBased*)pcOther.pc_unitation_problem->pcClone();
	else
		pc_unitation_problem = NULL;
	pc_parent = pcOther.pc_parent;
	d_multiplier = pcOther.d_multiplier;

	v_child_funcs_offsets = pcOther.v_child_funcs_offsets;
	v_genes = pcOther.v_genes;
}//CHiffNode::CHiffNode(CHiffNode  *pcOther)



double  CHiffNode::dGetFuncValue(int iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;
	
	i_unitation = iGetUnitation(iShift, piArgument,  iArgLen);
	if  (i_unitation  ==  -1)  return(0);

	double  d_value;
	d_value = pc_unitation_problem->dGetFuncValue(i_unitation);
	d_value *= d_multiplier;

	return(d_value);
}//CError  CHiffNode::eGetFuncValue(int iShift, int  *piArgument,  int  iArgLen,  double  *pdFuncVal)


int  CHiffNode::iGetBitValue(int iShift, int  *piArgument,  int  iArgLen)
{
	bool  b_zeroes, b_ones;

	b_zeroes = false;
	b_ones = false;

	int  i_gene_offset;

	i_gene_offset = 0;
	for  (int  ii = 0;  ii < v_genes.size(); ii++)
	{
		i_gene_offset = iShift + v_genes.at(ii);
		if  (i_gene_offset  <  0)  return(-1);
		if  (i_gene_offset  >=  iArgLen)  return(-1);

		if  (piArgument[i_gene_offset]  ==  0)  b_zeroes = true;
		if  (piArgument[i_gene_offset]  ==  1)  b_ones = true;

		if  ((b_zeroes == true)&&(b_ones == true))  return(-1);	
	}//for  (int  ii = 0;  ii < v_genes.size(); ii++)


	int  i_child_func_bit_value;
	for  (int  ii = 0;  ii < v_child_funcs_offsets.size(); ii++)
	{
		i_child_func_bit_value =  pc_parent->v_all_hiff_nodes.at(v_child_funcs_offsets.at(ii)).iGetBitValue(iShift, piArgument,  iArgLen);

		if  (i_child_func_bit_value  ==  -1)  return(-1);
		if  (i_child_func_bit_value  ==  0)  b_zeroes = true;
		if  (i_child_func_bit_value  ==  1)  b_ones = true;

		if  ((b_zeroes == true)&&(b_ones == true))  return(-1);		
	}//for  (int  ii = 0;  ii < v_child_funcs.size(); ii++)

	if  ((b_zeroes == true)&&(b_ones == true))  return(-1);	
	if  (b_zeroes == true)  return(0);
	return(1);
}//int  CHiffNode::iGetBitValue(int iShift, int  *piArgument,  int  iArgLen,  double  *pdFuncVal)




int  CHiffNode::iGetUnitation(int iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;
	bool  b_zeroes, b_ones;

	i_unitation = 0;


	int  i_gene_offset;

	i_gene_offset = 0;
	for  (int  ii = 0;  ii < v_genes.size(); ii++)
	{
		i_gene_offset = iShift + v_genes.at(ii);
		if  (i_gene_offset  <  0)  return(-1);
		if  (i_gene_offset  >=  iArgLen)  return(-1);

		if  (piArgument[i_gene_offset]  ==  1)  i_unitation++;
	}//for  (int  ii = 0;  ii < v_genes.size(); ii++)


	int  i_child_func_bit_value;
	for  (int  ii = 0;  ii < v_child_funcs_offsets.size(); ii++)
	{
		i_child_func_bit_value =  pc_parent->v_all_hiff_nodes.at(v_child_funcs_offsets.at(ii)).iGetBitValue(iShift, piArgument,  iArgLen);
		if  (i_child_func_bit_value  ==  1)  i_unitation++;
		if  (i_child_func_bit_value  ==  -1)  return(-1);	
	}//for  (int  ii = 0;  ii < v_child_funcs.size(); ii++)

	return(i_unitation);
}//int  CHiffNode::iGetUnitation(int iShift, int  *piArgument,  int  iArgLen,  double  *pdFuncVal)



int  CHiffNode::iGetBitLength()
{
	if  (pc_unitation_problem  !=  NULL)  return(pc_unitation_problem->iGetBitLength());

	return(0);
}//int  CHiffNode::iGetBitLength()



int  CHiffNode::iGetBitLengthTotal()
{
	int  i_result;

	i_result = 0;

	for  (int  ii = 0;  ii < v_child_funcs_offsets.size(); ii++)
		i_result += pc_parent->v_all_hiff_nodes.at(v_child_funcs_offsets.at(ii)).iGetBitLengthTotal();
	i_result += v_genes.size();

	return(i_result);
}//int  CHiffNode::iGetBitLengthTotal()


//---------------------------class  CDecFunc1Modal--------------------------------------

CDecFunc1Modal::CDecFunc1Modal()  :  CUnitationBased()
{

}//CDecFunc1Modal::CDecFunc1Modal()


CDecFunc1Modal::CDecFunc1Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr)  :  CUnitationBased(iBitLength,  dMaxVal,  dMinVal,  dLowerThanMaxScale, pcErr)
{
	*pcErr  =  eConfigure(iBitLength,  dMaxVal,  dMinVal, dLowerThanMaxScale);
}//CDecFunc1Modal::CDecFunc1Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  CError  *pcErr)


CDecFunc1Modal::~CDecFunc1Modal()
{

}//CDecFunc1Modal::~CDecFunc1Modal()



CError  CDecFunc1Modal::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale)
{
	CError  c_err;

	c_err  =  CUnitationBased::eConfigure(iBitLength,  dMaxVal,  dMinVal, dLowerThanMaxScale);

	d_min_max_dist_grid  =  (d_max_val  -  d_min_val)/(i_bit_length);
	d_min_max_dist_grid  =  d_min_max_dist_grid/dLowerThanMaxScale;

	return(c_err);
}//CError  CDecFunc2Modal::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal)




CError  CDecFunc1Modal::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	CString  s_setting_name,  s_buf;


	i_bit_length  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	
	d_max_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_min_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MIN_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MIN_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_lower_than_max_scaling  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)



	c_err  =  eConfigure(i_bit_length,  d_max_val,  d_min_val, d_lower_than_max_scaling);
		

	return(c_err);
}//CError  CDecFunc1Modal::eLoadSettings(FILE  *pfSettings)



double  CDecFunc1Modal::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument,  iArgLen);
	if  (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFunc1Modal::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)



double  CDecFunc1Modal::dGetFuncValue(int  iUnitation)
{
	int  i_mirror_uni;


	i_mirror_uni  =  i_bit_length  -  iUnitation;

	if  (i_mirror_uni  ==  0)  return(d_max_val);

	double  d_buf;
	d_buf  =  i_mirror_uni  -  1;

	double  d_func_val;
	d_func_val  =  d_min_val  +  d_buf  *  d_min_max_dist_grid;
	
	return(d_func_val);
}//double  CDecFunc1Modal::dGetFuncValue(int  iUnitation)



CCompProblem* CDecFunc1Modal::pcClone()
{
	CError  c_err;
	CDecFunc1Modal  *pc_result;

	pc_result = new CDecFunc1Modal(1, 1, 0, 1, &c_err);

	*pc_result = *this;

	//Tools::vShow(d_max_val);
	//Tools::vShow(pc_result->d_max_val);

	return(pc_result);
}//CCompProblem* CDecFunc1Modal::pcClone()





//---------------------------class  CDecFunc2Modal--------------------------------------

CDecFunc2Modal::CDecFunc2Modal()  :  CUnitationBased()
{

}//CDecFunc2Modal::CDecFunc2Modal()


CDecFunc2Modal::CDecFunc2Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale,  CError  *pcErr)  :  CUnitationBased(iBitLength,  dMaxVal,  dMinVal,  dLowerThanMaxScale,  pcErr)
{
	*pcErr  =  eConfigure(iBitLength,  dMaxVal,  dMinVal,  dLowerThanMaxScale);
}//CDecFunc2Modal::CDecFunc2Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  CError  *pcErr)


CDecFunc2Modal::~CDecFunc2Modal()
{

}//CDecFunc2Modal::~CDecFunc2Modal()


CError  CDecFunc2Modal::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale)
{
	CError  c_err;

	c_err  =  CUnitationBased::eConfigure(iBitLength,  dMaxVal,  dMinVal,  dLowerThanMaxScale);
	if  (c_err)  return(c_err);


	d_min_max_dist_grid  =  (d_max_val  -  d_min_val)/(i_bit_length/2);
	d_min_max_dist_grid  =  d_min_max_dist_grid/dLowerThanMaxScale;

	return(c_err);
}//CError  CDecFunc2Modal::eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal)



CError  CDecFunc2Modal::eLoadSettings(FILE  *pfSettings)
{
	CError  c_err;
	CString  s_setting_name,  s_buf;


	
	i_bit_length  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	
	d_max_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_min_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_MIN_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MIN_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_lower_than_max_scaling  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)



	c_err  =  eConfigure(i_bit_length,  d_max_val,  d_min_val, d_lower_than_max_scaling);
		

	return(c_err);
}//CError  CDecFunc2Modal::eLoadSettings(FILE  *pfSettings)




double  CDecFunc2Modal::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	int  i_unitation;

	i_unitation = i_get_unitation(iShift, piArgument,  iArgLen);
	if  (i_unitation < 0)  return(-1);

	return(dGetFuncValue(i_unitation));
}//double  CDecFunc2Modal::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)




double  CDecFunc2Modal::dGetFuncValue(int  iUnitation)
{
	int  i_mirror_uni;
	int  i_bit_length_div_2;

	i_bit_length_div_2  =  i_bit_length / 2;

	if  (iUnitation >  i_bit_length_div_2)
		i_mirror_uni  =  i_bit_length_div_2  -  (iUnitation  -  i_bit_length_div_2);
	else
		i_mirror_uni  =  iUnitation;

	if  (i_mirror_uni  ==  0)  return(d_max_val);

	double  d_buf;
	d_buf  =  i_mirror_uni  -  1;

	double  d_func_val;
	d_func_val  =  d_min_val  +  d_buf  *  d_min_max_dist_grid;

	return(d_func_val);

}//double  CDecFunc2Modal::dGetFuncValue(int  iUnitation)





CCompProblem* CDecFunc2Modal::pcClone()
{
	CError  c_err;
	CDecFunc2Modal  *pc_result;

	pc_result = new CDecFunc2Modal(1, 1, 0, 1, &c_err);

	*pc_result = *this;

	return(pc_result);
}//CCompProblem* CDecFunc1Modal::pcClone()



//---------------------------class  CKnapsackProblem--------------------------------------


CKnapsackProblem::CKnapsackProblem() : CConstraintBased()
{

}//CKnapsackProblem::CKnapsackProblem()


CKnapsackProblem::CKnapsackProblem(int iBitLength, double dMaxValue, double dCapacity, CError *pcErr) : CConstraintBased()
{
	*pcErr = eConfigure(iBitLength, dMaxValue, dCapacity);
}//CKnapsackProblem::CKnapsackProblem(int iBitLength, double dMaxValue, double dCapacity, CError *pcErr)


CError CKnapsackProblem::eConfigure(int iBitLength, double dMaxValue, double dCapacity)
{
	CError c_err;

	c_err = CConstraintBased::eConfigure(iBitLength, dMaxValue, 0, 0);
	if (c_err) return(c_err);

	if (dCapacity < 0)
	{
		c_err.vSetError("Capacity below 0");
		return(c_err);
	}//if (dCapacity < 0)

	v_constraints.clear();
	v_constraints.push_back(dCapacity);

	v_values.clear();
	v_values.resize(iBitLength, 0);

	v_weights.clear();
	v_weights.resize(iBitLength, 0);

	v_sorted_elements_indexes_by_profit.clear();
	v_sorted_elements_indexes_by_profit.resize(iBitLength, 0);

	for (int i = 0; i < v_sorted_elements_indexes_by_profit.size(); i++)
	{
		v_sorted_elements_indexes_by_profit.at(i) = i;
	}//for (int i = 0; i < v_sorted_elements_indexes_by_profit.size(); i++)

	return(c_err);
}//CError CKnapsackProblem::eConfigure(int iBitLength, double dMaxValue, double dCapacity)


CError CKnapsackProblem::eAddValue(int iValue, double dValueValue)
{
	CError c_err;

	if (iValue < 0)
	{
		CString  s_buf;

		s_buf.Format("Value index below zero: %d", iValue);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iValue < 0)

	if (dValueValue <= 0)
	{
		CString  s_buf;

		s_buf.Format("Value must be greater than zero: %f", dValueValue);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (dValueValue <= 0)

	if (iValue >= i_bit_length)
	{
		CString  s_buf;

		s_buf.Format("Value index too high: %d >= %d", iValue, i_bit_length);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iValue >= i_bit_length)

	v_values.at(iValue) = dValueValue;

	v_sort_elements_indexes_by_profit();

	return(c_err);
}//CError CKnapsackProblem::eAddValue(int iValue, double dValueValue)


CError CKnapsackProblem::eAddWeight(int iWeight, double dWeightValue)
{
	CError c_err;

	if (iWeight < 0)
	{
		CString  s_buf;

		s_buf.Format("Weight index below zero: %d", iWeight);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iWeight < 0)

	if (dWeightValue <= 0)
	{
		CString  s_buf;

		s_buf.Format("Weight must be greater than zero: %f", dWeightValue);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (dValueValue <= 0)

	if (iWeight >= i_bit_length)
	{
		CString  s_buf;

		s_buf.Format("Weight index too high: %d >= %d", iWeight, i_bit_length);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iWeight >= i_bit_length)

	v_weights.at(iWeight) = dWeightValue;

	v_sort_elements_indexes_by_profit();

	return(c_err);
}//CError CKnapsackProblem::eAddWeight(int iWeight, double dWeightValue)


double CKnapsackProblem::dGetFuncValue(int iShift, int *piArgument, int iArgLen)
{
	if (iShift < 0) return(-1);
	if (iShift + iGetBitLength() > iArgLen) return(-1);
	
	int *pi_copied_argument = new int[i_bit_length];

	int i_cur_pos;

	for (int i_pos = 0; i_pos < i_bit_length; i_pos++)
	{
		i_cur_pos = i_pos + iShift;
		pi_copied_argument[i_pos] = piArgument[i_cur_pos];
	}//for (int i_pos = 0; i_pos < i_bit_length; i_pos++)

	bRepairGreedy(0, pi_copied_argument, i_bit_length);

	//if (!bSatisfyConstraints(iShift, piArgument, iArgLen)) return(0);

	double d_value = 0;

	for (int i_pos = 0; i_pos < i_bit_length; i_pos++)
	{
		if (pi_copied_argument[i_pos] == 1) d_value += v_values.at(i_pos);
	}//for (int i_pos = 0; i_pos < i_bit_length; i_pos++)

	delete pi_copied_argument;

	return(d_value);
}//double CKnapsackProblem::dGetFuncValue(int iShift, int *piArgument, int iArgLen)



bool CKnapsackProblem::bRepairGreedy(int  iShift, int  *piArgument,  int  iArgLen)
{
	if (iShift < 0) return(false);
	if (iShift + iGetBitLength() > iArgLen) return(false);

	bool b_repaired = false;

	double pd_constraints_values[1];

	bGetConstraintsValues(iShift, piArgument, iArgLen, pd_constraints_values);

	double d_weight = pd_constraints_values[0];
	double d_capacity = v_constraints.at(0);

	int i_pos, i_cur_pos, *pi_val;

	for (int i = 0; d_weight > d_capacity; i++)
	{
		i_pos = v_sorted_elements_indexes_by_profit.at(i);
		i_cur_pos = i_pos + iShift;

		pi_val = piArgument + i_cur_pos;

		if (*pi_val == 1)
		{
			*pi_val = 0;
			d_weight -= v_weights.at(i_pos);
			b_repaired = true;
		}//if (*pi_val == 1)
	}//for (int i = 0; d_weight > d_capacity; i++)

	return b_repaired;
}//bool CKnapsackProblem::bRepairGreedy(int  iShift, int  *piArgument,  int  iArgLen)


bool CKnapsackProblem::bGetConstraintsValues(int  iShift, int  *piArgument,  int  iArgLen,  double  *pdValues)
{
	if (iShift < 0) return(false);
	if (iShift + iGetBitLength() > iArgLen) return(false);

	int i_cur_pos;

	double d_weight = 0;

	for (int i_pos = 0; i_pos < i_bit_length; i_pos++)
	{
		i_cur_pos = i_pos + iShift;
		if (piArgument[i_cur_pos] == 1) d_weight += v_weights.at(i_pos);
	}//for (int i_pos = 0; i_pos < i_bit_length; i_pos++)

	*(pdValues + 0) = d_weight;

	return(true);
}//void CKnapsackProblem::bGetConstraintsValues(int  iShift, int  *piArgument,  int  iArgLen,  double  *pdValues)



CError CKnapsackProblem::eLoadSettings(FILE *pfSettings)
{
	CError c_err;
	CString s_buf;
	CString s_values_for_value, s_values_for_weight, s_setting_name;

	int i_index;
	double d_value;

	i_bit_length = Tools::iReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != COMP_PROBLEM_PARAM_BIT_LENGTH)

	d_max_val = Tools::dReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != COMP_PROBLEM_PARAM_MAX_VAL)

	double d_capacity = Tools::dReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != KNAPSACK_PROBLEM_PARAM_CAPACITY)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", KNAPSACK_PROBLEM_PARAM_CAPACITY, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != KNAPSACK_PROBLEM_PARAM_CAPACITY)

	c_err = eConfigure(i_bit_length, d_max_val, d_capacity);
	if (c_err) return(c_err);

	s_values_for_value = Tools::sReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != KNAPSACK_PROBLEM_PARAM_VALUES)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", KNAPSACK_PROBLEM_PARAM_VALUES, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != KNAPSACK_PROBLEM_PARAM_VALUES)

	i_index = 0;
	int i_value = 0;
	while (i_index < s_values_for_value.GetLength())
	{
		d_value = Tools::dExtractFromString(s_values_for_value, i_index, &i_index);

		c_err = eAddValue(i_value, d_value);
		if (c_err) return(c_err);
		i_value++;
	}//while (i_index < s_values_for_value.GetLength())

	if (i_value < i_bit_length)
	{
		s_buf.Format("Not enough values defined. Elements: %d  Defined: %d", i_bit_length, i_value);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (i_value < i_bit_length)

	s_values_for_weight = Tools::sReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != KNAPSACK_PROBLEM_PARAM_WEIGHTS)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", KNAPSACK_PROBLEM_PARAM_WEIGHTS, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != KNAPSACK_PROBLEM_PARAM_WEIGHTS)

	i_index = 0;
	int i_weight = 0;
	while (i_index < s_values_for_weight.GetLength())
	{
		d_value = Tools::dExtractFromString(s_values_for_weight, i_index, &i_index);

		c_err = eAddWeight(i_weight, d_value);
		if (c_err) return(c_err);
		i_weight++;
	}//while (i_index < s_values_for_weight.GetLength())

	if (i_value < i_bit_length)
	{
		s_buf.Format("Not enough weights defined. Elements: %d  Defined: %d", i_bit_length, i_weight);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (i_value < i_bit_length)

	return(c_err);
}//CError CKnapsackProblem::eLoadSettings(FILE *pfSettings)



void  CKnapsackProblem::vSaveSettings(FILE  *pfSettings)
{
	CString s_buf;
	CString c_separator("\\\\");

	s_buf.Format("%d%s%s\n", i_bit_length, c_separator, COMP_PROBLEM_PARAM_BIT_LENGTH);
	fprintf(pfSettings, s_buf);

	s_buf.Format("%f%s%s\n", v_constraints.at(0), c_separator, KNAPSACK_PROBLEM_PARAM_CAPACITY);
	fprintf(pfSettings, s_buf);

	for (int i = 0; i < v_values.size(); i++)
	{
		if (i > 0) fprintf(pfSettings, " ");

		s_buf.Format("%f", v_values.at(i));
		fprintf(pfSettings, s_buf);
	}//for (int i = 0; i < v_values.size(); i++)

	s_buf.Format("%s%s\n", c_separator, KNAPSACK_PROBLEM_PARAM_VALUES);
	fprintf(pfSettings, s_buf);

	for (int i = 0; i < v_weights.size(); i++)
	{
		if (i > 0) fprintf(pfSettings, " ");

		s_buf.Format("%f", v_weights.at(i));
		fprintf(pfSettings, s_buf);
	}//for (int i = 0; i < v_weights.size(); i++)

	s_buf.Format("%s%s\n", c_separator, KNAPSACK_PROBLEM_PARAM_WEIGHTS);
	fprintf(pfSettings, s_buf);
}//void  CKnapsackProblem::vSaveSettings(FILE  *pfSettings)



CCompProblem* CKnapsackProblem::pcClone()
{
	CError  c_err;
	CKnapsackProblem  *pc_result;

	pc_result = new CKnapsackProblem(1, 0, 1, &c_err);
	(*pc_result) = (*this);

	return(pc_result);
}//CCompProblem* CKnapsackProblem::pcClone()



void CKnapsackProblem::v_sort_elements_indexes_by_profit()
{
	function<bool(int, int)> f_compare = [&](int iElementIndex0, int iElementIndex1) -> bool
	{
		double d_profit_0 = DBL_MAX;
		double d_profit_1 = DBL_MAX;

		if (v_weights.at(iElementIndex0) > 0) d_profit_0 = v_values.at(iElementIndex0) / v_weights.at(iElementIndex0);
		if (v_weights.at(iElementIndex1) > 0) d_profit_1 = v_values.at(iElementIndex1) / v_weights.at(iElementIndex1);

		return d_profit_0 < d_profit_1;
	};//function<bool(int, int)> f_compare = [&](int iElementIndex0, int iElementIndex1) -> bool

	sort(v_sorted_elements_indexes_by_profit.begin(), v_sorted_elements_indexes_by_profit.end(), f_compare);
}//void CKnapsackProblem::v_sort_elements_indexes_by_profit()



//---------------------------class  CMaxSatProblem--------------------------------------
CMaxSatProblem::CMaxSatProblem() : CCompProblem()
{

}//CMaxSatProblem::CMaxSatProblem() : CCompProblem()


CMaxSatProblem::CMaxSatProblem(int  iBitLength,  int  iClauseCount,  CError  *pcErr) : CCompProblem()
{
	*pcErr = eConfigure(iBitLength, iClauseCount);
}//CMaxSatProblem::CMaxSatProblem(int  iBitLength,  int  iClauseCount,  CError  *pcErr) : CCompProblem(iBitLength, 0, iBitLength, 0, pcErr)


CError CMaxSatProblem::eLoadSettings(FILE  *pfSettings)
{
	CError c_err;
	CString s_buf;
	CString s_values_for_variable, s_setting_name;

	int i_clause_count;
	int i_variable_count;
	int i_variable;
	int i_variable_value;
	int i_index;

	i_bit_length = Tools::iReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != COMP_PROBLEM_PARAM_BIT_LENGTH)

	i_clause_count = Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_CLAUSE_COUNT)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", MAX_SAT_PROBLEM_PARAM_CLAUSE_COUNT, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_CLAUSE_COUNT)

	c_err = eConfigure(i_bit_length, i_clause_count);
	if (c_err) return(c_err);

	for (int i = 0; i < i_clause_count; i++)
	{
		i_variable_count = Tools::iReadLine(pfSettings, &s_setting_name);
		if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_VARIABLE_COUNT)
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", MAX_SAT_PROBLEM_PARAM_VARIABLE_COUNT, s_setting_name);
			c_err.vSetError(s_buf);
			return(c_err);	
		}//if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_VARIABLE_COUNT)

		c_err = eAddClause(i, i_variable_count);
		if (c_err) return(c_err);

		s_values_for_variable = Tools::sReadLine(pfSettings, &s_setting_name);
		if (s_setting_name != MAX_SAT_PROBLEM_PARAM_VARIABLES)
		{
			s_buf.Format("Error at reading (%s) parameter read:(%s)", MAX_SAT_PROBLEM_PARAM_VARIABLES, s_setting_name);
			c_err.vSetError(s_buf);
			return(c_err);
		}//if (s_setting_name != MAX_SAT_PROBLEM_PARAM_VARIABLES)

		i_index = 0;
		i_variable = 0;
		while (i_index < s_values_for_variable.GetLength())
		{
			i_variable_value = Tools::iExtractFromString(s_values_for_variable, i_index, &i_index);

			c_err = eAddVariable(i, i_variable, i_variable_value);
			if (c_err) return(c_err);
			i_variable++;
		}//while (i_index < s_values_for_variable.GetLength())

		if (i_variable < i_variable_count)
		{
			s_buf.Format("Not enough variables defined. Elements: %d  Defined: %d", i_variable_count, i_variable);
			c_err.vSetError(s_buf);

			return(c_err);
		}//if (i_variable < i_variable_count)
	}//for (int i = 0; i < i_clause_count; i++)

	s_buf  =  Tools::sReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_CLAUSE_END)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", MAX_SAT_PROBLEM_PARAM_CLAUSE_END, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  MAX_SAT_PROBLEM_PARAM_CLAUSE_END)

	return(c_err);
}//CError CMaxSatProblem::eLoadSettings(FILE  *pfSettings)


void CMaxSatProblem::vSaveSettings(FILE  *pfSettings)
{
	CString s_buf;
	CString c_separator("\\\\");

	s_buf.Format("%d%s%s\n", i_bit_length, c_separator, COMP_PROBLEM_PARAM_BIT_LENGTH);
	fprintf(pfSettings, s_buf);

	s_buf.Format("%d%s%s\n", v_clauses.size(), c_separator, MAX_SAT_PROBLEM_PARAM_CLAUSE_COUNT);
	fprintf(pfSettings, s_buf);

	for (int i = 0; i < v_clauses.size(); i++)
	{
		s_buf.Format("%d%s%s\n", v_clauses.at(i).size(), c_separator, MAX_SAT_PROBLEM_PARAM_VARIABLE_COUNT);
		fprintf(pfSettings, s_buf);

		for (int ii = 0; ii < v_clauses.at(i).size(); ii++)
		{
			if (ii > 0) fprintf(pfSettings, " ");

			s_buf.Format("%d", v_clauses.at(i).at(ii));
			fprintf(pfSettings, s_buf);
		}//for (int ii = 0; ii < v_clauses.at(i).size(); ii++)

		s_buf.Format("%s%s\n", c_separator, MAX_SAT_PROBLEM_PARAM_VARIABLES);
		fprintf(pfSettings, s_buf);
	}//for (int i = 0; i < v_clauses.size(); i++)

	s_buf.Format("%s%s\n", c_separator, MAX_SAT_PROBLEM_PARAM_CLAUSE_END);
	fprintf(pfSettings, s_buf);
}//void CMaxSatProblem::vSaveSettings(FILE  *pfSettings)


CError CMaxSatProblem::eConfigure(int  iBitLength,  int  iClauseCount)
{
	CError c_err;

	c_err = CCompProblem::eConfigure(iBitLength, iBitLength, 0, 0);
	if (c_err) return(c_err);

	if  (iClauseCount  <  0)
	{
		c_err.vSetError("Clause count below 0");
		return(c_err);
	}//if  (iClauseCount  <  0)

	v_clauses.clear();
	v_clauses.resize(iClauseCount);

	d_max_val = iClauseCount;

	return(c_err);
}//CError CMaxSatProblem::eConfigure(int  iBitLength,  int  iClauseCount)


CError CMaxSatProblem::eAddClause(int  iClause, int  iVariableCount)
{
	CError c_err;

	if (iClause < 0)
	{
		CString  s_buf;

		s_buf.Format("Clause index below zero: %d", iClause);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iClause < 0)

	if (iClause >= v_clauses.size())
	{
		CString  s_buf;

		s_buf.Format("Clause index too high: %d >= %d", iClause, v_clauses.size());
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iClause >= v_clauses.size())

	if (iVariableCount <= 0)
	{
		CString  s_buf;

		s_buf.Format("Variable count must be greater than zero: %d", iVariableCount);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariableCount <= 0)

	v_clauses.at(iClause).clear();
	v_clauses.at(iClause).resize(iVariableCount, 0);

	return(c_err);
}//CError CMaxSatProblem::eAddClause(int  iClause, int  iVariableCount)


CError CMaxSatProblem::eAddVariable(int  iClause, int  iVariable, int  iVariableValue)
{
	CError c_err;

	if (iClause < 0)
	{
		CString  s_buf;

		s_buf.Format("Clause index below zero: %d", iClause);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iClause < 0)

	if (iClause >= v_clauses.size())
	{
		CString  s_buf;

		s_buf.Format("Clause index too high: %d >= %d", iClause, v_clauses.size());
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iClause >= v_clauses.size())

	if (iVariable < 0)
	{
		CString  s_buf;

		s_buf.Format("Variable index below zero: %d", iVariable);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariable < 0)

	if (iVariable > v_clauses.at(iClause).size())
	{
		CString  s_buf;

		s_buf.Format("Variable index too high: %d > %d", iClause, v_clauses.at(iClause).size());
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariable >= v_clauses.at(iClause).size())

	if (iVariableValue == 0)
	{
		CString  s_buf;

		s_buf.Format("Variable value equals zero: %d", iVariableValue);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariableValue == 0)

	if (iVariableValue < -i_bit_length)
	{
		CString  s_buf;

		s_buf.Format("Variable value too low: %d < %d", iVariableValue, -i_bit_length);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariableValue < -i_bit_length)

	if (iVariableValue > i_bit_length)
	{
		CString  s_buf;

		s_buf.Format("Variable value too high: %d > %d", iVariableValue, i_bit_length);
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (iVariableValue > i_bit_length)

	v_clauses.at(iClause).at(iVariable) = iVariableValue;

	return(c_err);
}//CError CMaxSatProblem::eAddVariable(int  iClause, int  iVariable, int  iVariableValue)


double  CMaxSatProblem::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	if (iShift < 0) return(-1);
	if (iShift + iGetBitLength() > iArgLen) return(-1);

	int i_cur_pos;
	bool b_negation;
	int i_desired_value;
	bool b_satisfied;

	double d_value = 0;

	for (size_t i = 0; i < v_clauses.size(); i++)
	{
		if (v_clauses.at(i).size() == 0) return (-1);

		b_satisfied = false;

		for (size_t ii = 0; ii < v_clauses.at(i).size() && !b_satisfied; ii++)
		{
			i_cur_pos = v_clauses.at(i).at(ii);

			if (i_cur_pos == 0) return (-1);

			b_negation = i_cur_pos < 0;
			
			if (b_negation) i_cur_pos = -i_cur_pos;
			i_cur_pos--;
			i_cur_pos += iShift;

			i_desired_value = b_negation ? 0 : 1;

			if (*(piArgument + i_cur_pos) == i_desired_value)
			{
				b_satisfied = true;
				d_value++;
			}//if (*(piArgument + i_cur_pos) == i_desired_value)
		}//for (size_t ii = 0; ii < v_clauses.at(i).size() && !b_satisfied; ii++)
	}//for (size_t i = 0; i < v_clauses.size(); i++)

	return(d_value);
}//double  CMaxSatProblem::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)


CCompProblem* CMaxSatProblem::pcClone()
{
	CError  c_err;
	CMaxSatProblem  *pc_result;

	pc_result = new CMaxSatProblem(1, 1, &c_err);
	(*pc_result) = (*this);

	return(pc_result);
}//CCompProblem* CMaxSatProblem::pcClone()



//---------------------------class  CReportGrid--------------------------------------

/*bool operator<(const CReportGrip*  const&pc_rep_g_1, const CReportGrip* const&pc_rep_g_2)
{
	return pc_rep_g_1->d_value < pc_rep_g_2->d_value;
}//bool operator<(const CReportGrid& pc_rep_g_1, const CReportGrid& pc_rep_g_2)


bool operator==(const CReportGrid* const&pc_rep_g_1, const CReportGrid* const&pc_rep_g_2)
{
	
	if  (pc_rep_g_1->d_value  !=  pc_rep_g_2->d_value)  return(false);
	if  (pc_rep_g_1->i_bit_len  !=  pc_rep_g_2->i_bit_len)  return(false);
	for  (int  ii = 0; ii < pc_rep_g_1->i_bit_len;  ii++)
		if  (pc_rep_g_1->pi_argument[ii]  !=  pc_rep_g_2->pi_argument[ii])  return(false);

	return(true);
}//bool operator==(const CReportGrid& pc_rep_g_1, const CReportGrid& pc_rep_g_2)*/


CReportGrid::CReportGrid()
{
	i_bit_len  =  0;
	d_value  =  0;

	pi_argument  =  NULL;
}//CReportGrid::CReportGrid()


CReportGrid::CReportGrid(int  iBitLength,  int  *piArgument,  double  dVal)
{
	i_bit_len  =  0;
	d_value  =  0;

	pi_argument  =  NULL;
	vSetData(iBitLength,  piArgument,  dVal);
}//CReportGrid::CReportGrid(int  iBitLength,  int  *piArgument,  double  dVal)


CReportGrid::~CReportGrid()
{
	if  (pi_argument  !=  NULL)  delete  []  pi_argument;
}//CReportGrid::~CReportGrid()



void  CReportGrid::vSetData(int  iBitLength,  int  *piArgument,  double  dVal)
{
	i_bit_len  =  iBitLength;
	d_value  =  dVal;

	if  (pi_argument  !=  NULL)  delete  []  pi_argument;
	pi_argument  =  NULL;

	pi_argument  =  new  int[i_bit_len];

	for  (int  ii = 0;  ii < i_bit_len;  ii++)
		pi_argument[ii]  =  piArgument[ii];


};//CReportGrid::vSetData(int  iBitLength,  int  *piArgument,  double  dVal)


void  CReportGrid::vPrint(FILE  *pfDest)
{
	for  (int  ii = 0;  ii < i_bit_len;  ii++)
		fprintf(pfDest, "%d", pi_argument[ii]);

	fprintf(pfDest, "\t");

	fprintf(pfDest, "%.8lf",  d_value);
};//void  CReportGrid::vPrint(FILE  *pfDest)











