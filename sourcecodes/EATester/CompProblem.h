#pragma once


//tools
#include  "atlstr.h"  //CString
#include  "error.h"
#include  <iostream>
#include  <vector>
#include  <functional>
#include  <algorithm>


//..and other tools...
#include  "util\tools.h"


using namespace std;


namespace ProblemTools
{
	class  CReportGrid;//predefinition
	class  CUnitationBased;//predefinition
	class  CHiffProblem;//predefinition
	class  CCompProblem;//predefinition

	#define  COMP_PROBLEM_PARAM_BIT_LENGTH  "bit length"
	#define  COMP_PROBLEM_PARAM_MAX_VAL	"max val"
	#define  COMP_PROBLEM_PARAM_MIN_VAL	"min val"
	#define  COMP_PROBLEM_PARAM_LOWER_THAN_MAX_SCALE	"scaling for lower than max values"
	#define  COMP_PROBLEM_PARAM_BIT_NUMBER_MULTIPLIER  "bit number multiplier"
	#define  COMP_PROBLEM_PARAM_UNITATION_NOISE_STEP  "unitation noise step"
	//#define  COMP_PROBLEM_PARAM_UNITATION_NOISE_SIZE  "unitation noise size"

	#define  HIFF_PROBLEM_PARAM_TREE_DEPTH	"HIFF tree depth"
	#define  HIFF_PROBLEM_PARAM_TREE_END	"HIFF tree END"
	#define  HIFF_PROBLEM_PARAM_MULTIPLIER	"HIFF multiplier"

	#define  KNAPSACK_PROBLEM_PARAM_CAPACITY  "capacity"
	#define  KNAPSACK_PROBLEM_PARAM_VALUES    "values"
	#define  KNAPSACK_PROBLEM_PARAM_WEIGHTS   "weights"

	#define  MAX_SAT_PROBLEM_PARAM_CLAUSE_COUNT    "clause count"
	#define  MAX_SAT_PROBLEM_PARAM_CLAUSE_END      "clause end"
	#define  MAX_SAT_PROBLEM_PARAM_VARIABLE_COUNT  "variable count"
	#define  MAX_SAT_PROBLEM_PARAM_VARIABLES	   "variables"

	#define  COMP_PROBLEM_DEF_END	  "PROBLEM DEFINITION END"
	#define  COMP_PROBLEM_MULTIPLIER  "func multiplier"
	#define  COMP_PROBLEM_TYPE        "function type"
	#define  DEC_FUNC_1_MODAL_NAME	  "1_MODAL"
	#define  DEC_FUNC_2_MODAL_NAME	  "2_MODAL"
	#define  DEC_FUNC_EASY_NAME		  "ONE_MAX"
	#define  DEC_FUNC_ANY_NAME		  "ANY_PROBLEM"
	#define  DEC_FUNC_ANY_SCALED_NAME "ANY_SCALED_PROBLEM"
	#define  DEC_FUNC_ANY_SCALED_NOISED_NAME "ANY_SCALED_NOISED_PROBLEM"
	#define  HIFF_PROBLEM_NAME		  "HIFF_PROBLEM"
	#define  KNAPSACK_PROBLEM_NAME    "KNAPSACK_PROBLEM"
	#define  MAX_SAT_PROBLEM_NAME	  "MAX_SAT_PROBLEM"



	CError  eLoadSingleUnitationProblem(FILE *pfSource, CUnitationBased **ppcUnitationProblem);
	CError  eLoadSingleCompProblem(FILE *pfSource, CCompProblem **ppcCompProblem);



	class  CConcatCompProblem
	{
	public:
		CConcatCompProblem();
		~CConcatCompProblem();

		CError  eGetFuncValue(int  *piArgument,  int  iArgLen,  double  *pdFuncVal);
		CError  eAddCompProblem(CCompProblem  *pcCompProblem);

		int  iGetProblemBitLength();
		CError  eLoad(FILE *pfSource);
		void  operator=(CConcatCompProblem  &cOther);

		CError  eCreateReportSeparateFuncSummary(FILE  *pfReport);
		CError  eCreateReportSeparateFunc(FILE  *pfReport);
		CError  eCreateReport(FILE  *pfReport);
		static  double  dGetFuncEval()  {return(d_func_evaluations);};
		static  void  vZeroEval()  {d_func_evaluations = 0;};

		double  dGetMaxFuncVal();
	
	private:
		CError  e_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  iBitLength,  int  iCurrentBit);

		CCompProblem  *pc_clone_comp_problem(CCompProblem  *pcCompProblem);

		vector  <CCompProblem  *>  v_COMP_PROBLEMs;	

		static  double  d_func_evaluations;
	};//class  CConcatCompProblem

	

	#define  COMP_PROBLEM_UNKNOWN  0
	class  CCompProblem
	{
	public:
		CCompProblem();
		CCompProblem(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr);
		virtual ~CCompProblem();

		virtual  CError  eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale);

		virtual  double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen) = 0;
		
		virtual  CError  eCreateReport(FILE  *pfReport);
		virtual  CError  eCreateReportSummary(FILE  *pfReport);
		virtual  CError  eCreateReportUnitation(FILE  *pfReport);

		double  dGetMaxFuncVal()  {return(d_max_val);};
		int  iGetBitLength()  {return(i_bit_length);};

		virtual  CError  eLoadSettings(FILE *pfSettings) = 0;
		virtual  CCompProblem *pcClone() = 0;

		virtual  int  iGetFuncType()  {return(COMP_PROBLEM_UNKNOWN);};

	protected:
		int  i_bit_length;
		double  d_max_val;
		double  d_min_val;
		double  d_min_max_dist_grid;//(d_max_val  -  d_min_val)/(i_bit_length/2)
		double  d_lower_than_max_scaling;//d_min_max_dist_grid =  d_min_max_dist_grid/d_lower_than_max_scaling

		virtual  CError  e_create_report(vector  <CReportGrid *>  *pvReportGrids, int  *piArgument, int  iCurrentBit, bool  bLinear = false, int iCurrentValForLinear = 0);
	};//class  CCompProblem



	#define  COMP_PROBLEM_UNITATION_BASED_UNKNOWN	1
	class  CUnitationBased : public CCompProblem
	{
	public:
		CUnitationBased();
		CUnitationBased(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr);

		virtual  double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		virtual  double  dGetFuncValue(int  iUnitation)  =  0;

		virtual  int  iGetFuncType()  {return(COMP_PROBLEM_UNITATION_BASED_UNKNOWN);};

	protected:
		int  i_get_unitation(int  iShift, int  *piArgument,  int  iArgLen);
	};//class  CUnitationBased



	#define  COMP_PROBLEM_CONSTRAINT_BASED_UNKNOWN 2
	class CConstraintBased : public CCompProblem
	{
	public:
		CConstraintBased();
		CConstraintBased(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr);

		int  iGetConstraintCount() {return(v_constraints.size());};

		double *pdGetConstraints() {return(&v_constraints.at(0));};
		virtual bool bGetConstraintsValues(int  iShift, int  *piArgument,  int  iArgLen,  double  *pdValues) = 0;

		bool bSatisfyConstraints(int  iShift, int  *piArgument,  int  iArgLen);

		virtual  int  iGetFuncType()  {return(COMP_PROBLEM_CONSTRAINT_BASED_UNKNOWN);};

	protected:
		vector<double> v_constraints;
	};//class CConstraintBased



	#define  DEC_FUNC_1_MODAL	3
	class  CDecFunc1Modal  :  public  CUnitationBased
	{
	public:
		CDecFunc1Modal();
		CDecFunc1Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr);
		~CDecFunc1Modal();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale);

		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType()  {return(DEC_FUNC_1_MODAL);};

	private:
		
	};//class  CDecFunc1Modal



	


	#define  DEC_FUNC_2_MODAL	4
	class  CDecFunc2Modal  :  public  CUnitationBased
	{
	public:
		CDecFunc2Modal();
		CDecFunc2Modal(int  iBitLength,  double  dMaxVal,  double  dMinVal,  double  dLowerThanMaxScale, CError  *pcErr);
		~CDecFunc2Modal();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength,  double  dMaxVal,  double  dLowerThanMaxScale, double  dMinVal);

		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType()  {return(DEC_FUNC_2_MODAL);};
	private:
	};//class  CDecFunc2Modal





	#define  DEC_FUNC_EASY	5
	class  CDecFuncEasy  :  public  CUnitationBased
	{
	public:
		CDecFuncEasy();
		CDecFuncEasy(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale, CError  *pcErr);
		~CDecFuncEasy();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale);

		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType()  {return(DEC_FUNC_EASY);};

	private:
		
	};//class  CDecFuncEasy



	#define  DEC_FUNC_ANY	6
	class  CDecFuncAny  :  public  CUnitationBased
	{
	public:
		CDecFuncAny();
		CDecFuncAny(int  iBitLength,  CError  *pcErr);
		~CDecFuncAny();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength,  double  dMaxVal,  double  dMinVal, double  dLowerThanMaxScale);
		CError  eAddUnitationLevel(int  iUnitation, double dUnitValue);

		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType()  {return(DEC_FUNC_ANY);};

	protected:

		vector<double>  v_unitation_values;
		
	};//class  CDecFuncAny




	class  CHiffNode
	{
	friend class CHiffProblem;
	public:
		CHiffNode();
		CHiffNode(CHiffProblem  *pcParent);
		CHiffNode(const CHiffNode  &pcOther);
		~CHiffNode();

		double  dGetFuncValue(int iShift, int  *piArgument,  int  iArgLen);
		int  iGetBitValue(int iShift, int  *piArgument,  int  iArgLen);
		int  iGetUnitation(int iShift, int  *piArgument,  int  iArgLen);

		int  iGetBitLength();
		int  iGetBitLengthTotal();

	private:
		CHiffProblem  *pc_parent;
		CUnitationBased  *pc_unitation_problem;
		double  d_multiplier;

		vector<int>  v_child_funcs_offsets;
		vector<int>  v_genes;
	};//class  CHiffNode




	#define  HIFF_PROBLEM	7
	class  CHiffProblem  :  public  CCompProblem
	{
	friend class CHiffNode;
	public:
		CHiffProblem();
		~CHiffProblem();

		CError  eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale);
		CError  eLoadSettings(FILE  *pfSettings);
	
		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);

		CCompProblem *pcClone();
		int  iGetFuncType()  {return(HIFF_PROBLEM);};

	private:

		CHiffNode  *pc_root;
		int  i_demanded_tree_level;
		vector<CHiffNode>  v_all_hiff_nodes;//can not be deleted by root, because single node may be a child of many other nodes...
		
	};//class  CHiffProblem



	#define  KNAPSACK_PROBLEM  8
	class  CKnapsackProblem : public CConstraintBased
	{
	public:
		CKnapsackProblem();
		CKnapsackProblem(int  iBitLength,  double  dMaxValue,  double  dCapacity,  CError  *pcErr);

		virtual  CError  eLoadSettings(FILE  *pfSettings);
		virtual  void  vSaveSettings(FILE  *pfSettings);

		CError  eConfigure(int  iBitLength,  double  dMaxValue,  double  dCapacity);
		CError  eAddValue(int  iValue,  double  dValueValue);
		CError  eAddWeight(int  iWeight,  double  dWeightValue);

		virtual  double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);

		bool  bRepairGreedy(int  iShift, int  *piArgument,  int  iArgLen);

		double  *pdGetValues()  {return(&v_values.at(0));};
		double  *pdGetWeights()  {return(&v_weights.at(0));};

		virtual bool bGetConstraintsValues(int  iShift, int  *piArgument,  int  iArgLen,  double  *pdValues);

		virtual  CCompProblem  *pcClone();
		virtual  int  iGetFuncType()  {return(KNAPSACK_PROBLEM);};

	private:
		void  v_sort_elements_indexes_by_profit();

		vector<double>  v_values;
		vector<double>  v_weights;
		vector<int>  v_sorted_elements_indexes_by_profit;
	};//class  CKnapsackProblem



	#define  MAX_SAT_PROBLEM  9
	class CMaxSatProblem : public CCompProblem
	{
	public:
		CMaxSatProblem();
		CMaxSatProblem(int  iBitLength,  int  iClauseCount,  CError  *pcErr);

		virtual  CError  eLoadSettings(FILE  *pfSettings);
		virtual  void  vSaveSettings(FILE  *pfSettings);

		CError  eConfigure(int  iBitLength,  int  iClauseCount);
		CError  eAddClause(int  iClause, int  iVariableCount);
		CError  eAddVariable(int  iClause, int  iVariable, int  iVariableValue);

		virtual  double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);

		virtual  CCompProblem  *pcClone();
		virtual  int  iGetFuncType()  {return(MAX_SAT_PROBLEM);};

	private:
		vector<vector<int>> v_clauses;
	};//class CMaxSatProblem



	#define  DEC_FUNC_ANY_SCALED	10
	class  CDecFuncAnyScaled : public  CDecFuncAny
	{
	public:
		CDecFuncAnyScaled();
		CDecFuncAnyScaled(int  iBitLength, CError  *pcErr);
		~CDecFuncAnyScaled();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale);
		
		double  dGetFuncValue(int  iShift, int  *piArgument, int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType() { return(DEC_FUNC_ANY_SCALED); };

	protected:

		int  i_bit_number_multiplier;

	};//class  CDecFuncAnyScaled




	#define  DEC_FUNC_ANY_SCALED_NOISED		11
	class  CDecFuncAnyScaledNoised : public  CDecFuncAnyScaled
	{
	public:
		CDecFuncAnyScaledNoised();
		CDecFuncAnyScaledNoised(int  iBitLength, CError  *pcErr);
		~CDecFuncAnyScaledNoised();

		CError  eLoadSettings(FILE  *pfSettings);
		CError  eConfigure(int  iBitLength, double  dMaxVal, double  dMinVal, double  dLowerThanMaxScale);

		double  dGetFuncValue(int  iShift, int  *piArgument, int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		CCompProblem *pcClone();
		int  iGetFuncType() { return(DEC_FUNC_ANY_SCALED_NOISED); };

	private:

		//double d_unitation_noise_size;
		int  i_unitation_noise_step_in_bits;

		vector  <double>  v_modifications;

	};//class  CDecFuncAnyScaled



	class  CReportGrid
	{
	public:
		CReportGrid();
		CReportGrid(int  iBitLength,  int  *piArgument,  double  dVal);
		~CReportGrid();

		void  vSetData(int  iBitLength,  int  *piArgument,  double  dVal);
		void  vPrint(FILE  *pfDest);

	//private:
		int  i_bit_len;
		int  *pi_argument;
		double  d_value;
	};//class  CReportGrid


	struct SReportGridSort
	{
		bool operator()(const CReportGrid* const& p1, const CReportGrid* const& p2)
		{
			if(!p1)
			return true;
			if(!p2)
			return false;
			return p1->d_value < p2->d_value;
		}
	};//struct SReportGridSort
};//namespace ProblemTools


/*namespace  std
{
	struct less<DeceptiveTools::CReportGrid*>
	{
		bool operator()(DeceptiveTools::CReportGrid const* p1, DeceptiveTools::CReportGrid const* p2)
		{
			if(!p1)
			return true;
			if(!p2)
			return false;
			return p1->d_value < p2->d_value;
		}
	};

	struct greater<DeceptiveTools::CReportGrid*>
	{
		bool operator()(DeceptiveTools::CReportGrid const* p1, DeceptiveTools::CReportGrid const* p2)
		{
			if(!p1)
			return false;
			if(!p2)
			return true;
			return p1->d_value > p2->d_value;
		}
	};
};//namespace  std*/





