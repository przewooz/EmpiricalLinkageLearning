#pragma once


//tools
#include  "atlstr.h"  //CString
#include  "Error.h"
#include  <iostream>
#include  <vector>
#include  <functional>
#include  <algorithm>

//for randomizer...
//#using <mscorlib.dll>
//using namespace System;



//time counter
#include  "util\timer.h"
//..and other tools...
#include  "util\tools.h"

#include "CompProblem.h"


using namespace std;
using namespace TimeCounters;


#define  SOLUTION_PRW_DEMAND_MARKER	"demand"
#define  SOLUTION_UNICAST_DEMAND_MARKER	"demand"
#define  SOLUTION_ANYCAST_DEMAND_MARKER	"anycastdemand"

#define  DEMAND_ANYCAST_EXT	"dea"
#define  DEMAND_UNICAST_EXT	"dem"

#define  DEMAND_TYPE_UNICAST	1
#define  DEMAND_TYPE_ANYCAST	2
#define  DEMAND_TYPE_MULTICAST	3



namespace  FiberNets
{
	class  CFiberNode;//predefinition
	class  CFiberConnection;//predefinition
	class  CFiberDemand;//predefinition
	class  CFiberPath;//predefinition
	class  CFiberSolution;//predefinition
	class  CFiberSolutionChannel;//predefinition
	class  CSteinerTreeSet;//predefinition
	class  CSSJNetwork;//predefinition
	class  CSSJConn;//predefinition
	class  CSearchShortestPathNode;//predefinition


	class  CFiberNet
	{
	friend class CFiberNode;
	friend class CFiberConnection;
	friend class CFiberDemand;
	friend class CFiberPath;
	friend class CFiberSolution;
	public:
		CFiberNet();
		~CFiberNet();


		double  dCompCost();
		CError  eResetNet();
		bool  bCheckIfReplicaNode(CFiberNode  *pcNode);
		bool  bCheckIfReplicaNode(int  iNodeId);
		bool  bCompareNodes(int  iNodeId0, int iNodeId1);

		
		CError  eCloneTopologyTo(CFiberNet  *pcOtherNetwork);
		CError  eLoadTopology(CString  sNetFile);
		CError  eLoadReplicas(CString  sRepFile);
		CError  eLoadDemands(CString  sUnicastFile, CString  sAnycastFile, CString sMulticastFile, int  iPathsToUse, int  iTreesToUse);
		
		CError  eLoadPaths(CString  sPathFile);
		CError  eLoadPathsSpec(CString  sUniSpecFile);
		CError  eLoadPathsReg(CString  sUniRegFile);
		CError  eLoadModulation(CString  sModFile);

		CError  eLoadTrees(CString  sTreeFile);
		CError  eLoadTreesSpec(CString  sTreeSpecFile);
		CError  eLoadTreesReg(CString  sTreeRegFile);

		
		CError  eSaveTopology(CString  sNetFile);
		CError  eSaveDemands(CString  sDemFile);
		CError  eSavePaths(CString  sPathFile);
		CError  eSaveTrees(CString  sTreeFile);
		CError  eSaveModulation(CString  sModFile);
		CError  eSaveReplicas(CString  sRepFile);

		CError  eSavePathsSpec(CString  sPathSpecFile);
		CError  eSavePathsReg(CString  sPathRegFile);
		CError  eSaveTreesSpec(CString  sTreeSpecFile);
		CError  eSaveTreesReg(CString  sTreeRegFile);


		int  iGetDemandNumber()  {return(i_demand_num);};
		CFiberDemand *pcGetDemands()  {return(pc_demands);};
		bool  bConstraintBroken()  {return(false);};

		CFiberNode  *pcGetNodes()  {return(pc_nodes);}
		int  iGetNodesNumber()  {return(i_nodes_num);}
		int  iGetConnsNumber()  {return(i_conns_num);}

		CFiberConnection  *pcGetConn(int  iOffset);
		CFiberConnection  ***pcGetConnMat()  {return(pc_conn_mat);};


		
		

	protected:

		CFiberNode *pc_nodes;
		int  i_nodes_num;

		CFiberNode **pc_replices;
		int  i_rep_num;

		CFiberConnection *pc_conns;
		int  i_conns_num;
		
		CFiberConnection  ***pc_conn_mat;

		CFiberDemand  *pc_demands;
		int  i_demand_num;
		int  i_uni_demands;
		int  i_any_demands;
		int  i_multi_demands;

		CFiberPath  *pc_paths;
		int  i_path_num;

		vector<int>  **pv_paths_between_nodes;


		CFiberPath  *pc_trees;
		int  i_tree_num;

		//CFiberPath  **pc_tree_paths;

		CError  e_load_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum);
		CError  e_load_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum);

		CError  e_save_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum);
		CError  e_save_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum);

	

		CError  e_load_uni_demands(FILE  *pfSource, int iMaxPathsToUseBetweenNodes);
		CError  e_load_any_demands(FILE  *pfSource);
		CError  e_load_multi_demands(FILE  *pfSource, int iMaxTreesToUseForNode);
	};//class  CFiberNet



	class  CFiberNode
	{
	friend class CFiberNet;
	friend class CFiberDemand;
	friend class CSteinerTree;
	public:
		CFiberNode();
		~CFiberNode();

		void  vSetNodeNumber(int  iNodeNumber);

		void  vAddInConn(CFiberConnection  *pcConnection);
		void  vAddOutConn(CFiberConnection  *pcConnection);


	private:
		CFiberConnection **pc_out_conns;
		CFiberConnection **pc_in_conns;

		int  i_node_number;

	};//class  CFiberNode


	class  CFiberConnection
	{
	friend class CFiberNet;
	friend class CFiberNode;
	friend class CFiberPath;
	public:
		CFiberConnection();
		~CFiberConnection();

		int iGetId()  {return(i_connection_id);}

		int  iGetMaxSlot();
		void  vReset();

		CError  eReserveSlots(int  iFirstSlot, int iSlotsNumber, bool  bErrorControl = true);
		int  iGetLowestStartSlotForChannel(int  iFirstSlot, int iSlotsNumber);
		int  iGetNumberOfFreeSlotsBelow(double  dMax);

		int  iGetNodeStart()  {return(i_fiber_node_start_offset);}
		int  iGetNodeEnd()  {return(i_fiber_node_end_offset);}

		CError    eConfig(int  iConectionId, int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, CFiberNet *pcParent, double  dConnectionCapacity);
		void  vSaveStatePRW(FILE  *pfDest);
		void  vSaveState(FILE  *pfDest, int iMaxSlotNumber);
		bool  bSlotInUse(int iSlotNumber);

		double  dGetCapacity()  {return(d_connection_capacity);}

		
	protected:
		int  i_connection_id;
		int  i_fiber_node_start_offset;
		int  i_fiber_node_end_offset;
		CFiberNet *pc_parent;
		double  d_connection_capacity;

		vector<int>  v_slots_in_use;//vector of pairs: (start-end)
	};//class  CFiberConnection



	class  CFiberDemand
	{
	friend class CFiberNet;
	friend class CFiberSolution;
	friend class CFiberSolutionChannel;
	public:
		CFiberDemand();
		~CFiberDemand();	

		vector  <CFiberPath  *>  *pvGetPossiblePaths()  {return(&v_possible_paths);};

		CError  eConfig
			(
			int  iDemandId, 
			int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, double  dDemandCapacity, 
			CFiberNet *pcParent,
			int  iDemandType,
			CFiberDemand *pcAnycastTwinDemand, int  iMaxPathTreesAvailable = -1
			);
		void  vSave(FILE  *pfDest);

		void  vSummarizePathFitness();

		int  iDemandType()  {return(i_demand_type);}
		bool  bDominatingTwin();

		CFiberDemand  *pcGetTwin()  {return(pc_anycast_twin_dfemand);}

		int  iGetId()  {return(i_demand_id);}

		int  iGetStartNode()  {return(i_fiber_node_start_offset);}
		int  iGetEndNode()  {return(i_fiber_node_end_offset);}
		
	protected:
		int  i_demand_id;
		int  i_fiber_node_start_offset;
		int  i_fiber_node_end_offset;
		CFiberNet *pc_parent;
		double d_path_fit_sum;

		vector  <CFiberPath  *>  v_possible_paths;

		CFiberDemand  *pc_anycast_twin_dfemand;
		int i_demand_type;
		double  d_demand_capacity;


		void  v_copy_offsets(vector  <int> *pvDest, vector<int>  *pvSource);
	};//class  CFiberConnection




	class  CFiberPath
	{
		friend class CFiberNet;
		friend class CFiberSolution;
	public:
		CFiberPath();
		~CFiberPath();	

		void  vSetActive(bool  bNewActive)  {b_active = bNewActive;}
		bool  bGetActive()  {return(b_active);}

		void  vSetTree(bool  bNewTree)  {b_tree = bNewTree;}
		bool  bGetTree()  {return(b_tree);}

		void  vConfigure(int  iPathId, CFiberNet *pcParent);
		CError  eLoad(CString  sPathLine, bool  bControlConnectionsOrder = true);
		CError  eLoadMod(CString  sModLine);

		void  vSetTransferMultiplier(double  dNewTransferMultiplier) {d_new_transfer_multiplier = dNewTransferMultiplier;};
		CError  eGetSlotsNumberForTransfer(CString  sLine);
		CError  eGetRegNumberForTransfer(CString  sLine);
		

		void  vSaveInGenotype(int  *piGenotype);
		void  vSave(FILE  *pfDest);
		void  vSaveWithLinks(FILE  *pfDest);
		void  vSaveMod(FILE  *pfDest);
		void  vSaveSpec(FILE  *pfDest);
		void  vSaveReg(FILE  *pfDest);
				

		int  iGetId()  {return(i_path_id);};

		bool  bIsConnPartOfPath(CFiberConnection  *pcFiberConn);

		CError  eFindAndSetPathSlots(int  *piFirstSlot, int  iSlotsNumber);
		CError  eFindPath(int  *piFirstSlot, int  iSlotsNumber);
		CError  eSetPath(int  iFirstSlot, int  iSlotsNumber);

		CError  eSetPathSlots(int  iFirstSlot, int  iSlotsNumber);
		CError  eSetPathCapacity(int  iFirstSlot, double  dDemandCapacity);
		int  iGetSlotsNumber(double  dDemandCapacity);

		void  vSetPathFitness(double  dSolutionValue);
		double  dGetPathFitness()  {return(d_path_fitness);}

		int  iGetStartNode()  {return(i_fiber_node_start_offset);}
		int  iGetEndNode()  {return(i_fiber_node_end_offset);}

		
	protected:
		int  i_get_first_connection_offset(vector<CFiberConnection *>  *pvConnections);
		int  i_get_linked_connection_offset(vector<CFiberConnection *>  *pvConnections, CFiberConnection *pcMotherConnection);



		bool  b_active;
		bool  b_tree;
		int  i_path_id;
		int  i_fiber_node_start_offset;
		int  i_fiber_node_end_offset;

		//modulation
		double  d_length;
		double  d_mp;
		double  d_delta_s;

		CFiberNet *pc_parent;

		vector<CFiberConnection *>  v_path_connections;

		double  d_new_transfer_multiplier;
		vector<int>  v_slot_number_for_transfer;//increases every 10 Gbps
		vector<int>  v_reg_number_for_transfer;//increases every 10 Gbps


		double  d_path_fitness;
	};//class  CFiberPath



	class  CFiberSolution
	{
	public:
		CFiberSolution(CFiberNet  *pcFibNet);
		CFiberSolution(CFiberSolution  *pcFibNet);
		~CFiberSolution();


		void  vSetPathFitness();
		//CFiberPath  *pcGetCleverPath(int  iDemandNumber, gcroot<Random*> pcRandomGen);
		//CFiberPath  *pcGetBestCleverPath(int  iStartNode, int iEndNode);
		//CFiberPath  *pcGetRandomPath(int  iDemandNumber, gcroot<Random*> pcRandomGen);


		CError  eLoadSolutionWalk(CString  sSource);
		CError  eLoadSolutionWalk(FILE  *pfSource);

		CError  eLoadSolutionAibin(CString  sSource);
		CError  eLoadSolutionAibin(FILE  *pfSource);

		CError  eLoadSolutionAFA(CString  sSource);
		CError  eLoadSolutionAFA(FILE  *pfSource);
		CError  eLoadSolution(CString  sSource);
		CError  eLoadSolution(FILE  *pfSource);
		CError  eLoadSolutionPRW(CString  sSource);
		CError  eLoadSolutionPRW(FILE  *pfSource);

		CError  eSaveSolution(CString  sSolDest);
		CError  eSaveSolution(FILE  *pfDest);


		double  dCompCost(bool  *pbConstraintBroken, bool  bNoPenalty = false);
		CFiberNet  *pcGetNet() {return(pc_fib_net);};
		CError  eFindAndSetChanel(int  iChanelOffset, CFiberPath  *pcPath);
		CError  eFindChanel(int  iChanelOffset, CFiberPath  *pcPath, int *piStartSlot, int  *piSlotsNumber);
		CError  eSetChanel(int  iChanelOffset, CFiberPath  *pcPath, int iStartSlot, int  iSlotsNumber);

		CError  eApplySolution();

		void  vCopyChannelsState(CFiberSolutionChannel  *pcChannelsDest);

		CFiberSolutionChannel  *pcGetChannels() {return(pc_channels);}

		static  double  dGetFuncEval()  {return(d_func_evaluations);};
		static  double  dGetCleverPathsEval()  {return(d_clever_paths_evaluations);};
		static  double  dGetCleverPathsInits()  {return(d_clever_paths_initializations);};
		static  double  dGetConstructionTime()  {return(d_construction_time);};
		static  void  vAddConstructionTime(double  dAdd)  {d_construction_time += dAdd;};
		
		static  void  vZeroEval()  {d_func_evaluations = 0; d_clever_paths_evaluations = 0; d_clever_paths_initializations = 0; d_construction_time = 0;};

		
	private:
		static  double  d_func_evaluations;
		static  double  d_clever_paths_evaluations;
		static  double  d_clever_paths_initializations;
		static  double  d_construction_time;
		double  d_last_raw_sol_value;

		CFiberNet  *pc_fib_net;
		CFiberSolutionChannel  *pc_channels;


		CError  e_load_anycast(FILE  *pfSource);
		CError  e_load_unicast(FILE  *pfSource);

		//CFiberPath  *pc_get_clever_path(vector<CFiberPath*>  *pvPossiblePaths, double dPathFitSumm, gcroot<Random*> pcRandomGen);
	
	};//class  CFiberSolution



	class  CFiberSolutionChannel
	{
	public:
		CFiberSolutionChannel();
		~CFiberSolutionChannel();


		CError  eSetChannel();
		CError  eFindAndSetChanel(CFiberPath  *pcChannelPath);
		CError  eFindChanel(CFiberPath  *pcChannelPath, int *piStartSlot, int  *piSlotsNumber);
		CError  eSetChanel(CFiberPath  *pcChannelPath, int iStartSlot, int  iSlotsNumber);

		void  vConfig(CFiberSolution  *pcSolution, CFiberDemand  *pcDemand);

		void  vSave(FILE  *pfDest);

		void  vCopy(CFiberSolutionChannel  *pcOther);
		bool operator==(CFiberSolutionChannel  &cOther);
		bool operator!=(CFiberSolutionChannel  &cOther);


		CFiberDemand  *pcGetDemand()  {return(pc_demand);};

		CFiberPath  *pcPath;
		int  iFirstSlotNumber;

	private:
		CFiberSolution  *pc_solution;
		CFiberDemand  *pc_demand;
	
	};//class  CFiberSolution




	#define  STEINER_TREE_NETWORK	"net file"
	#define  STEINER_TREE_TREE_DEMANDS	"dmc file"
	#define  STEINER_TREE_TREE_NUMBER	"tree number"
	#define  STEINER_TREE_OPTIMIZE	"intelligent tree build"
	
	class  CSteinerTree
	{
	friend class CSteinerTreeSet;
	public:
		CSteinerTree();
		CSteinerTree(CFiberNet  *pcBaseNetwork);
		~CSteinerTree();

		void  vSetBaseNetwork(CFiberNet  *pcBaseNetwork);

		void  vCopyTo(CSteinerTree  *pcOtherTree);


		CError  eLoadTreeDemands(CString  sDMCFile);
		CError  eLoadTreeDemands(FILE  *pfDMCFile);

		CError  eSaveTreeDemands(CString  sDMCFile);
		CError  eSaveTreeDemands(FILE  *pfDMCFile);

		
		CError  eLoadTree(CString  sTreeFile);
		CError  eLoadTree(FILE  *pfTreeFile, int  *piLoaded);

		CError  eSaveTree(CString  sTreeFile, bool  bDetails = false);
		CError  eSaveTree(FILE  *pfTreeFile);
		CError  eSaveTreeDetails(FILE  *pfTreeFile);
		CError  eSaveTreeStructure(FILE  *pfTreeFile);


		CError  eBuildTreeFromConns(int  *piConnsIncluded);
		void  vReadConns(int  *piConnsIncluded);
				

		double  dGetTreeCost();
		bool  bCheckTree();
		bool  bGenotypesTheSame(CSteinerTree  *pcOtherTree);
		CFiberNet  *pcGetNetwork()  {return(pc_base_network);}
		vector<int>  *pvGetDemandedNodes()  {return(&v_demanded_nodes);}
		int  iGetRoot()  {return(i_root);}


	private:
		void  v_zero_tree_links();
		void  v_zero_node_tool();

		void  v_save_nodes_state(FILE  *pfTreeFile);

		bool  b_check_tree_no_loops();
		bool  b_check_tree_nodes_covered();//RUN ONLY AFTER b_check_tree_no_loops!!!
		bool  b_check_tree_nodes_covered_and_no_loops();
	    bool  b_check_tree_all_connected();
		bool  b_check_tree_hanging_leaves();//RUN ONLY AFTER b_check_tree_all_connected!!!


		void  v_cut_leaves();
		bool  b_node_is_not_demanded(int  iNode);
		CError  e_build_one_node(int  iNodeOffset, int iNodeLevel, int  *piConnsIncluded, bool  bForceConnections);


		double  d_tree_cost;
		bool  b_tree_cost_actual;

		int  *pi_tree_links;
		int  *pi_nodes_tool;	

	protected:
		vector<int>  v_demanded_nodes;
		vector<int>  v_not_demanded_nodes;
		int  i_root;
		CFiberNet  *pc_base_network;
	};//class  CSteinerTree


	struct SSteinerTreeSort
	{
		bool operator()(CSteinerTree* p1, CSteinerTree* p2)
		{
			if(!p1)
			return true;
			if(!p2)
			return false;

			bool  b_res;

			return p1->dGetTreeCost() < p2->dGetTreeCost();
		}//bool operator()(CSteinerTree* p1, CSteinerTree* p2)
	};//struct SSteinerTreeSort


	//class  DeceptiveTools::CDecFunc;
	#define  DEC_STEINER_TREE_SET	120
	#define  STEINER_TREE_SET_MAX_PENALTY_FACTOR	1000000
	class  CSteinerTreeSet : public CSteinerTree, ProblemTools::CCompProblem //zeby moc ladowac wymogi dla drzewa rowniez w w CTreeSet bez przepisywania procedur
	{
	friend class CSSJNetwork;
	public:
		CSteinerTreeSet();
		~CSteinerTreeSet();

		void  vSetCapacity(int  iNewCapacity);

		
		bool  bTreeContained(CSteinerTree  *pcNewTree);
		bool  bAddTree(CSteinerTree  *pcNewTree);
		double  dGetTreeFitness(CSteinerTree  *pcTree, double  dDefaultFitValForTreesInSet);
		
		CError  eLoadTreeSet(CString  sTreeSetFile);
		CError  eLoadTreeSet(FILE  *pfTreeSetFile);

		CError  eSaveTreeSet(CString  sTreeFileSet);
		CError  eSaveTreeSet(FILE  *pfTreeSetFile);

		CError  eSaveTreeSetSummary(CString  sSummaryFile, CString  sTreeSetname);
		CError  eSaveTreeSetSummary(FILE  *pfSummaryFile, CString  sTreeSetname);

		double  dGetBestTree();
		double  dGetWorstTree();

		//from CDecFunc
		double  dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen);
		double  dGetFuncValue(int  iUnitation);

		bool  bCheckUpTreeAll(int  iShift, int  *piArgument,  int  iArgLen);
		bool  bCheckUpTreeLoops(int  iShift, int  *piArgument,  int  iArgLen);


		CError  eLoadSettings(FILE *pfSettings);
		int  iGetBitLength()  {return(pc_base_network->iGetConnsNumber());};
		int  iGetFuncType()  {return(DEC_STEINER_TREE_SET);};

		CCompProblem *pcClone();

		void  v_create_all_possible_r_graphs_of_size_step
			(
			vector<vector<int>>  *pvNodesToInclude,
			vector<int>  *pvNodesToChoose, vector<int>  *pvChosenNodes, int  iSize, int iCurPos, int iStep
			);
		void  v_create_all_possible_r_graphs_of_size(vector<vector<int>>  *pvNodesToInclude, vector<int>  *pvNodesToChoose, int  iSize);
		//void  vGenerateTreesByBrazilians(double  dMaxTime, System::Windows::Forms::ListBox*  listInfo, double *pdFfe, bool  bUseRealNodeCosts);
		//void  vGenerateTreesBySenOnSorensenJansens(double  dMaxTime, System::Windows::Forms::ListBox*  listInfo, double *pdFfe, bool  bUseRealNodeCosts);
		void  vGenerateTreesBySenOnSorensenJansensTest(bool  bUseRealNodeCosts);
		void  vGetShortestPathBetweenNodePair(int  iStartNode, int iEndNode, CSSJConn *pcResult, bool  bUseRealNodeCosts);
		
		vector <CSteinerTree  *>  *pvGetTreeSet()  {return(&v_tree_set);}

		double  dGetPenaltyFactor()  {return(d_penalty_factor);}
		void  vPenaltyFactorReset()  {d_penalty_factor = 1;}
		void  vPenaltyFactorSet(double  dNewPenaltyFactor)  {if  (dNewPenaltyFactor  >  STEINER_TREE_SET_MAX_PENALTY_FACTOR)  d_penalty_factor = STEINER_TREE_SET_MAX_PENALTY_FACTOR; else d_penalty_factor = dNewPenaltyFactor;}


	private:
		void  v_insert_spanning_tree(vector<CSSJNetwork *>  *pvSpanningTrees, CSSJNetwork *pcNewTree);
		void  v_add_and_order_by_cost(vector<CSearchShortestPathNode *>  *pvNodes, CSearchShortestPathNode *pcNodeToAdd);
		void  v_order_trees();
		vector <CSteinerTree  *>  v_tree_set;
		int  i_capacity;
		double  d_penalty_factor;
		int  i_optimize;

		vector <CSSJConn>  v_shortest_paths;

	
	};//class  CSteinerTreeSet


	class  CSearchShortestPathNode
	{
	public:
		CSearchShortestPathNode();
		~CSearchShortestPathNode();

		void  vCreateConn(CSSJConn  *pcResult);
		void  vCreateCandidates
			(
			vector<CSearchShortestPathNode *> *pvChosenNodes, vector<CSearchShortestPathNode *> *pvNewNodes, 
			CFiberConnection  ***pcConnMat, int iNodeNumber, 
			bool  bUseRealNodeCosts
			);


		bool  b_check_if_node_is_contained(vector<CSearchShortestPathNode *> *pvNodes, int iNodeOffset);


		int  iNode;
		double  dCost;
		CFiberConnection  *pc_get_to_node_connection;
		CSearchShortestPathNode *pc_parent;
		vector<CSearchShortestPathNode *> v_children;

	};//class  CSearchShortestPathNode


	class CSSJConn
	{
	public:
		CSSJConn() {};
		~CSSJConn() {};

		void  vSave(FILE  *pfDest);

		int iNodeStartOffset;
		int iNodeEndOffset;

		int iNodeStartReal;
		int iNodeEndReal;

		double  dCost;

		vector  <CFiberConnection *>  v_path;
	
	};//class CSSJConn



	class CSSJNode
	{
	public:
		CSSJNode() {};
		~CSSJNode() {};

		double  dGetCost(CSSJNode  *pcNode);
	
		int  i_node_offset;
		int  i_node_real_net_id;

		vector<CSSJConn  *>  v_conns;
	};//class CSSJNode



	class CSSJNetwork
	{
	public:
		CSSJNetwork();
		~CSSJNetwork();

		void  vCreateCompeleteGraph(vector<int>  *pvNodesToInclude, CSteinerTreeSet  *pcSTreeSet, int iRoot, bool  bUseRealNodeCosts);
		void  vCreateCompeleteGraph(CSteinerTreeSet  *pcStreeSet, int iRoot, bool  bUseRealNodeCosts);
		bool  bFindMinSpaningTreePrim(CSSJNetwork  *pcMinSpanTree);
		void  vFillGenotype(int  *piConnBasedGenotype);
		
		int iGetNodeCurrentOffset(int  iNodeOffset);
		CSSJNode  * CSSJNetwork::pcGetNode(int  iNodeOffset);
		bool  bDoIContainThisNode(int  iNodeOffset);
		void  vGenerateNetsToAnalyze(vector<CSSJNetwork *> *pvNetworksToAnalyze);
		void  vRemoveConnection(int iRealNodeStart, int iRealNodeEnd);

		double  dGetCost();

		void  vSaveGraph(FILE  *pfDest);

		CSSJNetwork  *pc_parent_network;

		CSSJNetwork  *pcClone();
	private:


		vector<CSSJNode>  v_nodes;
		vector<CSSJConn *>  v_conns;	
	};//class CSSJNetwork


  




	
}//namespace  FiberNets


