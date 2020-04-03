//#include "stdafx.h"

#include "FiberNet.h"
#include "CompProblem.h"


using  namespace  FiberNets;
using namespace TimeCounters;
using namespace std;
using namespace Tools;




//---------------------------class  CFiberNet--------------------------------------
CFiberNet::CFiberNet()
{
	pc_nodes = NULL;
	pc_conns = NULL;
	pc_conn_mat = NULL;
	pc_replices = NULL;
	pc_demands = NULL;
	pc_paths = NULL;
	pc_trees = NULL;
	pv_paths_between_nodes = NULL;

	
	i_nodes_num = 0;
	i_conns_num = 0;
	i_rep_num = 0;
	i_demand_num = 0;
	i_path_num = 0;
	i_tree_num = 0;
}//CFiberNet::CFiberNet()


CFiberNet::~CFiberNet()
{

	if  ( (i_nodes_num  >  0)&&(pc_conn_mat != NULL))
	{
		for  (int  ix = 0; ix < i_nodes_num; ix++)
		{
			delete  pc_conn_mat[ix];
			delete  pv_paths_between_nodes[ix];
		}//for  (int  ix = 0; ix < i_nodes_num; ix++)

		delete  pc_conn_mat;
		delete  pv_paths_between_nodes;
	}//if  ( (i_nodes_num  >  0)&&(pc_conn_mat != NULL))


	if  (pc_nodes != NULL)  delete  pc_nodes;
	if  (pc_conns != NULL)  delete  pc_conns;
	if  (pc_replices != NULL)  delete  pc_replices;
	if  (pc_demands != NULL)  delete  pc_demands;
	if  (pc_paths != NULL)  delete  pc_paths;
	
}//CFiberNet::~CFiberNet()




double  CFiberNet::dCompCost()
{
	double  d_result;

	d_result = 0;
	int  i_buf;
	for  (int  ii = 0; ii < i_conns_num; ii++)
	{
		pc_conns[ii].iGetMaxSlot();
		i_buf = pc_conns[ii].iGetMaxSlot();
		if  (i_buf  >  d_result)  d_result = i_buf;
	}//for  (int  ii = 0; ii < i_conns_num; ii++)

	return(d_result);
}//double  CFiberNet::dCompCost()





CError  CFiberNet::eResetNet()
{
	CError  c_err;

	for  (int  ii = 0; ii < i_conns_num; ii++)
		pc_conns[ii].vReset();

	return(c_err);
}//CError  CFiberNet::eResetNet()



CError  CFiberNet::eCloneTopologyTo(CFiberNet  *pcOtherNetwork)
{
	CError  c_err;

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "0", true);

	if  (pcOtherNetwork->pc_nodes  !=  NULL)  delete  pcOtherNetwork->pc_nodes;
	if  (pcOtherNetwork->pc_conns  !=  NULL)  delete  pcOtherNetwork->pc_conns;

	if  (pcOtherNetwork->pc_conn_mat  !=  NULL)
	{
		for  (int  ix = 0; ix < pcOtherNetwork->i_nodes_num; ix++)
		  delete pcOtherNetwork->pc_conn_mat[ix];
		delete pcOtherNetwork->pc_conn_mat;	
	}//if  (pcOtherNetwork->pc_conn_mat  !=  NULL)

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "1", true);

	if  (pcOtherNetwork->pv_paths_between_nodes  !=  NULL)
	{
		for  (int  ix = 0; ix < pcOtherNetwork->i_nodes_num; ix++)
		  delete pcOtherNetwork->pv_paths_between_nodes[ix];
		delete pcOtherNetwork->pv_paths_between_nodes;	
	}//if  (pcOtherNetwork->pc_conn_mat  !=  NULL)

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "2", true);


	pcOtherNetwork->i_nodes_num = i_nodes_num;
	pcOtherNetwork->i_conns_num = i_conns_num;

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "3", true);

	pcOtherNetwork->pc_nodes  = new CFiberNode[i_nodes_num];
	pcOtherNetwork->pc_conns = new CFiberConnection[i_conns_num];

	pcOtherNetwork->pc_conn_mat = new  CFiberConnection** [i_nodes_num];
	pcOtherNetwork->pv_paths_between_nodes = new vector<int>* [i_nodes_num];

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "4", true);

	for  (int  ix = 0; ix < i_nodes_num; ix++)
	{
		pcOtherNetwork->pc_conn_mat[ix] = new CFiberConnection* [i_nodes_num];
		pcOtherNetwork->pv_paths_between_nodes[ix] = new vector<int>[i_nodes_num];
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "5", true);

	for  (int  ix = 0; ix < i_nodes_num; ix++)
	{
		for  (int  iy = 0; iy < i_nodes_num; iy++)
		{
			pcOtherNetwork->pc_conn_mat[ix][iy] = NULL;
		}//for  (int  iy = 0; iy < i_nodes_num; iy++)
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "5", true);


	for  (int  ix = 0; ix < i_nodes_num; ix++)  pcOtherNetwork->pc_nodes[ix].vSetNodeNumber(ix);


	for  (int ii = 0; ii < i_conns_num; ii++)
	{
		pcOtherNetwork->pc_conns[ii].eConfig
			(
			pc_conns[ii].i_connection_id, pc_conns[ii].i_fiber_node_start_offset, pc_conns[ii].i_fiber_node_end_offset, pcOtherNetwork, 
			pc_conns[ii].d_connection_capacity
			);

		pcOtherNetwork->pc_conn_mat[pcOtherNetwork->pc_conns[ii].i_fiber_node_start_offset][pcOtherNetwork->pc_conns[ii].i_fiber_node_end_offset] = &(pcOtherNetwork->pc_conns[ii]);
		pcOtherNetwork->pc_nodes[pcOtherNetwork->pc_conns[ii].i_fiber_node_start_offset].vAddOutConn(&(pcOtherNetwork->pc_conns[ii]));
		pcOtherNetwork->pc_nodes[pcOtherNetwork->pc_conns[ii].i_fiber_node_end_offset].vAddInConn(&(pcOtherNetwork->pc_conns[ii]));
	}//for  (int ii = 0; ii < i_conns_num; ii++)

	//::vReportInFile("zzz_eCloneTopologyTo.txt", "6", true);


	return(c_err);
}//CError  CFiberNet::eCloneTopologyTo(CFiberNet  *pcOtherNetwork)


CError  CFiberNet::eLoadTopology(CString  sNetFile)
{
	CError  c_err;

	//Tools::vReportInFile("zzzz.txt",sNetFile);
	
	//c_err.vSetError("test");
	//return(c_err);

	FILE  *pf_source;

	pf_source = fopen(sNetFile, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No net file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	i_nodes_num = Tools::iReadLine(pf_source);
	i_conns_num = Tools::iReadLine(pf_source);

	if  (i_nodes_num  <  1)
	{
		c_err.vSetError("Nodes number must be greater than 0");
		return(c_err);
	}//if  (i_nodes_num  <  1)


	if  (i_conns_num  <  1)
	{
		c_err.vSetError("Connections number must be greater than 0");
		vShow(i_conns_num);
		return(c_err);
	}//if  (i_nodes_num  <  1)


	pc_nodes  = new CFiberNode[i_nodes_num];
	pc_conns = new CFiberConnection[i_conns_num];




	pc_conn_mat = new  CFiberConnection** [i_nodes_num];
	pv_paths_between_nodes = new vector<int>* [i_nodes_num];

	for  (int  ix = 0; ix < i_nodes_num; ix++)
	{
		pc_conn_mat[ix] = new CFiberConnection* [i_nodes_num];
		pv_paths_between_nodes[ix] = new vector<int>[i_nodes_num];
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)


	for  (int  ix = 0; ix < i_nodes_num; ix++)
	{
		for  (int  iy = 0; iy < i_nodes_num; iy++)
		{
			pc_conn_mat[ix][iy] = NULL;
		}//for  (int  iy = 0; iy < i_nodes_num; iy++)
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)


	for  (int  ix = 0; ix < i_nodes_num; ix++)  pc_nodes[ix].vSetNodeNumber(ix);



	CString  s_line;
	double  d_value;
	int  i_connection_id;
	int  i_index;
	int  i_buf;

	i_connection_id = 0;
	for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)
	{
		s_line  =  Tools::sReadLine(pf_source);
		//::MessageBox(NULL, s_line, s_line, MB_OK);
		
		i_index = 0;
		for  (int i_node_end = 0; i_node_end < i_nodes_num; i_node_end++)
		{
			/*i_buf = i_node_end;
			i_node_end = i_node_start;
			i_node_start = i_buf;//*/


			d_value  =  Tools::dExtractFromString(s_line, i_index, &i_index);

			if  (d_value  >  0)
			{
				c_err  =  pc_conns[i_connection_id].eConfig(i_connection_id, i_node_start, i_node_end, this, d_value);
				pc_conn_mat[i_node_start][i_node_end] = &(pc_conns[i_connection_id]);
				pc_nodes[i_node_start].vAddOutConn(&pc_conns[i_connection_id]);
				pc_nodes[i_node_end].vAddInConn(&pc_conns[i_connection_id]);

				i_connection_id++;

				if  (c_err)
				{
					fclose(pf_source);
					return(c_err);
				}//if  (c_err)
			}//if  (d_value  >  0)


			/*i_buf = i_node_end;
			i_node_end = i_node_start;
			i_node_start = i_buf;//*/
		}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)
		
	}//for  (int ii = 0; ii < i_nodes_num; ii++)



	fclose(pf_source);
	return(c_err);
}//CError  CFiberNet::eLoadTopology(CString  sNetFile)



CFiberConnection  *CFiberNet::pcGetConn(int  iOffset)
{
	if  (iOffset  <  0)  return(NULL);
	if  (iOffset  >=  i_conns_num)  return(NULL);

	return(&(pc_conns[iOffset]));
}//CFiberConnection  *CFiberNet::pcGetConn(int  iOffset)




bool  CFiberNet::bCheckIfReplicaNode(CFiberNode  *pcNode)
{
	for  (int  ii = 0; ii < i_rep_num; ii++)
	{
		if  (pcNode->i_node_number  ==  pc_replices[ii]->i_node_number)  return(true);
	}//for  (int  ii = 0; ii < i_rep_num; ii++)

	return(false);
}//bool  CFiberNet::bCheckIfReplicaNode(CFiberNode  *pcNode)


bool  CFiberNet::bCheckIfReplicaNode(int  iNodeId)
{
	for  (int ii = 0; ii < i_rep_num; ii++)
	{
		if  (pc_replices[ii]->i_node_number  ==  iNodeId)  return(true);	
	}//for  (int ii = 0; ii < i_rep_num; ii++)

	return(false);
}//bool  CFiberNet::bCheckIfReplicaNode(int  iNodeId)




bool  CFiberNet::bCompareNodes(int  iNodeId0, int iNodeId1)
{
	if  (iNodeId0 == iNodeId1)  return(true);

	if  (iNodeId0  ==  -1)
	{
		return(bCheckIfReplicaNode(iNodeId1));
	}//if  (iNodeId0  ==  -1)

	if  (iNodeId1  ==  -1)
	{
		return(bCheckIfReplicaNode(iNodeId0));
	}//if  (iNodeId1  ==  -1)

	return(false);
}//bool  CFiberNet::bCompareNodes(int  iNodeId0, int iNodeId1)





CError  CFiberNet::eLoadDemands(CString  sUnicastFile, CString  sAnycastFile, CString sMulticastFile, int  iPathsToUse, int  iTreesToUse)
{
	CError  c_err;
	CString  s_buf;
	FILE  *pf_unicast_source, *pf_anycast_source, *pf_multicast_source;

	pf_unicast_source = NULL;//musi byc bo to jest potem sprawdzane

	if  (sUnicastFile  !=  "")
	{
		pf_unicast_source = fopen(sUnicastFile, "r+");
		if  (pf_unicast_source  ==  NULL)
		{
			s_buf.Format("No unicast dem file defined 1 (%s)", sUnicastFile);
			c_err.vSetError(s_buf);
			return(c_err);
		}//if  (pf_source  ==  NULL)
	}//if  (sUnicastFile  !=  "")


	/*pf_anycast_source = fopen(sAnycastFile, "r+");
	if  (pf_unicast_source  ==  NULL)
	{
		fclose(pf_unicast_source);
		c_err.vSetError("No anycast dem file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)*/


	pf_multicast_source = fopen(sMulticastFile, "r+");
	if  (pf_multicast_source  ==  NULL)
	{
		c_err.vSetError("No multicast dem file defined 2");
		return(c_err);
	}//if  (pf_multicast_source  ==  NULL)
	
	if  (sUnicastFile  !=  "")
		i_uni_demands = Tools::iReadLine(pf_unicast_source);
	else
		i_uni_demands = 0;


	i_multi_demands = Tools::iReadLine(pf_multicast_source);


	//i_any_demands = Tools::iReadLine(pf_anycast_source);
	//i_any_demands *= 2;
	i_any_demands = 0;
	i_demand_num = i_uni_demands + i_any_demands + i_multi_demands;
	pc_demands = new CFiberDemand[i_demand_num];



	if  (sUnicastFile  !=  "")  c_err = e_load_uni_demands(pf_unicast_source, iPathsToUse);
	if  (c_err)
	{
		if  (pf_unicast_source  !=  NULL)  fclose(pf_unicast_source);
		fclose(pf_anycast_source);
		fclose(pf_multicast_source);
		return(c_err);	
	}//if  (c_err)

	/*c_err = e_load_any_demands(pf_anycast_source);
	if  (c_err)
	{
		if  (pf_unicast_source  !=  NULL)  fclose(pf_unicast_source);
		fclose(pf_anycast_source);
		fclose(pf_multicast_source);
		return(c_err);
	}//if  (c_err)*/

	c_err = e_load_multi_demands(pf_multicast_source, iTreesToUse);
	if  (c_err)
	{
		if  (pf_unicast_source  !=  NULL)  fclose(pf_unicast_source);
		fclose(pf_anycast_source);
		fclose(pf_multicast_source);
		return(c_err);	
	}//if  (c_err)

	
	
	

	if  (pf_unicast_source  !=  NULL)  fclose(pf_unicast_source);
	//fclose(pf_anycast_source);
	fclose(pf_multicast_source);
	return(c_err);
}//CError  CFiberNet::eLoadDemands(CString  sUnicastFile, CString  sAnycastFile)




CError  CFiberNet::e_load_any_demands(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;

	CString  s_line;
	int  i_client_node;
	double  d_downstream, d_upstream;

	int  i_index;
	int  i_demand_offset = 0;
	int  i_anycast_lines_to_read = i_any_demands / 2;
	for  (int i_demand_line = 0; i_demand_line < i_anycast_lines_to_read; i_demand_line++)
	{
		s_line  =  Tools::sReadLine(pfSource);

		i_index = 0;
		i_client_node = Tools::iExtractFromString(s_line, i_index, &i_index);

		d_downstream = Tools::iExtractFromString(s_line, i_index, &i_index);
		d_upstream  =  Tools::dExtractFromString(s_line, i_index, &i_index);

		if  (bCheckIfReplicaNode(&(pc_nodes[i_client_node]))  ==  true)
		{
			i_demand_num = i_demand_num - 2;
			i_any_demands = i_any_demands - 2;

			//s_buf.Format("%d  %d", i_client_node, pc_nodes[i_client_node].i_node_number);
			//vShow(s_buf);
		}//if  (bCheckIfReplicaNode(pc_nodes[i_client_node])  ==  true)
		else
		{
			pc_demands[i_uni_demands + i_demand_offset * 2].eConfig
				(
				i_uni_demands + i_demand_offset * 2, 
				-1, i_client_node, d_downstream, 
				this, DEMAND_TYPE_ANYCAST, &(pc_demands[i_uni_demands + i_demand_offset * 2 + 1])
				);


			pc_demands[i_uni_demands + i_demand_offset * 2 + 1].eConfig
				(
				i_uni_demands + i_demand_offset * 2 + 1, 
				i_client_node, -1, d_upstream, 
				this, DEMAND_TYPE_ANYCAST, &(pc_demands[i_uni_demands + i_demand_offset * 2])
				);

			i_demand_offset++;
		}//else  if  (bCheckIfReplicaNode(pc_nodes[i_client_node])  ==  true)

	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	return(c_err);
}//CError  CFiberNet::eLoadAnyDemands(CString  sDemFile)



CError  CFiberNet::e_load_uni_demands(FILE  *pfSource, int iMaxPathsToUseBetweenNodes)
{
	CError  c_err;


	CString  s_line;
	int  i_node_start, i_node_end;
	double  d_capacity;

	int  i_index;
	for  (int i_demand_offset = 0; i_demand_offset < i_uni_demands; i_demand_offset++)
	{
		s_line  =  Tools::sReadLine(pfSource);
		//Tools::vShow(s_line);

		i_index = 0;
		i_node_start = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_node_end = Tools::iExtractFromString(s_line, i_index, &i_index);

		d_capacity  =  Tools::dExtractFromString(s_line, i_index, &i_index);

		//Tools::vShow(i_node_start);
		//Tools::vShow(i_node_end);
		//Tools::vShow(d_capacity);

		c_err  = pc_demands[i_demand_offset].eConfig
			(
			i_demand_offset, 
			i_node_start, i_node_end, d_capacity, this,
			DEMAND_TYPE_UNICAST, NULL, iMaxPathsToUseBetweenNodes
			);
		if  (c_err)
		{
			fclose(pfSource);
			return(c_err);
		}//if  (c_err)
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	return(c_err);
}//CError  CFiberNet::e_load_unicast(CString  sDemFile)




CError  CFiberNet::e_load_multi_demands(FILE  *pfSource, int iMaxTreesToUseForNode)
{
	CError  c_err;


	CString  s_line;
	double  d_capacity;
	int  i_children_number;
	int  i_node_start;
	
	CString  s_report_buf;


	int  i_demand_offset;
	int  i_index;
	for  (int i_multi_demand = 0; i_multi_demand < i_multi_demands; i_multi_demand++)
	{
		s_line  =  Tools::sReadLine(pfSource);

		//::vShow(s_line);

		i_index = 0;

		d_capacity  =  Tools::dExtractFromString(s_line, i_index, &i_index);
		i_children_number = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_node_start = Tools::iExtractFromString(s_line, i_index, &i_index);

		//::vShow(d_capacity);
		//::vShow(i_children_number);
		//::vShow(i_node_start);

		

		i_demand_offset = i_multi_demand + i_uni_demands + i_any_demands * 2;

		s_report_buf.Format("%d  %d  %s\n", i_multi_demand, i_demand_offset, s_line);
		::Tools::vReportInFile("!!test.txt", s_report_buf);

		c_err  = pc_demands[i_demand_offset].eConfig
			(
			i_demand_offset, 
			i_node_start, -1, d_capacity, this,
			DEMAND_TYPE_MULTICAST, NULL,
			iMaxTreesToUseForNode
			);
		if  (c_err)
		{
			fclose(pfSource);
			return(c_err);
		}//if  (c_err)*/
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)

	return(c_err);
}//CError  CFiberNet::e_load_multi_demands(FILE  *pfSource)




CError  CFiberNet::eLoadReplicas(CString  sRepFile)
{
	CError  c_err;

	if  (pc_replices != NULL)  delete  pc_replices;
	pc_replices = NULL;


	FILE  *pf_source;
	pf_source = fopen(sRepFile, "r+");

	if  (pf_source  ==  NULL)
	{
		//NO ERROR HERE!!!
		//c_err.vSetError("No replices file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	i_rep_num = Tools::iReadLine(pf_source);
	pc_replices = new CFiberNode*[i_rep_num];


	CString  s_line;
	int  i_rep_node_offset;
	s_line  =  Tools::sReadLine(pf_source);
	int  i_last_offset = 0;

	for  (int ii = 0; ii < i_rep_num; ii++)
	{
		
		i_rep_node_offset = Tools::iExtractFromString(s_line, i_last_offset, &i_last_offset);

		if  (
			(i_rep_node_offset  >= 0)&&(i_rep_node_offset < i_nodes_num)
			)
		{
			pc_replices[ii] = &(pc_nodes[i_rep_node_offset]);
		}//if  (
		else
		{
			CString  s_buf;
			s_buf.Format("Wrong replica offset (%d)", i_rep_node_offset);
			c_err.vSetError(s_buf);
			fclose(pf_source);
			return(c_err);
		}//else  
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	fclose(pf_source);	

	return(c_err);
}//CError  CFiberNet::eLoadReplices(CString  sRepFile)


CError  CFiberNet::eLoadTrees(CString  sTreeFile)
{
	CError  c_err;

	if  (pc_paths != NULL)  
	{
		delete  pc_trees;
		pc_trees = NULL;
	}//if  (pc_paths != NULL)  
	i_tree_num = 0;


	FILE  *pf_source;
	pf_source = fopen(sTreeFile, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No tree file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	//first count lines;
	CString  s_line;
	int  i_lines;

	i_lines = 0;
	while (feof(pf_source)  ==  false)
	{
		s_line  =  Tools::sReadLine(pf_source);
		s_line.Replace(" ", "");
		if  (s_line.GetLength()  >  0)  i_lines++;
	}//while (foef(pf_source)  ==  false)
	fclose(pf_source);

	pf_source = fopen(sTreeFile, "r+");

	i_tree_num = i_lines;
	pc_trees = new CFiberPath[i_tree_num];

	

	for  (int i_tree_offset = 0; i_tree_offset < i_tree_num; i_tree_offset++)
	{
		s_line  =  Tools::sReadLine(pf_source);

		pc_trees[i_tree_offset].vConfigure(i_tree_offset, this);
		c_err = pc_trees[i_tree_offset].eLoad(s_line, false);
		pc_trees[i_tree_offset].vSetTree(true);
		if  (c_err)
		{
			fclose(pf_source);
			return(c_err);
		}//if  (c_err)

		/*pv_paths_between_nodes
			[pc_paths[i_path_offset].i_fiber_node_start_offset]
			[pc_paths[i_path_offset].i_fiber_node_end_offset]
			.push_back(i_path_offset);*/
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	fclose(pf_source);

	return(c_err);
}//CError  CFiberNet::eLoadTrees(CString  sTreeFile)


CError  CFiberNet::eLoadPaths(CString  sPathFile)
{
	CError  c_err;

	if  (pc_paths != NULL)  delete  pc_paths;


	FILE  *pf_source;
	pf_source = fopen(sPathFile, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No path file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	i_path_num = Tools::iReadLine(pf_source);
	pc_paths = new CFiberPath[i_path_num];


	CString  s_line;

	for  (int i_path_offset = 0; i_path_offset < i_path_num; i_path_offset++)
	{
		s_line  =  Tools::sReadLine(pf_source);

		pc_paths[i_path_offset].vConfigure(i_path_offset, this);
		c_err = pc_paths[i_path_offset].eLoad(s_line);
		if  (c_err)
		{
			fclose(pf_source);
			return(c_err);
		}//if  (c_err)

		pv_paths_between_nodes
			[pc_paths[i_path_offset].i_fiber_node_start_offset]
			[pc_paths[i_path_offset].i_fiber_node_end_offset]
			.push_back(i_path_offset);
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	fclose(pf_source);

	

	return(c_err);
}//CError  CFiberNet::eLoadPaths(CString  sPathFile)


CError  CFiberNet::e_load_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)
{
	CError  c_err;


	FILE  *pf_source;
	pf_source = fopen(sSpecFile, "r+");	

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No uni spec file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	CString  s_line;

	for  (int i_path_offset = 0; i_path_offset < iPathNum; i_path_offset++)
	{
		s_line  =  Tools::sReadLine(pf_source);

		c_err = pcPaths[i_path_offset].eGetSlotsNumberForTransfer(s_line);
		pcPaths[i_path_offset].vSetTransferMultiplier(10);
		if  (c_err)
		{
			fclose(pf_source);
			return(c_err);
		}//if  (c_err)

	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	fclose(pf_source);
	return(c_err);
}//CError  CFiberNet::e_load_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  *piPathNum)


CError  CFiberNet::eLoadPathsSpec(CString  sUniSpecFile)
{
	CError  c_err;

	c_err = e_load_spec(sUniSpecFile, pc_paths, i_path_num);

	return(c_err);
}//CError  CFiberNet::eLoadPathsSpec(CString  sUniSpecFile)


CError  CFiberNet::eLoadTreesSpec(CString  sTreeSpecFile)
{
	CError  c_err;

	c_err = e_load_spec(sTreeSpecFile, pc_trees, i_tree_num);

	return(c_err);
}//CError  CFiberNet::eLoadPathsSpec(CString  sUniSpecFile)


		

CError  CFiberNet::e_load_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)
{
	CError  c_err;


	FILE  *pf_source;
	pf_source = fopen(sSpecFile, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No uni reg file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	CString  s_line;

	for  (int i_path_offset = 0; i_path_offset < iPathNum; i_path_offset++)
	{
		s_line  =  Tools::sReadLine(pf_source);

		c_err = pcPaths[i_path_offset].eGetRegNumberForTransfer(s_line);
		if  (c_err)
		{
			fclose(pf_source);
			return(c_err);
		}//if  (c_err)

	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)


	fclose(pf_source);
	return(c_err);
}//CError  CFiberNet::e_load_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)



CError  CFiberNet::eLoadPathsReg(CString  sUniRegFile)
{
	CError  c_err;

	c_err = e_load_reg(sUniRegFile, pc_paths, i_path_num);

	return(c_err);
}//CError  CFiberNet::eLoadPathsReg(CString  sUniRegFile)



CError  CFiberNet::eLoadTreesReg(CString  sTreeRegFile)
{
	CError  c_err;

	c_err = e_load_reg(sTreeRegFile, pc_trees, i_tree_num);

	return(c_err);
}//CError  CFiberNet::eLoadTreesReg(CString  sTreeRegFile)




CError  CFiberNet::eLoadModulation(CString  sModFile)
{
	CError  c_err;

	CString  s_buf;

	FILE  *pf_source;
	pf_source = fopen(sModFile, "r+");

	if  (pf_source  ==  NULL)
	{
		s_buf.Format("No mod file defined (%s)", sModFile);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (pf_source  ==  NULL)

	int  i_mod_num;
	i_mod_num = Tools::iReadLine(pf_source);

	if  (i_mod_num  !=  i_path_num)
	{
		fclose(pf_source);
		CString  s_buf;
		s_buf.Format("different path and modulation number: %d != %d", i_path_num,  i_mod_num);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (i_mod_num  !=  i_path_num)


	CString  s_line;

	for  (int i_path_offset = 0; i_path_offset < i_path_num; i_path_offset++)
	{
		s_line  =  Tools::sReadLine(pf_source);

		pc_paths[i_path_offset].eLoadMod(s_line);
		if  (c_err)
		{
			fclose(pf_source);
			return(c_err);
		}//if  (c_err)
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)
	fclose(pf_source);

	return(c_err);
}//CError  CFiberNet::eLoadModulation(CString  sModFile)




CError  CFiberNet::eSaveModulation(CString  sModFile)
{
	CError  c_err;

	FILE  *pf_dest;
	pf_dest = fopen(sModFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No net file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)


	fprintf(pf_dest, "%d\n", i_path_num);

	

	for  (int i_path_offset = 0; i_path_offset < i_path_num; i_path_offset++)
	{
		pc_paths[i_path_offset].vSaveMod(pf_dest);
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)
	fclose(pf_dest);

	return(c_err);
}//CError  CFiberNet::eSaveModulation(CString  sModFile)




CError  CFiberNet::eSaveTopology(CString  sNetFile)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sNetFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No net file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	fprintf(pf_dest, "%d\n", i_nodes_num);
	fprintf(pf_dest, "%d\n", i_conns_num);


	for  (int  ix = 0; ix < i_nodes_num; ix++)
	{
		for  (int  iy = 0; iy < i_nodes_num; iy++)
		{
			if  (pc_conn_mat[ix][iy]  ==  NULL)
				fprintf(pf_dest, "0\t", i_nodes_num);
			else
				fprintf(pf_dest, "%d\t", (int) pc_conn_mat[ix][iy]->d_connection_capacity);
		}//for  (int  iy = 0; iy < i_nodes_num; iy++)
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)



	fclose(pf_dest);
	return(c_err);

}//void  CFiberNet::vSaveTopology(CString  sNetFile)



CError  CFiberNet::eSaveReplicas(CString  sRepFile)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sRepFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No rep file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	fprintf(pf_dest, "number: %d\n", i_rep_num);

	for  (int  ii = 0; ii < i_rep_num; ii++)
	{
		fprintf(pf_dest, "%d\n", pc_replices[ii]->i_node_number);
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	fprintf(pf_dest, "\n");



	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::eSaveReplicas(CString  sRepFile)




CError  CFiberNet::eSaveDemands(CString  sDemFile)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sDemFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No dem file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	fprintf(pf_dest, "ALL: %d\n", i_demand_num);
	fprintf(pf_dest, "Unicast: %d\n", i_uni_demands);
	fprintf(pf_dest, "Anycast: %d\n", i_any_demands);
	fprintf(pf_dest, "Multicast: %d\n", i_multi_demands);

	fprintf(pf_dest, "\n");

	for  (int  ii = 0; ii < i_uni_demands; ii++)
	{
		pc_demands[ii].vSave(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	fprintf(pf_dest, "\n");

	for  (int  ii = i_uni_demands; ii < i_uni_demands + i_any_demands; ii++)
	{
		pc_demands[ii].vSave(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)


	for  (int  ii = i_uni_demands + i_any_demands; ii < i_demand_num; ii++)
	{
		pc_demands[ii].vSave(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)



	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::eSaveDemands(CString  sDemFile)



CError  CFiberNet::eSaveTrees(CString  sTreeFile)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sTreeFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No tree file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	fprintf(pf_dest, "%d\n", this->i_tree_num);


	for  (int  ii = 0; ii < i_tree_num; ii++)
	{
		pc_trees[ii].vSave(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)



	fprintf(pf_dest, "\n");
	fprintf(pf_dest, "\n");
	/*for  (int  ii = 0; ii < pc_demands[5].pvGetPossiblePaths()->size(); ii++)
	{
		pc_demands[5].pvGetPossiblePaths()->at(ii)->vSaveWithLinks(pf_dest);
		fprintf(pf_dest, "\n");	
	}//for  (int  ii = 0; ii < pc_demands[5].pvGetPossiblePaths(); ii++)*/


	/*for  (int  ii = 0; ii < i_tree_num; ii++)
	{
		if  (pc_paths[ii].v_path_connections.size()  ==  5)
		{
			fprintf(pf_dest, "%d    root: %d\n",ii,  pc_trees[ii].i_fiber_node_start_offset);
			pc_trees[ii].vSaveWithLinks(pf_dest);
			fprintf(pf_dest, "\n");
			fprintf(pf_dest, "\n");
		}//if  (pc_paths[ii].v_path_connections.size()  ==  5)

	}//for  (int  ii = 0; ii < i_path_num; ii++)*/


	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::eSaveTrees(CString  sTreeFile)



CError  CFiberNet::eSavePaths(CString  sPathFile)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sPathFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No path file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	
	fprintf(pf_dest, "%d\n", i_path_num);


	for  (int  ii = 0; ii < i_path_num; ii++)
	{
		//if  ( (pc_paths[ii].i_fiber_node_start_offset == 18)&&(pc_paths[ii].i_fiber_node_end_offset == 25) )
		{
			//fprintf(pf_dest, "%d    %d->%d\n",ii,  pc_paths[ii].i_fiber_node_start_offset, pc_paths[ii].i_fiber_node_end_offset);
			pc_paths[ii].vSave(pf_dest);
			fprintf(pf_dest, "\n");
		}//if  ( (pc_paths[ii].i_fiber_node_start_offset == 18)&&(pc_paths[ii].i_fiber_node_end_offset == 25) )
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)



	fprintf(pf_dest, "\n");
	fprintf(pf_dest, "\n");
	/*for  (int  ii = 0; ii < pc_demands[5].pvGetPossiblePaths()->size(); ii++)
	{
		pc_demands[5].pvGetPossiblePaths()->at(ii)->vSaveWithLinks(pf_dest);
		fprintf(pf_dest, "\n");	
	}//for  (int  ii = 0; ii < pc_demands[5].pvGetPossiblePaths(); ii++)*/


	/*for  (int  ii = 0; ii < i_path_num; ii++)
	{
		if  (pc_paths[ii].v_path_connections.size()  ==  5)
		{
			fprintf(pf_dest, "%d    %d->%d\n",ii,  pc_paths[ii].i_fiber_node_start_offset, pc_paths[ii].i_fiber_node_end_offset);
			pc_paths[ii].vSaveWithLinks(pf_dest);
			fprintf(pf_dest, "\n");
			fprintf(pf_dest, "\n");
		}//if  (pc_paths[ii].v_path_connections.size()  ==  5)

	}//for  (int  ii = 0; ii < i_path_num; ii++)*/


	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::eSavePaths(CString  sPathFile)




CError  CFiberNet::e_save_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sSpecFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No path file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	


	for  (int  ii = 0; ii < iPathNum; ii++)
	{
		pcPaths[ii].vSaveSpec(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	fprintf(pf_dest, "\n");


	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::e_save_spec(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)




CError  CFiberNet::eSavePathsSpec(CString  sPathSpecFile)
{
	CError  c_err;
	
	c_err = e_save_spec(sPathSpecFile, pc_paths, i_path_num);

	return(c_err);
}//CError  CFiberNet::eSavePathsSpec(CString  sPathFile)


CError  CFiberNet::eSaveTreesSpec(CString  sTreeSpecFile)
{
	CError  c_err;
	
	c_err = e_save_spec(sTreeSpecFile, pc_trees, i_tree_num);

	return(c_err);
}//CError  CFiberNet::eSaveTreesSpec(CString  sTreeSpecFile)



CError  CFiberNet::e_save_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sSpecFile, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No reg file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)

	


	for  (int  ii = 0; ii < iPathNum; ii++)
	{
		pcPaths[ii].vSaveReg(pf_dest);
		fprintf(pf_dest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)

	fprintf(pf_dest, "\n");


	fclose(pf_dest);
	return(c_err);
}//CError  CFiberNet::e_save_reg(CString  sSpecFile, CFiberPath  *pcPaths, int  iPathNum)




CError  CFiberNet::eSaveTreesReg(CString  sTreeRegFile)
{
	CError  c_err;
	
	c_err = e_save_reg(sTreeRegFile, pc_paths, i_path_num);

	return(c_err);
}//CError  CFiberNet::eSaveTreesReg(CString  sTreeRegFile)




CError  CFiberNet::eSavePathsReg(CString  sPathRegFile)
{
	CError  c_err;
	
	c_err = e_save_reg(sPathRegFile, pc_paths, i_path_num);

	return(c_err);	
}//CError  CFiberNet::eSavePathsReg(CString  sPathRegFile)


//---------------------------class  CFiberNode--------------------------------------
CFiberNode::CFiberNode()
{
	pc_in_conns = NULL;
	pc_out_conns = NULL;
	i_node_number = 0;
}//CFiberNode::CFiberNode()


CFiberNode::~CFiberNode()
{
	if  (pc_in_conns != NULL)  delete  pc_in_conns;
	if  (pc_out_conns != NULL)  delete  pc_out_conns;
}//CFiberNode::~CFiberNode()


void  CFiberNode::vSetNodeNumber(int  iNodeNumber)
{
	if  (pc_in_conns != NULL)  delete  pc_in_conns;
	if  (pc_out_conns != NULL)  delete  pc_out_conns;

	i_node_number = iNodeNumber;

	pc_in_conns = new CFiberConnection*[i_node_number];
	pc_out_conns = new CFiberConnection*[i_node_number];

	for  (int  ii = 0; ii < i_node_number; ii++)  pc_in_conns[ii] = NULL;
	for  (int  ii = 0; ii < i_node_number; ii++)  pc_out_conns[ii] = NULL;
}//void  CFiberNode::vSetNodeNumber(int  iNodeNumber)



void  CFiberNode::vAddInConn(CFiberConnection  *pcConnection)
{
	if  (pcConnection->i_fiber_node_start_offset  <  i_node_number)
	{
		pc_in_conns[pcConnection->i_fiber_node_start_offset] = pcConnection;
	}//if  (pcConnection->i_fiber_node_start_offset  <  i_node_number)
}//void  CFiberNode::vAddInConn(CFiberConnection  *pcConnection)



void  CFiberNode::vAddOutConn(CFiberConnection  *pcConnection)
{
	if  (pcConnection->i_fiber_node_end_offset  <  i_node_number)
	{
		pc_out_conns[pcConnection->i_fiber_node_end_offset] = pcConnection;
	}//if  (pcConnection->i_fiber_node_start_offset  <  i_node_number)
}//void  CFiberNode::vAddOutConn(CFiberConnection  *pcConnection)






//---------------------------class  CFiberConnection--------------------------------------
CFiberConnection::CFiberConnection()
{
	i_connection_id = -1;
	i_fiber_node_start_offset = -1;
	i_fiber_node_end_offset = -1;
	pc_parent = NULL;
	d_connection_capacity = 0;
}//CFiberConnection::CFiberConnection(int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, CFiberNet *pcParent, double  dConnectionCapacity)



CFiberConnection::~CFiberConnection()
{

}//CFiberConnection::~CFiberConnection()



CError  CFiberConnection::eReserveSlots(int  iFirstSlot, int iSlotsNumber, bool  bErrorControl /*= true*/)
{
	CError  c_err;
	CString  s_buf;

	
	if  (bErrorControl  ==  true)
	{
		int  i_test;
		i_test = iGetLowestStartSlotForChannel(iFirstSlot, iSlotsNumber);
		if  (i_test  != iFirstSlot)
		{
			s_buf.Format("unable to settle channel in node connection %d at start slot: %d width: %d", i_connection_id, iFirstSlot, iSlotsNumber);			
			c_err.vSetError(s_buf);

			return(c_err);
		}//if  (i_test  != iFirstSlot)
	}//if  (bErrorControl  ==  true)

	bool  b_glue_with_after, b_glue_with_before;
	
	for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)
	{
		if  (v_slots_in_use.at(ii*2 + 1)  > iFirstSlot)
		{
			
			b_glue_with_after = false;
			b_glue_with_before = false;

			if  (iFirstSlot + iSlotsNumber - 1  ==  v_slots_in_use.at(ii*2) - 1)  b_glue_with_after = true;
			if  (ii > 0)
			{
				if  (iFirstSlot ==  v_slots_in_use.at((ii-1)*2 + 1) + 1)  b_glue_with_before = true;
			}//if  (ii > 0)


			if  ( (b_glue_with_before  ==  true)&&(b_glue_with_after ==  true) )
			{
				v_slots_in_use.erase(v_slots_in_use.begin() + ii*2 - 1);
				v_slots_in_use.erase(v_slots_in_use.begin() + ii*2 - 1);
				return(c_err);
			}//if  ( (b_glue_with_before  ==  true)&&(b_glue_with_after ==  true) )


			if  ( (b_glue_with_before  ==  true)&&(b_glue_with_after ==  false) )
			{
				v_slots_in_use.at((ii-1)*2 + 1) = iFirstSlot + iSlotsNumber - 1;
				return(c_err);
			}//if  ( (b_glue_with_before  ==  true)&&(b_glue_with_after ==  false) )


			if  ( (b_glue_with_before  ==  false)&&(b_glue_with_after ==  true) )
			{
				v_slots_in_use.at(ii*2) = iFirstSlot;
				return(c_err);
			}//if  ( (b_glue_with_before  ==  false)&&(b_glue_with_after ==  true) )


			if  ( (b_glue_with_before  ==  false)&&(b_glue_with_after ==  false) )
			{
				v_slots_in_use.insert(v_slots_in_use.begin() + ii*2,0);
				v_slots_in_use.insert(v_slots_in_use.begin() + ii*2,0);

				v_slots_in_use.at(ii * 2) = iFirstSlot;
				v_slots_in_use.at(ii * 2 + 1) = iFirstSlot + iSlotsNumber - 1;

				return(c_err);
			}//if  ( (b_glue_with_before  ==  false)&&(b_glue_with_after ==  false) )

		}//if  (v_slots_in_use.at(ii*2 + 1)  >= i_first_slot_allowed)
	}//for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)

	b_glue_with_before = false;
	if  (v_slots_in_use.size() > 0)
	{
		if  (iFirstSlot ==  v_slots_in_use.at(v_slots_in_use.size() - 1) + 1)  b_glue_with_before = true;
	}//if  (v_slots_in_use.size() > 0)


	if  (b_glue_with_before  ==  true)
	{
		v_slots_in_use.at(v_slots_in_use.size() - 1) = iFirstSlot + iSlotsNumber - 1;
		return(c_err);
	}//if  (b_glue_with_before  ==  true)
	else
	{
		v_slots_in_use.push_back(iFirstSlot);
		v_slots_in_use.push_back(iFirstSlot + iSlotsNumber - 1);
		return(c_err);
	}//else  if  (b_glue_with_before  ==  true)

	return(c_err);
}//CError  CFiberConnection::eReserveSlots(int  iFirstSlot, int iSlotsNumber)



int  CFiberConnection::iGetMaxSlot()
{
	if  (v_slots_in_use.size()  ==  0)  return(0);
	return(v_slots_in_use.at(v_slots_in_use.size() - 1));
}//int  CFiberConnection::iGetMaxSlot()




bool  CFiberConnection::bSlotInUse(int iSlotNumber)
{
	for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)
	{
		if  (
			(v_slots_in_use.at(ii * 2) <= iSlotNumber)&&
			iSlotNumber <= (v_slots_in_use.at(ii * 2 + 1))
			)
			return(true);
	}//for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)

	return(false);
}//bool  CFiberConnection::bSlotInUse(int iSlotNumber)



void  CFiberConnection::vSaveState(FILE  *pfDest, int iMaxSlotNumber)
{
	fprintf(pfDest, "edge %.3d    ", i_connection_id);

	for  (int  ii = 0; ii < iMaxSlotNumber; ii++)
	{
		if  (bSlotInUse(ii)  ==  true)
			fprintf(pfDest, "1   ");
		else
			fprintf(pfDest, "0   ");
	}//for  (int  ii = 0; ii < iMaxSlotNumber; ii++)
}//void  CFiberConnection::vSaveState(FILE  *pfDest, int iMaxSlotNumber)



void  CFiberConnection::vSaveStatePRW(FILE  *pfDest)
{
	fprintf(pfDest, "edge %.3d    ", i_connection_id);

	for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)
	{
		fprintf(pfDest, " (%d->%d) ", v_slots_in_use.at(ii*2), v_slots_in_use.at(ii*2 + 1));
	}//for  (int  ii = 0; ii < iMaxSlotNumber; ii++)
}//void  CFiberConnection::vSaveStatePRW(FILE  *pfDest)




int  CFiberConnection::iGetNumberOfFreeSlotsBelow(double  dMax)
{
	int  i_result = (int) dMax;
	int  i_buf;

	for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)
	{
		i_buf = v_slots_in_use.at(ii*2 + 1) - v_slots_in_use.at(ii*2);
		i_result = i_result - i_buf;
	}//for  (int  ii = 0; ii < v_slots_in_use.size(); ii++)

	return(i_result);
}//int  CFiberConnection::iGetNumberOfFreeSlotsBelow(double  dMax)




int  CFiberConnection::iGetLowestStartSlotForChannel(int  iFirstSlot, int iSlotsNumber)
{
	int  i_first_slot_allowed = iFirstSlot;


	for  (int  ii = 0; ii < v_slots_in_use.size()/2; ii++)
	{
		if  (v_slots_in_use.at(ii*2 + 1)  >= i_first_slot_allowed)
		{
			if  (v_slots_in_use.at(ii*2) <= i_first_slot_allowed + iSlotsNumber - 1)
			{
				i_first_slot_allowed = v_slots_in_use.at(ii*2 + 1) + 1;
			}//if  (v_slots_in_use.at(ii*2) <= i_first_slot_allowed + iSlotsNumber)
			else
			{
				return(i_first_slot_allowed);
			}//else  if  (v_slots_in_use.at(ii*2) <= i_first_slot_allowed + iSlotsNumber)

		}//if  (v_slots_in_use.at(ii*2 + 1)  >= i_first_slot_allowed)
	}//for  (int  ii = 0; ii < v_slots_in_use.size(); ii++)

	return(i_first_slot_allowed);
}//int  CFiberConnection::iGetLowestStartSlotForChannel(int  iFirstSlot, int iSlotsNumber)


void  CFiberConnection::vReset()
{
	v_slots_in_use.clear();
}//void  CFiberConnection::vReset()


CError  CFiberConnection::eConfig(int  iConectionId, int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, CFiberNet *pcParent, double  dConnectionCapacity)
{
	CError c_err;

	CString  s_buf;

	i_connection_id = iConectionId;
	i_fiber_node_start_offset = iFiberNodeStartOffset;
	i_fiber_node_end_offset = iFiberNodeEndOffset;
	pc_parent = pcParent;
	d_connection_capacity = dConnectionCapacity;

	if  (
		(i_fiber_node_start_offset  <  0)||(i_fiber_node_start_offset  >=  pc_parent->i_nodes_num)||
		(i_fiber_node_end_offset  <  0)||(i_fiber_node_end_offset  >=  pc_parent->i_nodes_num)
		)
	{
		s_buf.Format("%d connection - wrong start/end node", i_connection_id);	
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (	


	return(c_err);
}//void  CFiberConnection::vConfig(int  iConectionId, int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, CFiberNet *pcParent, double  dConnectionCapacity)





//---------------------------class  CFiberDemand--------------------------------------
CFiberDemand::CFiberDemand()
{
	i_demand_id = -1;
	i_fiber_node_start_offset = -1;
	i_fiber_node_end_offset = -1;
	pc_parent = NULL;
	d_demand_capacity = 0;
	i_demand_type = DEMAND_TYPE_UNICAST;
	pc_anycast_twin_dfemand = NULL;
	d_path_fit_sum = 0;
}//CFiberDemand::CFiberDemand()
		


CFiberDemand::~CFiberDemand()
{

}//CFiberDemand::~CFiberDemand()



void  CFiberDemand::v_copy_offsets(vector  <int> *pvDest, vector<int>  *pvSource)
{
	for  (int  ii  = 0; ii < pvSource->size(); ii++)
		pvDest->push_back(pvSource->at(ii));
}//void  CFiberDemand::v_copy_offsets(vector  <int> *pvDest, vector<int>  *pvSource)


CError  CFiberDemand::eConfig
	(
	int  iDemandId, 
	int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, double  dDemandCapacity, 
	CFiberNet *pcParent,
	int  iDemandType,
	CFiberDemand *pcAnycastTwinDemand,
	int  iMaxPathTreesAvailable
	)
{
	CError  c_err;

	CString  s_buf;

	i_demand_id = iDemandId;
	pc_parent = pcParent;
	d_demand_capacity = dDemandCapacity;
	pc_anycast_twin_dfemand = pcAnycastTwinDemand;
	i_demand_type = iDemandType;

	i_fiber_node_start_offset = iFiberNodeStartOffset;
	i_fiber_node_end_offset = iFiberNodeEndOffset;

	if  (i_demand_type  ==  DEMAND_TYPE_UNICAST)
	{
		if  (
			(i_fiber_node_start_offset  <  0)||(i_fiber_node_start_offset  >=  pc_parent->i_nodes_num)||
			(i_fiber_node_end_offset  <  0)||(i_fiber_node_end_offset  >=  pc_parent->i_nodes_num)
			)
		{
			s_buf.Format("%d demand - wrong start/end node", i_demand_id);	
			c_err.vSetError(s_buf);
			return(c_err);
		}//if  (	
	}//if  (i_demand_type  ==  DEMAND_TYPE_UNICAST)
	

	if  (i_demand_type  ==  DEMAND_TYPE_ANYCAST)
	{
		if  (
			(i_fiber_node_start_offset  <  -1)||(i_fiber_node_start_offset  >=  pc_parent->i_nodes_num)||
			(i_fiber_node_end_offset  <  -1)||(i_fiber_node_end_offset  >=  pc_parent->i_nodes_num)
			)
		{
			s_buf.Format("%d demand - wrong start/end node", i_demand_id);	
			c_err.vSetError(s_buf);
			return(c_err);
		}//if  (		
	}//if  (i_demand_type  ==  DEMAND_TYPE_ANYCAST)



	if  (i_demand_type  ==  DEMAND_TYPE_MULTICAST)
	{
		if  (
			(i_fiber_node_start_offset  <  0)||(i_fiber_node_start_offset  >=  pc_parent->i_nodes_num)
			)
		{
			s_buf.Format("%d demand - wrong tree root node", i_demand_id);	
			c_err.vSetError(s_buf);
			return(c_err);
		}//if  (		
	}//if  (i_demand_type  ==  DEMAND_TYPE_ANYCAST)





	v_possible_paths.clear();
	vector <int>  *pv_poss_paths_offsets;
	vector <int>  v_anycast_paths_buf;
	pv_poss_paths_offsets = NULL;

	if  (i_demand_type  ==  DEMAND_TYPE_UNICAST)
	{
		pv_poss_paths_offsets = &(pc_parent->pv_paths_between_nodes[i_fiber_node_start_offset][i_fiber_node_end_offset]);
	}//if  (i_demand_type  ==  DEMAND_TYPE_UNICAST)


	if  (i_demand_type  ==  DEMAND_TYPE_MULTICAST)
	{
		CFiberPath  *pc_poss_tree;

		//first find demand offset
		bool  b_found = false;
		for  (int  i_dem_offset = 0; (i_dem_offset < pc_parent->i_multi_demands)&&(b_found == false); i_dem_offset++)
		{
			if  (pc_parent->pc_demands[pc_parent->i_uni_demands + i_dem_offset].i_demand_id  ==  i_demand_id)
			{
				for  (int i_tree = 0; i_tree < 1000; i_tree++)
				{
					pc_poss_tree = &(pc_parent->pc_trees[i_dem_offset * 1000 + i_tree]);

					//if  (i_fiber_node_start_offset  ==  pc_poss_tree->iGetStartNode())
					{
						if  (
							(iMaxPathTreesAvailable  ==  -1)||
							(iMaxPathTreesAvailable  >  v_possible_paths.size())
							)
						{
							pc_poss_tree->vSetActive(true);
							v_possible_paths.push_back(pc_poss_tree);
						}//if  (
						
					}//if  (i_fiber_node_start_offset  ==  pc_poss_tree->iGetStartNode())
				
				}//for  (int i_tree = 0; i_tree < iMaxPathTreesAvailable; i_tree++)
			}//if  (pc_parent->pc_demands[pc_parent->i_uni_demands + i_dem_offset]  ==  this)
		}//for  (int  i_dem_offset = 0; (i_dem_offset < pc_parent->i_multi_demands)&&(b_found == false); i_dem_offset++)

		//for multicast we initialize exactly here and do nothing else
		/*for  (int  ii = 0; ii < pc_parent->i_tree_num; ii++)
		{
			pc_poss_tree = &(pc_parent->pc_trees[ii]);

			if  (i_fiber_node_start_offset  ==  pc_poss_tree->iGetStartNode())
			{
				if  (
					(iMaxPathTreesAvailable  ==  -1)||
					(iMaxPathTreesAvailable  >  v_possible_paths.size())
					)
				{
					pc_poss_tree->vSetActive(true);
					v_possible_paths.push_back(pc_poss_tree);
				}//if  (
				
			}//if  (i_fiber_node_start_offset  ==  pc_poss_tree->iGetStartNode())		
		}//for  (int  ii = 0; ii < pc_parent->i_tree_num; ii++)*/
	}//if  (i_demand_type  ==  DEMAND_TYPE_MULTICAST)


	if  (i_demand_type  ==  DEMAND_TYPE_ANYCAST)
	{
		if  (i_fiber_node_start_offset  ==  -1)
		{
			for  (int  ii = 0; ii < pc_parent->i_rep_num; ii++)
			{
				v_copy_offsets(&v_anycast_paths_buf, &(pc_parent->pv_paths_between_nodes[pc_parent->pc_replices[ii]->i_node_number][i_fiber_node_end_offset]));
			}//for  (int  ii = 0; ii < pc_parent->i_rep_num; ii++)
		}//if  (i_fiber_node_start_offset  ==  -1)


		if  (i_fiber_node_end_offset  ==  -1)
		{
			for  (int  ii = 0; ii < pc_parent->i_rep_num; ii++)
			{
				v_copy_offsets(&v_anycast_paths_buf, &(pc_parent->pv_paths_between_nodes[i_fiber_node_start_offset][pc_parent->pc_replices[ii]->i_node_number]));
			}//for  (int  ii = 0; ii < pc_parent->i_rep_num; ii++)
		}//if  (i_fiber_node_end_offset  ==  -1)
		
		pv_poss_paths_offsets = &v_anycast_paths_buf;	
	}//if  (i_demand_type  ==  DEMAND_TYPE_ANYCAST)


	if  (pv_poss_paths_offsets  !=  NULL)
	{
		CFiberPath  *pc_poss_path;
		for  (int  ii = 0; ii < pv_poss_paths_offsets->size(); ii++)
		{
			if  (
				(iMaxPathTreesAvailable  ==  -1)||
				(iMaxPathTreesAvailable  >  v_possible_paths.size())
				)
			{
				pc_poss_path = &(pc_parent->pc_paths[pv_poss_paths_offsets->at(ii)]);
				pc_poss_path->vSetActive(true);
				v_possible_paths.push_back(pc_poss_path);
			}//if  (
		}//for  (int  ii = 0; ii < pv_poss_paths_offsets->size(); ii++)
	}//if  (pv_poss_paths_offsets  !=  NULL)

	
	

	return(c_err);
}//CError  CFiberDemand::eConfig(int  iDemandId, int  iFiberNodeStartOffset, int  iFiberNodeEndOffset, double  dDemandCapacity, CFiberNet *pcParent)


bool  CFiberDemand::bDominatingTwin()
{
	if  (pc_anycast_twin_dfemand  ==  NULL)  return(true);

	if  (d_demand_capacity == pc_anycast_twin_dfemand->d_demand_capacity)
	{
		if  (i_demand_id  >  pc_anycast_twin_dfemand->i_demand_id)
			return(true);
		else
			return(false);
	}//if  (d_demand_capacity == pc_anycast_twin_dfemand->d_demand_capacity)

	if  (d_demand_capacity > pc_anycast_twin_dfemand->d_demand_capacity)
		return(true);
	else
		return(false);

	return(false);
}//bool  CFiberDemand::bDominatingTwin()



void  CFiberDemand::vSummarizePathFitness()
{
	d_path_fit_sum = 0;

	for  (int  ii = 0; ii < v_possible_paths.size(); ii++)
	{
		d_path_fit_sum +=  v_possible_paths.at(ii)->dGetPathFitness();
	}//for  (int  ii = 0; ii < v_possible_paths.size(); ii++)
}//void  CFiberDemand::vSummarizePathFitness()



void  CFiberDemand::vSave(FILE  *pfDest)
{
	if  (i_demand_type ==  DEMAND_TYPE_UNICAST)  fprintf(pfDest, "unicastdemand ");
	if  (i_demand_type ==  DEMAND_TYPE_ANYCAST)  fprintf(pfDest, "anycastdemand ");
	if  (i_demand_type ==  DEMAND_TYPE_MULTICAST)  fprintf(pfDest, "multicastdemand ");



	fprintf(pfDest, "%d\t", i_fiber_node_start_offset);
	fprintf(pfDest, "%d\t", i_fiber_node_end_offset);
	fprintf(pfDest, "%d\t", (int) d_demand_capacity);
}//void  CFiberDemand::vSave(FILE  *pfDest)



//---------------------------class  CFiberPath--------------------------------------
CFiberPath::CFiberPath()
{
	i_path_id = -1;
	i_fiber_node_start_offset = -1;
	i_fiber_node_end_offset = -1;

	b_active = false;
	b_tree = false;

	d_length = 0;
	d_mp = 1;
	d_delta_s = 6.25;

	d_new_transfer_multiplier = 0;

	pc_parent = NULL;
}//CFiberPath::CFiberPath()



CFiberPath::~CFiberPath()
{

}//CFiberPath::~CFiberPath()




void  CFiberPath::vConfigure(int  iPathId, CFiberNet *pcParent)
{
	i_path_id = iPathId;
	pc_parent = pcParent;
}//void  CFiberPath::vConfigure(int  iPathId, CFiberNet *pcParent)



bool  CFiberPath::bIsConnPartOfPath(CFiberConnection  *pcFiberConn)
{
	for  (int  ii = 0; ii < v_path_connections.size(); ii++)
	{
		if  (v_path_connections.at(ii)  ==  pcFiberConn)  return(true);
	}//for  (int  ii = 0; ii < v_path_connections.size(); ii++)

	return(false);
}//bool  bIsConnPartOfPath(CFiberConnection  *pcFiberConn)



void  CFiberPath::vSaveWithLinks(FILE  *pfDest)
{
	for  (int  ii = 0; ii < v_path_connections.size(); ii++)
	{
		fprintf(pfDest, "%d ", v_path_connections.at(ii)->i_connection_id);
	}//for  (int  ii = 0; ii < v_path_connections.size(); ii++)
}//void  CFiberPath::vSaveWithLinks(FILE  *pfDest)


void  CFiberPath::vSave(FILE  *pfDest)
{
	for  (int  ii = 0; ii < pc_parent->i_conns_num; ii++)
	{
		if  (ii  >  0)  fprintf(pfDest, " ");

		if  (bIsConnPartOfPath(&(pc_parent->pc_conns[ii]))  ==  true)
			fprintf(pfDest, "1");
		else
			fprintf(pfDest, "0");
	}//for  (int  ii = 0; ii < pc_parent->i_conns_num; ii++)
}//void  CFiberPath::vSave(FILE  *pfDest)



void  CFiberPath::vSaveInGenotype(int  *piGenotype)
{
	for  (int  ii = 0; ii < pc_parent->i_conns_num; ii++)
	{
		if  (bIsConnPartOfPath(&(pc_parent->pc_conns[ii]))  ==  true)  piGenotype[ii] = 1;			
	}//for  (int  ii = 0; ii < pc_parent->i_conns_num; ii++)
}//void  CFiberPath::vSave(FILE  *pfDest)



void  CFiberPath::vSaveMod(FILE  *pfDest)
{
	fprintf(pfDest, "%d\t%d\n", (int) d_length, (int) d_mp);
}//void  CFiberPath::vSaveMod(FILE  *pfDest)



void  CFiberPath::vSaveSpec(FILE  *pfDest)
{
	for  (int  ii = 0; ii < v_slot_number_for_transfer.size(); ii++)
	{
		fprintf(pfDest, "%d\t", v_slot_number_for_transfer.at(ii));
	}//for  (int  ii = 0; ii < v_slot_number_for_transfer.size(); ii++)

}//void  CFiberPath::vSaveSpec(FILE  *pfDest)


void  CFiberPath::vSaveReg(FILE  *pfDest)
{
	for  (int  ii = 0; ii < v_reg_number_for_transfer.size(); ii++)
	{
		fprintf(pfDest, "%d\t", v_reg_number_for_transfer.at(ii));
	}//for  (int  ii = 0; ii < v_slot_number_for_transfer.size(); ii++)
}//void  CFiberPath::vSaveReg(FILE  *pfDest)



CError  CFiberPath::eLoadMod(CString  sModLine)
{
	CError  c_err;

	int  i_index;
	i_index = 0;
	d_length = Tools::dExtractFromString(sModLine, i_index, &i_index);
	d_mp = Tools::dExtractFromString(sModLine, i_index, &i_index);

	return(c_err);
}//CError  CFiberPath::eLoadMod(CString  sModLine)


int  CFiberPath::iGetSlotsNumber(double  dDemandCapacity)
{
	if  (b_tree  ==  true)
	{
		if  (v_slot_number_for_transfer.size()  !=  1)
		{
			CString  s_buf;
			s_buf.Format("Number of slots in the spec file for tree different than 1 (%d)", v_slot_number_for_transfer.size());

			CError  c_err;
			c_err.vSetError(s_buf);

			return(0);
		}//if  (v_slot_number_for_transfer.size()  !=  1)

		return(v_slot_number_for_transfer.at(0));
	}//if  (b_tree  ==  true)

	double  d_slot_number_offset;
	d_slot_number_offset = dDemandCapacity / d_new_transfer_multiplier;

	int  i_slot_number_offset;
	i_slot_number_offset = (int) round(d_slot_number_offset);//::Math::Round(d_slot_number_offset);
	i_slot_number_offset--;//we start from zero

	if  (i_slot_number_offset  >=  v_slot_number_for_transfer.size())
	{
		CString  s_buf;
		s_buf.Format
			(
			"Number of slots in the spec file for tree too small (%d positions). Required for %.2lf transfer: %d slots",
			v_slot_number_for_transfer.size(), dDemandCapacity,
			i_slot_number_offset
			);

		CError  c_err;
		c_err.vSetError(s_buf);

		return(0);
	}//if  (v_slot_number_for_transfer.size()  !=  1)

	return(v_slot_number_for_transfer.at(i_slot_number_offset));
	


	//old slot computation
	double  d_slot_number;

	d_slot_number = dDemandCapacity / (d_mp * 4.0 * d_delta_s);
	int  i_slot_number;
	i_slot_number = (int) d_slot_number;
	if  (i_slot_number  <  d_slot_number)  i_slot_number++;
	i_slot_number *= 2;

	return(i_slot_number);
}//int  CFiberPath::iGetSlotsNumber(double  dDemandCapacity)


CError  CFiberPath::eSetPathSlots(int  iFirstSlot, int  iSlotsNumber)
{
	CError  c_err;

	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		c_err  = v_path_connections.at(i_conn)->eReserveSlots(iFirstSlot, iSlotsNumber);
		if  (c_err)  return(c_err);
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	return(c_err);
}//CError  CFiberPath::eSetPathModules(int  iModulesNumber)


CError  CFiberPath::eFindPath(int  *piFirstSlot, int  iSlotsNumber)
{
	CError  c_err;

	int  i_current_start_slot = 0;
	int  i_proposed_start_slot;


	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		i_proposed_start_slot  = v_path_connections.at(i_conn)->iGetLowestStartSlotForChannel(i_current_start_slot, iSlotsNumber);
		if  (i_proposed_start_slot  !=  i_current_start_slot)
		{
			i_conn = -1;
			i_current_start_slot = i_proposed_start_slot;
		}//if  (i_proposed_start_slot  !=  i_current_start_slot)
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	*piFirstSlot = i_current_start_slot;

	return(c_err);
}//CError  CFiberPath::eFindPath(int  *piFirstSlot, int  iSlotsNumber)



CError  CFiberPath::eSetPath(int  iFirstSlot, int  iSlotsNumber)
{
	CError  c_err;

	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		c_err  = v_path_connections.at(i_conn)->eReserveSlots(iFirstSlot, iSlotsNumber);
		if  (c_err)  
		{
			c_err.vSetError("v_path_connections.at(i_conn)->eSetPath");
			return(c_err);
		}//if  (c_err)  
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	return(c_err);
}//CError  CFiberPath::eSetPath(int  iFirstSlot, int  iSlotsNumber)



CError  CFiberPath::eFindAndSetPathSlots(int  *piFirstSlot, int  iSlotsNumber)
{
	CError  c_err;

	int  i_current_start_slot = 0;
	int  i_proposed_start_slot;


	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		i_proposed_start_slot  = v_path_connections.at(i_conn)->iGetLowestStartSlotForChannel(i_current_start_slot, iSlotsNumber);
		if  (i_proposed_start_slot  !=  i_current_start_slot)
		{
			i_conn = -1;
			i_current_start_slot = i_proposed_start_slot;
		}//if  (i_proposed_start_slot  !=  i_current_start_slot)
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	*piFirstSlot = i_current_start_slot;

	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		c_err  = v_path_connections.at(i_conn)->eReserveSlots(i_current_start_slot, iSlotsNumber);
		if  (c_err)  
		{
			c_err.vSetError("v_path_connections.at(i_conn)->eReserveSlots");
			return(c_err);
		}//if  (c_err)  
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	return(c_err);
}//CError  CFiberPath::eFindAndSetPathSlots(int  *piFirstSlot, int  iSlotsNumber)



CError  CFiberPath::eSetPathCapacity(int  iFirstSlot, double  dDemandCapacity)
{
	CError  c_err;

	int  i_modules = iGetSlotsNumber(dDemandCapacity);
	c_err  = eSetPathSlots(iFirstSlot, i_modules);

	return(c_err);
}//CError  CFiberPath::eSetPathCapacity(double  dDemandCapacity)



void  CFiberPath::vSetPathFitness(double  dSolutionValue)
{
	//if the path is not active we dont comute its fitness
	if  (b_active  ==  false)
	{
		d_path_fitness = 0;
		return;
	}//if  (b_active  ==  false)

	int  i_min_slots;
	int  i_buf;

	i_min_slots = 0;
	for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)
	{
		i_buf  =  v_path_connections.at(i_conn)->iGetNumberOfFreeSlotsBelow(dSolutionValue);
		if  (i_conn == 0)  i_min_slots = i_buf;
		if  (i_buf  <  i_min_slots)  i_min_slots = i_buf;		
	}//for  (int  i_conn = 0; i_conn < v_path_connections.size(); i_conn++)

	d_path_fitness = i_min_slots;
	if  (d_path_fitness  <  1)  d_path_fitness  =  1;
}//void  CFiberPath::vSetPathWeight(double  dSolutionValue)




CError  CFiberPath::eLoad(CString  sPathLine, bool  bControlConnectionsOrder)
{
	CError  c_err;

	int  i_index;
	i_fiber_node_start_offset = -1;
	i_fiber_node_start_offset = -1;
	v_path_connections.clear();
	i_index = 0;
	int  i_conn_in_path;
	vector<CFiberConnection *>  v_path_connections_buf;

	
	for  (int i_conn_offset = 0; i_conn_offset < pc_parent->i_conns_num; i_conn_offset++)
	{
		
		i_conn_in_path = Tools::iExtractFromString(sPathLine, i_index, &i_index);
		if  ( (i_conn_offset  <  0)||(pc_parent->i_conns_num <= i_conn_offset) )
		{
			CString  s_buf;
			s_buf.Format("Connection %d offset in path %d outside borders", i_conn_offset, i_path_id);
			c_err.vSetError(s_buf);
			return(c_err);
		}//if  ( (i_conn_offset  <  0)||(pc_parent->i_conns_num <= i_conn_offset) )

		if  (i_conn_in_path == 1)  v_path_connections_buf.push_back(&(pc_parent->pc_conns[i_conn_offset]));
		
	}//for  (int i_node_start = 0; i_node_start < i_nodes_num; i_node_start++)



	int  i_conn_offset;
	CFiberConnection  *pc_conn_buf;

	i_conn_offset = i_get_first_connection_offset(&v_path_connections_buf);
	if  (i_conn_offset  < 0)
	{
		CString  s_err;
		s_err.Format("First connection not found. Path: %d", i_path_id);
		c_err.vSetError(s_err);
		return(c_err);
	}//if  (i_conn_offset  < 0)

	pc_conn_buf = v_path_connections_buf.at(i_conn_offset);
	v_path_connections_buf.erase(v_path_connections_buf.begin() + i_conn_offset);	
	v_path_connections.push_back(pc_conn_buf);

	if  (bControlConnectionsOrder  ==  true)
	{
		while  ( (v_path_connections_buf.size()  >  0)&&(i_conn_offset >= 0) )
		{
			i_conn_offset = i_get_linked_connection_offset(&v_path_connections_buf, pc_conn_buf);
			if  (i_conn_offset  < 0)
			{
				CString  s_err;
				s_err.Format("Next connection not found. Path: %d", i_path_id);
				c_err.vSetError(s_err);
				return(c_err);
			}//if  (i_conn_offset  < 0)

			pc_conn_buf = v_path_connections_buf.at(i_conn_offset);
			v_path_connections_buf.erase(v_path_connections_buf.begin() + i_conn_offset);	
			v_path_connections.push_back(pc_conn_buf);
		}//while  ( (v_path_connections_buf.size()  >  0)&&(i_conn_offset > 0) )

		i_fiber_node_start_offset = v_path_connections.at(0)->i_fiber_node_start_offset;
		i_fiber_node_end_offset = v_path_connections.at(v_path_connections.size() - 1)->i_fiber_node_end_offset;
	}//if  (bControlConnectionsOrder  ==  true)
	else
	{
		while  (v_path_connections_buf.size()  >  0)
		{
			v_path_connections.push_back(v_path_connections_buf.at(0));
			v_path_connections_buf.erase(v_path_connections_buf.begin());			
		}//while  (v_path_connections_buf.size()  >  0)

		i_fiber_node_start_offset = v_path_connections.at(0)->i_fiber_node_start_offset;
		i_fiber_node_end_offset = -1;
	}//else  if  (bControlConnectionsOrder  ==  true)

	


	return(c_err);
}//CError  CFiberPath::eLoad(CString  sPathLine)



CError  CFiberPath::eGetSlotsNumberForTransfer(CString  sLine)
{
	CError  c_err;

	v_slot_number_for_transfer.clear();
	int  i_index = 0;
	int  i_buf;
	while  (i_index  <  sLine.GetLength())
	{
		i_buf = Tools::iExtractFromString(sLine, i_index, &i_index);
		v_slot_number_for_transfer.push_back(i_buf);
	}//while  (i_index  <  sLine.GetLength())

	
	return(c_err);
}//CError  CFiberPath::eGetSlotsNumberForTransfer(CString  sLine)



CError  CFiberPath::eGetRegNumberForTransfer(CString  sLine)
{
	CError  c_err;

	v_reg_number_for_transfer.clear();
	int  i_index = 0;
	int  i_buf;
	while  (i_index  <  sLine.GetLength())
	{
		i_buf = Tools::iExtractFromString(sLine, i_index, &i_index);
		v_reg_number_for_transfer.push_back(i_buf);
	}//while  (i_index  <  sLine.GetLength())

	
	return(c_err);
}//CError  CFiberPath::eGetRegNumberForTransfer(CString  sLine)




int  CFiberPath::i_get_first_connection_offset(vector<CFiberConnection *>  *pvConnections)
{
	bool  b_prev_conn_found;


	for  (int  i_checked = 0; i_checked < pvConnections->size(); i_checked++)
	{
		b_prev_conn_found = false;
		
		for  (int  ii = 0; (ii < pvConnections->size())&&(b_prev_conn_found == false); ii++)
		{
			if  (
				pvConnections->at(i_checked)->i_fiber_node_start_offset
				==
				pvConnections->at(ii)->i_fiber_node_end_offset
				)
				b_prev_conn_found = true;
		}//for  (int  ii = 0; (ii < pvConnections->size())&&(b_prev_conn_found == true); ii++)
	
		if  (b_prev_conn_found == false)  return(i_checked);
	}//for  (int  i_checked = 0; i_checked < pvConnections->size(); i_checked++)


	return(-1);
}//int  CFiberPath::i_get_first_connection_offset(vector<CFiberConnection *>  vConnections)



int  CFiberPath::i_get_linked_connection_offset(vector<CFiberConnection *>  *pvConnections, CFiberConnection *pcMotherConnection)
{
	for  (int  ii = 0; ii < pvConnections->size(); ii++)
	{
		if  (
			pcMotherConnection->i_fiber_node_end_offset
			==
			pvConnections->at(ii)->i_fiber_node_start_offset
			)
			return(ii);
	}//for  (int  ii = 0; (ii < pvConnections->size())&&(b_prev_conn_found == true); ii++)

	return(-1);
}//int  CFiberPath::i_get_linked_connection_offset(vector<CFiberConnection *>  vConnections, CFiberConnection *pcMotherConnection)





//---------------------------class  CFiberSolution--------------------------------------
CFiberSolution::CFiberSolution(CFiberNet  *pcFibNet)
{
	pc_channels = NULL;
	pc_fib_net = pcFibNet;

	if  (pc_fib_net->i_demand_num  <= 0)  return;
	pc_channels = new CFiberSolutionChannel[pc_fib_net->i_demand_num];

	for  (int  ii = 0; ii < pc_fib_net->i_demand_num; ii++)
		pc_channels[ii].vConfig(this, &(pc_fib_net->pc_demands[ii]));

	d_last_raw_sol_value = 1;
}//CFiberSolution::CFiberSolution(CFiberNet  *pcFibNet)



CFiberSolution::CFiberSolution(CFiberSolution  *pcFibSol)
{
	pc_channels = NULL;
	pc_fib_net = pcFibSol->pc_fib_net;

	if  (pc_fib_net->i_demand_num  <= 0)  return;
	pc_channels = new CFiberSolutionChannel[pc_fib_net->i_demand_num];

	for  (int  ii = 0; ii < pc_fib_net->i_demand_num; ii++)
	{
		pc_channels[ii].vConfig(this, &(pc_fib_net->pc_demands[ii]));
		pc_channels[ii].iFirstSlotNumber = pcFibSol->pc_channels[ii].iFirstSlotNumber;
		pc_channels[ii].pcPath = pcFibSol->pc_channels[ii].pcPath;
	}//for  (int  ii = 0; ii < pc_fib_net->i_demand_num; ii++)

}//CFiberSolution::CFiberSolution(CFiberSolution  *pcFibNet)



CFiberSolution::~CFiberSolution()
{
	if  (pc_channels  !=  NULL)  delete  pc_channels;
}//CFiberSolution::~CFiberSolution()



void  CFiberSolution::vCopyChannelsState(CFiberSolutionChannel  *pcChannelsDest)
{
	for  (int ii = 0; ii < pc_fib_net->iGetDemandNumber(); ii++)
	{
		pc_channels[ii].vCopy(&(pcChannelsDest[ii]));
	}//for  (int ii = 0; ii < pc_fib_net->iGetDemandNumber(); ii++)

}//void  CFiberSolution::vCopyChannelsState(CFiberSolutionChannel  *pcChannelsDest)




CError  CFiberSolution::eFindChanel(int  iChanelOffset, CFiberPath  *pcPath, int *piStartSlot, int  *piSlotsNumber)
{
	CError  c_err;

	if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )
	{
		CString  s_buf;

		s_buf.Format("CFiberSolution::eFindAndSetChanel -> wrong chanel offset: %d allowed: %d", iChanelOffset, pc_fib_net->iGetDemandNumber());
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )

	c_err = pc_channels[iChanelOffset].eFindChanel(pcPath, piStartSlot, piSlotsNumber);

	return(c_err);
}//CError  CFiberSolution::eFindAndSetChanel(int  iChanelOffset, CFiberPath  *pcPath)




CError  CFiberSolution::eSetChanel(int  iChanelOffset, CFiberPath  *pcPath, int iStartSlot, int  iSlotsNumber)
{
	CError  c_err;

	if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )
	{
		CString  s_buf;

		s_buf.Format("CFiberSolution::eFindAndSetChanel -> wrong chanel offset: %d allowed: %d", iChanelOffset, pc_fib_net->iGetDemandNumber());
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )

	c_err = pc_channels[iChanelOffset].eSetChanel(pcPath, iStartSlot, iSlotsNumber);

	return(c_err);
}//CError  CFiberSolution::eFindAndSetChanel(int  iChanelOffset, CFiberPath  *pcPath)





CError  CFiberSolution::eFindAndSetChanel(int  iChanelOffset, CFiberPath  *pcPath)
{
	CError  c_err;

	if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )
	{
		CString  s_buf;

		s_buf.Format("CFiberSolution::eFindAndSetChanel -> wrong chanel offset: %d allowed: %d", iChanelOffset, pc_fib_net->iGetDemandNumber());
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  ( (iChanelOffset  <  0)||(iChanelOffset  >=  pc_fib_net->iGetDemandNumber()) )

	c_err = pc_channels[iChanelOffset].eFindAndSetChanel(pcPath);

	return(c_err);
}//CError  CFiberSolution::eFindAndSetChanel(int  iChanelOffset, CFiberPath  *pcPath)


void  CFiberSolution::vSetPathFitness()
{
	d_clever_paths_initializations++;

	CFiberPath  *pc_all_paths;

	pc_all_paths = pc_fib_net->pc_paths;
	for  (int  ii = 0; ii < pc_fib_net->i_path_num; ii++)
	{
		pc_all_paths[ii].vSetPathFitness(d_last_raw_sol_value);
	}//for  (int  ii = 0; ii < pc_fib_net->i_path_num; ii++)


	int  i_dem_number = pc_fib_net->iGetDemandNumber();
	for  (int  ii = 0; ii < i_dem_number; ii++)
	{
		pc_fib_net->pcGetDemands()[ii].vSummarizePathFitness();
	}//for  (int  ii = 0; ii < i_dem_number; ii++)

}//void  CFiberSolution::vSetPathFitness()



//CFiberPath  *CFiberSolution::pcGetRandomPath(int  iDemandNumber, gcroot<Random*> pcRandomGen)
//{
//	vector  <CFiberPath  *>  *pv_possible_paths;
//	pv_possible_paths = (pc_fib_net->pcGetDemands()[iDemandNumber]).pvGetPossiblePaths();
//
//	if  (pv_possible_paths->size()  <=  0)  return(NULL);
//
//	int i_chosen_path; 
//	i_chosen_path = pcRandomGen->Next(pv_possible_paths->size());
//	if  (i_chosen_path  >=  pv_possible_paths->size())  i_chosen_path  =  pv_possible_paths->size() - 1;
//
//	return(pv_possible_paths->at(i_chosen_path));
//}//CFiberPath  *CFiberSolution::pcGetRandomPath(int  iDemandNumber, gcroot<Random*> pcRandomGen)


//CFiberPath  *CFiberSolution::pcGetCleverPath(int  iDemandNumber, gcroot<Random*> pcRandomGen)
//{
//	d_clever_paths_evaluations++;
//
//
//	double  d_fit_sum;
//
//	vector<CFiberPath*>  *pv_possible_paths = NULL;
//	pv_possible_paths = (pc_fib_net->pcGetDemands()[iDemandNumber]).pvGetPossiblePaths();
//	d_fit_sum = (pc_fib_net->pcGetDemands()[iDemandNumber]).d_path_fit_sum;
//
//	if  (d_fit_sum < 1)  return(pcGetRandomPath(iDemandNumber, pcRandomGen));
//
//	double  d_random;
//	d_random = pcRandomGen->NextDouble();
//	d_random *= d_fit_sum;
//
//	return(pc_get_clever_path(pv_possible_paths, d_fit_sum, pcRandomGen));
//}//CFiberPath  *CFiberSolution::pcGetCleverPath(int  iDemandNumber, gcroot<Random*> pcRandomGen)




//CFiberPath  *CFiberSolution::pcGetBestCleverPath(int  iStartNode, int iEndNode)
//{
//	double  d_current_fit;
//	double  d_best_fit;
//	int  i_best_fit_offset;
//
//	vector<int>  *pv_possible_paths_offsets = NULL;
//	pv_possible_paths_offsets = &(pc_fib_net->pv_paths_between_nodes[iStartNode][iEndNode]);
//
//	CFiberPath  *pc_all_paths = NULL;
//	pc_all_paths = pc_fib_net->pc_paths;
//
//	d_best_fit = pc_all_paths[pv_possible_paths_offsets->at(0)].dGetPathFitness();
//	i_best_fit_offset = pv_possible_paths_offsets->at(0);
//	for  (int  ii = 1; ii < pv_possible_paths_offsets->size(); ii++)
//	{
//		d_current_fit = pc_all_paths[pv_possible_paths_offsets->at(ii)].dGetPathFitness();
//		if  (d_best_fit  <  d_current_fit)
//		{
//			d_best_fit = d_current_fit;
//			i_best_fit_offset = pv_possible_paths_offsets->at(ii);
//		}//if  (d_best_fit  <  d_fit_sum)
//	}//for  (int  ii = 1; ii < pv_possible_paths->size(); ii++)
//		
//	return(&(pc_all_paths[i_best_fit_offset]));
//}//CFiberPath  *CFiberSolution::pcGetBestCleverPath(int  iStartNode, int iEndNode)





//CFiberPath  *CFiberSolution::pc_get_clever_path(vector<CFiberPath*>  *pvPossiblePaths, double dPathFitSumm, gcroot<Random*> pcRandomGen)
//{
//	double  d_random;
//	d_random = pcRandomGen->NextDouble();
//	d_random *= dPathFitSumm;
//
//	double  d_one_after_one_sum;
//	d_one_after_one_sum = 0;
//	for  (int  ii = 0; ii < pvPossiblePaths->size(); ii++)
//	{
//		d_one_after_one_sum += pvPossiblePaths->at(ii)->dGetPathFitness();
//		if  (d_random < d_one_after_one_sum)  return(pvPossiblePaths->at(ii));
//	}//for  (int  ii = 0; ii < pv_possible_paths->size(); ii++)
//
//	return(pvPossiblePaths->at(pvPossiblePaths->size() - 1));
//}//CFiberPath  *CFiberSolution::pc_get_clever_path(vector<CFiberPath*>  pcPossiblePaths, double dPathFitSumm, gcroot<Random*> pcRandomGen)




CError  CFiberSolution::eApplySolution()
{
	CError  c_err;

	pc_fib_net->eResetNet();

	for  (int  ii = 0; ii < pc_fib_net->i_demand_num; ii++)
	{
		c_err = pc_channels[ii].eSetChannel();
		if  (c_err)  return(c_err);
	}//for  (int  ii = 0; ii < pc_fib_net->i_demand_num; ii++)

	return(c_err);
}//CError  CFiberSolution::eApplySolution()



CError  CFiberSolution::eLoadSolution(CString  sSource)
{
	CError  c_err;

	FILE  *pf_source;
	pf_source = fopen(sSource, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No solf file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err = eLoadSolution(pf_source);
	fclose(pf_source);

	if  (c_err)  return(c_err);
	c_err = eApplySolution();

	return(c_err);
}//CError  CFiberSolution::eLoadSolution(CString  sSource)



CError  CFiberSolution::eLoadSolutionWalk(CString  sSource)
{
	CError  c_err;

	FILE  *pf_source;
	pf_source = fopen(sSource, "r+");

	if  (pf_source  ==  NULL)
	{
		CString  s_buf;

		s_buf.Format("No solf file defined (%s)", sSource);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err = eLoadSolutionWalk(pf_source);
	fclose(pf_source);


	if  (c_err)  return(c_err);
	::vReportInFile("ai_solution.txt", "", true);
	c_err = eApplySolution();
	eSaveSolution("zz_filled_what_i_can.txt");


	return(c_err);
}//CError  CFiberSolution::eLoadSolutionWalk(CString  sSource)




CError  CFiberSolution::eLoadSolutionAibin(CString  sSource)
{
	CError  c_err;

	FILE  *pf_source;
	pf_source = fopen(sSource, "r+");

	if  (pf_source  ==  NULL)
	{
		CString  s_buf;

		s_buf.Format("No solf file defined (%s)", sSource);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err = eLoadSolutionAibin(pf_source);
	fclose(pf_source);


	if  (c_err)  return(c_err);
	::vReportInFile("ai_solution.txt", "", true);
	c_err = eApplySolution();
	eSaveSolution("zz_filled_what_i_can.txt");


	return(c_err);
}//CError  CFiberSolution::eLoadSolutionAibin(CString  sSource)

		


CError  CFiberSolution::eLoadSolutionAFA(CString  sSource)
{
	CError  c_err;

	FILE  *pf_source;
	pf_source = fopen(sSource, "r+");

	if  (pf_source  ==  NULL)
	{
		CString  s_buf;

		s_buf.Format("No solf file defined (%s)", sSource);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err = eLoadSolutionAFA(pf_source);
	fclose(pf_source);

	if  (c_err)  return(c_err);
	c_err = eApplySolution();

	return(c_err);
}//CError  CFiberSolution::eLoadSolution(CString  sSource)





CError  CFiberSolution::eLoadSolutionPRW(CString  sSource)
{
	CError  c_err;

	FILE  *pf_source;
	pf_source = fopen(sSource, "r+");

	if  (pf_source  ==  NULL)
	{
		c_err.vSetError("No prw solf file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err = eLoadSolutionPRW(pf_source);
	
	fclose(pf_source);
	return(c_err);
}//CError  CFiberSolution::eLoadSolutionPRW(CString  sSource)


CError  CFiberSolution::eLoadSolutionPRW(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;
    CString  s_line  = "";
	
	int  i_demand_offset = 0;

	int  i_marker;
	int  i_index;

	int  i_start_node, i_end_node;
	int  i_path_number;
	int  i_start_slot_number;

	CString  s_marker = SOLUTION_PRW_DEMAND_MARKER;


	while  ( (i_demand_offset  <  pc_fib_net->i_demand_num)&&(!feof(pfSource)) )
	{
		s_line  =  Tools::sReadLine(pfSource);

		i_marker = s_line.Find(SOLUTION_PRW_DEMAND_MARKER);

		if  (i_marker  ==  0)
		{
			i_index = s_marker.GetLength();

			i_start_node = Tools::iExtractFromString(s_line, i_index, &i_index);//the first number is unimportant
			i_start_node = Tools::iExtractFromString(s_line, i_index, &i_index);
			i_end_node = Tools::iExtractFromString(s_line, i_index, &i_index);
			i_start_slot_number = Tools::iExtractFromString(s_line, i_index, &i_index);

			i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);//this is width in slots
			i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);

			if  ( (i_path_number < 0)||(i_path_number >= pc_fib_net->i_path_num) )
			{
				s_buf.Format("Wrong path offset (%d) at demand %d", i_path_number, i_demand_offset);
				c_err.vSetError(s_buf);
				return(c_err);
			}//if  ( (i_path_number < 0)||(i_path_number >= pc_fib_net->i_path_num) )

			pc_channels[i_demand_offset].pcPath = &(pc_fib_net->pc_paths[i_path_number]);
			pc_channels[i_demand_offset].iFirstSlotNumber = i_start_slot_number;

			i_demand_offset++;
		}//if  (i_marker  ==  0)

	}//while  ( (i_demands_read  <  pc_fib_net->i_demand_num)&&(feof(pfSource)) )


	
	return(c_err);
}//CError  CFiberSolution::eLoadSolution(FILE  *pfSource)



CError  CFiberSolution::eLoadSolutionWalk(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf, s_slot_info, s_buf2;
	CString  s_line, s_setting_name;

	int  i_index;
	int  i_offset, i_path, i_first_slot, i_slot_number, i_slot_number_mine;
	int  i_counter = 0;


	s_buf.Format("Dem number: %d", pc_fib_net->iGetDemandNumber());
	::vReportInFile("test.txt", s_buf, true);
	s_buf.Format("uni number: %d", pc_fib_net->i_uni_demands);
	::vReportInFile("test.txt", s_buf);
	s_buf.Format("multi number: %d", pc_fib_net->i_multi_demands);
	::vReportInFile("test.txt", s_buf);

	while  (!feof(pfSource)  &&  (i_counter < pc_fib_net->i_demand_num))
	{
		i_counter++;

		i_index = 0;
		s_line  =  Tools::sReadLine(pfSource, &s_setting_name);
		i_offset = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_path = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_slot_number = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_first_slot = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_first_slot--;
		i_slot_number_mine = 
			pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->iGetSlotsNumber
					(
					pc_fib_net->pcGetDemands()[i_offset].d_demand_capacity
					);

		s_slot_info = "";
		for  (int  ii = 0; ii < pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.size(); ii++)
		{
			s_buf2.Format("%d ", pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.at(ii));
			s_slot_info += s_buf2;
		}//for  (int  ii = 0; ii < pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.size(); ii++)

		s_buf.Format("%d %d %d (pos paths: %d) slice_num: %d (mine:%d)", 
			i_offset, i_path, i_first_slot, 
			pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->size(), i_slot_number, i_slot_number_mine
			);
		s_buf += s_slot_info;

		if  (i_slot_number !=  i_slot_number_mine)
			s_buf = "diffr ->" + s_buf;

		::vReportInFile("test.txt", s_buf);

		pc_channels[i_offset].pcPath = pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path);
		pc_channels[i_offset].iFirstSlotNumber = i_first_slot;
	
	}//while  (!feof(pfSource)  &&  (i_counter < pc_fib_net->i_demand_num))

	return(c_err);
}//CError  CFiberSolution::eLoadSolutionWalk(FILE  *pfSource)


CError  CFiberSolution::eLoadSolutionAibin(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf, s_slot_info, s_buf2;
	CString  s_line, s_setting_name;

	int  i_index;
	int  i_offset, i_path, i_first_slot, i_slot_number, i_slot_number_mine;
	int  i_counter = 0;


	s_buf.Format("Dem number: %d", pc_fib_net->iGetDemandNumber());
	::vReportInFile("test.txt", s_buf, true);
	s_buf.Format("uni number: %d", pc_fib_net->i_uni_demands);
	::vReportInFile("test.txt", s_buf);
	s_buf.Format("multi number: %d", pc_fib_net->i_multi_demands);
	::vReportInFile("test.txt", s_buf);

	while  (!feof(pfSource)  &&  (i_counter < pc_fib_net->i_demand_num))
	{
		i_counter++;

		i_index = 0;
		s_line  =  Tools::sReadLine(pfSource, &s_setting_name);
		i_offset = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_path = Tools::iExtractFromString(s_line, i_index, &i_index);
		//i_slot_number = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_first_slot = Tools::iExtractFromString(s_line, i_index, &i_index);
		//i_first_slot--;
		i_slot_number_mine = 
			pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->iGetSlotsNumber
					(
					pc_fib_net->pcGetDemands()[i_offset].d_demand_capacity
					);

		s_slot_info = "";
		for  (int  ii = 0; ii < pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.size(); ii++)
		{
			s_buf2.Format("%d ", pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.at(ii));
			s_slot_info += s_buf2;
		}//for  (int  ii = 0; ii < pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path)->v_slot_number_for_transfer.size(); ii++)

		s_buf.Format("%d %d %d (pos paths: %d) (mine slot num:%d)", 
			i_offset, i_path, i_first_slot, 
			pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->size(), i_slot_number_mine
			);
		s_buf += s_slot_info;

		//if  (i_slot_number !=  i_slot_number_mine)
		//	s_buf = "diffr ->" + s_buf;

		::vReportInFile("test.txt", s_buf);

		pc_channels[i_offset].pcPath = pc_fib_net->pcGetDemands()[i_offset].pvGetPossiblePaths()->at(i_path);
		pc_channels[i_offset].iFirstSlotNumber = i_first_slot;
	
	}//while  (!feof(pfSource)  &&  (i_counter < pc_fib_net->i_demand_num))

	return(c_err);
}//CError  CFiberSolution::eLoadSolutionAibin(FILE  *pfSource)




CError  CFiberSolution::eLoadSolutionAFA(FILE  *pfSource)
{
	CError  c_err;

	CString  s_line, s_setting_name;
	CString  s_buf;

	int  i_index;
	int  i_number, i_order, i_route_global, i_route_local, i_first_slot, i_modulation, i_replica, i_slots;
	int  i_counter = 0;

	int  i_replica_shift = 0;

	while  (!feof(pfSource)  &&  (i_counter < pc_fib_net->i_demand_num))
	{
		i_index = 0;
		s_line  =  Tools::sReadLine(pfSource, &s_setting_name);

		i_number = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_order = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_route_global = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_route_local = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_first_slot = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_modulation = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_replica = Tools::iExtractFromString(s_line, i_index, &i_index);
		i_slots = Tools::iExtractFromString(s_line, i_index, &i_index);


		pc_channels[i_number].pcPath = &(pc_fib_net->pc_paths[i_route_global]);
		pc_channels[i_number].iFirstSlotNumber = i_first_slot;

		i_counter++;

		if  (
			(pc_fib_net->bCompareNodes(pc_channels[i_number].pcGetDemand()->iGetStartNode(), pc_channels[i_number].pcPath->iGetStartNode())  == false)
			||
			(pc_fib_net->bCompareNodes(pc_channels[i_number].pcGetDemand()->iGetEndNode(),  pc_channels[i_number].pcPath->iGetEndNode())  == false)
			)
		{
			s_buf.Format
				(
				"[%d %d]   (%d %d) == (%d %d)", 
				i_number, i_counter,
				pc_channels[i_number].pcGetDemand()->iGetStartNode(), pc_channels[i_number].pcGetDemand()->iGetEndNode(),
				pc_channels[i_number].pcPath->iGetStartNode(), pc_channels[i_number].pcPath->iGetEndNode()
				);

			::MessageBox(NULL, s_buf, s_line, MB_OK);

			s_buf.Format("%d  %d  %d  %d  %d  %d  %d  %d", i_number, i_order, i_route_global, i_route_local, i_first_slot, i_modulation, i_replica, i_slots);
			::MessageBox(NULL, s_buf, s_line, MB_OK);
		}//if  ((
		
	}//while  (!feof(pfSource))
	
	return(c_err);
}//CError  CFiberSolution::eLoadSolutionAFA(FILE  *pfSource)




CError  CFiberSolution::eLoadSolution(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;
    CString  s_line  = "";
	
	int  i_demand_offset = 0;

	int  i_ind_uni, i_ind_any;
	int  i_index;

	int  i_path_number;
	int  i_start_slot_number;

	CString  s_uni_marker = SOLUTION_UNICAST_DEMAND_MARKER;
	CString  s_any_marker = SOLUTION_ANYCAST_DEMAND_MARKER;

	while  ( (i_demand_offset  <  pc_fib_net->i_demand_num)&&(!feof(pfSource)) )
	{
		s_line  =  Tools::sReadLine(pfSource);

		i_ind_uni = s_line.Find(SOLUTION_UNICAST_DEMAND_MARKER);
		i_ind_any = s_line.Find(SOLUTION_ANYCAST_DEMAND_MARKER);

		if  ( (i_ind_uni  ==  0)||(i_ind_any  ==  0) )
		{
			if  (i_ind_uni  ==  0)
			{
				i_index = s_uni_marker.GetLength();

				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);//the first number is unimportant
				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);
				i_start_slot_number = Tools::iExtractFromString(s_line, i_index, &i_index);


				if  (pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->size()  <= i_path_number)
				{
					CString  s_buf;

					s_buf.Format
						(
						"Demand: %d(%d -> %d) possible paths: %d demanded path number: %d", 
						i_demand_offset, pc_fib_net->pc_demands[i_demand_offset].i_fiber_node_start_offset, pc_fib_net->pc_demands[i_demand_offset].i_fiber_node_end_offset,
						pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->size(), i_path_number
						);
					c_err.vSetError(s_buf);
					return(c_err);
				}//if  (pv_possible_path_offsets->size()  >= i_path_number)

				pc_channels[i_demand_offset].pcPath = pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->at(i_path_number);
				pc_channels[i_demand_offset].iFirstSlotNumber = i_start_slot_number;


			}//if  (i_ind_uni  ==  0)
			else
			{
				i_index = s_any_marker.GetLength();

				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);//the first number is unimportant
				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);
				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);
				i_path_number = Tools::iExtractFromString(s_line, i_index, &i_index);
				i_start_slot_number = Tools::iExtractFromString(s_line, i_index, &i_index);


				if  (pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->size()  <= i_path_number)
				{
					CString  s_buf;

					s_buf.Format
						(
						"Demand: %d(%d -> %d) possible paths: %d demanded path number: %d", 
						i_demand_offset, pc_fib_net->pc_demands[i_demand_offset].i_fiber_node_start_offset, pc_fib_net->pc_demands[i_demand_offset].i_fiber_node_end_offset,
						pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->size(), i_path_number
						);
					c_err.vSetError(s_buf);
					return(c_err);
				}//if  (pv_possible_path_offsets->size()  >= i_path_number)

				pc_channels[i_demand_offset].pcPath = pc_fib_net->pc_demands[i_demand_offset].pvGetPossiblePaths()->at(i_path_number);
				pc_channels[i_demand_offset].iFirstSlotNumber = i_start_slot_number;
			
			}//else  if  (i_ind_uni  ==  0)

			i_demand_offset++;
		}//if  ( (i_ind_uni  ==  0)||(i_ind_any  ==  0) )
	}//while  ( (i_demands_read  <  pc_fib_net->i_demand_num)&&(feof(pfSource)) )


	
	return(c_err);
}//CError  CFiberSolution::eLoadSolution(FILE  *pfSource)



CError  CFiberSolution::eSaveSolution(FILE  *pfDest)
{
	CError  c_err;


	for  (int  i_dem = 0; i_dem < pc_fib_net->i_demand_num; i_dem++)
	{
		pc_channels[i_dem].vSave(pfDest);
		fprintf(pfDest, "\n");
	}//for  (int  ix = 0; ix < i_nodes_num; ix++)


	//c_err = eApplySolution();
	//if  (c_err)  return(c_err);

	fprintf(pfDest, "\n");
	fprintf(pfDest, "\n");

	for  (int  i_con = 0; i_con < pc_fib_net->i_conns_num; i_con++)
	{
		pc_fib_net->pc_conns[i_con].vSaveStatePRW(pfDest);
		fprintf(pfDest, "\n");
	}//for  (int  i_con = 0; i_con < pc_fib_net->i_conns_num; i_con++)


	bool  b_con_broken;
	fprintf(pfDest, "\n\n RESULT VALUE: %.2lf",  pc_fib_net->dCompCost() /*dCompCost(&b_con_broken)*/);

	return(c_err);
}//CError  CFiberSolution::eSaveSolution(FILE  *pfDest)


CError  CFiberSolution::eSaveSolution(CString  sSolDest)
{
	CError  c_err;
	FILE  *pf_dest;

	pf_dest = fopen(sSolDest, "w+");

	if  (pf_dest  ==  NULL)
	{
		c_err.vSetError("No sol file defined");
		return(c_err);
	}//if  (pf_source  ==  NULL)


	c_err  =  eSaveSolution(pf_dest);


	fclose(pf_dest);
	return(c_err);
}//CError  CFiberSolution::eSaveSolution(CString  sSource)






double  CFiberSolution::dCompCost (bool  *pbConstraintBroken, bool  bNoPenalty /*= false*/)
{
	double  d_result = 0;
	d_func_evaluations++;

	/*CError  c_err;
	c_err = eApplySolution();
	if  (c_err)  
	{
		c_err.csMessage += "//CFiberSolution::dCompCost";
		c_err.vShowWindow();
		d_result = -1;
		return(d_result);
	}//if  (c_err)  */

	

	*pbConstraintBroken  =  pc_fib_net->bConstraintBroken();
	d_result = pc_fib_net->dCompCost();
	//vShow("result");
	//vShow(d_result);

	if  ( (bNoPenalty  ==  false)&&(*pbConstraintBroken  ==  true) )
	{
		double  d_penalty = 0;
		d_result += d_penalty;
	}//if  ( (bNoPenalty  ==  false)&&(*pbConstraintBroken  ==  true) )

	d_last_raw_sol_value = d_result;

	//vShow("result2");
	//vShow(d_result);

	return(d_result);
}//double  CFiberSolution::dCompCost()

double CFiberSolution::d_clever_paths_initializations = 0;
double CFiberSolution::d_func_evaluations = 0;


//---------------------------class  CFiberSolutionChannel--------------------------------------
CFiberSolutionChannel::CFiberSolutionChannel()
{
	pc_solution = NULL;

	pcPath = NULL;
	iFirstSlotNumber  = -1;
}//CFiberSolution::CFiberSolution(CFiberNet  *pcFibNet)




CFiberSolutionChannel::~CFiberSolutionChannel()
{

}//CFiberSolutionChannel::~CFiberSolutionChannel()



CError  CFiberSolutionChannel::eSetChannel()
{
	CError  c_err;

	CString  s_buf;

	if  (pcPath  ==  NULL)
	{
		s_buf.Format("No path for the demand %d", pc_demand->i_demand_id);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (pcPath  ==  NULL)


	if  (iFirstSlotNumber < 0)
	{
		s_buf.Format("No slot defined for the demand %d", pc_demand->i_demand_id);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (iSlotNumber < 0)



	int  i_slots_number;
	i_slots_number = pcPath->iGetSlotsNumber(pc_demand->d_demand_capacity);
	c_err  =  pcPath->eSetPathSlots(iFirstSlotNumber, i_slots_number);

	if  (c_err)
	{
		s_buf.Format("  //demand %d  first slot: %d + %d", pc_demand->i_demand_id,  iFirstSlotNumber, i_slots_number);
		c_err.vSetError(c_err.sGetMessage() + s_buf);

		::vReportInFile("ai_solution.txt", c_err.sGetMessage());
	}//if  (c_err)



	return(c_err);
}//CError  CFiberSolutionChannel::eSetChannel()



CError  CFiberSolutionChannel::eFindChanel(CFiberPath  *pcChannelPath, int *piStartSlot, int  *piSlotsNumber)
{
	CError  c_err;

	*piSlotsNumber = pcChannelPath->iGetSlotsNumber(pc_demand->d_demand_capacity);	
	c_err  =  pcChannelPath->eFindPath(piStartSlot, *piSlotsNumber);

	return(c_err);
}//CError  CFiberSolutionChannel::eFindChanel(CFiberPath  *pcChannelPath)


CError  CFiberSolutionChannel::eSetChanel(CFiberPath  *pcChannelPath, int iStartSlot, int  iSlotsNumber)
{
	CError  c_err;

	iFirstSlotNumber = iStartSlot;
	pcPath =  pcChannelPath;
	c_err  =  pcChannelPath->eSetPath(iStartSlot, iSlotsNumber);

	return(c_err);
}//CError  CFiberSolutionChannel::eFindChanel(CFiberPath  *pcChannelPath)






CError  CFiberSolutionChannel::eFindAndSetChanel(CFiberPath  *pcChannelPath)
{
	CError  c_err;

	pcPath = pcChannelPath;
	
	int  i_slots_number;
	i_slots_number = pcPath->iGetSlotsNumber(pc_demand->d_demand_capacity);	
	c_err  =  pcPath->eFindAndSetPathSlots(&iFirstSlotNumber, i_slots_number);

	return(c_err);
}//CError  CFiberSolutionChannel::eFindAndSetChanel




void  CFiberSolutionChannel::vConfig(CFiberSolution  *pcSolution, CFiberDemand  *pcDemand)
{
	pc_solution = pcSolution;
	pc_demand = pcDemand;
}//void  CFiberSolutionChannel::vConfig(CFiberSolution  *pcSolution)


void  CFiberSolutionChannel::vCopy(CFiberSolutionChannel  *pcOther)
{
	pcOther->iFirstSlotNumber = iFirstSlotNumber;
	pcOther->pc_demand = pc_demand;
	pcOther->pc_solution = pc_solution;
	pcOther->pcPath = pcPath;
}//void  CFiberSolutionChannel::vCopy(CFiberSolutionChannel  *pcOther)


bool CFiberSolutionChannel::operator==(CFiberSolutionChannel  &cOther)
{
	if  (
		(iFirstSlotNumber == cOther.iFirstSlotNumber)&&
		(pc_demand == cOther.pc_demand)&&
		(pcPath == cOther.pcPath)
		)
		return(true);

	return(false);
}//bool CFiberSolutionChannel::operator==(CFiberSolutionChannel  &cOther)



bool CFiberSolutionChannel::operator!=(CFiberSolutionChannel  &cOther)
{
	if  (
		(iFirstSlotNumber != cOther.iFirstSlotNumber)||
		(pc_demand != cOther.pc_demand)||
		(pcPath != cOther.pcPath)
		)
		return(true);

	return(false);
}//bool CFiberSolutionChannel::operator!=(CFiberSolutionChannel  &cOther)



void  CFiberSolutionChannel::vSave(FILE  *pfDest)
{
	CString  s_buf;

	if  (pc_demand->i_demand_type  ==  DEMAND_TYPE_UNICAST)
		s_buf.Format("uni demand  %d  %d  %d  %d  %d %d", pc_demand->i_demand_id, pc_demand->i_fiber_node_start_offset, pc_demand->i_fiber_node_end_offset, iFirstSlotNumber, pcPath->iGetSlotsNumber(pc_demand->d_demand_capacity), pcPath->iGetId());

	if  (pc_demand->i_demand_type  ==  DEMAND_TYPE_MULTICAST)
		s_buf.Format("multi demand  %d  %d  %d  %d %d", pc_demand->i_demand_id, pc_demand->i_fiber_node_start_offset, iFirstSlotNumber, pcPath->iGetSlotsNumber(pc_demand->d_demand_capacity), pcPath->iGetId());
	
	if  (pc_demand->i_demand_type  ==  DEMAND_TYPE_ANYCAST)
		s_buf.Format("anycast demand  %d  %d(%d)  %d(%d)  %d  %d %d", pc_demand->i_demand_id, pc_demand->i_fiber_node_start_offset, pcPath->iGetStartNode(),  pc_demand->i_fiber_node_end_offset, pcPath->iGetEndNode(), iFirstSlotNumber, pcPath->iGetSlotsNumber(pc_demand->d_demand_capacity), pcPath->iGetId());
	fprintf(pfDest, s_buf);
}//void  CFiberSolutionChannel::vSave(FILE  *pfDest)






//---------------------------class  CSteinerTree--------------------------------------------------
CSteinerTree::CSteinerTree()
{
	i_root = -1;
	pc_base_network = NULL;

	pi_tree_links = NULL;
	pi_nodes_tool = NULL;

	b_tree_cost_actual = false;
}//CSteinerTree::CSteinerTree()


CSteinerTree::CSteinerTree(CFiberNet  *pcBaseNetwork)
{
	i_root = -1;
	vSetBaseNetwork(pcBaseNetwork);

	b_tree_cost_actual = false;
}//CSteinerTree::CSteinerTree(CFiberNet  *pcBaseNetwork)



CSteinerTree::~CSteinerTree()
{
	//vShow("CSteinerTree::~CSteinerTree()");
	if  (pi_tree_links != NULL)  delete  pi_tree_links;
	if  (pi_nodes_tool != NULL)  delete  pi_nodes_tool;
}//CSteinerTree::~CSteinerTree()


void  CSteinerTree::vSetBaseNetwork(CFiberNet  *pcBaseNetwork)
{
	b_tree_cost_actual = false;

	if  (pi_tree_links != NULL)  delete  pi_tree_links;
	if  (pi_nodes_tool != NULL)  delete  pi_nodes_tool;

	pc_base_network = pcBaseNetwork;

	if  (pc_base_network  ==  NULL)
	{
		pi_tree_links = NULL;
		pi_nodes_tool = NULL;

		return;
	}//if  (pc_base_network  ==  NULL)

	pi_tree_links = new int[pc_base_network->iGetConnsNumber()];
	pi_nodes_tool = new int[pc_base_network->iGetNodesNumber()];
}//void  CSteinerTree::vSetBaseNetwork(CFiberNet  *pcBaseNetwork)



void  CSteinerTree::vCopyTo(CSteinerTree  *pcOtherTree)
{
	pcOtherTree->v_demanded_nodes = v_demanded_nodes;
	pcOtherTree->v_not_demanded_nodes = v_not_demanded_nodes;
	pcOtherTree->i_root = i_root;
	pcOtherTree->b_tree_cost_actual = b_tree_cost_actual;
	pcOtherTree->d_tree_cost = d_tree_cost;

	if  (pcOtherTree->pc_base_network  !=  pc_base_network)  pcOtherTree->vSetBaseNetwork(pc_base_network);

	if  (pc_base_network  ==  NULL)  return;

	for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
		pcOtherTree->pi_tree_links[ii] = pi_tree_links[ii];
}//void  CSteinerTree::vCopyTo(CSteinerTree  *pcOtherTree)



CError  CSteinerTree::eSaveTreeDemands(CString  sDMCFile)
{
	CError  c_err;
	CString  s_buf;

	if  (sDMCFile  ==  "")
	{
		c_err.vSetError("no dmc file for save defined");
		return(c_err);	
	}//if  (sDMCFile  ==  "")


	FILE  *pf_dmc_file;
	pf_dmc_file = fopen(sDMCFile, "w+");

	if  (pf_dmc_file  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sDMCFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_dmc_file  ==  NULL)

	c_err = eSaveTreeDemands(pf_dmc_file);
	

	fclose(pf_dmc_file);
	return(c_err);
}//CError  CSteinerTree::eSaveTreeDemands(CString  sDMCFile)



CError  CSteinerTree::eSaveTreeDemands(FILE  *pfDMCFile)
{
	CError  c_err;
	CString  s_buf;

	CString  s_line;

	s_line.Format("%d %d  ", v_demanded_nodes.size(), i_root);

	CString  s_children;
	s_children = "";

	for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)
	{
		s_buf.Format("%d ", v_demanded_nodes.at(ii));
		s_children += s_buf;			
	}//for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)

	s_line += s_children + "\n";
	fprintf(pfDMCFile, s_line);

	return(c_err);
}//CError  CSteinerTree::eSaveTreeDemands(FILE  *pfDMCFile)




CError  CSteinerTree::eLoadTreeDemands(CString  sDMCFile)
{
	CError  c_err;
	CString  s_buf;

	b_tree_cost_actual = false;

	if  (sDMCFile  ==  "")
	{
		c_err.vSetError("no dmc file");
		return(c_err);	
	}//if  (sDMCFile  ==  "")


	FILE  *pf_dmc_file;
	pf_dmc_file = fopen(sDMCFile, "r");

	if  (pf_dmc_file  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sDMCFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_dmc_file  ==  NULL)

	c_err = eLoadTreeDemands(pf_dmc_file);
	

	fclose(pf_dmc_file);
	return(c_err);
}//CError  CSteinerTree::eLoadTreeDemands(CString  sDMCFile)



CError  CSteinerTree::eLoadTreeDemands(FILE  *pfDMCFile)
{
	CError  c_err;
	CString  s_buf;
	CString  s_line;

	b_tree_cost_actual = false;

	s_line  =  Tools::sReadLine(pfDMCFile);

	int  i_index = 0;
	int  i_children_number;

	i_children_number = Tools::iExtractFromString(s_line, i_index, &i_index);
	i_root = Tools::iExtractFromString(s_line, i_index, &i_index);

	if  (i_children_number  < 0)
	{
		s_buf.Format("Children number: %d", i_children_number);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (i_children_number  < 0)

	if  (i_root  < 0)
	{
		s_buf.Format("Root offset: %d", i_root);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (i_children_number  < 0)


	int  i_child_offset;
	for  (int  i_child = 0; i_child < i_children_number; i_child++)
	{
		i_child_offset = Tools::iExtractFromString(s_line, i_index, &i_index);
		v_demanded_nodes.push_back(i_child_offset);
	}//for  (int  i_child = 0; i_child < i_children_number; i_child++)

	bool  b_demanded_node;
	for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	{
		b_demanded_node = false;
		
		for  (int i_offset = 0; (i_offset < v_demanded_nodes.size())&&(b_demanded_node == false) ; i_offset++)
		{
			if  (v_demanded_nodes.at(i_offset) == i_node)  b_demanded_node = true;
		}//for  (int i_offset = 0; i_offset < v_demanded_nodes.size(); i_offset++)


		if  (b_demanded_node == false)
		{
			v_not_demanded_nodes.push_back(i_node);
			//::vShow(i_node);
		}//if  (b_demanded_node == false)
	}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)



	return(c_err);
}//CError  CSteinerTree::eLoadTreeDemands(FILE  *pfDMCFile)




CError  CSteinerTree::eLoadTree(CString  sTreeFile)
{
	CError  c_err;
	CString  s_buf;

	b_tree_cost_actual = false;

	if  (sTreeFile  ==  "")
	{
		c_err.vSetError("no tree file");
		return(c_err);	
	}//if  (sDMCFile  ==  "")


	FILE  *pf_tree;
	pf_tree = fopen(sTreeFile, "r");

	if  (pf_tree  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sTreeFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_dmc_file  ==  NULL)

	int  i_loaded;
	c_err = eLoadTree(pf_tree, &i_loaded);

	fclose(pf_tree);

	if  (i_loaded  ==  0)  c_err.vSetError("Tree not loaded");

	return(c_err);
}//CError  CSteinerTree::eLoadTree(CString  sTreeFile)



CError  CSteinerTree::eLoadTree(FILE  *pfTreeFile, int  *piLoaded)
{
	CError  c_err;
	CString  s_buf;
	CString  s_line;

	b_tree_cost_actual = false;

	*piLoaded = 1;
	s_line  =  Tools::sReadLine(pfTreeFile);

	s_buf = s_line;
	s_buf.Replace(" ", "");
	if  (s_buf  ==  "")
	{
		*piLoaded = 0;
		return(c_err);	
	}//if  (s_buf  ==  "")

	
	int  i_index = 0;
	int  i_buf;

	v_zero_tree_links();

	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		i_buf = Tools::iExtractFromString(s_line, i_index, &i_index);

		if  (i_buf > 0) pi_tree_links[i_link] = 1;
	}//for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)


	return(c_err);
}//CError  CSteinerTree::eLoadTree(FILE  *pfTreeFile)




CError  CSteinerTree::eSaveTree(CString  sTreeFile, bool  bDetails)
{
	CError  c_err;
	CString  s_buf;

	if  (sTreeFile  ==  "")
	{
		c_err.vSetError("no tree file for save defined");
		return(c_err);	
	}//if  (sTreeFile  ==  "")


	FILE  *pf_tree;
	pf_tree = fopen(sTreeFile, "w+");

	if  (pf_tree  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sTreeFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_dmc_file  ==  NULL)

	if  (bDetails == false)
		c_err = eSaveTree(pf_tree);
	else
		c_err = eSaveTreeDetails(pf_tree);

	

	fclose(pf_tree);
	return(c_err);
}//CError  CSteinerTree::eSaveTree(CString  sTreeFile)



CError  CSteinerTree::eSaveTree(FILE  *pfTreeFile)
{
	CError  c_err;
	CString  s_buf, s_line_to_save;

	s_line_to_save = "";

	if  (b_check_tree_no_loops() == true)
		s_line_to_save += "1\t";
	else
	{
		//musze zestawic wszystkie, bo b_check_tree_no_loops() konczy po pierwszym bledzie
		v_zero_node_tool();
		for  (int  i_conn = 0; i_conn < pc_base_network->iGetConnsNumber(); i_conn++)
		{
			if  (pi_tree_links[i_conn] > 0)
			{
				(pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()])++;
				//if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()]  >  1)  return(false);	
			}//if  (pi_tree_links[i_conn] > 0)
		}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)

		s_line_to_save += "0\t";
	}//else  
	
	if  (b_check_tree_nodes_covered() == true)
		s_line_to_save += "1\t";
	else
		s_line_to_save += "0\t";

	if  (b_check_tree_all_connected() == true)
		s_line_to_save += "1\t";
	else
		s_line_to_save += "0\t";


	if  (b_check_tree_hanging_leaves() == true)
		s_line_to_save += "1\t";
	else
		s_line_to_save += "0\t";
	


	s_line_to_save += "|\t";
	//::vShow(pc_base_network->iGetConnsNumber());
	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		
		s_buf.Format("%d\t", pi_tree_links[i_link]);
		
		s_line_to_save += s_buf;
	}//for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)

	s_line_to_save += "\n";

	fprintf(pfTreeFile, s_line_to_save);

	return(c_err);
}//CError  CSteinerTree::eSaveTree(FILE  *pfTreeFile)



void  CSteinerTree::v_save_nodes_state(FILE  *pfTreeFile)
{
	CString  s_line, s_enumeration_line, s_buf;

	fprintf(pfTreeFile, "\n\nNODES STATE:\n");
	s_line = "";
	s_enumeration_line = "";
	for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	{
		s_buf.Format("%d\t", i_node);
		s_enumeration_line += s_buf;

		s_buf.Format("%d\t", pi_nodes_tool[i_node]);
		s_line += s_buf;	
	}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	s_enumeration_line += "\n";
	s_line += "\n\n";
	fprintf(pfTreeFile, s_enumeration_line);
	fprintf(pfTreeFile, s_line);
}//void  CSteinerTree::v_save_nodes_state(FILE  *pfTreeFile)



CError  CSteinerTree::eSaveTreeDetails(FILE  *pfTreeFile)
{
	CError  c_err;
	CString  s_line, s_enumeration_line, s_buf;

	c_err = eSaveTree(pfTreeFile);
	if  (c_err)  return(c_err);


	//loops check
	v_zero_node_tool();

	pi_nodes_tool[i_root] = 1;
	for  (int  i_conn = 0; i_conn < pc_base_network->iGetConnsNumber(); i_conn++)
	{
		if  (pi_tree_links[i_conn] > 0)
		{
			(pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()])++;
			//if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()]  >  1)  return(false);	
		}//if  (pi_tree_links[i_conn] > 0)
	}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)



	v_save_nodes_state(pfTreeFile);
	

	bool  b_error;
	b_error = false;
	fprintf(pfTreeFile, "\n\nLOOPS ERRORS: ");
	for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	{
		if  (pi_nodes_tool[i_node]  >  1)
		{
			b_error = true;
			s_line.Format("\nNode: %d Links leading to it:\n", i_node);
			fprintf(pfTreeFile, s_line);

			for  (int  i_conn = 0; i_conn < pc_base_network->iGetConnsNumber(); i_conn++)
			{
				if  (pi_tree_links[i_conn] > 0)
				{
					if  (pc_base_network->pcGetConn(i_conn)->iGetNodeEnd() == i_node)
					{
						s_line.Format("Connection: %d  (%d->%d)\n", i_conn, pc_base_network->pcGetConn(i_conn)->iGetNodeStart(), pc_base_network->pcGetConn(i_conn)->iGetNodeEnd());
						fprintf(pfTreeFile, s_line);
					}//if  (pc_base_network->pcGetConn(i_conn)->iGetNodeEnd() == i_node)
				}//if  (pi_tree_links[i_conn] > 0)
			}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
			
		}//if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()]  >  1)	
	}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)


	if  (b_error == false)
		fprintf(pfTreeFile, "NONE\n\n");
	else
		fprintf(pfTreeFile, "\n");
		


	//missing demanded nodes
	s_line = "";
	b_error = false;
	fprintf(pfTreeFile, "\n\nDEMANDED NODES MISSING: ");
	for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)
	{
		if  (pi_nodes_tool[v_demanded_nodes.at(ii)]  ==  0)
		{
			b_error = true;
			s_buf.Format("%d\t", v_demanded_nodes.at(ii));
			s_line += s_buf;
		}//if  (pi_nodes_tool[v_demanded_nodes.at(ii)]  ==  0)
	}//for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)

	if  (b_error == false)
		fprintf(pfTreeFile, "NONE\n\n");
	else
	{
		s_line += "\n\n";
		fprintf(pfTreeFile, s_line);
	}//else  if  (b_error == false)


	//unconnected links
	s_line = "";
	b_error = false;
	fprintf(pfTreeFile, "\n\nUNCONNECTED LINKS: ");
	for  (int  i_conn = 0; i_conn < pc_base_network->iGetConnsNumber(); i_conn++)
	{
		if  (pi_tree_links[i_conn] > 0)
		{
			if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeStart()] ==  0)
			{
				s_line.Format("\nConnection: %d  (%d->%d)", i_conn, pc_base_network->pcGetConn(i_conn)->iGetNodeStart(), pc_base_network->pcGetConn(i_conn)->iGetNodeEnd());
				fprintf(pfTreeFile, s_line);

				b_error = true;
			}//if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeStart()] ==  0)
			//if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()]  >  1)  return(false);	
		}//if  (pi_tree_links[i_conn] > 0)
	}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
	

	if  (b_error == false)
		fprintf(pfTreeFile, "NONE\n\n");
	else
		fprintf(pfTreeFile, "\n\n");



	fprintf(pfTreeFile, "\n\n");

	//unnecessary leaves
	b_check_tree_all_connected();//to build the proper nodes state
	c_err = eSaveTreeDemands(pfTreeFile);
	if  (c_err)  return(c_err);

	v_save_nodes_state(pfTreeFile);

	

	s_line = "";
	b_error = false;
	fprintf(pfTreeFile, "\n\nUNNECESSARY LEAVES: ");

	bool  b_node_found;
	for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	{
		if  (pi_nodes_tool[i_node] ==  1)
		{
			b_node_found = false;
			for  (int  ii = 0; (ii < v_demanded_nodes.size())&&(b_node_found == false) ; ii++)
			{
				if  (v_demanded_nodes.at(ii) == i_node)  b_node_found = true;
			}//for  (int  ii = 0; (ii < v_demanded_nodes.size())&&(b_node_found == false) ; ii++)

			if  (b_node_found == false)
			{
				b_error = true;
				s_buf.Format("%d\t", i_node);
				s_line += s_buf;
			}//if  (b_node_found == false)
		}//if  (pi_nodes_tool[i_node] ==  1)	
	}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)

	if  (b_error == false)
		fprintf(pfTreeFile, "NONE\n\n");
	else
	{
		s_line += "\n\n";
		fprintf(pfTreeFile, s_line);
	}//if  (b_error == false)




	c_err = eSaveTreeStructure(pfTreeFile);
	if  (c_err)  return(c_err);


	return(c_err);
}//CError  CSteinerTree::eSaveTreeDetails(FILE  *pfTreeFile)




CError  CSteinerTree::eSaveTreeStructure(FILE  *pfTreeFile)
{
	CError  c_err;

	CString  s_buf;

	fprintf(pfTreeFile, "\n\n TREE STRUCTURE: \n\n");


	v_zero_node_tool();
	pi_nodes_tool[i_root] = 1;

	vector<int>  v_links;
	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		if  (pi_tree_links[i_link])  v_links.push_back(i_link);	
	}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)


	fprintf(pfTreeFile, "Root: %d\n", i_root);
	

	CString  s_up_line, s_med_line, s_down_line;
	int  i_prev_links;
	i_prev_links = v_links.size() + 1;
	vector<int>  v_links_to_save;
	while  (v_links.size()  <  i_prev_links)
	{
		i_prev_links  =  v_links.size();
		v_links_to_save.clear();
		
		for  (int  i_link_offset = 0; i_link_offset < v_links.size(); i_link_offset++)
		{
			if  (
				pi_nodes_tool[pc_base_network->pcGetConn(v_links.at(i_link_offset))->iGetNodeStart()] > 0
				)
			{
				v_links_to_save.push_back(v_links.at(i_link_offset));

				v_links.erase(v_links.begin() + i_link_offset);
				i_link_offset--;
			}//if  (
		}//for  (int  i_link_offset = 0; i_link_offset < v_links.size(); i_link_offset++)


		for (int  i_link_to_save = 0; i_link_to_save < v_links_to_save.size(); i_link_to_save++)
		{
			pi_nodes_tool[pc_base_network->pcGetConn(v_links_to_save.at(i_link_to_save))->iGetNodeStart()] = 2;
			pi_nodes_tool[pc_base_network->pcGetConn(v_links_to_save.at(i_link_to_save))->iGetNodeEnd()] = 1;
		}//for (int  i_link_to_save = 0; i_link_to_save < v_links_to_save.size(); i_link_to_save++)


		s_up_line = "";
		s_med_line = "";
		s_down_line = "";		
		for (int  i_link_to_save = 0; i_link_to_save < v_links_to_save.size(); i_link_to_save++)
		{
			s_buf.Format("[%.2d] ", pc_base_network->pcGetConn(v_links_to_save.at(i_link_to_save))->iGetNodeStart());
			s_up_line += s_buf;

			s_buf.Format("L:%.2d ", v_links_to_save.at(i_link_to_save));
			s_med_line += s_buf;

			s_buf.Format("[%.2d] ", pc_base_network->pcGetConn(v_links_to_save.at(i_link_to_save))->iGetNodeEnd());
			s_down_line += s_buf;		
		}//for (int  i_link_to_save = 0; i_link_to_save < v_links_to_save.size(); i_link_to_save++)

		fprintf(pfTreeFile, "%s\n", s_up_line);
		fprintf(pfTreeFile, "%s\n", s_med_line);
		fprintf(pfTreeFile, "%s\n", s_down_line);
		fprintf(pfTreeFile, "\n");
	}//while  (v_links.size()  <  i_prev_links)

	fprintf(pfTreeFile, "\n\nNOT USED TREE LINKS:\n");

	for (int  i_link_to_save = 0; i_link_to_save < v_links.size(); i_link_to_save++)
	{
		s_buf.Format("[%.2d] ", pc_base_network->pcGetConn(v_links.at(i_link_to_save))->iGetNodeStart());
		s_up_line += s_buf;

		s_buf.Format("L:%.2d ", v_links.at(i_link_to_save));
		s_med_line += s_buf;

		s_buf.Format("[%.2d] ", pc_base_network->pcGetConn(v_links.at(i_link_to_save))->iGetNodeEnd());
		s_down_line += s_buf;		
	}//for (int  i_link_to_save = 0; i_link_to_save < v_links.size(); i_link_to_save++)

	fprintf(pfTreeFile, "%s\n", s_up_line);
	fprintf(pfTreeFile, "%s\n", s_med_line);
	fprintf(pfTreeFile, "%s\n", s_down_line);
	fprintf(pfTreeFile, "\n");


	return(c_err);
}//CError  CSteinerTree::eSaveTreeStructure(FILE  *pfTreeFile)




bool  CSteinerTree::bGenotypesTheSame(CSteinerTree  *pcOtherTree)
{
	if  (pc_base_network->iGetConnsNumber()  !=  pcOtherTree->pc_base_network->iGetConnsNumber())  return(false);

	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		if  (pi_tree_links[i_link]  !=  pcOtherTree->pi_tree_links[i_link])  return(false);
	}//for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)

	return(true);
}//bool  CSteinerTree::bGenotypesTheSame(CSteinerTree  *pcOtherTree)



CError  CSteinerTree::e_build_one_node(int  iNodeOffset, int iNodeLevel, int  *piConnsIncluded, bool  bForceConnections)
{
	CError  c_err;

	pi_nodes_tool[iNodeOffset] = iNodeLevel;

	bool  b_add_link_to_tree;
	CFiberConnection  *pc_link;
	//first add what can be added...
	for  (int  i_end_node = 0; i_end_node < pc_base_network->iGetNodesNumber();  i_end_node++)
	{
		pc_link = pc_base_network->pcGetConnMat()[iNodeOffset][i_end_node];
		if  (pc_link  !=  NULL)
		{
			if  (pi_nodes_tool[i_end_node] == 0)
			{
				if  (pi_tree_links[pc_link->iGetId()] == 0)//to jest nadmiarowe dla "if  (pi_nodes_tool[i_end_node] == 0)", ale powinno byc
				{
					b_add_link_to_tree = bForceConnections;

					//if we do not force then we check the connections included by genotype
					if  (b_add_link_to_tree  ==  false)
					{
						if  (piConnsIncluded[pc_link->iGetId()] == 1)  b_add_link_to_tree = true;				
					}//if  (b_add_link_to_tree  ==  false)


					if  (b_add_link_to_tree  ==  true)
					{
						pi_tree_links[pc_link->iGetId()] = 1;
						e_build_one_node(i_end_node, iNodeLevel-1, piConnsIncluded, false);
					}//if  (b_add_link_to_tree  ==  true)
				}//if  (pi_tree_links[pc_link->iGetId()] == 0) 						
			}//if  (pi_nodes_tool[i_end_node] == 0)
		}//if  (pc_link  !=  NULL)
	}//for  (int  i_out_link = 0; i_out_link < pc_base_network->iGetNodesNumber();  i_out_link++)


	return(c_err);
}//CError  CSteinerTree::e_build_one_node(int  iNodeOffset, int  *piConnsIncluded, bool  bForceConnections)



void  CSteinerTree::v_cut_leaves()
{
	CString  s_buf;

	CFiberConnection  *pc_link;
	int  i_link_id;
	bool  b_node_cut_off;

	b_node_cut_off = true;

	while (b_node_cut_off == true)
	{
		b_node_cut_off = false;

		b_check_tree_all_connected();
		for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
		{
			if  (pi_nodes_tool[i_node] ==  1)
			{
				if  (b_node_is_not_demanded(i_node) == true)
				{
					for  (int  i_start_node = 0; i_start_node < pc_base_network->iGetNodesNumber(); i_start_node++)
					{
						pc_link = pc_base_network->pcGetConnMat()[i_start_node][i_node];

						if  (pc_link  !=  NULL)
						{
							i_link_id  =  pc_link->iGetId();
							//s_buf.Format("%d\n",i_link_id);
							pi_tree_links[i_link_id] = 0;
							//::vReportInFile("zzz_cut.txt", s_buf);
							b_node_cut_off = true;
						}//if  (pc_link  !=  NULL)
					}//for  (int  i_link_off_in_node = 0; i_link_off_in_node < pc_base_network->iGetNodesNumber(); i_link_off_in_node++)			
				}//if  (b_node_is_not_demanded(i_node) == true)
			}//if  (pi_nodes_tool[i_node] ==  1)	
		}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	}//while (b_node_cut_off == true)

}//void  CSteinerTree::v_cut_leaves()



bool  CSteinerTree::b_node_is_not_demanded(int  iNode)
{
	//::vShow("b_node_is_not_demanded(int  iNode)");
	for  (int  ii = 0; ii < v_not_demanded_nodes.size(); ii++)
	{
		//::vShow(v_not_demanded_nodes.at(ii));
		if  (v_not_demanded_nodes.at(ii)  ==  iNode)  return(true);
	}//for  (int  ii = 0; ii < v_not_demanded_nodes.size(); ii++)

	return(false);
}//bool  CSteinerTree::b_node_is_not_demanded(int  iNode)




void  CSteinerTree::vReadConns(int  *piConnsIncluded)
{
	v_zero_tree_links();
	v_zero_node_tool();

	for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
		pi_tree_links[ii] = piConnsIncluded[ii];
}//void  CSteinerTree::vReadConns(int  *piConnsIncluded)



CError  CSteinerTree::eBuildTreeFromConns(int  *piConnsIncluded)
{
	CError  c_err;

	v_zero_tree_links();
	v_zero_node_tool();


	int  i_current_level = pc_base_network->iGetNodesNumber() + 1;

	c_err = e_build_one_node(i_root, i_current_level, piConnsIncluded, false);
	c_err = e_build_one_node(i_root, i_current_level, piConnsIncluded, true);
	

	for  (
		int  i_current_level = pc_base_network->iGetNodesNumber();
		i_current_level > 0;
		i_current_level--
		)
	{
		for  (int i_node = 0; i_node < pc_base_network->iGetNodesNumber();  i_node++)
		{
			if  (pi_nodes_tool[i_node] == i_current_level)  c_err = e_build_one_node(i_node, i_current_level, piConnsIncluded, false);
		}//for  (int i_node = 0; i_node < pc_base_network->iGetNodesNumber();  i_node++)


		for  (int i_node = 0; i_node < pc_base_network->iGetNodesNumber();  i_node++)
		{
			if  (pi_nodes_tool[i_node] == i_current_level)  c_err = e_build_one_node(i_node, i_current_level, piConnsIncluded, true);
		}//for  (int i_node = 0; i_node < pc_base_network->iGetNodesNumber();  i_node++)

	}//for  (


	v_cut_leaves();


	return(c_err);
}//CError  CSteinerTree::eBuildTreeFromConns(int  *piConnsIncluded)


double  CSteinerTree::dGetTreeCost()
{
	if  (b_tree_cost_actual ==  true)  return(d_tree_cost);

	if  (i_root < 0)  
	{
		b_tree_cost_actual =  false;
		return(-1);
	}//if  (i_root < 0)  


	d_tree_cost = 0;
	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		if  (pi_tree_links[i_link]  >  0)
			d_tree_cost += pc_base_network->pcGetConn(i_link)->dGetCapacity();
	
	}//for  (int  ii = 0; ii < v_tree_links.size(); ii++)
	
	b_tree_cost_actual = true;

	return(d_tree_cost);
}//double  CSteinerTree::dGetTreeCost()


void  CSteinerTree::v_zero_tree_links()
{
	for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
		pi_tree_links[ii] = 0;
	b_tree_cost_actual =  false;
}//CError  CSteinerTree::v_zero_tree_links()


void  CSteinerTree::v_zero_node_tool()
{
	for  (int  ii = 0; ii < pc_base_network->iGetNodesNumber(); ii++)
		pi_nodes_tool[ii] = 0;
}//void  CSteinerTree::v_zero_node_tool()


bool  CSteinerTree::bCheckTree()
{
	if  (b_check_tree_nodes_covered_and_no_loops() ==  false)  return(false);
	if  (b_check_tree_all_connected() ==  false)  return(false);
	if  (b_check_tree_hanging_leaves() ==  false)  return(false);

	return(true);
};//bool  CSteinerTree::bCheckTree()



bool  CSteinerTree::b_check_tree_no_loops()
{
	v_zero_node_tool();

	pi_nodes_tool[i_root] = 1;
	for  (int  i_conn = 0; i_conn < pc_base_network->iGetConnsNumber(); i_conn++)
	{
		if  (pi_tree_links[i_conn] > 0)
		{
			(pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()])++;
			if  (pi_nodes_tool[pc_base_network->pcGetConn(i_conn)->iGetNodeEnd()]  >  1)  return(false);	
		}//if  (pi_tree_links[i_conn] > 0)
	}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
	return(true);
};//bool  CSteinerTree::b_check_tree_no_loops()



//RUN ONLY AFTER b_check_tree_no_loops!!!
bool  CSteinerTree::b_check_tree_nodes_covered()
{
	for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)
	{
		if  (pi_nodes_tool[v_demanded_nodes.at(ii)]  ==  0)  return(false);	
	}//for  (int  ii = 0; ii < v_demanded_nodes.size(); ii++)

	return(true);
}//bool  CSteinerTree::b_check_tree_nodes_covered()



bool  CSteinerTree::b_check_tree_nodes_covered_and_no_loops()
{
	v_zero_node_tool();

	if  (b_check_tree_no_loops() == false) return(false);
	if  (b_check_tree_nodes_covered() == false) return(false);

	return(true);
};//bool  CSteinerTree::b_check_tree_nodes_covered()




bool  CSteinerTree::b_check_tree_all_connected()
{
	v_zero_node_tool();

	pi_nodes_tool[i_root] = 1;

	vector<int>  v_links;
	for  (int  i_link = 0; i_link < pc_base_network->iGetConnsNumber(); i_link++)
	{
		if  (pi_tree_links[i_link])  v_links.push_back(i_link);	
	}//for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)

	int  i_prev_links;

	i_prev_links = v_links.size() + 1;

	while  (v_links.size()  <  i_prev_links)
	{
		i_prev_links  =  v_links.size();
		
		for  (int  i_link_offset = 0; i_link_offset < v_links.size(); i_link_offset++)
		{
			if  (
				pi_nodes_tool[pc_base_network->pcGetConn(v_links.at(i_link_offset))->iGetNodeStart()] > 0
				)
			{
				pi_nodes_tool[pc_base_network->pcGetConn(v_links.at(i_link_offset))->iGetNodeStart()] = 2;
				pi_nodes_tool[pc_base_network->pcGetConn(v_links.at(i_link_offset))->iGetNodeEnd()] = 1;

				v_links.erase(v_links.begin() + i_link_offset);
				i_link_offset--;
			}//if  (
		}//for  (int  i_link_offset = 0; i_link_offset < v_links.size(); i_link_offset++)
	}//while  (v_links.size()  <  i_prev_links)

	if  (v_links.size() > 0)  return(false);


	return(true);
};//bool  CSteinerTree::b_check_tree_all_connected()


//RUN ONLY AFTER b_check_tree_all_connected!!!
bool  CSteinerTree::b_check_tree_hanging_leaves()
{
	bool  b_node_found;
	for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)
	{
		if  (pi_nodes_tool[i_node] ==  1)
		{
			b_node_found = false;
			for  (int  ii = 0; (ii < v_demanded_nodes.size())&&(b_node_found == false) ; ii++)
			{
				if  (v_demanded_nodes.at(ii) == i_node)  b_node_found = true;
			}//for  (int  ii = 0; (ii < v_demanded_nodes.size())&&(b_node_found == false) ; ii++)

			if  (b_node_found == false) return(false);
		}//if  (pi_nodes_tool[i_node] ==  1)	
	}//for  (int  i_node = 0; i_node < pc_base_network->iGetNodesNumber(); i_node++)

	return(true);
}//CSteinerTree::b_check_tree_hanging_leaves()





//---------------------------class  CSteinerTreeSet--------------------------------------------------

CSteinerTreeSet::CSteinerTreeSet()
{
	d_min_val = 0;
	d_max_val = DBL_MAX;
	d_lower_than_max_scaling = 1;
	i_capacity = 1;

	d_penalty_factor = 1;
}//CSteinerTreeSet::CSteinerTreeSet()



CSteinerTreeSet::~CSteinerTreeSet()
{
	for  (int  ii = 0; ii < v_tree_set.size(); ii++)
		delete  v_tree_set.at(ii);

	if  (pc_base_network != NULL)  delete  pc_base_network;
}//CSteinerTreeSet::~CSteinerTreeSet()


void  CSteinerTreeSet::vSetCapacity(int  iNewCapacity)
{
	i_capacity = iNewCapacity;

	while( (v_tree_set.size()  >  i_capacity)&&(v_tree_set.size() > 0) )
	{
		delete  v_tree_set.at(v_tree_set.size() - 1);
		v_tree_set.erase(v_tree_set.begin() + v_tree_set.size() - 1);
	}//while( (v_tree_set.size()  >  i_capacity)&&(v_tree_set.size() > 0) )

}//void  CSteinerTreeSet::vSetCapacity(int  iNewCapacity)



void  CSteinerTreeSet::v_order_trees()
{
	sort(v_tree_set.begin(), v_tree_set.end(), SSteinerTreeSort());
}//void  CSteinerTreeSet::v_order_trees()


bool  CSteinerTreeSet::bTreeContained(CSteinerTree  *pcNewTree)
{
	for  (int  ii = 0; ii < v_tree_set.size(); ii++)
	{
		if  (v_tree_set.at(ii)->bGenotypesTheSame(pcNewTree) == true)  return(true);
	}//for  (int  ii = 0; ii < v_tree_set.size(); ii++)

	return(false);
}//bool  CSteinerTreeSet::bTreeContained(CSteinerTree  *pcNewTree)




bool  CSteinerTreeSet::bAddTree(CSteinerTree  *pcNewTree)
{
	//first check if we can add
	if  (bTreeContained(pcNewTree)  ==  true)  return(false);

	bool  b_can_add;
	b_can_add = false;

	if  (v_tree_set.size() < i_capacity)
	{
		b_can_add = true;
	}//if  (v_tree_set.size() < i_capacity)
	else
	{
		if  (v_tree_set.at(v_tree_set.size() -1)->dGetTreeCost()  >  pcNewTree->dGetTreeCost())
		{
			b_can_add = true;
			delete  v_tree_set.at(v_tree_set.size() -1);
			v_tree_set.erase(v_tree_set.begin() + v_tree_set.size() -1);		
		}//if  (v_tree_set.at(v_tree_set.size() -1)->dGetTreeCost()  >  pcNewTree->dGetTreeCost())
	}//else  if  (v_tree_set.size() < i_capacity)

	if  (b_can_add ==  false)  return(false);

	CSteinerTree  *pc_copy;
	pc_copy = new CSteinerTree();

	pcNewTree->vCopyTo(pc_copy);
	v_tree_set.push_back(pc_copy);

	v_order_trees();
}//bool  CSteinerTreeSet::bAddTree(CSteinerTree  *pcNewTree)




bool  CSteinerTreeSet::bCheckUpTreeAll(int  iShift, int  *piArgument,  int  iArgLen)
{
	CSteinerTree  c_buf;
	vCopyTo(&c_buf);

	c_buf.vReadConns(piArgument + iShift);
	
	if  (b_check_tree_no_loops() == false)  return(false);
	if  (b_check_tree_nodes_covered() == false)  return(false);
	if  (b_check_tree_all_connected() == false)  return(false);
	if  (b_check_tree_hanging_leaves() == false)  return(false);

	return(true);
}//bool  CSteinerTreeSet::bCheckUpTree(int  iShift, int  *piArgument,  int  iArgLen)




bool  CSteinerTreeSet::bCheckUpTreeLoops(int  iShift, int  *piArgument,  int  iArgLen)
{
	CSteinerTree  c_buf;
	vCopyTo(&c_buf);

	c_buf.vReadConns(piArgument + iShift);
		
	if  (c_buf.b_check_tree_no_loops() == false)  return(false);

	return(true);
}//bool  CSteinerTreeSet::bCheckUpTree(int  iShift, int  *piArgument,  int  iArgLen)




double  CSteinerTreeSet::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
{
	//CString  s_buf;

	//s_buf.Format("iShift: %d   iArgLen: %d", iShift, iArgLen);
	//::vReportInFile("zzz_dGetFuncValue.txt", s_buf);

	CSteinerTree  c_buf;
	vCopyTo(&c_buf);

	CError  c_err;

	if  (i_optimize > 0)
	{
		c_err = c_buf.eBuildTreeFromConns(piArgument + iShift);
		//if  (c_err)  c_err.vShowWindow();
	}//if  (i_optimize > 0)
	else
	{
		c_buf.vReadConns(piArgument + iShift);
		
		if  (c_buf.b_check_tree_no_loops() == false)  return(0);
		if  (c_buf.b_check_tree_nodes_covered() == false)  return(0);
		if  (c_buf.b_check_tree_all_connected() == false)  return(0);
		if  (c_buf.b_check_tree_hanging_leaves() == false)  return(0);	
	}//else  if  (i_optimize > 0)




	double  d_def_fintes;
	d_def_fintes = 0;
	if  (v_tree_set.size() > 0)
	{
		d_def_fintes = 	v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost();
	}//if  (v_tree_set.size() > 0)

	double  d_result;
	//d_result = c_buf.dGetTreeCost();
	d_result = dGetTreeFitness(&c_buf, d_def_fintes * d_penalty_factor);
	d_result = 1.0 / (d_result + 1);


	CString  s_before;
	bool  b_add_this_tree = false;
	if  (v_tree_set.size()  < i_capacity)  
	{
		b_add_this_tree = true;
		for  (int  ii = 0; ii < v_tree_set.size(); ii++)
		{
			if  (v_tree_set.at(ii)->bGenotypesTheSame(&c_buf) == true)  return(d_result);
			//if  (v_tree_set.at(ii)->dGetTreeCost()  >  c_buf.dGetTreeCost())  b_add_this_tree = true;
		}//for  (int  ii = 0; ii < v_tree_set.size(); ii++)


		if  (b_add_this_tree  ==  true)
		{
			CSteinerTree  *pc_new_tree;
			pc_new_tree = new CSteinerTree();
			c_buf.vCopyTo(pc_new_tree);

			v_tree_set.push_back(pc_new_tree);	
			b_add_this_tree = true;
			s_before = "empty";
		}//if  (b_add_this_tree  ==  true)
		
	}//if  (v_tree_set.size()  <  i_capacity)  
	else
	{
		for  (int  ii = 0; ii < v_tree_set.size(); ii++)
		{
			if  (v_tree_set.at(ii)->bGenotypesTheSame(&c_buf) == true)  return(d_result);
			if  (v_tree_set.at(ii)->dGetTreeCost()  >  c_buf.dGetTreeCost())  b_add_this_tree = true;
		}//for  (int  ii = 0; ii < v_tree_set.size(); ii++)

		s_before.Format("%.2lf", v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost());
		if  (b_add_this_tree  ==  true)
			c_buf.vCopyTo(v_tree_set.at(v_tree_set.size() - 1));
	}//else  if  (v_tree_set.size()  <  i_capacity)  

	if  ( b_add_this_tree  ==  true )
	{
		CString  s_buf;
		vPenaltyFactorReset();

		//s_buf.Format("before: %s in: %.2lf  set:%.2lf", (LPCSTR) s_before, c_buf.dGetTreeCost(), v_tree_set.at(0)->dGetTreeCost());
		//::Tools::vReportInFile("zzz_fit_test.txt", s_buf);
	}//if  ( (d_result == 1013)&&(b_add_this_tree  ==  true) )*/




	if  (b_add_this_tree  ==  true)  v_order_trees();


	return(d_result);
}//double  CSteinerTreeSet::dGetFuncValue(int  iShift, int  *piArgument,  int  iArgLen)
		

double  CSteinerTreeSet::dGetFuncValue(int  iUnitation)
{
	//::vReportInFile("zzz_dGetFuncValue.txt", "dGetFuncValue(int  iUnitation)");
	//return(0);
	throw -1;
}//double  CSteinerTreeSet::dGetFuncValue(int  iUnitation)



double  CSteinerTreeSet::dGetTreeFitness(CSteinerTree  *pcTree, double  dDefaultFitValForTreesInSet)
{
	if  (bTreeContained(pcTree)  ==  true)  return(dDefaultFitValForTreesInSet);

	return(pcTree->dGetTreeCost());
}//double  CSteinerTreeSet::dGetTreeFitness(CSteinerTree  *pcTree, double  dDefaultFitValForTreesInSet)




CError  CSteinerTreeSet::eLoadTreeSet(CString  sTreeSetFile)
{
	CError  c_err;
	CString  s_buf;

	if  (sTreeSetFile  ==  "")
	{
		c_err.vSetError("no tree set file");
		return(c_err);	
	}//if  (sTreeSetFile  ==  "")


	FILE  *pf_tree_set;
	pf_tree_set = fopen(sTreeSetFile, "r");

	if  (pf_tree_set  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sTreeSetFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_dmc_file  ==  NULL)

	c_err = eLoadTreeSet(pf_tree_set);
	
	v_order_trees();

	fclose(pf_tree_set);
	return(c_err);
}//CError  CSteinerTreeSet::eLoadTreeSet(CString  sTreeSetFile)



CError  CSteinerTreeSet::eLoadTreeSet(FILE  *pfTreeSetFile)
{
	CError  c_err;


	int  i_loaded;
	CSteinerTree  *pc_tree_buf;
	while  (!feof(pfTreeSetFile))
	{
		pc_tree_buf = new CSteinerTree();
		vCopyTo(pc_tree_buf);

		pc_tree_buf->eLoadTree(pfTreeSetFile, &i_loaded);

		if  ( (i_loaded == 1)&&(v_tree_set.size() < i_capacity) )
			v_tree_set.push_back(pc_tree_buf);
		else
			delete  pc_tree_buf;	
	}//while  (!feof(pfTreeSetFile))

	return(c_err);
};//CError  CSteinerTreeSet::eLoadTreeSet(FILE  *pfTreeSetFile)


double  CSteinerTreeSet::dGetBestTree()
{
	if  (v_tree_set.size() == 0)  return(0);
	return(v_tree_set.at(0)->dGetTreeCost());
}//double  CSteinerTreeSet::dGetBestTree()


double  CSteinerTreeSet::dGetWorstTree()
{
	if  (v_tree_set.size() == 0)  return(0);
	return(v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost());
}//double  CSteinerTreeSet::dGetWorstTree()




CError  CSteinerTreeSet::eSaveTreeSetSummary(CString  sSummaryFile, CString  sTreeSetname)
{
	CError  c_err;
	CString  s_buf;

	if  (sSummaryFile  ==  "")
	{
		c_err.vSetError("no tree summary file for save defined");
		return(c_err);	
	}//if  (sSummaryFile  ==  "")

	FILE  *pf_tree_set_summary;
	pf_tree_set_summary = fopen(sSummaryFile, "a");

	if  (pf_tree_set_summary  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sSummaryFile);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_tree_set  ==  NULL)

	c_err = eSaveTreeSetSummary(pf_tree_set_summary, sTreeSetname);
	

	fclose(pf_tree_set_summary);
	return(c_err);
}//CError  CSteinerTreeSet::eSaveTreeSetSummary(CString  sSummaryFile)



CError  CSteinerTreeSet::eSaveTreeSetSummary(FILE  *pfSummaryFile, CString  sTreeSetname)
{
	CError  c_err;
	CString  s_report_line;
	CString  s_buf;
	bool  b_tree_ok;
	int  i_error = 0;

	s_buf = "";
	if  (v_tree_set.size() > 0)
	{
		s_buf.Format("%.0lf\t %.0lf\t", v_tree_set.at(0)->dGetTreeCost(), v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost());
	}//if  (v_tree_set.size() > 0)

	if  (i_capacity != v_tree_set.size()) i_error += 1;
	
	b_tree_ok = true;
	for  (int  ii = 0; (ii < v_tree_set.size())&&(b_tree_ok == true); ii++)
	{
		if  (v_tree_set.at(ii)->bCheckTree() == false)  b_tree_ok = false;
	}//for  (int  ii = 0; ii < v_tree_set.size(); ii++)

	if  (b_tree_ok == false) i_error += 10;

	s_report_line.Format("%s\t%d\t%d\t%s\n", sTreeSetname, i_error, v_tree_set.size(), s_buf);
	fprintf(pfSummaryFile, s_report_line);


	return(c_err);
}//CError  CSteinerTreeSet::eSaveTreeSetSummary(FILE  *pfSummaryFile)



CError  CSteinerTreeSet::eSaveTreeSet(CString  sTreeFileSet)
{
	CError  c_err;
	CString  s_buf;

	if  (sTreeFileSet  ==  "")
	{
		c_err.vSetError("no tree file for save defined");
		return(c_err);	
	}//if  (sTreeFileSet  ==  "")


	FILE  *pf_tree_set;
	pf_tree_set = fopen(sTreeFileSet, "w+");

	if  (pf_tree_set  ==  NULL)
	{
		s_buf.Format("Couldnt open file: '%s'", sTreeFileSet);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (pf_tree_set  ==  NULL)

	c_err = eSaveTreeSet(pf_tree_set);
	

	fclose(pf_tree_set);
	return(c_err);
}//CError  CSteinerTreeSet::eSaveTreeSet(CString  sTreeFileSet)



CError  CSteinerTreeSet::eSaveTreeSet(FILE  *pfTreeFileSet)
{
	CError  c_err;
	CString  s_buf;

	s_buf.Format("tree number: %d\n", v_tree_set.size());
	fprintf(pfTreeFileSet, s_buf);

	for  (int  ii = 0; ii < v_tree_set.size(); ii++)
	{
		fprintf(pfTreeFileSet, "%.0lf\t", v_tree_set.at(ii)->dGetTreeCost());
		c_err = v_tree_set.at(ii)->eSaveTree(pfTreeFileSet);
		if  (c_err)  return(c_err);


		if  (v_tree_set.at(ii)->bCheckTree() == false)
		{
			//v_tree_set.at(ii)->b_tree_cost_actual = false;
			s_buf.Format("%.0lf.txt", v_tree_set.at(ii)->dGetTreeCost());
			c_err = v_tree_set.at(ii)->eSaveTree(s_buf, true);

			if  (c_err)  return(c_err);
		}//if  (v_tree_set.at(ii)->bCheckTree() == false)


	}//for  (int  ii = 0; ii < v_tree_set.size(); ii++)

	return(c_err);
};//CError  CSteinerTreeSet::eSaveTreeSet(FILE  *pfTreeFileSet)



void  CSteinerTreeSet::vGetShortestPathBetweenNodePair(int  iStartNode, int iEndNode, CSSJConn *pcResult, bool  bUseRealNodeCosts)
{
	for  (int  ii = 0;  ii < v_shortest_paths.size(); ii++)
	{
		if  (
			(v_shortest_paths.at(ii).iNodeStartReal == iStartNode)&&
			(v_shortest_paths.at(ii).iNodeEndReal== iEndNode)
			)
		{
			pcResult->dCost = v_shortest_paths.at(ii).dCost;
			pcResult->v_path = v_shortest_paths.at(ii).v_path;
			return;
		}//if  (
	}//for  (int  ii = 0;  ii < v_shortest_paths.size(); ii++)
	

	CFiberConnection  ***pc_conn_mat;

	pc_conn_mat = pc_base_network->pcGetConnMat();
	vector<CSearchShortestPathNode *>  v_chosen_nodes;
	vector<CSearchShortestPathNode *>  v_candidates;
	vector<CSearchShortestPathNode *>  v_candidates_new;

	
	CSearchShortestPathNode  c_root;
	c_root.iNode = iStartNode;

	v_candidates.push_back(&c_root);

	CSearchShortestPathNode  *pc_chosen_node;
	bool  b_found;
	b_found = false;

	while  ((b_found  ==  false)  && (v_candidates.size() > 0) )
	{
		pc_chosen_node = v_candidates.at(0);
		v_candidates.erase(v_candidates.begin());

		if  (pc_chosen_node->iNode == iEndNode)
		{
			if  (iEndNode == 4)
			{
				int  ig = 0;
				ig++;
			}//if  (iEndNode == 4)
			
			pc_chosen_node->vCreateConn(pcResult);
			//pcResult->iNodeStart = iStartNode;
			//pcResult->iNodeEnd = iEndNode;

			v_shortest_paths.push_back(*pcResult);
			
			return;
		}//if  (pc_chosen_node->iNode == iEndNode)

		pc_chosen_node->vCreateCandidates(&v_chosen_nodes, &v_candidates_new, pc_conn_mat, pc_base_network->iGetNodesNumber(), bUseRealNodeCosts);
		v_chosen_nodes.push_back(pc_chosen_node);

		for  (int  ii = 0; ii < v_candidates_new.size(); ii++)
			v_add_and_order_by_cost(&v_candidates, v_candidates_new.at(ii));
		v_candidates_new.clear();
	}//while  ((b_found  ==  false)  && (v_candidates.size() > 0) )

}//void  CSteinerTreeSet::vGetShortestPathBetweenNodePair()


void  CSteinerTreeSet::v_add_and_order_by_cost(vector<CSearchShortestPathNode *>  *pvNodes, CSearchShortestPathNode *pcNodeToAdd)
{
	if  (pvNodes->size() == 0)
	{
		pvNodes->push_back(pcNodeToAdd);
		return;
	}//if  (pvNodes->size() == 0)

	for  (int  ii = 0; ii < pvNodes->size(); ii++)
	{
		if  (pvNodes->at(ii)->dCost > pcNodeToAdd->dCost)
		{
			pvNodes->insert(pvNodes->begin() + ii, pcNodeToAdd);
			return;
		}//if  (pvNodes->at(ii)->dCost > pcNodeToAdd->dCost)
	}//for  (int  ii = 0; ii < pvNodes->size(); ii++)

	pvNodes->push_back(pcNodeToAdd);
}//void  CSteinerTreeSet::v_add_and_order_by_cost(vector<CSearchShortestPathNode *>  *pvNodes, CSearchShortestPathNode *pcNodeToAdd)





void  CSteinerTreeSet::v_insert_spanning_tree(vector<CSSJNetwork *>  *pvSpanningTrees, CSSJNetwork *pcNewTree)
{
	if  (pvSpanningTrees->size() == 0)
		pvSpanningTrees->push_back(pcNewTree);
	else
	{
		for  (int  ii = 0; ii < pvSpanningTrees->size(); ii++)
		{
			if  (pvSpanningTrees->at(ii)->dGetCost() > pcNewTree->dGetCost())
			{
				pvSpanningTrees->insert(pvSpanningTrees->begin() + ii, pcNewTree);
				return;
			}//if  (pvSpanningTrees->at(ii)->dGetCost() > pcNewTree->dGetCost())
		}//for  (int  ii = 0; ii < pvSpanningTrees->size(); ii++)

		pvSpanningTrees->push_back(pcNewTree);
	}//else  if  (v_spanning_trees.size() == 0)

}//void  CSteinerTreeSet::v_insert_spanning_tree(vector<CSSJNetwork *>  *pvSpanningTrees, CSSJNetwork *pcNewTree)



void  CSteinerTreeSet::v_create_all_possible_r_graphs_of_size_step
	(
	vector<vector<int>>  *pvNodesToInclude,
	vector<int>  *pvNodesToChoose, vector<int>  *pvChosenNodes, int  iSize, int iCurPos, int iStep
	)
{
	if  (iSize == pvChosenNodes->size())
	{
		pvNodesToInclude->push_back(*pvChosenNodes);
		return;
	}//if  (iSize == pvChosenNodes->size())

	for  (int  ii = iCurPos; ii < pvNodesToChoose->size() /*- (iSize - 1 - iStep)*/; ii++)
	{
		if  (pvNodesToChoose->at(ii)  !=  i_root)
		{
			pvChosenNodes->push_back(pvNodesToChoose->at(ii));
			v_create_all_possible_r_graphs_of_size_step(pvNodesToInclude, pvNodesToChoose, pvChosenNodes, iSize, ii + 1, iStep + 1);
			pvChosenNodes->pop_back();
		}//if  (pvNodesToChoose->at(ii)  !=  i_root)
	}//for  (int  ii = 0; ii < vNodesToChoose.size() - (iSize - 1); ii++)

}//void  CSteinerTreeSet::v_create_all_possible_r_graphs_of_size_step(vector<vector<int>>  *pvNodesToInclude, vector<int>  *pvNodesToChoose, int  iSize, int iCurPos, int iStep)



void  CSteinerTreeSet::v_create_all_possible_r_graphs_of_size(vector<vector<int>>  *pvNodesToInclude, vector<int>  *pvNodesToChoose, int  iSize)
{
	vector<int>  v_chosen_nodes;

	v_create_all_possible_r_graphs_of_size_step(pvNodesToInclude, pvNodesToChoose, &v_chosen_nodes, iSize, 0, 0);
	
}//void  CSteinerTreeSet::v_create_all_possible_r_graphs_of_size(vector<vector<int>>  *pvNodesToInclude, int  iSize)




//void  CSteinerTreeSet::vGenerateTreesByBrazilians(double  dMaxTime, System::Windows::Forms::ListBox*  listInfo, double *pdFfe, bool  bUseRealNodeCosts)
//{
//	bool  b_debug = false;
//
//	CString  s_buf;
//	int  i_no_tree_increase, i_trees_before_new;
//	CTimeCounter  c_time_counter;
//	c_time_counter.vSetStartNow();
//	double  d_time_passed;
//
//	CString  s_debug_graphs, s_debug_trees;
//	s_debug_graphs = "C:\\Projekty\\projekty_pwr\\VmGA_Steiner_trees\\release\\zz_graph.txt";
//	s_debug_trees = "C:\\Projekty\\projekty_pwr\\VmGA_Steiner_trees\\release\\zz_span_tree.txt";
//
//
//	int  *pi_conn_based_genotype;
//	pi_conn_based_genotype = new int[pc_base_network->iGetConnsNumber()];
//
//	vector<vector<int>>  v_nodes_to_include;
//
//	i_no_tree_increase = 0;
//	for  (int i_size = 0; i_size <= v_demanded_nodes.size() - 1 + 1; i_size++)
//	//for  (int i_size = 2; i_size <= 2; i_size++)
//	{
//		s_debug_graphs.Format("C:\\Projekty\\projekty_pwr\\VmGA_Steiner_trees\\release\\zz_graph_%.2d.txt", i_size);
//		s_debug_trees.Format("C:\\Projekty\\projekty_pwr\\VmGA_Steiner_trees\\release\\zz_span_tree_%.2d.txt", i_size);
//
//		if  (i_size <= v_not_demanded_nodes.size())
//		{
//			v_nodes_to_include.clear();
//			v_create_all_possible_r_graphs_of_size(&v_nodes_to_include, &v_not_demanded_nodes, i_size);
//
//			
//
//			vector<CSSJNetwork *>  v_nets_to_analyze;
//			CSSJNetwork  *pc_ssj_net_buf;
//
//			FILE  *pf_graph_test;
//			if  (b_debug)  pf_graph_test = fopen(s_debug_graphs, "w+");
//			for  (int  i_graph = 0; i_graph < v_nodes_to_include.size(); i_graph++)
//			{
//				if  (v_nodes_to_include.at(i_graph).size() > 0)
//				{
//					if  (v_nodes_to_include.at(i_graph).at(0) == 4)
//					{
//						int ig = 0;
//						ig++;
//					}//if  (v_nodes_to_include.at(i_graph).at(0) == 4)
//				}//if  (v_nodes_to_include.at(i_graph).size() > 0)
//					
//				
//
//				pc_ssj_net_buf = new CSSJNetwork();
//				pc_ssj_net_buf->vCreateCompeleteGraph(&(v_nodes_to_include.at(i_graph)), this, i_root, bUseRealNodeCosts);
//				v_nets_to_analyze.push_back(pc_ssj_net_buf);
//
//				if  (b_debug)
//				{
//					fprintf(pf_graph_test, "ADDITIONAL NODES: ");
//					for  (int  ii = 0; ii < v_nodes_to_include.at(i_graph).size(); ii++)
//						fprintf(pf_graph_test, "%d ", v_nodes_to_include.at(i_graph).at(ii));
//					fprintf(pf_graph_test, "\n");
//					
//
//					fprintf(pf_graph_test, "GRAPH: %d\n", i_graph);
//					pc_ssj_net_buf->vSaveGraph(pf_graph_test);
//					fprintf(pf_graph_test, "\n\n");
//				}//if  (b_debug)
//				
//			}//for  (int  i_graph = 0; i_graph < v_nodes_to_include.size(); i_graph++)
//			if  (b_debug)  fclose(pf_graph_test);
//
//
//			CSSJNetwork  *pc_cur_tree;
//			vector<CSSJNetwork *>  v_spanning_trees;
//
//
//			FILE  *pf_span_tree_test;
//			if  (b_debug)  pf_span_tree_test = fopen(s_debug_trees, "w+");
//			if  (b_debug)  fclose(pf_span_tree_test);
//			
//			int  i_span_tree_counter = 0;
//			int  i_graphs_counter = v_nets_to_analyze.size();
//
//			while  ( (v_nets_to_analyze.size() > 0)||(v_spanning_trees.size() > 0) )
//			{
//				while  (v_nets_to_analyze.size() > 0)
//				{
//					c_time_counter.bGetTimePassed(&d_time_passed);
//					if  (dMaxTime  <  d_time_passed)  return;
//
//					pc_ssj_net_buf = v_nets_to_analyze.at(0);
//					v_nets_to_analyze.erase(v_nets_to_analyze.begin());
//					
//
//					pc_cur_tree = new CSSJNetwork();
//					if  (pc_ssj_net_buf->bFindMinSpaningTreePrim(pc_cur_tree) == true)
//					{
//						if  (b_debug)  
//						{
//							pf_span_tree_test = fopen(s_debug_trees, "a");
//							fprintf(pf_span_tree_test, "SPAN TREE: %d\n", i_span_tree_counter);
//							pc_cur_tree->vSaveGraph(pf_span_tree_test);
//							fprintf(pf_span_tree_test, "\n\n");
//							fclose(pf_span_tree_test);
//						}//if  (b_debug)  
//
//						i_span_tree_counter++;
//						pc_cur_tree->pc_parent_network = pc_ssj_net_buf;
//						//now insert new spanning tree to spanning tree list so we know from which network we shall generate next candidates
//						v_insert_spanning_tree(&v_spanning_trees, pc_cur_tree);
//
//
//						for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
//							pi_conn_based_genotype[ii] = 0;
//						pc_cur_tree->vFillGenotype(pi_conn_based_genotype);
//
//						i_trees_before_new = v_tree_set.size();
//						dGetFuncValue(0,pi_conn_based_genotype,  pc_base_network->iGetConnsNumber());
//						(*pdFfe)++;
//
//						if  (i_trees_before_new == v_tree_set.size())  
//							i_no_tree_increase++;
//						else
//						{
//							if  (listInfo  !=  NULL)
//							{
//								c_time_counter.bGetTimePassed(&d_time_passed);
//								s_buf.Format
//									(
//									"[tree add] Size: %d  graphs:%d [time passed:%.4lf] ComputedShortestPaths: %d  Trees: %d  [%.2lf  -  %.2lf] AllSpanTrees: %d AllGraphs: %d",
//									i_size, v_nodes_to_include.size(),
//									d_time_passed,
//									v_shortest_paths.size(),
//									v_tree_set.size(),
//									v_tree_set.at(0)->dGetTreeCost(), v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost(),
//									i_span_tree_counter, i_graphs_counter
//									);
//								listInfo->Items->Add((String *) s_buf);
//								listInfo->SelectedIndex  =  listInfo->Items->Count - 1;
//								listInfo->Refresh();
//
//								System::Windows::Forms::Application::DoEvents();
//							}//if  (listInfo  !=  NULL)
//						}//else  if  (i_trees_before_new == v_tree_set.size())  
//					}//if  (pc_ssj_net_buf->bFindMinSpaningTreePrim(pc_cur_tree) == true)
//					else
//						delete  pc_cur_tree;
//				}//while  (v_nets_to_analyze.size() > 0)
//
//
//				if  (v_spanning_trees.size() > 0)
//				{
//					pc_cur_tree = v_spanning_trees.at(0);
//					v_spanning_trees.erase(v_spanning_trees.begin());
//
//					pc_cur_tree->vGenerateNetsToAnalyze(&v_nets_to_analyze);
//
//
//					if  (b_debug)  
//					{
//						FILE  *pf_graph_test;
//						pf_graph_test = fopen(s_debug_graphs, "a");
//						for  (int  i_graph = 0; i_graph < v_nets_to_analyze.size(); i_graph++)
//						{
//							fprintf(pf_graph_test, "GRAPH: %d\n", i_graphs_counter + i_graph);
//							v_nets_to_analyze.at(i_graph)->vSaveGraph(pf_graph_test);
//							fprintf(pf_graph_test, "\n\n");
//						}//for  (int  i_graph = 0; i_graph < v_nodes_to_include.size(); i_graph++)
//						fclose(pf_graph_test);//*/
//					}//if  (b_debug)  
//
//
//					i_graphs_counter += v_nets_to_analyze.size();
//					delete  pc_cur_tree;
//				}//if  (v_spanning_trees.size() > 0)
//			}//while  ( (v_nets_to_analyze.size() > 0)||(v_spanning_trees.size() > 0) )
//
//
//			if  (listInfo  !=  NULL)
//			{
//				c_time_counter.bGetTimePassed(&d_time_passed);
//				s_buf.Format
//					(
//					"Size: %d  graphs:%d [time passed:%.4lf] NoTreeIncrease: %d ComputedShortestPaths: %d  Trees: %d  [%.2lf  -  %.2lf]",
//					i_size, v_nodes_to_include.size(),
//					d_time_passed,
//					i_no_tree_increase, v_shortest_paths.size(),
//					v_tree_set.size(),
//					v_tree_set.at(0)->dGetTreeCost(), v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost()
//					);
//				listInfo->Items->Add((String *) s_buf);
//				listInfo->SelectedIndex  =  listInfo->Items->Count - 1;
//				listInfo->Refresh();
//
//				System::Windows::Forms::Application::DoEvents();
//			}//if  (listInfo  !=  NULL)
//		
//		}//if  (i_size <= v_not_demanded_nodes.size())	
//	}//for  (int i_size = 0; i_size < v_demanded_nodes.size() - 2; i_size++)
//
//	
//	delete  pi_conn_based_genotype;
//}//void  CSteinerTreeSet::vGenerateTreesByBrazilians(double  dMaxTime, System::Windows::Forms::ListBox*  listInfo, double *pdFfe)
//
//
//
//
//void  CSteinerTreeSet::vGenerateTreesBySenOnSorensenJansens(double  dMaxTime, System::Windows::Forms::ListBox*  listInfo, double *pdFfe, bool  bUseRealNodeCosts)
//{
//	CString  s_buf;
//
//	//first create a complete network of demanded nodes
//	CSSJNetwork  *pc_ssj_net_buf;
//	CSSJNetwork  *pc_cur_tree;
//	vector<CSSJNetwork *>  v_networks_to_analyze;
//	vector<CSSJNetwork *>  v_spanning_trees;
//	double  d_time_passed;
//
//	int  *pi_conn_based_genotype;
//	pi_conn_based_genotype = new int[pc_base_network->iGetConnsNumber()];
//	
//	pc_ssj_net_buf = new CSSJNetwork();
//	pc_ssj_net_buf->vCreateCompeleteGraph(this, i_root, bUseRealNodeCosts);
//	v_networks_to_analyze.push_back(pc_ssj_net_buf);
//
//	CTimeCounter  c_time_counter;
//	c_time_counter.vSetStartNow();
//
//	int  i_trees_before_new, i_no_tree_increase;
//
//	i_no_tree_increase = 0;
//	while  
//		( 
//		( (v_spanning_trees.size() > 0)||(v_networks_to_analyze.size() > 0) )
//		&&
//		(v_tree_set.size() < i_capacity) 
//		&&
//		(d_time_passed < dMaxTime)
//		)
//	{
//		if  (v_networks_to_analyze.size()  >  0)
//		{
//			pc_ssj_net_buf = v_networks_to_analyze.at(0);
//			v_networks_to_analyze.erase(v_networks_to_analyze.begin());
//
//			pc_cur_tree = new CSSJNetwork();
//			
//			if  (pc_ssj_net_buf->bFindMinSpaningTreePrim(pc_cur_tree) == true)
//			{
//				for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
//					pi_conn_based_genotype[ii] = 0;
//				pc_cur_tree->vFillGenotype(pi_conn_based_genotype);
//
//				i_trees_before_new = v_tree_set.size();
//				dGetFuncValue(0,pi_conn_based_genotype,  pc_base_network->iGetConnsNumber());
//				(*pdFfe)++;
//
//				if  (i_trees_before_new == v_tree_set.size())  
//				{
//					i_no_tree_increase++;
//				}//if  (i_trees_before_new == v_tree_set.size())  
//				else
//				{
//					pc_cur_tree->pc_parent_network = pc_ssj_net_buf;
//					//now insert new spanning tree to spanning tree list so we know from which network we shall generate next candidates
//					v_insert_spanning_tree(&v_spanning_trees, pc_cur_tree);				
//				}//else  if  (i_trees_before_new == v_tree_set.size())  
//
//								
//			}//if  (c_ssj_net.bFindMinSpaningTreePrim(&c_ssj_tree) == true)
//			else
//				delete  pc_cur_tree;
//		}//if  (v_networks_to_analyze.size()  >  0)
//		else
//		{
//			if  (v_spanning_trees.size() > 0)
//			{
//				pc_cur_tree = v_spanning_trees.at(0);
//				v_spanning_trees.erase(v_spanning_trees.begin());
//
//				pc_cur_tree->vGenerateNetsToAnalyze(&v_networks_to_analyze);
//				delete  pc_cur_tree;
//			}//if  (v_spanning_trees.size() > 0)
//
//		}//else  if  (v_networks_to_analyze.size()  >  0)	
//
//		if  (listInfo  !=  NULL)
//		{
//			c_time_counter.bGetTimePassed(&d_time_passed);
//			s_buf.Format
//				(
//				"Trees: %d  [%.2lf  -  %.2lf] NetsToAnalyze: %d  SpanningTrees: %d  NoTreeIncrease: %d [time passed:%.4lf]",
//				v_tree_set.size(),
//				v_tree_set.at(0)->dGetTreeCost(), v_tree_set.at(v_tree_set.size() - 1)->dGetTreeCost(),
//				v_networks_to_analyze.size(), v_spanning_trees.size(),
//				i_no_tree_increase,
//				d_time_passed
//				);
//			listInfo->Items->Add((String *) s_buf);
//			listInfo->SelectedIndex  =  listInfo->Items->Count - 1;
//			listInfo->Refresh();
//
//			System::Windows::Forms::Application::DoEvents();
//		}//if  (listInfo  !=  NULL)
//	}//while  ( (v_networks_to_analyze.size() > 0)&&(v_tree_set.size() < i_capacity) )
//
//
//	
//	
//	delete  pi_conn_based_genotype;
//}//void  CSteinerTreeSet::vGenerateTreesBySenOnSorensenJansens()




void  CSteinerTreeSet::vGenerateTreesBySenOnSorensenJansensTest(bool  bUseRealNodeCosts)
{
	//first create a complete network of demanded nodes
	CSSJNetwork  c_ssj_net;
	CSSJNetwork  c_ssj_tree;

	c_ssj_net.vCreateCompeleteGraph(this, i_root, bUseRealNodeCosts);

	FILE  *pf_test;
	pf_test = fopen("z_complete.txt", "w+");
	c_ssj_net.vSaveGraph(pf_test);
	fclose(pf_test);
	
	c_ssj_net.bFindMinSpaningTreePrim(&c_ssj_tree);


	pf_test = fopen("z_MinSpanTree.txt", "w+");
	c_ssj_tree.vSaveGraph(pf_test);
	fclose(pf_test);


	int  *pi_conn_based_genotype;
	pi_conn_based_genotype = new int[pc_base_network->iGetConnsNumber()];

	for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
		pi_conn_based_genotype[ii] = 0;
	c_ssj_tree.vFillGenotype(pi_conn_based_genotype);

	double  d_val;
	d_val = dGetFuncValue(0,pi_conn_based_genotype,  pc_base_network->iGetConnsNumber());

	if  (v_tree_set.size() > 0)  v_tree_set.at(0)->eSaveTree("z_tree.txt", true);

	pf_test = fopen("z_tree_summary.txt", "w+");

	CString  s_line, s_buf;
	s_line = "";
	for  (int  ii = 0; ii < pc_base_network->iGetConnsNumber(); ii++)
	{
		s_buf.Format("%d\t", pi_conn_based_genotype[ii]);
		s_line += s_buf;
	}//for  (int  ii = 0; ii < i_bit_length; ii++)

	fprintf(pf_test, "%s\n\n", s_line);
	

	eSaveTreeSet(pf_test);
	fclose(pf_test);
	

	delete  pi_conn_based_genotype;
	


}//void  CSteinerTreeSet::vGenerateTreesBySenOnSorensenJansens()



ProblemTools::CCompProblem *CSteinerTreeSet::pcClone()
{
	CSteinerTreeSet  *pc_result;

	pc_result = new CSteinerTreeSet();


	pc_result->d_max_val = d_max_val;
	pc_result->d_min_val = 0;
	pc_result->d_lower_than_max_scaling = 1;

	pc_result->pc_base_network = new CFiberNet();
	pc_base_network->eCloneTopologyTo(pc_result->pc_base_network);
	pc_result->vSetBaseNetwork(pc_result->pc_base_network);//TO JEST METODA ZE CSteinerTree, ktore nie jest wlascicielem topologii, wiec nie kasuje topologii, Wlascicielem topologii jest tylko CSteinerTreeSet

	//pc_base_network->eSaveTopology("test_net.txt");
	//pc_result->pc_base_network->eSaveTopology("test_net_cloned.txt");
	//::vShow("Topology clone!");

	

	pc_result->i_root = i_root;

	pc_result->v_demanded_nodes  =  v_demanded_nodes;
	pc_result->v_not_demanded_nodes  =  v_not_demanded_nodes;
	
	pc_result->vSetCapacity(i_capacity);
	i_bit_length  =  pc_base_network->iGetConnsNumber();

	pc_result->i_optimize = i_optimize;


	return(pc_result);
}//CDecFunc *CSteinerTreeSet::pcClone()




CError  CSteinerTreeSet::eLoadSettings(FILE *pfSettings)
{
	CError  c_err;
	CString  s_setting_name,  s_buf;


	CString  s_topology_file, s_tree_demands_file;
	int  i_capacity_buf;


	

	s_topology_file  =  Tools::sReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  STEINER_TREE_NETWORK)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", STEINER_TREE_NETWORK, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  STEINER_TREE_NETWORK)



	if  (pc_base_network  != NULL)  delete  pc_base_network;
	pc_base_network = new CFiberNet();
	c_err = pc_base_network->eLoadTopology(s_topology_file);
	if  (c_err)  return(c_err);
	vSetBaseNetwork(pc_base_network);//TO JEST METODA ZE CSteinerTree, ktore nie jest wlascicielem topologii, wiec nie kasuje topologii, Wlascicielem topologii jest tylko CSteinerTreeSet

	

	s_tree_demands_file  =  Tools::sReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  STEINER_TREE_TREE_DEMANDS)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", STEINER_TREE_TREE_DEMANDS, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  STEINER_TREE_TREE_DEMANDS)
	
	c_err = this->eLoadTreeDemands(s_tree_demands_file);
	if  (c_err)  return(c_err);



	i_capacity_buf  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  STEINER_TREE_TREE_NUMBER)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", STEINER_TREE_TREE_NUMBER, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  STEINER_TREE_TREE_NUMBER)
	vSetCapacity(i_capacity_buf);


	i_optimize  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  STEINER_TREE_OPTIMIZE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", STEINER_TREE_OPTIMIZE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  STEINER_TREE_TREE_NUMBER)*/


	i_bit_length  =  pc_base_network->iGetConnsNumber();

	/*i_bit_length  =  Tools::iReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  DEC_PARAM_DEC_BIT_LEN)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", DEC_PARAM_DEC_BIT_LEN, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	
	d_max_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  DEC_PARAM_DEC_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", DEC_PARAM_DEC_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_min_val  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  DEC_PARAM_DEC_MIN_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", DEC_PARAM_DEC_MIN_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)


	d_lower_than_max_scaling  =  Tools::dReadLine(pfSettings, &s_setting_name);
	if  (s_setting_name  !=  DEC_PARAM_DEC_LOWER_THAN_MAX_SCALE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", DEC_PARAM_DEC_LOWER_THAN_MAX_SCALE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);	
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)



	c_err  =  eConfigure(i_bit_length,  d_max_val,  d_min_val, d_lower_than_max_scaling);*/
		

	return(c_err);
};//CError  CSteinerTreeSet::eLoadSettings(FILE *pfSettings)



//---------------------------class  CSearchShortestPathNode--------------------------------------------------
CSearchShortestPathNode::CSearchShortestPathNode()
{
	iNode = -1;
	dCost = 0;
	pc_get_to_node_connection = NULL;
	pc_parent = NULL;
}//CSearchShortestPathNode::CSearchShortestPathNode()


CSearchShortestPathNode::~CSearchShortestPathNode()
{
	for  (int  ii = 0; ii < v_children.size(); ii++)
		delete  v_children.at(ii);

}//CSearchShortestPathNode::~CSearchShortestPathNode()


void  CSearchShortestPathNode::vCreateConn(CSSJConn  *pcResult)
{
	pcResult->v_path.clear();
	//pcResult->iNodeEnd = iNode;
	pcResult->dCost = dCost;
	if  (pc_get_to_node_connection != NULL)  pcResult->v_path.push_back(pc_get_to_node_connection);

	CSearchShortestPathNode  *pc_node;
	pc_node = pc_parent;
	
	while  (pc_node !=  NULL)
	{
		if  (pc_node->pc_get_to_node_connection != NULL)  pcResult->v_path.push_back(pc_node->pc_get_to_node_connection);
		//pcResult->iNodeStart = pc_node->iNode;
		pc_node = pc_node->pc_parent;		
	}//while  (pc_node !=  NULL)
}//void  CSearchShortestPathNode::vCreateConn(CSSJConn  *pcResult)



void  CSearchShortestPathNode::vCreateCandidates
	(
	vector<CSearchShortestPathNode *> *pvChosenNodes, vector<CSearchShortestPathNode *> *pvNewNodes, 
	CFiberConnection  ***pcConnMat, int iNodeNumber, 
	bool  bUseRealNodeCosts
	)
{
	CSearchShortestPathNode *pc_candidate_node;

	for  (int  ii = 0; ii < iNodeNumber; ii++)
	{
		if  (pcConnMat[iNode][ii]  !=  NULL)
		{
			if  (pcConnMat[iNode][ii]->dGetCapacity()  >  0)
			{
				if  (
					(b_check_if_node_is_contained(pvChosenNodes, ii)  ==  false)&&
					(b_check_if_node_is_contained(pvNewNodes, ii)  ==  false)
					)
				{
					pc_candidate_node = new CSearchShortestPathNode();

					if  (bUseRealNodeCosts  ==  true)
						pc_candidate_node->dCost = this->dCost + pcConnMat[iNode][ii]->dGetCapacity();
					else
						pc_candidate_node->dCost = this->dCost + 1;

					pc_candidate_node->iNode = ii;
					pc_candidate_node->pc_get_to_node_connection = pcConnMat[iNode][ii];
					pc_candidate_node->pc_parent = this;

					v_children.push_back(pc_candidate_node);
					pvNewNodes->push_back(pc_candidate_node);
				}//if  (			
			}//if  (pcConnMat[iNode][ii]  >  0)
		}//if  (pcConnMat[iNode][ii]  !=  NULL)
	}//for  (int  ii = 0; ii < iNodeNumber; ii++)


}//void  CSearchShortestPathNode::vCreateCandidates(vector<CSearchShortestPathNode *> *pvChosenNodes, vector<CSearchShortestPathNode *> *pvNewNodes, CFiberConnection  ***pcConnMat)



bool  CSearchShortestPathNode::b_check_if_node_is_contained(vector<CSearchShortestPathNode *> *pvNodes, int iNodeOffset)
{
	for  (int ii = 0; ii < pvNodes->size(); ii++)
	{
		if  (pvNodes->at(ii)->iNode  ==  iNodeOffset)  return(true);
	}//for  (int ii = 0; ii < pvNodes->size(); ii++)

	return(false);
}//void  CSearchShortestPathNode::v_check_if_node_is_contained(vector<CSearchShortestPathNode *> *pvNodes, int iNodeOffset)


//---------------------------class  CSSJNode--------------------------------------------------

double  CSSJNode::dGetCost(CSSJNode  *pcNode)
{
	for  (int i_conn = 0; i_conn < v_conns.size(); i_conn++)
	{
		if  (
				(
				(v_conns.at(i_conn)->iNodeStartOffset  ==  i_node_offset)&&
				(v_conns.at(i_conn)->iNodeEndOffset  ==  pcNode->i_node_offset)
				)
			)
		{
			return(v_conns.at(i_conn)->dCost);
		}//if  (
	}//for  (int i_conn = 0; i_conn < v_conns.size(); ii++)

	return(0);
}//double  CSSJNode::dGetCost(CSSJNode  *pcNode)

//---------------------------class  CSSJNetwork--------------------------------------------------
CSSJNetwork::CSSJNetwork()
{
	pc_parent_network = NULL;
};//CSSJNetwork::CSSJNetwork()


CSSJNetwork::~CSSJNetwork()
{
	for  (int  ii = 0;  ii < v_conns.size();  ii++)
		delete  v_conns.at(ii);
	
	if  (pc_parent_network  !=  NULL)  delete  pc_parent_network;
};//CSSJNetwork::~CSSJNetwork()




void  CSSJNetwork::vCreateCompeleteGraph(vector<int>  *pvNodesToInclude, CSteinerTreeSet  *pcSTreeSet, int iRoot, bool  bUseRealNodeCosts)
{
	CSSJNode  c_node_buf;

	c_node_buf.i_node_real_net_id = iRoot;	
	c_node_buf.i_node_offset = 0;
	v_nodes.push_back(c_node_buf);

	for  (int  ii = 0; ii <  pcSTreeSet->v_demanded_nodes.size(); ii++)
	{
		c_node_buf.i_node_real_net_id = pcSTreeSet->v_demanded_nodes.at(ii);	
		c_node_buf.i_node_offset = v_nodes.size();
		v_nodes.push_back(c_node_buf);
	}//for  (int  ii = 0; ii <  v_demanded_nodes.size(); ii++)


	for  (int  ii = 0; ii <  pvNodesToInclude->size(); ii++)
	{
		if  (pvNodesToInclude->at(ii) != iRoot)
		{
			c_node_buf.i_node_real_net_id = pvNodesToInclude->at(ii);
			c_node_buf.i_node_offset = v_nodes.size();
			v_nodes.push_back(c_node_buf);
		}//if  (pvNodesToInclude->at(ii) != iRoot)
	}//for  (int  ii = 0; ii <  v_demanded_nodes.size(); ii++)




	CSSJConn  *pc_conn_buf;
	for  (int  i_node_start = 0; i_node_start < v_nodes.size(); i_node_start++)
	{
		for  (int  i_node_end = 0; i_node_end < v_nodes.size(); i_node_end++)
		{
			if  (i_node_start != i_node_end)
			{
				pc_conn_buf = new CSSJConn();

				pc_conn_buf->iNodeStartOffset = i_node_start;
				pc_conn_buf->iNodeEndOffset = i_node_end;

				pc_conn_buf->iNodeStartReal = v_nodes.at(pc_conn_buf->iNodeStartOffset).i_node_real_net_id;
				pc_conn_buf->iNodeEndReal = v_nodes.at(pc_conn_buf->iNodeEndOffset).i_node_real_net_id;

				pcSTreeSet->vGetShortestPathBetweenNodePair
					(
					v_nodes.at(pc_conn_buf->iNodeStartOffset).i_node_real_net_id, 
					v_nodes.at(pc_conn_buf->iNodeEndOffset).i_node_real_net_id, 
					pc_conn_buf,
					bUseRealNodeCosts
					);
				v_conns.push_back(pc_conn_buf);

				v_nodes.at(i_node_start).v_conns.push_back(pc_conn_buf);
				//v_nodes.at(i_node_end).v_conns.push_back(pc_conn_buf);
			}//if  (i_node_start != i_node_end)
		}//for  (int  i_node_end = i_node_start + 1; i_node_end < v_full_network.size() - 1; i_node_end++)
	}//for  (int  ii = 0; ii < v_full_network.size(); ii++)
}//void  CSSJNetwork::vCreateCompeleteGraph(vector<int>  *pvNodesToInclude, CSteinerTreeSet  *pcSTreeSet, int iRoot)




void  CSSJNetwork::vCreateCompeleteGraph(CSteinerTreeSet  *pcSTreeSet, int iRoot, bool  bUseRealNodeCosts)
{
	CSSJNode  c_node_buf;

	c_node_buf.i_node_real_net_id = iRoot;	
	c_node_buf.i_node_offset = 0;
	v_nodes.push_back(c_node_buf);

	for  (int  ii = 0; ii <  pcSTreeSet->v_demanded_nodes.size(); ii++)
	{
		c_node_buf.i_node_real_net_id = pcSTreeSet->v_demanded_nodes.at(ii);	
		c_node_buf.i_node_offset = ii + 1;
		v_nodes.push_back(c_node_buf);
	}//for  (int  ii = 0; ii <  v_demanded_nodes.size(); ii++)


	CSSJConn  *pc_conn_buf;
	for  (int  i_node_start = 0; i_node_start < v_nodes.size(); i_node_start++)
	{
		for  (int  i_node_end = 0; i_node_end < v_nodes.size(); i_node_end++)
		{
			if  (i_node_start != i_node_end)
			{
				pc_conn_buf = new CSSJConn();

				pc_conn_buf->iNodeStartOffset = i_node_start;
				pc_conn_buf->iNodeEndOffset = i_node_end;

				pc_conn_buf->iNodeStartReal = v_nodes.at(pc_conn_buf->iNodeStartOffset).i_node_real_net_id;
				pc_conn_buf->iNodeEndReal = v_nodes.at(pc_conn_buf->iNodeEndOffset).i_node_real_net_id;

				pcSTreeSet->vGetShortestPathBetweenNodePair
					(
					v_nodes.at(pc_conn_buf->iNodeStartOffset).i_node_real_net_id, 
					v_nodes.at(pc_conn_buf->iNodeEndOffset).i_node_real_net_id, 
					pc_conn_buf,
					bUseRealNodeCosts
					);
				v_conns.push_back(pc_conn_buf);

				v_nodes.at(i_node_start).v_conns.push_back(pc_conn_buf);
				//v_nodes.at(i_node_end).v_conns.push_back(pc_conn_buf);
			}//if  (i_node_start != i_node_end)
		}//for  (int  i_node_end = i_node_start + 1; i_node_end < v_full_network.size() - 1; i_node_end++)
	}//for  (int  ii = 0; ii < v_full_network.size(); ii++)

};//void  CSSJNetwork::vCreateCompeleteGraph(vector<int>  *pvDemandedNodes)



bool  CSSJNetwork::bFindMinSpaningTreePrim(CSSJNetwork  *pcMinSpanTree)
{
	CSSJNode  *pc_node_buf;
	CSSJConn  *pc_conn_to_add, *pc_candidate_conn, *pc_new_conn;
	int  i_node_to_check;

	pcMinSpanTree->v_nodes.push_back(v_nodes.at(0));
	pcMinSpanTree->v_nodes.at(0).v_conns.clear();
	pc_conn_to_add = NULL;

	while  ( pcMinSpanTree->v_nodes.size() < this->v_nodes.size())
	{
		pc_conn_to_add = NULL;

		for  (int i_node = 0; i_node < pcMinSpanTree->v_nodes.size(); i_node++)
		{
			for  (int  i_conn = 0; i_conn <  v_nodes.at(pcMinSpanTree->v_nodes.at(i_node).i_node_offset).v_conns.size(); i_conn++)
			{
				pc_candidate_conn = NULL;
				i_node_to_check = v_nodes.at(pcMinSpanTree->v_nodes.at(i_node).i_node_offset).v_conns.at(i_conn)->iNodeEndOffset;

				if  (pcMinSpanTree->bDoIContainThisNode(i_node_to_check)  ==  false)  pc_candidate_conn = v_nodes.at(pcMinSpanTree->v_nodes.at(i_node).i_node_offset).v_conns.at(i_conn);

								
				if  (pc_candidate_conn != NULL)
				{

					if  (pc_conn_to_add == NULL)
						pc_conn_to_add = pc_candidate_conn;
					else
					{
						if  (pc_conn_to_add->dCost > pc_candidate_conn->dCost)  pc_conn_to_add = pc_candidate_conn;
					}//else  if  (pc_conn_to_add == NULL)
				}//if  (pc_candidate_conn != NULL)			
			}//for  (int  i_conn = 0; i_conn <  pcMinSpanTree->v_nodes.at(i_node).v_conns.size(); i_conn++)		
		}//for  (int i_node = 0; i_node < pcMinSpanTree->v_nodes.size(); i_node++)


		if  (pc_conn_to_add != NULL)  
		{
			if  (
				(pcMinSpanTree->bDoIContainThisNode(pc_conn_to_add->iNodeEndOffset)  ==  false)  
				)
			{
				pc_new_conn = new CSSJConn();
				*pc_new_conn = *pc_conn_to_add;

				pc_node_buf = pcMinSpanTree->pcGetNode(pc_conn_to_add->iNodeStartOffset);
				pc_node_buf->v_conns.push_back(pc_new_conn);

				pcMinSpanTree->v_nodes.push_back(v_nodes.at(pc_conn_to_add->iNodeEndOffset));
				pcMinSpanTree->v_nodes.at(pcMinSpanTree->v_nodes.size() - 1).v_conns.clear();
				pcMinSpanTree->v_conns.push_back(pc_new_conn);
			}//if  (
			else
			{
				::vShow("Error at min span tree creation!");
			}//else  if  (			

		}//if  (pc_conn_to_add != NULL)  
		else
			return(false);

	
	}//while  (pcMinSpanTree->v_nodes.size() < this->v_nodes.size())

	CSSJConn  *pc_conn_buf;
	//now make offset in min spanning tree straight...
	for  (int  i_conn = 0; i_conn < pcMinSpanTree->v_conns.size(); i_conn++)
	{
		pc_conn_buf = pcMinSpanTree->v_conns.at(i_conn);

		pc_conn_buf->iNodeStartOffset = pcMinSpanTree->iGetNodeCurrentOffset(pc_conn_buf->iNodeStartOffset);
		pc_conn_buf->iNodeEndOffset = pcMinSpanTree->iGetNodeCurrentOffset(pc_conn_buf->iNodeEndOffset);
	}//for  (int  i_conn = 0; i_conn < pcMinSpanTree->v_nodes.at(i_node).v_conns.size(); i_conn++)
	


	for  (int  i_node = 0; i_node < pcMinSpanTree->v_nodes.size(); i_node++)
	{
		pcMinSpanTree->v_nodes.at(i_node).i_node_offset = i_node;		
	}//for  (int  i_node = 0; i_node < pcMinSpanTree->v_nodes.size(); i_node++)

	return(true);
}//void  CSSJNetwork::vFindMinSpaningTreePrim(CSSJNetwork  *pcMinSpanTree)


void  CSSJNetwork::vFillGenotype(int  *piConnBasedGenotype)
{
	for  (int  i_conn = 0; i_conn < v_conns.size(); i_conn++)
	{
		for  (int  i_fiber_link = 0; i_fiber_link < v_conns.at(i_conn)->v_path.size(); i_fiber_link++)
		{
			piConnBasedGenotype[v_conns.at(i_conn)->v_path.at(i_fiber_link)->iGetId()] = 1;
		}//for  (int  i_fiber_link = 0; i_fiber_link < v_conns.at(i_conn)->v_path.size(); i_fiber_link++)
	}//for  (int  ii = 0; ii < v_conns.size(); ii++)
}//void  CSSJNetwork::vFillGenotype(int  *piConnBasedGenotype)


int CSSJNetwork::iGetNodeCurrentOffset(int  iNodeOffset)
{
	for  (int ii = 0; ii < v_nodes.size(); ii++)
	{
		if  (v_nodes.at(ii).i_node_offset == iNodeOffset)  return(ii);
	}//for  (int ii = 0; ii < v_nodes.size(); ii++)

	return(-1);
}//CSSJNode  * CSSJNetwork::pcGetNode(int  iNodeOffset)


CSSJNode* CSSJNetwork::pcGetNode(int  iNodeOffset)
{
	for  (int ii = 0; ii < v_nodes.size(); ii++)
	{
		if  (v_nodes.at(ii).i_node_offset == iNodeOffset)  return(&(v_nodes.at(ii)));
	}//for  (int ii = 0; ii < v_nodes.size(); ii++)

	return(NULL);
}//CSSJNode  * CSSJNetwork::pcGetNode(int  iNodeOffset)



bool  CSSJNetwork::bDoIContainThisNode(int  iNodeOffset)
{
	if  (pcGetNode(iNodeOffset)  !=  NULL)  return(true);

	return(false);
}//bool  CSSJNetwork::bDoIContainThisNode(int  iNode)


double  CSSJNetwork::dGetCost()
{
	double  d_cost;

	d_cost = 0;
	for  (int  ii = 0; ii < v_conns.size(); ii++)
	{
		d_cost	+=  v_conns.at(ii)->dCost;
	}//for  (int  ii = 0; ii < v_conns.size(); ii++)

	return(d_cost);
}//double  CSSJNetwork::dGetCost()




void  CSSJNetwork::vGenerateNetsToAnalyze(vector<CSSJNetwork *> *pvNetworksToAnalyze)
{
	if  (pc_parent_network == NULL)  return;

	CSSJNetwork  *pc_new_net;

	for  (int  i_conn = 0; i_conn < v_conns.size(); i_conn++)
	{
		pc_new_net = pc_parent_network->pcClone();

		pc_new_net->vRemoveConnection
			(
			v_nodes.at(v_conns.at(i_conn)->iNodeStartOffset).i_node_real_net_id,
			v_nodes.at(v_conns.at(i_conn)->iNodeEndOffset).i_node_real_net_id
			);
		pvNetworksToAnalyze->push_back(pc_new_net);
	}//for  (int  i_conn = 0; i_conn < v_conns.size(); i_conn++)

}//void  CSSJNetwork::vGenerateNetsToAnalyze(vector<CSSJNetwork *> *pvNetworksToAnalyze)



void  CSSJNetwork::vRemoveConnection(int iRealNodeStart, int iRealNodeEnd)
{
	CSSJConn  *pc_conn;

	for  (int  i_conn = 0; i_conn < v_conns.size(); i_conn++)
	{
		if  (
			(v_nodes.at(v_conns.at(i_conn)->iNodeStartOffset).i_node_real_net_id == iRealNodeStart)&&
			(v_nodes.at(v_conns.at(i_conn)->iNodeEndOffset).i_node_real_net_id == iRealNodeEnd)
			)
		{
			pc_conn = v_conns.at(i_conn);
			v_conns.erase(v_conns.begin() + i_conn);

			for  (int i_node = 0; i_node < v_nodes.size(); i_node++)
			{
				for  (int  i_node_conn = 0; i_node_conn < v_nodes.at(i_node).v_conns.size(); i_node_conn++)
				{
					if  (v_nodes.at(i_node).v_conns.at(i_node_conn) == pc_conn)
					{
						v_nodes.at(i_node).v_conns.erase(v_nodes.at(i_node).v_conns.begin() + i_node_conn);
						delete  pc_conn;					
					}//if  (v_nodes.at(i_node).v_conns.at(i_node_conn) == pcConnToDelete)
				}//for  (int  i_node_conn = 0; i_node_conn < v_nodes.at(i_node).v_conns.size(); i_node_conn++)
			}//for  (int i_node = 0; i_node < v_nodes.size(); i_node++)
		}//	if  (v_conns.at(i_conn) == pc_conn_to_delete)
	}//for  (int  i_conn = 0; i_conn < v_conns.size(); i_conn++)
}//void  CSSJNetwork::vRemoveConnection(CSSJConn  *pcConnToDelete)




CSSJNetwork  *CSSJNetwork::pcClone()
{
	CSSJNetwork  *pc_result;
	pc_result = new CSSJNetwork();

	CSSJConn  *pc_conn;

	
	for  (int  i_node = 0; i_node  <  v_nodes.size(); i_node++)
	{
		pc_result->v_nodes.push_back(v_nodes.at(i_node));
		pc_result->v_nodes.at(i_node).v_conns.clear();

		for  (int  i_conn = 0; i_conn  <  v_nodes.at(i_node).v_conns.size(); i_conn++)
		{
			pc_conn = new CSSJConn();
			*pc_conn = *(v_nodes.at(i_node).v_conns.at(i_conn));

			pc_result->v_nodes.at(i_node).v_conns.push_back(pc_conn);
			pc_result->v_conns.push_back(pc_conn);
		}//for  (int  i_conn = 0; i_conn  <  v_nodes.at(i_node).v_conns.size(); i_conn++)
	}//for  (int  i_node = 0; i_node  <  v_nodes.size(); ii++)

	return(pc_result);
}//CSSJNetwork  *CSSJNetwork::pcClone()




void  CSSJNetwork::vSaveGraph(FILE  *pfDest)
{
	CString  s_line, s_buf;

	//fisrt create headers
	s_line = "* \t";
	for  (int  ii = 0; ii < v_nodes.size(); ii++)
	{
		s_buf.Format("%d(%d) \t", v_nodes.at(ii).i_node_offset, v_nodes.at(ii).i_node_real_net_id);
		s_line += s_buf;
	}//for  (int  ii = 0; ii < v_nodes.size(); ii++)
	fprintf(pfDest, "%s\n", s_line);

	
	double  d_cost;
	for  (int  i_start_node = 0; i_start_node < v_nodes.size(); i_start_node++)
	{
		s_line.Format("%d(%d) \t", v_nodes.at(i_start_node).i_node_offset, v_nodes.at(i_start_node).i_node_real_net_id);
		for  (int  i_end_node = 0; i_end_node < v_nodes.size(); i_end_node++)
		{
			d_cost = v_nodes.at(i_start_node).dGetCost(&v_nodes.at(i_end_node));

			s_buf.Format("%.0lf \t", d_cost);
			s_line += s_buf;		
		}//for  (int  i_end_node = 0; i_end_node < v_nodes.size(); i_end_node++)
		fprintf(pfDest, "%s\n", s_line);
	}//for  (int  i_start_node = 0; i_start_node < v_nodes.size(); i_start_node++)


	fprintf(pfDest, "\nPaths:\n", s_line);
	for  (int ii = 0; ii < v_conns.size();  ii++)
	{
		s_buf.Format
			(
			"%d(%d)->%d(%d): ", 
			v_nodes.at(v_conns.at(ii)->iNodeStartOffset).i_node_offset, v_nodes.at(v_conns.at(ii)->iNodeStartOffset).i_node_real_net_id,
			v_nodes.at(v_conns.at(ii)->iNodeEndOffset).i_node_offset, v_nodes.at(v_conns.at(ii)->iNodeEndOffset).i_node_real_net_id
			);
		v_conns.at(ii)->vSave(pfDest);
	}//for  (int ii = 0; ii < v_conns.size();  ii++)

	
}//void  CSSJNetwork::vSaveGraph(FILE  *pfDest)



void  CSSJConn::vSave(FILE  *pfDest)
{
	CString  s_buf;

	fprintf(pfDest, "[");
	for  (int ii = v_path.size()-1; ii >= 0; ii--)
	{
		s_buf.Format("|%d -(%.0lf)-> %d|", v_path.at(ii)->iGetNodeStart(), v_path.at(ii)->dGetCapacity(), v_path.at(ii)->iGetNodeEnd());
		fprintf(pfDest, "%s", s_buf);
	}//for  (int ii = 0; ii < v_path.size(); ii++)
	fprintf(pfDest, "]");
};//void  CSSJConn::vSaveGraph(FILE  *pfDest)


