#ifndef LOG_H
#define LOG_H

#include "Error.h"

#include <atlstr.h>
#include <cstdint>
#include <ostream>
#include <vector>

using namespace std;

class CLog
{
public:
	CLog() { i_log_system_to_show = 0; }
	~CLog();

	void vClear() { v_messages.clear(); };
	
	void vPrintLine(CString sMessage, bool bEcho = false, int iMessOffset = 0);
	void vPrintEmptyLine(bool bEcho = false) { vPrintLine("", bEcho); };

	CError eSave(ostream *psOutput, int  iLogOffset = 0);
	CError eSave(CString sFilePath, int  iLogOffset = 0);

	bool bIsClear() { return v_messages.size() == 0; };

	void  vLogSystemToShow(int iLogSystemToShow) { i_log_system_to_show = iLogSystemToShow; }

private:
	CError e_save_single_mess_system(ostream *psOutput, int  iLogOffset);

	static uint32_t iERROR_PARENT_CLOG;

	static uint32_t iERROR_CODE_LOG_CANNOT_WRITE_TO_STREAM;
	
	int i_log_system_to_show;
	CString s_log_file_path;
	vector<vector<CString>> v_messages;
};//class CLog

#endif//LOG_H