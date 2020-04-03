#include "Error.h"
#include "StringUtils.h"
#include "System.h"

#include <Windows.h>

#include "../DSMGA2/global.h"

#ifdef DEBUG

//#include <vld.h>

#endif

void main(int iArgCount, char **ppcArgValues)
{
	try
	{
		CError c_error;
		
		CString  s_path(*(ppcArgValues));
		CString  s_dir;
		s_dir = StringUtils::sGetDirFromPath(s_path);
		eLoadZKeyFile(s_dir);

		
		CSystem c_system(iArgCount, ppcArgValues, &c_error);

		if (!c_error)
		{
			c_error = c_system.eRun();
		}//if (!c_error)
	}//try
	catch (exception  *p_ex)
	{
		CError c_error;
	
		CString  s_buf;
		s_buf.Format("General error: %s", p_ex->what());
		c_error.vSetError(s_buf);
	}//catch (exception  *p_ex)

}//void main()