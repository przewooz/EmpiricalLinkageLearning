#include "zkey.h"



uint32_t ZKey::iERROR_PARENT_CLOG = CError::iADD_ERROR_PARENT("ZKey");
uint32_t ZKey::iERROR_CODE_ZKEY_FILE_NOT_FOUND = CError::iADD_ERROR("iERROR_CODE_ZKEY_FILE_NOT_FOUND");

CError  ZKey::eLoadZKey(CString  sDirectory)
{
	CError  c_err(iERROR_PARENT_CLOG);

	//MessageBox(NULL, sDirectory, sDirectory, MB_OK);

	FILE *fp = fopen("zobristkey", "rb");
	if (fp == NULL)
	{
		c_err.vSetError(iERROR_CODE_ZKEY_FILE_NOT_FOUND);
		return(c_err);	
	}//if (fp == NULL)

	fread(keys, sizeof(unsigned long), 1000, fp);
	fclose(fp);

	return(c_err);
}//CError  ZKey::eLoadZKey(CString  sDirectory)