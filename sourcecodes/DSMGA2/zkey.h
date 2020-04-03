
#ifndef __ZKEY__
#define __ZKEY__

#include <cstdio>
#include <cstdlib>

#include <Windows.h>
#include <atlstr.h>

#include "../EATester/Error.h"


using namespace std;

class ZKey {

public:
	static uint32_t iERROR_PARENT_CLOG;
	static uint32_t iERROR_CODE_ZKEY_FILE_NOT_FOUND;


    ZKey() {}

	CError  eLoadZKey(CString  sDirectory);

    unsigned long operator[](int i) const {return keys[i]; }

protected:

    unsigned long keys[1000];

};

#endif
