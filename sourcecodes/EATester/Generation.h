#ifndef GENERATION_H
#define GENERATION_H

#include "Error.h"

#include <cstdint>
#include <istream>
#include <vector>

using namespace std;

template <class TGenotype>
class CGeneration
{
public:
	static uint32_t iERROR_PARENT_CGeneration;

	CGeneration();

	virtual ~CGeneration();

	virtual CError eConfigure(istream *psSettings) { return CError(iERROR_PARENT_CGeneration); };

	virtual TGenotype *pcGenerate() = 0;
	virtual TGenotype *pcGenerateEmpty() = 0;
};//class CGeneration

#endif//GENERATION_H