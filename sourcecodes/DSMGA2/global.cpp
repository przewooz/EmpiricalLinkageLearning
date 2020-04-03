/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cfloat>
#include "myrand.h"
#include "statistics.h"
#include "doublelinkedlistarray.h"
#include "zkey.h"
#include "chromosome.h"
#include "sat.h"

#include "../EATester/BinaryCoding.h"
#include "../EATester/Evaluation.h"

int maxMemory = 0;

bool GHC = true;
bool SELECTION = true;
bool CACHE = false;
bool SHOW_BISECTION = true;

char outputFilename[100];
Chromosome::Function Chromosome::function;
int Chromosome::nfe;
int Chromosome::lsnfe;
int Chromosome::hitnfe;
bool Chromosome::hit;
bool Chromosome::twoEdge;
unordered_map<unsigned long, double> Chromosome::cache;
CEvaluation<CBinaryCoding>* Chromosome::injectedEvaluation;
CDSMGA2* Chromosome::pcParent;

ZKey zKey;
MyRand myRand;
BitwiseDistance myBD;
SPINinstance mySpinGlassParams;
NKWAProblem nkwa;
SATinstance mySAT;



CError  eLoadZKeyFile(CString  sDirectory)
{
	CError  c_err(CError::iERROR_PARENT_LOAD_ZKEY_FILE);
	c_err = zKey.eLoadZKey(sDirectory);
	return(c_err);
}//void  vLoadZKeyFile(CString  sDirectory)

void outputErrMsg(const char *errMsg) {
    printf("%s\n", errMsg);
    exit(1);
}

int pow2(int x) {
    return (1 << x);
}

