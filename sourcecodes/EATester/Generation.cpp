#include "Generation.h"

#include "BinaryCoding.h"

template <class TGenotype>
uint32_t CGeneration<TGenotype>::iERROR_PARENT_CGeneration = CError::iADD_ERROR_PARENT("CGeneration");

template <class TGenotype>
CGeneration<TGenotype>::CGeneration()
{

}//CGeneration<TGenotype>::CGeneration()

template <class TGenotype>
CGeneration<TGenotype>::~CGeneration()
{
}//CGeneration<TGenotype>::~CGeneration()

template class CGeneration<CBinaryCoding>;