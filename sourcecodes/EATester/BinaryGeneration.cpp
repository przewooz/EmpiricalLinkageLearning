#include "BinaryGeneration.h"

#include "CommandParam.h"
#include "FloatCommandParam.h"
#include "RandUtils.h"

#include <vector>

using namespace std;

CBinaryGeneration::CBinaryGeneration(uint16_t iNumberOfBits)
{
	i_number_of_bits = iNumberOfBits;
}//CBinaryGeneration::CBinaryGeneration(uint16_t iNumberOfBits)

CBinaryGeneration::CBinaryGeneration(CBinaryCoding *pcSample)
	: CBinaryGeneration(pcSample->iGetNumberOfBits())
{

}//CBinaryGeneration::CBinaryGeneration(CBinaryCoding *pcSample)