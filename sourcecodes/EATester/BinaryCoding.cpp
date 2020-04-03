#include "BinaryCoding.h"

using namespace std;

CBinaryCoding::CBinaryCoding(uint16_t iNumberOfBits)
{
	i_number_of_bits = iNumberOfBits;
	pi_bits = new int32_t[iNumberOfBits];
	b_own_bits = true;
}//CBinaryCoding::CBinaryCoding(uint16_t iNumberOfBits, int32_t *piBits)

CBinaryCoding::CBinaryCoding(uint16_t iNumberOfBits, int32_t *piBits, bool bOwnBits)
{
	i_number_of_bits = iNumberOfBits;
	pi_bits = piBits;
	b_own_bits = bOwnBits;
}//CBinaryCoding::CBinaryCoding(uint16_t iNumberOfBits, int32_t *piBits, bool bOwnBits)

CBinaryCoding::CBinaryCoding(CBinaryCoding *pcOther)
{
	i_number_of_bits = pcOther->i_number_of_bits;

	pi_bits = new int32_t[pcOther->i_number_of_bits];
	
	for (uint16_t i = 0; i < pcOther->i_number_of_bits; i++)
	{
		*(pi_bits + i) = *(pcOther->pi_bits + i);
	}//for (uint16_t i = 0; i < pcOther->i_number_of_bits; i++)

	b_own_bits = true;
}//CBinaryCoding::CBinaryCoding(CBinaryCoding *pcOther)

CBinaryCoding::~CBinaryCoding()
{
	if (b_own_bits)
	{
		delete pi_bits;
	}//if (b_own_bits)
}//CBinaryCoding::~CBinaryCoding()

void CBinaryCoding::vMerge(CBinaryCoding *pcOther)
{
	//TODO
}//void CBinaryCoding::vMerge(CBinaryCoding *pcOther)

double CBinaryCoding::dGetUnitation()
{
	uint16_t i_number_of_ones = 0;

	for (uint16_t i = 0; i < i_number_of_bits; i++)
	{
		if (*(pi_bits + i) == 1)
		{
			i_number_of_ones++;
		}//if (*(pi_bits + i) == 1)
	}//for (uint16_t i = 0; i < i_number_of_bits; i++)

	return (double)i_number_of_ones / (double)i_number_of_bits;
}//double CBinaryCoding::dGetUnitation()

void CBinaryCoding::vSetBits(int32_t *piBits, bool bOwnBits)
{
	if (b_own_bits)
	{
		delete pi_bits;
	}//if (b_own_bits)

	pi_bits = piBits;
	b_own_bits = bOwnBits;
}//void CBinaryCoding::vSetBits(int32_t *piBits, bool bOwnBits)

CString CBinaryCoding::sToString(uint16_t iSpaceFrequency)
{
	CString s_result;

	for (uint16_t i = 0; i < i_number_of_bits; i++)
	{
		if (i > 0 && i % iSpaceFrequency == 0)
		{
			s_result.Append(" ");
		}//if (i > 0 && i % iSpaceFrequency == 0)

		s_result.AppendFormat("%d", *(pi_bits + i));
	}//for (uint16_t i = 0; i < i_number_of_bits; i++)

	return s_result;
}//CString CBinaryCoding::sToString(uint16_t iSpaceFrequency)

bool CBinaryCoding::operator==(CBinaryCoding &cOther)
{
	bool b_equal = true;

	for (uint16_t i = 0; i < i_number_of_bits && b_equal; i++)
	{
		b_equal &= *(pi_bits + i) == *(cOther.pi_bits + i);
	}//for (uint16_t i = 0; i < i_number_of_bits && b_equal; i++)

	return b_equal;
}//bool CBinaryCoding::operator==(CBinaryCoding &cOther)

bool CBinaryCoding::operator!=(CBinaryCoding &cOther)
{
	return !(*this == cOther);
}//bool CBinaryCoding::operator!=(CBinaryCoding &cOther)

ostream& operator<<(ostream &sOutput, CBinaryCoding *pcBinaryCoding)
{
	sOutput << pcBinaryCoding->sToString();

	return sOutput;
}//ostream& operator<<(ostream &sOutput, CBinaryCoding &cBinaryCoding)