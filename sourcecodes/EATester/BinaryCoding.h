#ifndef BINARY_CODING_H
#define BINARY_CODING_H

#include <atlstr.h>
#include <cstdint>
#include <ostream>

using namespace std;

class CBinaryCoding
{
public:
	CBinaryCoding(uint16_t iNumberOfBits);
	CBinaryCoding(uint16_t iNumberOfBits, int32_t *piBits, bool bOwnBits = false);
	CBinaryCoding(CBinaryCoding *pcOther);
	
	~CBinaryCoding();

	void vMerge(CBinaryCoding *pcOther);

	void vRepair() {};

	bool bHasSameParams(CBinaryCoding *pcOther) { return i_number_of_bits == pcOther->i_number_of_bits; };

	uint16_t iGetNumberOfBits() { return i_number_of_bits; };
	int32_t *piGetBits() { return pi_bits; };

	double dGetUnitation();

	void vSetNumberOfBits(uint16_t iNumberOfBits) { i_number_of_bits = iNumberOfBits; };
	void vSetBits(int32_t *piBits, bool bOwnBits = false);

	CString sToString() { return sToString(i_number_of_bits); };
	CString sToString(uint16_t iSpaceFrequency);

	bool operator==(CBinaryCoding &cOther);
	bool operator!=(CBinaryCoding &cOther);

	friend ostream& operator<<(ostream &sOutput, CBinaryCoding *pcBinaryCoding);

private:
	uint16_t i_number_of_bits;
	int32_t *pi_bits;

	bool b_own_bits;
};//class CBinaryCoding

#endif//BINARY_CODING_H