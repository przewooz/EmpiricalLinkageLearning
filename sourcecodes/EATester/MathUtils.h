#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cfloat>
#include <cstdint>
#include <functional>

using namespace std;

namespace MathUtils
{
	double dComputeAngle(double *pdValues0, double *pdValues1, uint32_t iLength);
	double dComputeDotProduct(double *pdValues0, double *pdValues1, uint32_t iLength);
	double dComputeSecondNorm(double *pdValues, uint32_t iLength);
	double dComputeSquareDistance(double *pdValues0, double *pdValues1, uint32_t iLength, double dMaxValue = DBL_MAX);
	
	double dComputeEntropy(uint32_t *piCounts, uint32_t iLength);
	double dComputeEntropy(uint32_t *piCounts, uint32_t iLength, uint32_t iTotal);

	double dMaxValue(double *pdValues, uint32_t iLength);
	double dSum(double *pdValues, uint32_t iLength);

	void vNormalize(double *pcValues, uint32_t iLength);
}//namespace MathUtils

#endif//MATH_UTILS_H