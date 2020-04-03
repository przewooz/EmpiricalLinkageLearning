#include "MathUtils.h"

#include "RandUtils.h"

#include <cmath>

double MathUtils::dComputeAngle(double *pdValues0, double *pdValues1, uint32_t iLength)
{
	double d_cosinus_value = dComputeDotProduct(pdValues0, pdValues1, iLength) / dComputeSecondNorm(pdValues0, iLength) / dComputeSecondNorm(pdValues1, iLength);
	
	if (d_cosinus_value > 1 && round(d_cosinus_value) == 1)
	{
		d_cosinus_value = 1;
	}//if (d_cosinus_value > 1 && round(d_cosinus_value) == 1)
	
	return acos(d_cosinus_value);
}//double MathUtils::dComputeAngle(double *pdValues0, double *pdValues1, uint32_t iLength)

double MathUtils::dComputeDotProduct(double *pdValues0, double *pdValues1, uint32_t iLength)
{
	double d_dot_product = 0;

	for (uint32_t i = 0; i < iLength; i++)
	{
		d_dot_product += *(pdValues0 + i) * *(pdValues1 + i);
	}//for (uint32_t i = 0; i < iLength; i++)

	return d_dot_product;
}//double MathUtils::dComputeDotProduct(double *pdValues0, double *pdValues1, uint32_t iLength)

double MathUtils::dComputeSecondNorm(double *pdValues, uint32_t iLength)
{
	return sqrt(dComputeDotProduct(pdValues, pdValues, iLength));
}//double MathUtils::dComputeSecondNorm(double *pdValues, uint32_t iLength)

double MathUtils::dComputeSquareDistance(double* pdValues0, double* pdValues1, uint32_t iLength, double dMaxValue)
{
	double d_square_distance = 0;

	for (uint32_t i = 0; i < iLength && d_square_distance < dMaxValue; i++)
	{
		d_square_distance += (*(pdValues0 + i) - *(pdValues1 + i)) * (*(pdValues0 + i) - *(pdValues1 + i));
	}//for (uint32_t i = 0; i < iLength && d_square_distance < dMaxValue; i++)

	if (d_square_distance > dMaxValue)
	{
		d_square_distance = dMaxValue;
	}//if (d_square_distance > dMaxValue)

	return d_square_distance;
}//double MathUtils::dComputeSquareDistance(double* pdValues0, double* pdValues1, uint32_t iLength, double dMaxValue)

double MathUtils::dComputeEntropy(uint32_t *piCounts, uint32_t iLength)
{
	uint32_t i_total = 0;

	for (uint32_t i = 0; i < iLength; i++)
	{
		i_total += *(piCounts + i);
	}//for (uint32_t i = 0; i < iLength; i++)

	return dComputeEntropy(piCounts, iLength, i_total);
}//double MathUtils::dComputeEntropy(uint32_t *piCounts, uint32_t iLength)

double MathUtils::dComputeEntropy(uint32_t *piCounts, uint32_t iLength, uint32_t iTotal)
{
	double d_entropy = 0;

	if (iTotal > 0)
	{
		double d_total = (double)iTotal;

		double d_probability;

		for (uint32_t i = 0; i < iLength; i++)
		{
			if (*(piCounts + i) > 0)
			{
				d_probability = *(piCounts + i) / d_total;
				d_entropy -= d_probability * log(d_probability);
			}//if (*(piCounts + i) > 0)
		}//for (uint32_t i = 0; i < iLength; i++)
	}//if (iTotal > 0)

	return d_entropy;
}//double MathUtils::dComputeEntropy(uint32_t *piCounts, uint32_t iLength, uint32_t iTotal)

double MathUtils::dMaxValue(double *pdValues, uint32_t iLength)
{
	double d_max_value = -DBL_MAX;

	for (uint32_t i = 0; i < iLength; i++)
	{
		if (*(pdValues + i) > d_max_value)
		{
			d_max_value = *(pdValues + i);
		}//if (*(pdValues + i) > d_max_value)
	}//for (uint32_t i = 0; i < iLength; i++)

	return d_max_value;
}//double MathUtils::dMaxValue(double *pdValues, uint32_t iLength)

double MathUtils::dSum(double *pdValues, uint32_t iLength)
{
	double d_sum = 0;

	for (uint32_t i = 0; i < iLength; i++)
	{
		d_sum += *(pdValues + i);
	}//for (uint32_t i = 0; i < iLength; i++)

	return d_sum;
}//double MathUtils::dSum(double *pdValues, uint32_t iLength)

void MathUtils::vNormalize(double *pdValues, uint32_t iLength)
{
	double d_second_norm = dComputeSecondNorm(pdValues, iLength);

	if (d_second_norm > 0)
	{
		for (uint32_t i = 0; i < iLength; i++)
		{
			*(pdValues + i) /= d_second_norm;
		}//for (uint32_t i = 0; i < iLength; i++)
	}//if (d_sum > 0)
}//void MathUtils::vNormalize(double *pdValues, uint32_t iLength)