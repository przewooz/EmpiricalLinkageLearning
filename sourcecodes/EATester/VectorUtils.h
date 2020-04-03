#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <vector>

using namespace std;

namespace VectorUtils
{
	template <typename T>
	void vDeleteElementsAndClear(vector<T*> *pvVector)
	{
		for (size_t i = 0; i < pvVector->size(); i++)
		{
			delete pvVector->at(i);
		}//for (size_t i = 0; i < pvVector->size(); i++)

		pvVector->clear();
	}//void vDeleteElementsAndClear(vector<T*> *pvVector)
}//namespace VectorUtils

#endif//VECTOR_UTILS_H