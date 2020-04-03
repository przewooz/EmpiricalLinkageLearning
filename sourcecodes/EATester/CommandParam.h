#ifndef COMMAND_PARAM_H
#define COMMAND_PARAM_H

#include "Error.h"

#include <atlstr.h>
#include <istream>

using namespace std;

class CCommandParam
{
public:
	CCommandParam(CString sName, bool bIsObligatory, bool bHasDefaultValue);
	
	virtual ~CCommandParam();

	CString sGetName() { return s_name; };
	bool bIsObligatory() { return b_is_obligatory; };
	bool bHasDefaultValue() { return b_has_default_value; };
	bool bHasValue() { return b_is_obligatory || b_has_default_value || b_has_value_been_set; };

	CError eSetValue(istream *psSettings);
	virtual CError eSetValue(CString sValue) = 0;
	virtual void vSetDefaultValue() = 0;

protected:
	void v_set_value() { b_has_value_been_set = true; };

private:
	CString s_name;
	bool b_is_obligatory;
	bool b_has_default_value;
	bool b_has_value_been_set;
};//class CCommandParam


class CBoolCommandParam : public CCommandParam
{
public:
	CBoolCommandParam(CString sName, bool bIsObligatory = true);
	CBoolCommandParam(CString sName, bool bDefaultValue, bool bIsObligatory);

	bool bGetValue() { return b_value; };
	bool bGetValue(istream *psSettings, CError *pcError);

	using CCommandParam::eSetValue;
	virtual CError eSetValue(CString sValue);
	virtual void vSetDefaultValue() { b_value = b_default_value; };

private:
	void v_init(bool bDefaultValue) { b_default_value = bDefaultValue; };

	bool b_value;
	bool b_default_value;
};//class CBoolCommandParam : public CCommandParam

#endif//COMMAND_PARAM_H