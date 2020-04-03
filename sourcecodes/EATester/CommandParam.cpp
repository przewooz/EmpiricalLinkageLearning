#include "CommandParam.h"

#include <cstdlib>
#include <string>

CCommandParam::CCommandParam(CString sName, bool bIsObligatory, bool bHasDefaultValue)
{
	s_name = sName;
	b_is_obligatory = bIsObligatory;
	b_has_default_value = bHasDefaultValue;
	b_has_value_been_set = false;
}//CCommandParam::CCommandParam(CString sName, bool bIsObligatory, bool bHasDefaultValue)

CCommandParam::~CCommandParam()
{

}//CCommandParam::~CCommandParam()

CError CCommandParam::eSetValue(istream *psSettings)
{
	CError c_error;
	string s_line;
	
	streamoff i_settings_position_before_read = psSettings->tellg();

	if (getline(*psSettings, s_line))
	{
		CString s_cstring_line(s_line.c_str());

		CString s_param_value; 
		CString s_param_name; 

		int i_line_separator_position = s_cstring_line.Find("\\\\");

		if (i_line_separator_position != -1)
		{
			s_param_value = s_cstring_line.Mid(0, i_line_separator_position);
			s_param_name = s_cstring_line.Mid(i_line_separator_position + 2);
		}//if (i_line_separator_position != -1)

		if (s_name != s_param_name)
		{
			if (!bIsObligatory())
			{
				if (bHasDefaultValue())
				{
					vSetDefaultValue();
				}//if (bHasDefaultValue())

				psSettings->clear();
				psSettings->seekg(i_settings_position_before_read);
			}//if (!bIsObligatory())
			else
			{
				CString s_error_message;
				s_error_message.AppendFormat("line: %s; expected: %s; found: %s", s_cstring_line, s_name, s_param_name);

				c_error.vSetError(CError::iERROR_CODE_SYSTEM_MISMATCHING_ARGUMENT, s_error_message);
			}//else if (!bIsObligatory())
		}//if (s_name != s_param_name)
		else
		{
			c_error = eSetValue(s_param_value);
		}//else if (s_name != s_param_name)
	}//if (getline(*psSettings, s_line))
	else
	{
		if (bIsObligatory())
		{
			c_error.vSetError(CError::iERROR_CODE_SYSTEM_SETTINGS_CANNOT_READ_LINE, s_name);
		}//if (bIsObligatory())
	}//else if (getline(*psSettings, s_line))

	return c_error;
}//CError CCommandParam::eSetValue(istream *psSettings)


CBoolCommandParam::CBoolCommandParam(CString sName, bool bIsObligatory)
	: CCommandParam(sName, bIsObligatory, false)
{
	v_init(b_default_value);
}//CBoolCommandParam::CBoolCommandParam(CString sName, bool bIsObligatory)

CBoolCommandParam::CBoolCommandParam(CString sName, bool bDefaultValue, bool bIsObligatory)
	: CCommandParam(sName, bIsObligatory, true)
{
	v_init(bDefaultValue);
}//CBoolCommandParam::CBoolCommandParam(CString sName, bool bDefaultValue, bool bIsObligatory)

bool CBoolCommandParam::bGetValue(istream *psSettings, CError *pcError)
{
	*pcError = eSetValue(psSettings);
	return b_value;
}//bool CBoolCommandParam::bGetValue(istream *psSettings, CError *pcError)

CError CBoolCommandParam::eSetValue(CString sValue)
{
	CError c_error;

	char *c_parsing_end_point;

	long l_value = strtol(sValue, &c_parsing_end_point, 10);

	if (*c_parsing_end_point == '\0' && c_parsing_end_point != sValue)
	{
		if (l_value == 0 || l_value == 1)
		{
			b_value = (bool)l_value;
			v_set_value();
		}//if (l_value == 0 || l_value == 1)
		else
		{
			CString s_error_message;
			s_error_message.Format("%s must be 0 or 1", sGetName());
			
			c_error.vSetError(CError::iERROR_CODE_SYSTEM_OUT_OF_BOUND_ARGUMENT, s_error_message);
		}//else if (l_value == 0 || l_value == 1)
	}//if (*c_parsing_end_point == '\0' && c_parsing_end_point != sValue)
	else
	{
		c_error.vSetError(CError::iERROR_CODE_SYSTEM_CANNOT_PARSE_ARGUMENT, sGetName());
	}//else if (*c_parsing_end_point == '\0' && c_parsing_end_point != sValue)

	return c_error;
}//CError CBoolCommandParam::eSetValue(CString sValue)