#include "global_options.h"

string trim(const string s)
{
	string trimmed = s;
	string::size_type pos = trimmed.find_last_not_of(' ');
	if(pos != string::npos)
	{
		if (trimmed.length()!=pos+1)//if there are trailing whitespaces erase them
			trimmed.erase(pos+1);
		pos = trimmed.find_first_not_of(' ');
		if(pos!=0) //if there are leading whitespaces erase them
			trimmed.erase(0, pos);
	}
	else trimmed="";
	return trimmed;
}

