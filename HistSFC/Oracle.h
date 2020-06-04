#pragma once
#include "occi.h"
#include "BaseStruct.h"
#include <iomanip>

using namespace oracle::occi;
using namespace std;

/* --------------------------------------------------------------------------------------------- *
* CONSTANTS
* --------------------------------------------------------------------------------------------- */

inline ConnParam orclconn(){
	ConnParam conn = {"orcl", "c##haicheng", "123456"};
	return conn;
}
