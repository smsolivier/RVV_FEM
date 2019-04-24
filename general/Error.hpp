#pragma once 

#include <string> 

namespace error  
{

/// runs libunwind 
void backtrace(int mode=0, int depth=-1); 

} // end namespace error 
