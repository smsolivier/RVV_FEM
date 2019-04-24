#include "Error.hpp"
// #include "Array.hpp"

#ifdef USE_UNWIND 
#define UNW_NAME_LEN 512
#define UNW_LOCAL_ONLY
#include <libunwind.h>
#include <cxxabi.h>
#if defined(__APPLE__) || defined(__linux__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <dlfcn.h>
#endif
#endif

#include <iostream> 
#include "Array.hpp"

using namespace fem; 

namespace error
{

void backtrace(int mode, int depth) {
#if defined USE_UNWIND && !defined NDEBUG
	std::cout << "\nBacktrace: " << std::endl; 
	char name[UNW_NAME_LEN];
	unw_cursor_t cursor;
	unw_context_t uc;
	unw_word_t ip, offp;

	int err = unw_getcontext(&uc);
	err = err ? err : unw_init_local(&cursor, &uc);

	Array<unw_word_t> addrs;
	while (unw_step(&cursor) > 0 && addrs.GetSize() != depth) {
		err = err ? err : unw_get_proc_name(&cursor, name, UNW_NAME_LEN, &offp);
		err = err ? err : unw_get_reg(&cursor, UNW_REG_IP, &ip);
		if (err) { break; }
		char *name_p = name;
		int demangle_status;

		// __cxa_demangle is not standard, but works with GCC, Intel, PGI, Clang
		char *name_demangle =
			abi::__cxa_demangle(name, NULL, NULL, &demangle_status);
		if (demangle_status == 0) // use mangled name if something goes wrong
		{
			name_p = name_demangle;
		}

		std::cout << "\t" << addrs.GetSize() << ") [0x" << std::hex << ip - 1 << std::dec
			<< "]: " << name_p << std::endl;
		addrs.Append(ip - 1);

		if (demangle_status == 0) {
			free(name_demangle);
		}
	}
#if defined(__APPLE__) || defined(__linux__)
   if (addrs.GetSize() > 0 && (mode & 1))
   {
      std::cout << "\nLookup backtrace source lines:";
      const char *fname = NULL;
      for (int i = 0; i < addrs.GetSize(); i++)
      {
         Dl_info info;
         err = !dladdr((void*)addrs[i], &info);
         if (err)
         {
            fname = "<exe>";
         }
         else if (fname != info.dli_fname)
         {
            fname = info.dli_fname;
            std::cout << '\n';
#ifdef __linux__
            std::cout << "\taddr2line -C -e " << fname;
#else
            std::cout << "\tatos -o " << fname << " -l "
                      << (err ? 0 : info.dli_fbase);
#endif
         }
         std::cout << " 0x" << std::hex << addrs[i] << std::dec;
      }
      std::cout << '\n';
   }
#endif
   std::cout << "\n"; 
#endif // MFEM_USE_LIBUNWIND
}

} // end namespace error 