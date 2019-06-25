/*
 *  Print message to stderr and exit
 */
#include "HW3.h"

void Fatal(const char* format , ...)
{
   va_list args;
   va_start(args,format);
   vfprintf(stderr,format,args);
   va_end(args);
   exit(1);
}
