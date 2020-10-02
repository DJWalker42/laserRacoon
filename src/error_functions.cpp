//Headers
//System headers
#include <sys/types.h>  /* Type definitions used by many programs */
#include <sys/stat.h>
#include <stdarg.h>		/* va_* macros and types */
#include <stdio.h>      /* Standard I/O functions */
#include <stdlib.h>     /* EXIT_SUCCESS and EXIT_FAILURE constants */
#include <unistd.h>     /* Prototypes for many system calls (e.g. _exit) */
#include <errno.h>      /* Declares errno and defines error constants */
#include <string.h>     /* Commonly used string-handling functions */
#include <time.h>		/* time functions, e.g. localtime*/

//project headers
#include "error_functions.h"  	/* Declares our error-handling functions */
#include "ename.c.inc"          /* Defines ename and MAX_ENAME */
            /* - to create: ./Build_ename.sh > ./src/ename.c.inc */


//anonymous namespace == static function
namespace {

void printPreamble(FILE *stream, const char *label) {
    fprintf(stream, "%s", label);
    char s[20];
    time_t t = time(0); //time now
    struct tm *lt = localtime(&t);
    strftime(s, 20, "%F %T", lt);// format 'YYYY-MM-DD HH:MM:SS'
    fprintf(stream, "%s", s);
    fprintf(stream, " ");
}

}

namespace phys {

#ifdef __GNUC__                 /* Prevent 'gcc -Wall' complaining  */
__attribute__ ((__noreturn__))  /* if we call this function as last */
#endif                          /* statement in a non-void function */
static void
terminate(bool useExit3)
{
    char *s;

    /* Dump core if EF_DUMPCORE environment variable is defined and
       is a nonempty string; otherwise call exit(3) or _exit(2),
       depending on the value of 'useExit3'.
   */

    s = getenv("EF_DUMPCORE");

    if (s != NULL && *s != '\0')
        abort();
    else if (useExit3)
        exit(EXIT_FAILURE);
    else
        _exit(EXIT_FAILURE);
}

/* Diagnose 'errno' error by:

      * outputting a string containing the error name (if available
        in 'ename' array) corresponding to the value in 'err', along
        with the corresponding error message from strerror(), and

      * outputting the caller-supplied error message specified in
        'format' and 'ap'.
*/

static void
outputError(bool useErr, int err, bool flushStdout,
        const char *format, va_list ap)
{
#define BUF_SIZE 500
    char buf[BUF_SIZE], userMsg[BUF_SIZE], errText[BUF_SIZE];

    vsnprintf(userMsg, BUF_SIZE, format, ap);

    if (useErr)
        snprintf(errText, BUF_SIZE, " [%s %s]",
                (err > 0 && err <= MAX_ENAME) ?
                ename[err] : "?UNKNOWN? ", strerror(err));
    else
        snprintf(errText, BUF_SIZE, ":");

    snprintf(buf, BUF_SIZE, "ERROR%s %s\n", errText, userMsg);

    if (flushStdout)
        fflush(stdout);       /* Flush any pending stdout */
    fputs(buf, stderr);
    fflush(stderr);           /* In case stderr is not line-buffered */
}


/* Display message on stdout and return to caller -
 * adds a timestamp to the front of the message, format to be no more than 999 characters */
void
infoMsg(const char *format, ...)
{
    va_list argList;

    printPreamble(stdout, "I ");
    va_start(argList, format);
    vfprintf(stdout, format, argList);
    va_end(argList);

    if (format[strlen(format) - 1] != '\n') {
        fprintf(stderr, "\n");
    }

    fflush(stdout);           /* In case stdout is not line-buffered */
}

/* Display message on stderr and return to caller - doesn't add newline at end of message  */
void
warnMsg(const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "W ");

    int savedErrno;

    savedErrno = errno;       /* In case we change it here */

    fflush(stdout);           /* Flush any pending stdout */

    va_start(argList, format);
    vfprintf(stderr, format, argList);
    va_end(argList);

    if (format[strlen(format) - 1] != '\n') {
        fprintf(stderr, "\n");
    }

    fflush(stderr);           /* In case stderr is not line-buffered */

    errno = savedErrno;
}


/* Display error message excluding errno diagnostic and return to caller
*/

void
errMsg(const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "E ");

    int savedErrno;

    savedErrno = errno;       /* In case we change it here */

    va_start(argList, format);
    outputError(false, errno, true, format, argList);
    va_end(argList);

    errno = savedErrno;
}

/* Display error message including 'errno' diagnostic, and
   terminate the process
*/

void
errExit(const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "E ");

    va_start(argList, format);
    outputError(true, errno, true, format, argList);
    va_end(argList);

    terminate(true);
}

/* Display error message including 'errno' diagnostic, and
   terminate the process by calling _exit().

   The relationship between this function and errExit() is analogous
   to that between _exit(2) and exit(3): unlike errExit(), this
   function does not flush stdout and calls _exit(2) to terminate the
   process (rather than exit(3), which would cause exit handlers to be
   invoked).

   These differences make this function especially useful in a library
   function that creates a child process that must then terminate
   because of an error: the child must terminate without flushing
   stdio buffers that were partially filled by the caller and without
   invoking exit handlers that were established by the caller.
*/

void
err_exit(const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "E ");

    va_start(argList, format);
    outputError(true, errno, false, format, argList);
    va_end(argList);

    terminate(false);
}

/* The following function does the same as errExit(), but expects
   the error number in 'errnum'
*/

void
errExitEN(int errnum, const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "E ");

    va_start(argList, format);
    outputError(true, errnum, true, format, argList);
    va_end(argList);

    terminate(true);
}

/* Print an error message (without an 'errno' diagnostic) */

void
fatal(const char *format, ...)
{
    va_list argList;

    printPreamble(stderr, "F ");

    va_start(argList, format);
    outputError(false, 0, true, format, argList);
    va_end(argList);

    terminate(true);
}

/* Print a command usage error message and terminate the process */

void
usageErr(const char *format, ...)
{
    va_list argList;

    fflush(stdout);           /* Flush any pending stdout */

    fprintf(stderr, "Usage: ");
    va_start(argList, format);
    vfprintf(stderr, format, argList);
    va_end(argList);

    fflush(stderr);           /* In case stderr is not line-buffered */
    exit(EXIT_FAILURE);
}

} //namespace
