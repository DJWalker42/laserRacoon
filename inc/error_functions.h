#ifndef LASER_RACOON_ERROR_FUNCTIONS_H
#define LASER_RACOON_ERROR_FUNCTIONS_H

/* Error diagnostic routines - all similar to printf(format_string, ...) */

namespace phys {
/*
 * @brief: prints the provided message on stdout and returns to caller (alt. to printf)
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 */
void infoMsg(const char *format, ...);


/*
 * @brief: prints the warning message on stderr and returns to caller
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 */
void warnMsg(const char *format, ...);

/*
 * @brief: prints an error message including errno diagnostic on stderr and returns
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 */
void errMsg(const char *format, ...);

#ifdef __GNUC__

    /* This macro stops 'gcc -Wall' complaining that "control reaches
       end of non-void function" if we use the following functions to
       terminate main() or some other non-void function. */

#define NORETURN __attribute__ ((__noreturn__))
#else
#define NORETURN
#endif

/*
 * @brief: same as errMsg but also terminates the program using exit() or abort() (see notes)
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 * @notes: if environment variable EF_DUMPCORE is defined with a non-empty string
 * terminates by calling abort() that produces a core dump file.
 */
void errExit(const char *format, ...) NORETURN ;


/*
 * @brief: same as errExit but does not flush stdout before printing message and calls _exit() instead of exit()
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 * @notes: The relationship between this function and errExit() is analogous to that between _exit(2) and exit(3):
 * unlike errExit(), this function does not flush stdout and calls _exit(2) to terminate the process (rather than
 * exit(3), which would cause exit handlers to be invoked). If EF_DUMPCORE is defined calls abort() regardless.
 */
void err_exit(const char *format, ...) NORETURN ;

/*
 * @brief: similar to errExit but prints information on the errnum supplied rather than using errno (see notes)
 * @param: "errnum" the error number for which we want information
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 * @notes: use with programs that use POSIX threads (in threaded programs errno is macro that expands to a function call)
 */
void errExitEN(int errnum, const char *format, ...) NORETURN ;

/*
 * @brief: diagnose general errors then terminates (c.f. errExit())
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 * @notes: use if errors produced from library functions that do not set errno
 */
void fatal(const char *format, ...) NORETURN ;

/*
 * @brief: general function to diagnose errors in command-line arguments
 * @param: "format" character string containing the output format of the message
 * @param: varargs contains all the required variables to print
 * @notes: prints "Usage: " followed by the formatted output
 */
void usageErr(const char *format, ...) NORETURN ;

} //namespace

#endif //header guard
