#ifndef MACROS_H_
#define MACROS_H_

#ifndef NDEBUG
#define PRINT(x) std::cout << (x)
#define PRINTLN(x) std::cout << (x) << std::endl
#define ERROR(x) std::cerr << (x)
#define ERRORLN(x) std::cerr << (x) << std::endl
#else
/*! Prints a message */
/* print */

#define PRINT(x)
#define PRINTLN(x)
#define ERROR(x)
#define ERRORLN(x)

#endif // !NDEBUG
#endif // !MACROS_H_

