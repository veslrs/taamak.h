#ifndef INCLUDE_TAAMAK_H
#define INCLUDE_TAAMAK_H

#ifndef TMKDEF
#ifdef TAAMAK_STATIC
#define TMKDEF static
#else
#ifdef __cplusplus
#define TMKDEF extern "C"
#else
#define TMKDEF extern
#endif
#endif
#endif

TMKDEF int add(const int a, const int b);

#endif // INCLUDE_TAAMAK_H

#ifdef TAAMAK_IMPLEMENTATION

#ifdef _WIN32
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#endif
#endif

TMKDEF int add(const int a, const int b) {
        int c;
        c = a + b;

        return c;
}

#endif
