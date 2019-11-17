#ifndef _H_TESTS
#define _H_TESTS

#undef _CCODE_START
#undef _CCODE_END
#ifdef __cplusplus
#define _CCODE_START extern "C" {
#define _CCODE_END }
#else
#define _CCODE_START
#define _CCODE_END
#endif

#include "./vector.h"
#include "./matrix.h"
#include "./environment.h"
#include "./mathfuncs.h"
#include "./mpdata.h"

_CCODE_START

int test_vector_funcs();
int test_matrix_funcs();
int test_sounding();
void test_advection();
void test_advection_2d();
void test_droplet_advc();
void test_droplet_coal();

_CCODE_END
#endif