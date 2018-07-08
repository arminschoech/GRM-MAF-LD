#ifndef RANDOM_H_
#define RANDOM_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern gsl_rng* ran;
// now all files including "random.h" can use object "ran"
// object "ran" will be declared and defined in "main.cpp"
















#endif