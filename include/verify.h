#ifndef _VERIFY_H
#define _VERIFY_H

#include <iostream>
#include "serial.h"

bool isResultCorrect(const matrix_t<double> *const L,
                     const matrix_t<double> *const b,
                     const double *const X);

#endif // _VERIFY_H