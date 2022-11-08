/**
 * @brief Collection of miscellaneous C++ utils
 */
#include "spu_utils.hpp"

#include <cmath>

/**
 * @brief Returns area under ROC curve for given distributions
 *
 * @param value       floating point number to round
 * @param n_sigfigs   number of significant figures
 *
 * @returns rounded floating point number
 */
float round_sigfigs(float value, int n_sigfigs) {
  //credit to https://stackoverflow.com/questions/13094224/a-c-routine-to-round-a-float-to-n-significant-digits
  if (value == 0.0) return 0.0;
  float factor = pow(10.0, n_sigfigs - ceil(log10(fabs(value))));
  return roundf(value * factor) / factor;
}

