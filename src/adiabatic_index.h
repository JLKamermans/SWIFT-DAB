/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_ADIABATIC_INDEX_H
#define SWIFT_ADIABATIC_INDEX_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "const.h"
#include "debug.h"
#include "inline.h"

/* First define some constants */
#if defined(HYDRO_GAMMA_5_3)

#define hydro_gamma 1.66666666666666667f
#define hydro_gamma_minus_one 0.66666666666666667f
#define hydro_one_over_gamma_minus_one 1.5f

#elif defined(HYDRO_GAMMA_4_3)

#define hydro_gamma 1.33333333333333333f
#define hydro_gamma_minus_one 0.33333333333333333f
#define hydro_one_over_gamma_minus_one 3.f

#elif defined(HYDRO_GAMMA_2_1)

#define hydro_gamma 2.f
#define hydro_gamma_minus_one 1.f
#define hydro_one_over_gamma_minus_one 1.f

#endif

/**
 * @brief Returns the argument to the power given by the adiabatic index
 *
 * Computes \f$x^\gamma\f$.
 */
__attribute__((always_inline)) INLINE static float pow_gamma(float x) {

#if defined(HYDRO_GAMMA_5_3)

  const float cbrt = cbrtf(x); /* x^(1/3) */
  return cbrt * cbrt * x;      /* x^(5/3) */

#elif defined(HYDRO_GAMMA_4_3)

  return cbrtf(x) * x; /* x^(4/3) */

#elif defined(HYDRO_GAMMA_2_1)

  return x * x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns the argument to the power given by the adiabatic index minus
 * one
 *
 * Computes \f$x^{(\gamma-1)}\f$.
 */
__attribute__((always_inline)) INLINE static float pow_gamma_minus_one(
    float x) {

#if defined(HYDRO_GAMMA_5_3)

  const float cbrt = cbrtf(x); /* x^(1/3) */
  return cbrt * cbrt;          /* x^(2/3) */

#elif defined(HYDRO_GAMMA_4_3)

  return cbrtf(x); /* x^(1/3) */

#elif defined(HYDRO_GAMMA_2_1)

  return x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns one over the argument to the power given by the adiabatic
 * index minus one
 *
 * Computes \f$x^{-(\gamma-1)}\f$.
 */
__attribute__((always_inline)) INLINE static float pow_minus_gamma_minus_one(
    float x) {

#if defined(HYDRO_GAMMA_5_3)

  const float inv_cbrt = 1.f / cbrtf(x); /* x^(-1/3) */
  return inv_cbrt * inv_cbrt;            /* x^(-2/3) */

#elif defined(HYDRO_GAMMA_4_3)

  return 1.f / cbrtf(x); /* x^(-1/3) */

#elif defined(HYDRO_GAMMA_2_1)

  return 1.f / x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

#endif /* SWIFT_ADIABATIC_INDEX_H */