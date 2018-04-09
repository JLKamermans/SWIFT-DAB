/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int num_tests = 100;
const double eps = 0.02;

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(double a, double b, const char *s) {
  if (fabs(a - b) / fabs(a + b) > 1e-6 && fabs(a - b) > 1.e-6)
    error("Values are inconsistent: %12.15e %12.15e (%s)!", a, b, s);
}

/* Definitions of the potential and force that match
   exactly the theory document */
double S(double x) { return good_approx_exp(x) / (1. + good_approx_exp(x)); }

double S_prime(double x) {
  return good_approx_exp(x) /
         ((1. + good_approx_exp(x)) * (1. + good_approx_exp(x)));
}

double potential(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double pot;
  if (u > 1.)
    pot = -mass / r;
  else
    pot = -mass *
          (-3. * u * u * u * u * u * u * u + 15. * u * u * u * u * u * u -
           28. * u * u * u * u * u + 21. * u * u * u * u - 7. * u * u + 3.) /
          H;

  return pot * (2. - 2. * S(2. * x));
}

double acceleration(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double acc;
  if (u > 1.)
    acc = -mass / (r * r * r);
  else
    acc = -mass * (21. * u * u * u * u * u - 90. * u * u * u * u +
                   140. * u * u * u - 84. * u * u + 14.) /
          (H * H * H);

  return r * acc * (4. * x * S_prime(2 * x) - 2. * S(2. * x) + 2.);
}

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Initialise a few things to get us going */
  struct engine e;
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.time_base = 1e-10;

  struct space s;
  s.width[0] = 0.2;
  s.width[1] = 0.2;
  s.width[2] = 0.2;
  e.s = &s;

  struct gravity_props props;
  props.a_smooth = 1.25;
  props.epsilon_cur = eps;
  e.gravity_properties = &props;

  struct runner r;
  bzero(&r, sizeof(struct runner));
  r.e = &e;

  const double rlr = props.a_smooth * s.width[0];

  /* Init the cache for gravity interaction */
  gravity_cache_init(&r.ci_gravity_cache, num_tests * 2);

  /* Let's create one cell with a massive particle and a bunch of test particles
   */
  struct cell c;
  bzero(&c, sizeof(struct cell));
  c.width[0] = 1.;
  c.width[1] = 1.;
  c.width[2] = 1.;
  c.loc[0] = 0.;
  c.loc[1] = 0.;
  c.loc[2] = 0.;
  c.gcount = 1 + num_tests;
  c.ti_old_gpart = 8;
  c.ti_gravity_end_min = 8;
  c.ti_gravity_end_max = 8;

  if (posix_memalign((void **)&c.gparts, gpart_align,
                     c.gcount * sizeof(struct gpart)) != 0)
    error("Impossible to allocate memory for the gparts.");
  bzero(c.gparts, c.gcount * sizeof(struct gpart));

  /* Create the massive particle */
  c.gparts[0].x[0] = 0.;
  c.gparts[0].x[1] = 0.5;
  c.gparts[0].x[2] = 0.5;
  c.gparts[0].mass = 1.;
  c.gparts[0].time_bin = 1;
  c.gparts[0].type = swift_type_dark_matter;
  c.gparts[0].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
  c.gparts[0].ti_drift = 8;
#endif

  /* Create the mass-less particles */
  for (int n = 1; n < num_tests + 1; ++n) {

    struct gpart *gp = &c.gparts[n];

    gp->x[0] = n / ((double)num_tests);
    gp->x[1] = 0.5;
    gp->x[2] = 0.5;
    gp->mass = 0.;
    gp->time_bin = 1;
    gp->type = swift_type_dark_matter;
    gp->id_or_neg_offset = n + 1;
#ifdef SWIFT_DEBUG_CHECKS
    gp->ti_drift = 8;
#endif
  }

  /* Now compute the forces */
  runner_doself_grav_pp_truncated(&r, &c);

  /* Verify everything */
  for (int n = 1; n < num_tests + 1; ++n) {
    const struct gpart *gp = &c.gparts[n];

    const double epsilon = gravity_get_softening(gp, &props);

    double pot_true = potential(c.gparts[0].mass, gp->x[0], epsilon, rlr);
    double acc_true = acceleration(c.gparts[0].mass, gp->x[0], epsilon, rlr);

    // message("x=%e f=%e f_true=%e", gp->x[0], gp->a_grav[0], acc_true);

    check_value(gp->potential, pot_true, "potential");
    check_value(gp->a_grav[0], acc_true, "acceleration");
  }

  free(c.gparts);
  return 0;
}