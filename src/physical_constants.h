/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PHYSICAL_CONSTANTS_H
#define SWIFT_PHYSICAL_CONSTANTS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "units.h"

/* physical constants in in defined programme units */
struct phys_const {

  /* Newton's gravitationl constant */
  double const_newton_G;

  /* Speed of light in vacuum */
  double const_speed_light_c;

  /* Planck's constant */
  double const_planck_h;

  /* Planck's reduced constant */
  double const_planck_hbar;

  /* Boltzmann's constant */
  double const_boltzmann_k;

  /* Thomson cross-section */
  double const_thomson_cross_section;

  /* Charge of the electron  */
  double const_electron_charge;

  /* Electron-Volt */
  double const_electron_volt;

  /* Mass of the electron */
  double const_electron_mass;

  /* Mass of the proton */
  double const_proton_mass;

  /* (Tropical) Year */
  double const_year;

  /* Astronomical unit */
  double const_astronomical_unit;

  /* Parsec */
  double const_parsec;

  /* Light-year */
  double const_light_year;

  /* Mass of the Sun */
  double const_solar_mass;

  /* Mass of the Earth */
  double const_earth_mass;
};

void phys_const_init(struct UnitSystem* us, struct phys_const* internal_const);

void phys_const_print(struct phys_const* internal_const);

#endif /* SWIFT_PHYSICAL_CONSTANTS_H */