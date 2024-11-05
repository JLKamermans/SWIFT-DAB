/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DAB_BLACK_HOLES_PROPERTIES_H
#define SWIFT_DAB_BLACK_HOLES_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "exp10.h"
#include "hydro_properties.h"

/* Includes. */
#include <string.h>

/**
 * @brief Modes of energy injection for AGN feedback
 */
enum AGN_feedback_models {
  AGN_random_ngb_model,       /*< Random neighbour model for AGN feedback */
  AGN_isotropic_model,        /*< Isotropic model of AGN feedback */
  AGN_minimum_distance_model, /*< Minimum-distance model of AGN feedback */
  AGN_minimum_density_model   /*< Minimum-density model of AGN feedback */
};

/**
 * @brief Properties of black holes in the DAB model.
 */
struct black_holes_props {

  /* ----- Basic neighbour search properties ------ */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* ----- Initialisation properties  ------ */

  /*! Mass of a BH seed at creation time */
  float subgrid_seed_mass;

  /*! Should we use the subgrid mass specified in ICs? */
  int use_subgrid_mass_from_ics;

  /*! Should we enforce positive subgrid masses initially? */
  int with_subgrid_mass_check;

  /* ---- Properties of the repositioning model --- */

  /*! Maximal mass of BH to reposition */
  float max_reposition_mass;

  /*! Maximal distance to reposition, in units of softening length */
  float max_reposition_distance_ratio;

  /*! Switch to enable a relative velocity limit for particles to which the
   * black holes can reposition */
  int with_reposition_velocity_threshold;

  /*! Maximal velocity offset of particles to which the black hole can
   * reposition, in units of the ambient sound speed of the black hole */
  float max_reposition_velocity_ratio;

  /*! Minimum value of the velocity repositioning threshold */
  float min_reposition_velocity_threshold;

  /*! Switch to enable repositioning at fixed (maximum) speed */
  int set_reposition_speed;

  /*! Normalisation factor for repositioning velocity */
  float reposition_coefficient_upsilon;

  /*! Reference black hole mass for repositioning scaling */
  float reposition_reference_mass;

  /*! Repositioning velocity scaling with black hole mass */
  float reposition_exponent_mass;

  /*! Reference gas density for repositioning scaling */
  float reposition_reference_n_H;

  /*! Repositioning velocity scaling with gas density */
  float reposition_exponent_n_H;

  /* ---- Properties of the merger model ---------- */

  /*! Mass ratio above which a merger is considered 'minor' */
  float minor_merger_threshold;

  /*! Mass ratio above which a merger is considered 'major' */
  float major_merger_threshold;

  /*! Type of merger threshold (0: standard, 1: improved) */
  int merger_threshold_type;

  /*! Maximal distance over which BHs merge, in units of softening length */
  float max_merging_distance_ratio;

  /* ---- Black hole time-step properties ---------- */

  /*! Minimum allowed time-step of BH (internal units) */
  float time_step_min;

  /* ---- Common conversion factors --------------- */

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /*! Conversion factor from physical density to n_H [cgs] */
  float rho_to_n_cgs;

  /*! Conversion factor from internal mass to solar masses */
  float mass_to_solar_mass;

  /*! Allow the movement of a black hole, or not*/
  float allow_movement_black_hole;
};

/**
 * @brief Initialise the black hole properties from the parameter file.
 *
 * For the basic black holes neighbour finding properties we use the
 * defaults from the hydro scheme if the users did not provide specific
 * values.
 *
 * @param bp The #black_holes_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void black_holes_props_init(struct black_holes_props *bp,
                                          const struct phys_const *phys_const,
                                          const struct unit_system *us,
                                          struct swift_params *params,
                                          const struct hydro_props *hydro_props,
                                          const struct cosmology *cosmo) {

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);


  /* Kernel properties */
  bp->eta_neighbours = parser_get_opt_param_float(
      params, "BlackHoles:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  bp->h_tolerance = parser_get_opt_param_float(params, "BlackHoles:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  bp->target_neighbours = pow_dimension(bp->eta_neighbours) * kernel_norm;
  const float delta_eta = bp->eta_neighbours * (1.f + bp->h_tolerance);
  bp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(bp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  bp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BlackHoles:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Initialisation properties  ---------------------------- */

  bp->subgrid_seed_mass =
      parser_get_param_float(params, "DAB:subgrid_seed_mass_Msun");

  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics =
      parser_get_opt_param_int(params, "DAB:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics)
    bp->with_subgrid_mass_check =
        parser_get_opt_param_int(params, "DAB:with_subgrid_mass_check", 1);

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "DAB:max_reposition_mass") *
      phys_const->const_solar_mass;
  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "DAB:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "DAB:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "DAB:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "DAB:min_reposition_velocity_threshold");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "DAB:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "DAB:reposition_coefficient_upsilon");

    /* Prevent the user from making silly wishes */
    if (bp->reposition_coefficient_upsilon <= 0)
      error(
          "reposition_coefficient_upsilon must be positive, not %f "
          "km/s/M_sun.",
          bp->reposition_coefficient_upsilon);

    /* Convert from km/s to internal units */
    bp->reposition_coefficient_upsilon *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    /* Scaling parameters with BH mass and gas density */
    bp->reposition_reference_mass =
        parser_get_param_float(params, "DAB:reposition_reference_mass") *
        phys_const->const_solar_mass;
    bp->reposition_exponent_mass = parser_get_opt_param_float(
        params, "DAB:reposition_exponent_mass", 2.0);
    bp->reposition_reference_n_H =
        parser_get_param_float(params, "DAB:reposition_reference_n_H");
    bp->reposition_exponent_n_H = parser_get_opt_param_float(
        params, "DAB:reposition_exponent_n_H", 1.0);
  }

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "DAB:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "DAB:threshold_major_merger");

  bp->merger_threshold_type =
      parser_get_param_int(params, "DAB:merger_threshold_type");

  bp->max_merging_distance_ratio =
      parser_get_param_float(params, "DAB:merger_max_distance_ratio");

  /* ---- Black hole time-step properties ------------------ */

  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 60. * 60.;

  const double time_step_min_Myr = parser_get_opt_param_float(
      params, "DAB:minimum_timestep_Myr", FLT_MAX);

  bp->time_step_min = time_step_min_Myr * Myr_in_cgs /
                      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Common conversion factors ----------------------------- */

  /* Calculate conversion factor from rho to n_H.
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  bp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Conversion factor for internal mass to M_solar */
  bp->mass_to_solar_mass = 1. / phys_const->const_solar_mass;

  bp->allow_movement_black_hole =
      parser_get_param_float(params, "DAB:allow_movement_black_hole");
}

/**
 * @brief Write a black_holes_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_dump(
    const struct black_holes_props *props, FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct black_holes_props), 1,
                       stream, "black_holes props", "black holes props");
}

/**
 * @brief Restore a black_holes_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_restore(
    const struct black_holes_props *props, FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct black_holes_props), 1,
                      stream, NULL, "black holes props");
}

#endif /* SWIFT_DAB_BLACK_HOLES_PROPERTIES_H */
