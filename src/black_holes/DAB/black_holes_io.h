/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *		 2022 Jasper Leonora Kamermans(leon.kamermans@student.uva.nl)	
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
#ifndef SWIFT_DAB_BLACK_HOLES_IO_H
#define SWIFT_DAB_BLACK_HOLES_IO_H

#include "adiabatic_index.h"
#include "black_holes_part.h"
#include "black_holes_properties.h"
#include "io_properties.h"
#include "kick.h"

/**
 * @brief Specifies which b-particle fields to read from a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void black_holes_read_particles(struct bpart* bparts,
                                              struct io_props* list,
                                              int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, bparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, bparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                bparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, bparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, bparts, h);
  list[5] = io_make_input_field("SubgridMasses", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_MASS, bparts, subgrid_mass);
  list[6] = io_make_input_field("MassGrowthGyr", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_MASS, bparts, mass_growth_gyr);
  list[7] = io_make_input_field("MaxMass", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_MASS, bparts, max_mass);
}

INLINE static void convert_bpart_pos(const struct engine* e,
                                     const struct bpart* bp, double* ret) {

  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(bp->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(bp->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(bp->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = bp->x[0];
    ret[1] = bp->x[1];
    ret[2] = bp->x[2];
  }
}

INLINE static void convert_bpart_vel(const struct engine* e,
                                     const struct bpart* bp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const float dt_kick_grav_mesh = e->dt_kick_grav_mesh_for_io;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, bp->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, bp->time_bin);

  /* Get time-step since the last kick */
  const float dt_kick_grav =
      kick_get_grav_kick_dt(ti_beg, ti_current, time_base, with_cosmology,
                            cosmo) -
      kick_get_grav_kick_dt(ti_beg, (ti_beg + ti_end) / 2, time_base,
                            with_cosmology, cosmo);

  /* Extrapolate the velocites to the current time */
  const struct gpart* gp = bp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Extrapolate the velocites to the current time (mesh forces) */
  ret[0] += gp->a_grav_mesh[0] * dt_kick_grav_mesh;
  ret[1] += gp->a_grav_mesh[1] * dt_kick_grav_mesh;
  ret[2] += gp->a_grav_mesh[2] * dt_kick_grav_mesh;

  /* Conversion from internal to physical units */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}



/**
 * @brief Specifies which b-particle fields to write to a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void black_holes_write_particles(const struct bpart* bparts,
                                               struct io_props* list,
                                               int* num_fields,
                                               int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 15;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_bpart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, bparts,
      convert_bpart_pos, "Co-moving position of the particles");

  list[1] = io_make_output_field_convert_bpart(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts, convert_bpart_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");

  list[2] =
      io_make_output_field("DynamicalMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, mass, "Dynamical masses of the particles");

  list[3] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, id, "Unique ID of the particles");

  list[4] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, bparts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[5] = io_make_output_field("SubgridMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 bparts, subgrid_mass,
                                 "Subgrid masses of the particles");

  if (with_cosmology) {
    list[6] = io_make_output_field(
        "FormationScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
        formation_scale_factor, "Scale-factors at which the BHs were formed");
  } else {
    list[6] = io_make_output_field("FormationTimes", FLOAT, 1, UNIT_CONV_TIME,
                                   0.f, bparts, formation_time,
                                   "Times at which the BHs were formed");
  }

  list[7] = io_make_output_field(
      "CumulativeNumberOfSeeds", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      cumulative_number_seeds,
      "Total number of BH seeds that have merged into this black hole");

  list[7] =
      io_make_output_field("NumberOfMergers", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, number_of_mergers,
                           "Number of mergers the black holes went through. "
                           "This does not include the number of mergers "
                           "accumulated by any merged black hole.");

  if (with_cosmology) {
    list[8] = io_make_output_field(
        "LastMinorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_minor_merger_scale_factor,
        "Scale-factors at which the black holes last had a minor merger.");
  } else {
    list[8] = io_make_output_field(
        "LastMinorMergerTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_minor_merger_time,
        "Times at which the black holes last had a minor merger.");
  }

  if (with_cosmology) {
    list[9] = io_make_output_field(
        "LastMajorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_major_merger_scale_factor,
        "Scale-factors at which the black holes last had a major merger.");
  } else {
    list[9] = io_make_output_field(
        "LastMajorMergerTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_major_merger_time,
        "Times at which the black holes last had a major merger.");
  }

  list[10] =
      io_make_output_field("TimeBins", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
                           time_bin, "Time-bins of the particles");

  list[11] = io_make_output_field(
      "NumberOfRepositions", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_repositions,
      "Number of repositioning events the black holes went through. This does "
      "not include the number of reposition events accumulated by any merged "
      "black holes.");

  list[12] = io_make_output_field(
      "NumberOfRepositionAttempts", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_reposition_attempts,
      "Number of time steps in which the black holes had an eligible particle "
      "to reposition to. They may or may not have ended up moving there, "
      "depending on their subgrid mass and on whether these particles were at "
      "a lower or higher potential than the black holes themselves. It does "
      "not include attempted repositioning events accumulated by any merged "
      "black holes.");

  list[13] = io_make_output_field(
      "NumberOfTimeSteps", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_time_steps,
      "Total number of time steps at which the black holes were active.");

  list[14] = io_make_output_field(
      "LastRepositionVelocities", FLOAT, 1, UNIT_CONV_SPEED, 0.f, bparts,
      last_repos_vel,
      "Physical speeds at which the black holes repositioned most recently. "
      "This is 0 for black holes that have never repositioned, or if the "
      "simulation has been run without prescribed repositioning speed.");

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES

  list += *num_fields;
  *num_fields += 4;

  list[0] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_density);
  list[1] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_force);
  list[2] = io_make_output_field("Ids_ngb_density", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_density);
  list[3] = io_make_output_field("Ids_ngb_force", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_force);
#endif
}

#endif /* SWIFT_DAB_BLACK_HOLES_IO_H */
