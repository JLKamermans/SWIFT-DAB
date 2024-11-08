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
#ifndef SWIFT_DAB_BLACK_HOLE_PART_H
#define SWIFT_DAB_BLACK_HOLE_PART_H

/*! The total number of rays used in AGN feedback */
#define dab_blackhole_number_of_rays FEEDBACK_NR_RAYS_AGN

#include "black_holes_struct.h"
#include "particle_splitting_struct.h"
#include "rays_struct.h"
#include "timeline.h"

/**
 * @brief Particle fields for the black hole particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct bpart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v[3];
  
  /*! Black hole mass */
  float mass;

  /*! The mass at which the BH should stop growing */
  float max_mass;

  /*! Mass growth per giga year */
  float mass_growth_gyr;

  /*! Time of last growth*/
  double time_last_growth;

  /*! Black hole mass at the start of each step, prior to any nibbling */
  float mass_at_start_of_step;

  /* Particle cutoff radius. */
  float h;

  /*! Particle time bin */
  timebin_t time_bin;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Union for the formation time and formation scale factor */
  union {

    /*! Formation time */
    float formation_time;

    /*! Formation scale factor */
    float formation_scale_factor;
  };

  /*! Subgrid mass of the black hole */
  float subgrid_mass;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Number of seeds in this BH (i.e. itself + the merged ones) */
  int cumulative_number_seeds;

  /*! Total number of BH merger events (i.e. not including all progenies) */
  int number_of_mergers;

  /*! Total number of times the black hole has been repositioned (excluding
   * repositionings of merged-in black holes) */
  int number_of_repositions;

  /*! Total number of times a black hole attempted repositioning (including
   * cases where it was aborted because the black hole was already at a
   * lower potential than all eligible neighbours) */
  int number_of_reposition_attempts;

  /* Velocity of most recent reposition jump */
  float last_repos_vel;

  /*! Total number of time steps in which the black hole was active. */
  int number_of_time_steps;

  /*! Total (physical) angular momentum accumulated by swallowing particles */
  float swallowed_angular_momentum[3];

  /*! Union for the last minor merger point in time */
  union {

    /*! Last time the BH had a a high Eddington fraction */
    float last_minor_merger_time;

    /*! Last scale factor the BH had a a high Eddington fraction */
    float last_minor_merger_scale_factor;
  };

  /*! Union for the last major merger point in time */
  union {

    /*! Last time the BH had a a high Eddington fraction */
    float last_major_merger_time;

    /*! Last scale factor the BH had a a high Eddington fraction */
    float last_major_merger_scale_factor;
  };

  struct {

    /*! Gravitational potential copied from the #gpart. */
    float potential;

    /*! Value of the minimum potential across all neighbours. */
    float min_potential;

    /*! Delta position to apply after the reposition procedure */
    double delta_x[3];

  } reposition;

  /*! Splitting structure */
  struct particle_splitting_data split_data;

  /*! Black holes merger information (e.g. merging ID) */
  struct black_holes_bpart_data merger_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DAB_BLACK_HOLE_PART_H */
