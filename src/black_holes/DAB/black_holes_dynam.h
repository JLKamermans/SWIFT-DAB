/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) Jasper Leonora Kamermans (leon.kamermans@student.uva.nl)
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
#ifndef SWIFT_DAB_BLACK_HOLES_DYNAM_H
#define SWIFT_DAB_BLACK_HOLES_DYNAM_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "cosmology.h"
#include "dimension.h"
#include "gravity.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"
#include "rays.h"

/* Standard includes */
#include <float.h>
#include <math.h>

// Ma et al: https://arxiv.org/abs/2208.12275

#ifdef KEES

/**
 * @brief This function will compute the dynamical friction upon the BH by surrounding DM particles
 *
 * 
 * TODO: Find out if I need to multiply by G here, or if the code internally works with G=1
 * 
 *
 * @param bp The particle to act upon (pointer)
 * @param gp The collection of DM particles (array of pointers)
 * @param const_G Newtons G constant in internal units
 * 
 */
__attribute__((always_inline)) INLINE static void black_holes_dynam_fric(
    struct bpart* bp) {

    //todo: Find some sort of radius beyond which we won't calculate the DF
    // Preferably, we should have already selected which particles we will calculate the DF before passing them to this function.

    // For every DM particle in the array of pointers, execute this task:
    for (int i = 0; i < sizeof(gp)/sizeof(gp[0]); i++){
        
        // x_vec is the distance between the BH and the particle, x_sq its length squared
        float x_vec[3] = {gp[i]->x[0] - bp->x[0],gp[i]->x[1] - bp->x[1],gp[i]->x[2] - bp->x[2]};
        float x_sq = x_vec[0]*x_vec[0] + x_vec[1]*x_vec[1] + x_vec[2]*x_vec[2];

        // if the particle is too far away, continue with the next particle.
        if(x_sq > 10.) continue;


        // v_vec is the velocity in the restframe of the BH, v_sq its length squared
        float v_vec[3] = {gp[i]->v_full[0] - bp->v[0],gp[i]->v_full[1] - bp->v[1],gp[i]->v_full[2] - bp->v[2]};
        float v_sq = v_vec[0]*v_vec[0] + v_vec[1]*v_vec[1] + v_vec[2]*v_vec[2];

        // Normalised v_vec
        float v_vec_norm[3] = {v_vec[0]/sqrt(v_sq),v_vec[1]/sqrt(v_sq),v_vec[2]/sqrt(v_sq)};


        //define some stuff (ask Camila) (plummer softening as in equation (7) of Ma et al.)
        float impact_parameter = 1.f;
        float plummer_softening = 1.f;

        //determine the value of the kernel (which accounts for the force softening):
        float kernel = 1.f;
        float q = sqrt(x_sq)/(2.8*plummer_softening);
        if (q < 0.5){
            kernel = 32./5.*q*q*q - 192./5.*q*q*q*q*q + 32.*q*q*q*q*q*q;
        }
        else if(q < 1){
            kernel = -1./15. + 64./3.*q*q*q -48.*q*q*q*q + 192./5. *q*q*q*q*q
                     + 192./5.*q*q*q*q*q - 32./3.*q*q*q*q*q*q;
        }

        // alpha as defined in equation (9) of Ma et al.
        float alpha = impact_parameter*v_sq/(const_G*bp->mass);


        // The dot product between x_vec and normalised v_vec
        float r_dot_v = x_vec[0] * v_vec_norm[0] + x_vec[1] * v_vec_norm[1] + x_vec[2] * v_vec_norm[2];

        // The b vector is x_vec - (x_vec * norm(v_vec))v_vec, and the direction perpendicular to v, b_sq its length squared
        float b_vec[3] = {x_vec[0] - r_dot_v*v_vec_norm[0],
                          x_vec[1] - r_dot_v*v_vec_norm[1],
                          x_vec[2] - r_dot_v*v_vec_norm[2]};
        float b_sq = b_vec[0]*b_vec[0] + b_vec[1]*b_vec[1] + b_vec[2]*b_vec[2];

        // This term just occurs everywhere, so better to compute it just once instead of 6 times
        float temp_term = impact_parameter/(1 + alpha*alpha)/sqrt(x_sq)*const_G*gp[i]->mass/x_sq*kernel;

        // This is the perpendicular acceleration due to dynamical friction (equation (8) of Ma et al.)
        float a_perpendicular[3] = {temp_term*b_vec[0]/sqrt(b_sq),
                                    temp_term*b_vec[1]/sqrt(b_sq),
                                    temp_term*b_vec[2]/sqrt(b_sq)};

        // main acceleration (equation (9) of Ma et al.)
        float a_dynam[3] = {alpha*temp_term*v_vec_norm[0],
                            alpha*temp_term*v_vec_norm[1],
                            alpha*temp_term*v_vec_norm[2]};


        // update the dynamical friction acceleration of the BH
        bp->gpart->a_dynam[0] += a_perpendicular[0] + a_dynam[0];
        bp->gpart->a_dynam[1] += a_perpendicular[1] + a_dynam[1];
        bp->gpart->a_dynam[2] += a_perpendicular[2] + a_dynam[2];

        // add the backreaction on the gparticle, so momentum is conserved.
        gp[i]->a_dynam[0] = -(bp->mass/gp[i]->mass)*(a_perpendicular[0] + a_dynam[0]);
        gp[i]->a_dynam[1] = -(bp->mass/gp[i]->mass)*(a_perpendicular[1] + a_dynam[1]);
        gp[i]->a_dynam[2] = -(bp->mass/gp[i]->mass)*(a_perpendicular[2] + a_dynam[2]);
    }
}
#endif /* KEES */
#endif /* SWIFT_DAB_BLACK_HOLES_DYNAM_H */
