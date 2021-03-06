/*!
 * \file SU2_DOT.hpp
 * \brief Headers of the main subroutines of the code SU2_DOT.
 *        The subroutines and functions are in the <i>SU2_DOT.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
//#include "../../SU2_CFD/include/solver_structure.hpp"  //JRH 06202017

using namespace std;


/*!
 * \brief Projection of the surface sensitivity using finite differences (FD).
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement class of the problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, ofstream& Gradient_file);

/*!
 * \brief Projection of the surface sensitivity using algorithmic differentiation (AD).
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement class of the problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, ofstream& Gradient_file);

/*!
 * \brief Prints the gradient information to a file.
 * \param[in] Gradient - The gradient data.
 * \param[in] config - Definition of the particular problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file);

void GetBetaFimlGrad(CGeometry *geometry, CConfig* config);
