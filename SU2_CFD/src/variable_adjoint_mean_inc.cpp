/*!
 * \file variable_adjoint_mean.cpp
 * \brief Definition of the solution fields.
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

#include "../include/variable_structure.hpp"

CAdjIncEulerVariable::CAdjIncEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  Psi = NULL;
  ForceProj_Vector = NULL;
  ObjFuncSource = NULL;
  IntBoundary_Jump = NULL;
  // mskim
  AxiAuxVar = NULL;
  Grad_AxiAuxVar = NULL;

}

CAdjIncEulerVariable::CAdjIncEulerVariable(su2double val_psirho, su2double *val_phi, su2double val_psie, unsigned short val_nDim,
                                     unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
  Psi = NULL;
  ForceProj_Vector = NULL;
  ObjFuncSource = NULL;
  IntBoundary_Jump = NULL;
  // mskim
  AxiAuxVar = NULL;
  Grad_AxiAuxVar = NULL;


  /*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Only for residual smoothing (multigrid) ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
  
  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED)
    Undivided_Laplacian = new su2double [nVar];
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_UPWIND) {
    Limiter = new su2double [nVar];
    Solution_Max = new su2double [nVar];
    Solution_Min = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Limiter[iVar] = 0.0;
      Solution_Max[iVar] = 0.0;
      Solution_Min[iVar] = 0.0;
    }
  }
  
  /*--- Allocate and initialize solution ---*/
  Solution[0] = 0.0;   Solution_Old[0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[iDim+1] = 0.0;
    Solution_Old[iDim+1] = 0.0;
  }

  
  /*--- Allocate and initialize solution for dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0] = 0.0;
    Solution_time_n1[0] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_time_n[iDim+1] = 0.0;
      Solution_time_n1[iDim+1] = 0.0;
    }
  }
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
  Grad_AuxVar = new su2double [nDim];
// mskim
  nAuxVar = 1;
  Grad_AxiAuxVar = new su2double* [nAuxVar];
  for (iVar = 0; iVar < nAuxVar; iVar++) AxiAuxVar[iVar] = 0.0;
  
  Grad_AxiAuxVar = new su2double* [nAuxVar];
  for (iVar = 0; iVar < nAuxVar; iVar++) {
    Grad_AxiAuxVar[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Grad_AxiAuxVar[iVar][iDim] = 0.0;
  }
// mskim-end
  
  /*--- Allocate and initialize projection vector for wall boundary condition ---*/
  ForceProj_Vector = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    ForceProj_Vector[iDim] = 0.0;
  
  /*--- Allocate and initialize interior boundary jump vector for near field boundary condition ---*/
  IntBoundary_Jump = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    IntBoundary_Jump[iVar] = 0.0;
  
  
}

CAdjIncEulerVariable::CAdjIncEulerVariable(su2double *val_solution, unsigned short val_nDim,
                                     unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
  Psi = NULL;
  ForceProj_Vector = NULL;
  ObjFuncSource = NULL;
  IntBoundary_Jump = NULL;
  
  /*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Only for residual smoothing (multigrid) ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
  
  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED)
    Undivided_Laplacian = new su2double [nVar];
  
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_UPWIND) {
    Limiter = new su2double [nVar];
    Solution_Max = new su2double [nVar];
    Solution_Min = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Limiter[iVar] = 0.0;
      Solution_Max[iVar] = 0.0;
      Solution_Min[iVar] = 0.0;
    }
  }
  
  /*--- Solution initialization ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }
  
  /*--- Allocate and initializate solution for dual time strategy ---*/
  if (dual_time) {
    Solution_time_n = new su2double [nVar];
    Solution_time_n1 = new su2double [nVar];
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
  Grad_AuxVar = new su2double [nDim];
// mskim  
  nAuxVar = 1;
  Grad_AxiAuxVar = new su2double* [nAuxVar];
  for (iVar = 0; iVar < nAuxVar; iVar++) AxiAuxVar[iVar] = 0.0;
  
  Grad_AxiAuxVar = new su2double* [nAuxVar];
  for (iVar = 0; iVar < nAuxVar; iVar++) {
    Grad_AxiAuxVar[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Grad_AxiAuxVar[iVar][iDim] = 0.0;
  }
// mskim-end

  /*--- Allocate and initializate projection vector for wall boundary condition ---*/
  ForceProj_Vector = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    ForceProj_Vector[iDim] = 0.0;
  
  /*--- Allocate and initializate interior boundary jump vector for near field boundary condition ---*/
  IntBoundary_Jump = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    IntBoundary_Jump[iVar] = 0.0;
  
  
}

CAdjIncEulerVariable::~CAdjIncEulerVariable(void) {
    unsigned short iVar;

  if (Psi               != NULL) delete [] Psi;
  if (ForceProj_Vector  != NULL) delete [] ForceProj_Vector;
  if (ObjFuncSource     != NULL) delete [] ObjFuncSource;
  if (IntBoundary_Jump  != NULL) delete [] IntBoundary_Jump;
  // mskim
  if (AxiAuxVar       != NULL) delete [] AxiAuxVar;
  nAuxVar = 1;
  if (Grad_AxiAuxVar  != NULL) {
    for (iVar = 0; iVar < nAuxVar; iVar++)
      if (Grad_AxiAuxVar != NULL) delete [] Grad_AxiAuxVar[iVar];
        delete [] Grad_AxiAuxVar;
  }

}

bool CAdjIncEulerVariable::SetPrimVar(su2double SharpEdge_Distance, bool check, CConfig *config) {
  unsigned short iVar;
  bool check_dens = false, RightVol = true;
  
  su2double adj_limit = config->GetAdjointLimit();
  
  check_dens = (fabs(Solution[0]) > adj_limit);
  
  /*--- Check that the adjoint solution is bounded ---*/
  
  if (check_dens) {
    
    /*--- Copy the old solution ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    RightVol = false;
    
  }
  
  return RightVol;
  
}


CAdjIncNSVariable::CAdjIncNSVariable(void) : CAdjIncEulerVariable() { }

CAdjIncNSVariable::CAdjIncNSVariable(su2double *val_solution, unsigned short val_nDim,
                               unsigned short val_nvar, CConfig *config) : CAdjIncEulerVariable(val_solution, val_nDim, val_nvar, config) {
  
}

CAdjIncNSVariable::CAdjIncNSVariable(su2double val_psirho, su2double *val_phi, su2double val_psie,
                               unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CAdjIncEulerVariable(val_psirho, val_phi, val_psie, val_nDim, val_nvar, config) {

}

CAdjIncNSVariable::~CAdjIncNSVariable(void) { }
