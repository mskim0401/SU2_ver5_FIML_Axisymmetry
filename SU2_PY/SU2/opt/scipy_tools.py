#!/usr/bin/env python 

## \file scipy_tools.py
#  \brief tools for interfacing with scipy
#  \author T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy

from .. import eval as su2eval
from numpy import array, zeros
from numpy.linalg import norm


# -------------------------------------------------------------------
#  Scipy SLSQP
# -------------------------------------------------------------------

def scipy_slsqp(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_slsqp(project,x0=[],xb=[],its=100,accu=1e-10)
    
        Runs the Scipy implementation of SLSQP with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH-ACTUALLY USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_slsqp

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []
    
    # function handles
    func           = obj_f
    f_eqcons       = con_ceq
    f_ieqcons      = con_cieq 
    
    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE': 
        fprime         = None
        fprime_eqcons  = None
        fprime_ieqcons = None
    else:
        fprime         = obj_df
        fprime_eqcons  = con_dceq
        fprime_ieqcons = con_dcieq        
    
    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    if project.config.NPOIN == 0 :
            dv_size = project.config['DEFINITION_DV']['SIZE']
            n_dv = sum( dv_size)
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] =x0[k]/dv_scl;
            k = k + 1

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = []
    for this_obj in obj.keys():
        obj_scale = obj_scale + [obj[this_obj]['SCALE']]
    
    # Only scale the accuracy for single-objective problems: 
    if len(obj.keys())==1:
        accu = accu*obj_scale[0]

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Sequential Least SQuares Programming (SLSQP) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(len(dv_size)) + ' ( ' + str(n_dv) + ' ) \n' )
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Run Optimizer
    outputs = fmin_slsqp( x0             = x0             ,
                          func           = func           , 
                          f_eqcons       = f_eqcons       , 
                          f_ieqcons      = f_ieqcons      ,
                          fprime         = fprime         ,
                          fprime_eqcons  = fprime_eqcons  , 
                          fprime_ieqcons = fprime_ieqcons , 
                          args           = (project,)     , 
                          bounds         = xb             ,
                          iter           = its            ,
                          iprint         = 2              ,
                          full_output    = True           ,
                          acc            = accu           ,
                          epsilon        = eps            )
    
    # Done
    return outputs
    
# -------------------------------------------------------------------
#  Scipy CG
# -------------------------------------------------------------------

def scipy_cg(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_cg(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of CG with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_cg

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        dv_size = project.config['DEFINITION_DV']['SIZE']
        n_dv = sum( dv_size)
    else :
        n_dv = int(project.config.NPOIN)
        
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv

    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    if int(project.config.NPOIN) == 0 :
        x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Conjugate gradient (CG) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = fmin_cg( x0             = x0             ,
                       f              = func           ,
                       fprime         = fprime         ,
                       args           = (project,)     ,
                       gtol           = accu           ,
                       epsilon        = eps            ,
                       maxiter        = its            ,
                       full_output    = True           ,
                       disp           = True           ,
                       retall         = True           )


    # Done
    return outputs

# -------------------------------------------------------------------
#  Scipy BFGS
# -------------------------------------------------------------------

def scipy_bfgs(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_bfgs

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    if project.config['NPOIN'] == 0 :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]
    else :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        for i in range(len(x0)) :
            x0[i] = x0[i]/dv_scales[0]
        
        #x0 = [ x0[i]/dv_scales[0] for i in enumerate(x0)]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    if project.config.NPOIN == 0 :
        accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Broyden-Fletcher-Goldfarb-Shanno (BFGS) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = fmin_bfgs( x0             = x0             ,
                         f              = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         gtol           = accu           ,
                         epsilon        = eps            ,
                         maxiter        = its            ,
                         full_output    = True           ,
                         disp           = True           ,
                         retall         = True           )

    # Done
    return outputs

# -------------------------------------------------------------------
#  Jon's BFGS - Modified SciPy BFGS
# -------------------------------------------------------------------

def jon_bfgs(project,x0=None,xb=None,its=999999,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from optim_custom import jon_bfgs

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    if project.config['NPOIN'] == 0 :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]
    else :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        for i in range(len(x0)) :
            x0[i] = x0[i]/dv_scales[0]
        
        #x0 = [ x0[i]/dv_scales[0] for i in enumerate(x0)]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    if project.config.NPOIN == 0 :
        accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Broyden-Fletcher-Goldfarb-Shanno (BFGS) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = jon_bfgs( x0             = x0             ,
                         f              = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         gtol           = accu           ,
                         epsilon        = eps            ,
                         maxiter        = its            ,
                         full_output    = True           ,
                         disp           = True           ,
                         retall         = True           )

    # Done
    return outputs

# -------------------------------------------------------------------
#  Jon's Steepest Descent - Modified SciPy BFGS
# -------------------------------------------------------------------

def jon_steep(project,x0=None,xb=None,its=999999,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from optim_custom import jon_steep

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    if project.config['NPOIN'] == 0 :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]
    else :
        dv_scales = project.config['DEFINITION_DV']['SCALE']
        for i in range(len(x0)) :
            x0[i] = x0[i]/dv_scales[0]
        
        #x0 = [ x0[i]/dv_scales[0] for i in enumerate(x0)]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    l_rate = float(project.config['LEARN_RATE'])
    
    obj_scale = obj[obj.keys()[0]]['SCALE']
    if project.config.NPOIN == 0 :
        accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Steepest Descent (STEEPJ) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Setting learning rate to ' + str(l_rate)+'\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')
    
    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = jon_steep( x0             = x0             ,
                         f              = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         gtol           = accu           ,
                         epsilon        = eps            ,
                         maxiter        = its            ,
                         full_output    = True           ,
                         disp           = True           ,
                         retall         = True           ,
                         learn_rate     = l_rate)

    # Done
    return outputs

# -------------------------------------------------------------------
#  Scipy L-BFGS-G
# -------------------------------------------------------------------

def scipy_bfgsg(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_l_bfgs_b

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    #dv_scales = project.config['DEFINITION_DV']['SCALE']
    #x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    if project.config.NPOIN == 0 :
        accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Broyden-Fletcher-Goldfarb-Shanno (BFGS) parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    #maxls option not available on Deepthought2 for some reason...
    outputs = fmin_l_bfgs_b( x0         = x0             ,
                         func           = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         pgtol          = accu           ,
                         bounds         = xb             ,
                         epsilon        = eps            ,
                         factr          = 10.0           ,
                         maxiter        = its            ,
                         iprint         = 99             ,
                         m              = 50             )

    # Done
    return outputs

# -------------------------------------------------------------------
#  Scipy Steepest Descent (Truncated Newton Conjugate Gradient)
# -------------------------------------------------------------------

def scipy_steep(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_tnc

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    #dv_scales = project.config['DEFINITION_DV']['SCALE']
    #x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Steepest Descent Parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = fmin_tnc( x0              = x0             ,
                         func           = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         pgtol          = accu           ,
                         bounds         = xb             ,
                         epsilon        = eps            ,
                         disp           = 5              ,
                         maxCGit        = 0) #maxCGit = 0 for steepest descent
    

    # Done
    return outputs

# -------------------------------------------------------------------
#  Scipy Steepest Descent (Truncated Newton Conjugate Gradient)
# -------------------------------------------------------------------

def scipy_tnc(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_bfgs(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of BFGS with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_tnc

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f

    # gradient handles
    if project.config.get('GRADIENT_METHOD','NONE') == 'NONE':
        fprime         = None
    else:
        fprime         = obj_df

    # number of design variables
    if project.config.NPOIN == 0 :
        n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    else :
        n_dv = int(project.config.NPOIN)
            
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    

    # prescale x0
    # commented out below for testing... JRH
    #dv_scales = project.config['DEFINITION_DV']['SCALE']
    #x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Truncated Newton Conjugate Gradient (TNC) Parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = fmin_tnc( x0              = x0             ,
                         func           = func           ,
                         fprime         = fprime         ,
                         args           = (project,)     ,
                         pgtol          = accu           ,
                         bounds         = xb             ,
                         epsilon        = eps            ,
                         disp           = 5              ,
                         maxCGit        = -1) #maxCGit = 0 for steepest descent
    

    # Done
    return outputs

def scipy_powell(project,x0=None,xb=None,its=100,accu=1e-10,grads=False):
    """ result = scipy_powell(project,x0=[],xb=[],its=100,accu=1e-10)

        Runs the Scipy implementation of Powell's method with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds JRH - NOT USED!!
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_powell

    # handle input cases
    if x0 is None: x0 = []

    # function handles
    func           = obj_f

    # number of design variables
    n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv

    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Powells method parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(n_dv) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')

    # Evaluate the objective function (only 1st iteration)
    obj_f(x0,project)

    # Run Optimizer
    outputs = fmin_powell( x0             = x0             ,
                           func           = func           ,
                           args           = (project,)     ,
                           ftol           = accu           ,
                           maxiter        = its            ,
                           full_output    = True           ,
                           disp           = True           ,
                           retall         = True           )

    # Done
    return outputs
 
def obj_f(x,project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """
        
    obj_list = project.obj_f(x)
    obj = 0
    for this_obj in obj_list:
        obj = obj+this_obj
    
    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)
        
        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         df(x), list[nobj x dim]
        scipy_slsqp: df(x), ndarray[dim]
    """    
    
    dobj_list = project.obj_df(x)
    dobj=[0.0]*len(dobj_list[0])
    
    for this_dobj in dobj_list:
        idv=0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv]+this_dv_dobj;
            idv+=1
    dobj = array( dobj )
    
    return dobj

def con_ceq(x,project):
    """ cons = con_ceq(x,project)
        
        Equality Constraint Functions
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         ceq(x) = 0.0, list[nceq]
        scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """
    
    cons = project.con_ceq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
        
    return cons

def con_dceq(x,project):
    """ dcons = con_dceq(x,project)
        
        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dceq(x), list[nceq x dim]
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    
    dcons = project.con_dceq(x)

    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return dcons

def con_cieq(x,project):
    """ cons = con_cieq(x,project)
        
        Inequality Constraints
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         cieq(x) < 0.0, list[ncieq]
        scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """
    
    cons = project.con_cieq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
    
    return -cons
    
def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)
        
        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dcieq(x), list[ncieq x dim]
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    
    dcons = project.con_dcieq(x)
    
    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return -dcons
