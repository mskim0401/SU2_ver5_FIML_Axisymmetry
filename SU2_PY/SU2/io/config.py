#!/usr/bin/env python

## \file config.py
#  \brief python package for config 
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import numpy as np
from ..util import bunch, ordered_bunch, switch
from .tools import *
from config_options import *

try:
    from collections import OrderedDict
except ImportError:
    from ..util.ordered_dict import OrderedDict

inf = 1.0e20


# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class Config(ordered_bunch):
    """ config = SU2.io.Config(filename="")
        
        Starts a config class, an extension of 
        ordered_bunch()
        
        use 1: initialize by reading config file
            config = SU2.io.Config('filename')
        use 2: initialize from dictionary or bunch
            config = SU2.io.Config(param_dict)
        use 3: initialize empty
            config = SU2.io.Config()
        
        Parameters can be accessed by item or attribute
        ie: config['MESH_FILENAME'] or config.MESH_FILENAME
        
        Methods:
            read()       - read from a config file
            write()      - write to a config file (requires existing file)
            dump()       - dump a raw config file
            unpack_dvs() - unpack a design vector 
            diff()       - returns the difference from another config
            dist()       - computes the distance from another config
    """    

    _filename = 'config.cfg'
    
    def __init__(self,*args,**kwarg):
        
        # look for filename in inputs
        if args and isinstance(args[0],str):
            filename = args[0]
            args = args[1:]
        elif kwarg.has_key('filename'):
            filename = kwarg['filename']
            del kwarg['filename']
        else:
            filename = ''
        
        # initialize ordered bunch
        super(Config,self).__init__(*args,**kwarg)
        
        # read config if it exists
        if filename:
            try:
                self.read(filename)
            except IOError:
                print 'Could not find config file: %s' % filename
	    except:
		print 'Unexpected error: ',sys.exc_info()[0]
		raise
        
        self._filename = filename
    
    def read(self,filename):
        """ reads from a config file """
        konfig = read_config(filename)
        self.update(konfig)
        
    def write(self,filename=''):
        """ updates an existing config file """
        if not filename: filename = self._filename
        assert os.path.exists(filename) , 'must write over an existing config file'
        write_config(filename,self)
        
    def dump(self,filename=''):
        """ dumps all items in the config bunch, without comments """
        if not filename: filename = self._filename
        dump_config(filename,self)
    
    def __getattr__(self,k):
        try:
            return super(Config,self).__getattr__(k)
        except AttributeError:
            raise AttributeError , 'Config parameter not found'
        
    def __getitem__(self,k):
        try:
            return super(Config,self).__getitem__(k)
        except KeyError:
            raise KeyError , 'Config parameter not found: %s' % k

    def unpack_dvs(self,dv_new,dv_old=None):
        """ updates config with design variable vectors
            will scale according to each DEFINITION_DV scale parameter
            
            Modifies:
                DV_KIND
                DV_MARKER
                DV_PARAM
                DV_VALUE_OLD
                DV_VALUE_NEW
            
            Inputs:
                dv_new - list or array of new dv values
                dv_old - optional, list or array of old dv values, defaults to zeros
                         
        """
        
        dv_new = copy.deepcopy(dv_new)
        dv_old = copy.deepcopy(dv_old)
        
        # handle unpacking cases
        def_dv = self['DEFINITION_DV']
	if int(self.NPOIN) == 0:
	    n_dv   = sum(def_dv['SIZE'])
	else:
	    if self.KIND_TRAIN_NN == 'WEIGHTS':
		n_dv = len(dv_new)
	    else: 
		n_dv   = int(self.NPOIN)
	if not dv_old: dv_old = [0.0]*n_dv
	#if self.NPOIN == 0:
	#assert len(dv_new) == len(dv_old) , 'unexpected design vector length'
        
        # handle param
        param_dv = self['DV_PARAM']

        # apply scale
        dv_scales = def_dv['SCALE']

        k = 0
	if int(self.NPOIN) == 0 :
	    for i, dv_scl in enumerate(dv_scales):
		for j in range(def_dv['SIZE'][i]):
		    dv_new[k] = dv_new[k]*dv_scl;
		    dv_old[k] = dv_old[k]*dv_scl;
		    k = k + 1
        else:
	    for i in range(len(dv_new)) :
		dv_new[i] = dv_new[i]*dv_scales[0]
		dv_old[i] = dv_old[i]*dv_scales[0]
	    
        # Change the parameters of the design variables

        self['DV_KIND'] = def_dv['KIND']
        param_dv['PARAM'] = def_dv['PARAM']
        param_dv['FFDTAG'] = def_dv['FFDTAG']
        param_dv['SIZE']   = def_dv['SIZE']
	
	if int(self.NPOIN) == 0 :
	    self.update({ 'DV_MARKER'        : def_dv['MARKER'][0] ,
	                  'DV_VALUE_OLD'     : dv_old              ,
	                  'DV_VALUE_NEW'     : dv_new              })
	else :
	    self.update({ 'DV_MARKER'        : def_dv['MARKER'][0] ,
	                  'DV_VALUE_OLD'     : dv_old              ,
	                  'DV_VALUE_NEW'     : dv_new              })		
	    
    def __eq__(self,konfig):
        return super(Config,self).__eq__(konfig)
    def __ne__(self,konfig):
        return super(Config,self).__ne__(konfig)
    
    
    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key,value in self.iteritems():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)    
    
    def diff(self,konfig):
        """ compares self to another config
            
            Inputs: 
                konfig - a second config
            
            Outputs:
                config_diff - a config containing only the differing 
                    keys, each with values of a list of the different 
                    config values.
                for example: 
                config_diff.MATH_PROBLEM = ['DIRECT','CONTINUOUS_ADJOINT']
                
        """
        
        keys = set([])
        keys.update( self.keys() )
        keys.update( konfig.keys() )
        
        konfig_diff = Config()
        
        for key in keys:
            value1 = self.get(key,None)
            value2 = konfig.get(key,None)
            if not value1 == value2:
                konfig_diff[key] = [value1,value2]
        
        return konfig_diff
    
    def dist(self,konfig,keys_check='ALL'):
        """ calculates a distance to another config
            
            Inputs: 
                konfig     - a second config
                keys_check - optional, a list of keys to check
            
            Outputs:
                distance   - a float
                
            Currently only works for DV_VALUE_NEW and DV_VALUE_OLD
            Returns a large value otherwise
                
        """        

        konfig_diff = self.diff(konfig)
        
        if keys_check == 'ALL':
            keys_check = konfig_diff.keys()
    
        distance = 0.0
        
        for key in keys_check:
            if konfig_diff.has_key(key):
                
                val1 = konfig_diff[key][0]
                val2 = konfig_diff[key][1]
                
                if key in ['DV_VALUE_NEW',
                           'DV_VALUE_OLD']:
                    val1 = np.array( val1 )
                    val2 = np.array( val2 )
                    this_diff = np.sqrt( np.sum( (val1-val2)**2 ) )
                
                else:
                    print 'Warning, unexpected config difference'
                    this_diff = inf
                    
                distance += this_diff
            
            #: if key different
        #: for each keys_check
        
        return distance
    
    def __repr__(self):
        #return '<Config> %s' % self._filename
        return self.__str__()
    
    def __str__(self):
        output = 'Config: %s' % self._filename
        for k,v in self.iteritems():
            output +=  '\n    %s= %s' % (k,v)
        return output
#: class Config







# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_config(filename):
    """ reads a config file """
      
    # initialize output dictionary
    data_dict = OrderedDict()
    
    input_file = open(filename)
    
    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:
            break
        
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equals sign
        line = line.split("=",1)
        this_param = line[0].strip()
        this_value = line[1].strip()
        
        assert not data_dict.has_key(this_param) , ('Config file has multiple specifications of %s' % this_param )
        for case in switch(this_param):
            
            # comma delimited lists of strings with or without paren's
            if case("MARKER_EULER")      : pass
            if case("MARKER_FAR")        : pass
            if case("MARKER_PLOTTING")   : pass
            if case("MARKER_MONITORING") : pass
            if case("MARKER_SYM")        : pass
            if case("DV_KIND")           : 
                # remove white space
                this_value = ''.join(this_value.split())   
                # remove parens
                this_value = this_value.strip('()')
                # split by comma
                data_dict[this_param] = this_value.split(",")
                break
            
            # semicolon delimited lists of comma delimited lists of floats
            if case("DV_PARAM"):
                # remove white space
                info_General = ''.join(this_value.split())
                # split by semicolon
                info_General = info_General.split(';')
                # build list of dv params, convert string to float
                dv_Parameters = []
                dv_FFDTag     = []
                dv_Size       = []

                for this_dvParam in info_General:
                    this_dvParam = this_dvParam.strip('()')
                    this_dvParam = this_dvParam.split(",")
                    this_dvSize  = 1

                    # if FFD change the first element to work with numbers and float(x)
                    if data_dict["DV_KIND"][0] in ['FFD_SETTING','FFD_ANGLE_OF_ATTACK','FFD_CONTROL_POINT','FFD_NACELLE','FFD_GULL','FFD_TWIST_2D','FFD_TWIST','FFD_ROTATION','FFD_CAMBER','FFD_THICKNESS','FFD_CONTROL_POINT_2D','FFD_CAMBER_2D','FFD_THICKNESS_2D']:
                        this_dvFFDTag = this_dvParam[0]
                        this_dvParam[0] = '0'
                    else:
                        this_dvFFDTag = []

                    if not data_dict["DV_KIND"][0] in ['NO_DEFORMATION']:
                        this_dvParam = [ float(x) for x in this_dvParam ]

                    if data_dict["DV_KIND"][0] in ['FFD_CONTROL_POINT_2D']:
                        if this_dvParam[3] == 0 and this_dvParam[4] == 0:
                            this_dvSize = 2

                    if data_dict["DV_KIND"][0]in ['FFD_CONTROL_POINT']:
                        if this_dvParam[4] == 0 and this_dvParam[5] == 0 and this_dvParam[6] == 0:
                            this_dvSize = 3

                    dv_FFDTag     = dv_FFDTag     + [this_dvFFDTag]
                    dv_Parameters = dv_Parameters + [this_dvParam]
                    dv_Size       = dv_Size       + [this_dvSize]
            
            # store in a dictionary
                dv_Definitions = { 'FFDTAG' : dv_FFDTag     ,
                                   'PARAM'  : dv_Parameters ,
                                   'SIZE'   : dv_Size}

                data_dict[this_param] = dv_Definitions
                break
            
            # comma delimited lists of floats
            if case("DV_VALUE_OLD")    : pass
            if case("DV_VALUE_NEW")    : pass
            if case("DV_VALUE")        :           
                # remove white space
                this_value = ''.join(this_value.split())                
                # split by comma, map to float, store in dictionary
                data_dict[this_param] = map(float,this_value.split(","))
                break              

            # float parameters
            if case("MACH_NUMBER")            : pass
            if case("AoA")                    : pass
            if case("FIN_DIFF_STEP")          : pass
            if case("CFL_NUMBER")             : pass
            if case("HB_PERIOD")    : pass
            if case("WRT_SOL_FREQ")           :
                data_dict[this_param] = float(this_value)
                break   
            
            # int parameters
            if case("NUMBER_PART")            : pass
            if case("AVAILABLE_PROC")         : pass
            if case("EXT_ITER")               : pass
            if case("TIME_INSTANCES")         : pass
            if case("UNST_ADJOINT_ITER")      : pass
            if case("ITER_AVERAGE_OBJ")       : pass
            if case("ADAPT_CYCLES")           :
                data_dict[this_param] = int(this_value)
                break                
            
            
            # unitary design variable definition
            if case("DEFINITION_DV"):
                # remove white space
                this_value = ''.join(this_value.split())                
                # split into unitary definitions
                info_Unitary = this_value.split(";")
                # process each Design Variable
                dv_Kind       = []
                dv_Scale      = []
                dv_Markers    = []
                dv_FFDTag     = []
                dv_Parameters = []
                dv_Size       = []

                for this_General in info_Unitary:
                    if not this_General: continue
                    # split each unitary definition into one general definition
                    info_General = this_General.strip("()").split("|") # check for needed strip()?
                    # split information for dv Kinds
                    info_Kind    = info_General[0].split(",")
                    # pull processed dv values
                    this_dvKind       = get_dvKind( int( info_Kind[0] ) )     
                    this_dvScale      = float( info_Kind[1] )
                    this_dvMarkers    = info_General[1].split(",")
                    this_dvSize       = 1

                    if this_dvKind=='MACH_NUMBER' or this_dvKind=='AOA':
                        this_dvParameters = []
                    else:
                        this_dvParameters = info_General[2].split(",")
                        # if FFD change the first element to work with numbers and float(x), save also the tag
                        if this_dvKind in ['FFD_SETTING','FFD_ANGLE_OF_ATTACK','FFD_CONTROL_POINT','FFD_NACELLE','FFD_GULL','FFD_TWIST','FFD_TWIST_2D','FFD_TWIST_ANGLE','FFD_ROTATION','FFD_CAMBER','FFD_THICKNESS','FFD_CONTROL_POINT_2D','FFD_CAMBER_2D','FFD_THICKNESS_2D']:
                          this_dvFFDTag = this_dvParameters[0]
                          this_dvParameters[0] = '0'
                        else:
                          this_dvFFDTag = []
                        
                        this_dvParameters = [ float(x) for x in this_dvParameters ]

                        if this_dvKind in ['FFD_CONTROL_POINT_2D']:
                            if this_dvParameters[3] == 0 and this_dvParameters[4] == 0:
                                this_dvSize = 2

                        if this_dvKind in ['FFD_CONTROL_POINT']:
                            if this_dvParameters[4] == 0 and this_dvParameters[5] == 0 and this_dvParameters[6] == 0:
                                this_dvSize = 3

                    # add to lists
                    dv_Kind       = dv_Kind       + [this_dvKind]
                    dv_Scale      = dv_Scale      + [this_dvScale]
                    dv_Markers    = dv_Markers    + [this_dvMarkers]
                    dv_FFDTag     = dv_FFDTag     + [this_dvFFDTag]
                    dv_Parameters = dv_Parameters + [this_dvParameters]
                    dv_Size       = dv_Size       + [this_dvSize]
                # store in a dictionary
                dv_Definitions = { 'KIND'   : dv_Kind       ,
                                   'SCALE'  : dv_Scale      ,
                                   'MARKER' : dv_Markers    ,
                                   'FFDTAG' : dv_FFDTag     ,
                                   'PARAM'  : dv_Parameters ,
                                   'SIZE'   : dv_Size}

                # save to output dictionary
                data_dict[this_param] = dv_Definitions
                break  
            
            # unitary objective definition
            if case('OPT_OBJECTIVE'):
                # remove white space
                this_value = ''.join(this_value.split())
                #split by + 
                this_def={}
                this_value = this_value.split(";")
                
                for  this_obj in this_value:       
                    # split by scale
                    this_obj = this_obj.split("*")
                    this_name  = this_obj[0]
                    this_scale = 1.0
                    if len(this_obj) > 1:
                        this_scale = float( this_obj[1] )
                    this_def.update({ this_name : {'SCALE':this_scale} })
                # save to output dictionary
                data_dict[this_param] = this_def
                break
            
            # unitary constraint definition
            if case('OPT_CONSTRAINT'):
                # remove white space
                this_value = ''.join(this_value.split())                    
                # check for none case
                if this_value == 'NONE':
                    data_dict[this_param] = {'EQUALITY':OrderedDict(), 'INEQUALITY':OrderedDict()}
                    break                    
                # split definitions
                this_value = this_value.split(';')
                this_def = OrderedDict()
                for this_con in this_value:
                    if not this_con: continue # if no definition
                    # defaults
                    this_obj = 'NONE'
                    this_sgn = '='
                    this_scl = 1.0
                    this_val = 0.0
                    # split scale if present
                    this_con = this_con.split('*')
                    if len(this_con) > 1:
                        this_scl = float( this_con[1] )
                    this_con = this_con[0]
                    # find sign
                    for this_sgn in ['<','>','=']:
                        if this_sgn in this_con: break
                    # split sign, store objective and value
                    this_con = this_con.strip('()').split(this_sgn)
                    assert len(this_con) == 2 , 'incorrect constraint definition'
                    this_obj = this_con[0]
                    this_val = float( this_con[1] )
                    # store in dictionary
                    this_def[this_obj] = { 'SIGN'  : this_sgn ,
                                           'VALUE' : this_val ,
                                           'SCALE' : this_scl  }
                #: for each constraint definition
                # sort constraints by type
                this_sort = { 'EQUALITY'   : OrderedDict() ,
                              'INEQUALITY' : OrderedDict()  }
                for key,value in this_def.iteritems():
                    if value['SIGN'] == '=':
                        this_sort['EQUALITY'][key]   = value
                    else:
                        this_sort['INEQUALITY'][key] = value
                #: for each definition                
                # save to output dictionary
                data_dict[this_param] = this_sort
                break
            
            # otherwise
            # string parameters
            if case():
                data_dict[this_param] = this_value
                break              
            
            #: if case DEFINITION_DV
                        
        #: for case
        
    #: for line

    #hack - twl
    if not data_dict.has_key('DV_VALUE_NEW'):
        data_dict['DV_VALUE_NEW'] = [0]
    if not data_dict.has_key('DV_VALUE_OLD'):
        data_dict['DV_VALUE_OLD'] = [0]
    if not data_dict.has_key('OPT_ITERATIONS'):
        data_dict['OPT_ITERATIONS'] = 100
    if not data_dict.has_key('OPT_ACCURACY'):
        data_dict['OPT_ACCURACY'] = 1e-10
    if not data_dict.has_key('OPT_BOUND_UPPER'):
        data_dict['OPT_BOUND_UPPER'] = 1e10
    if not data_dict.has_key('OPT_BOUND_LOWER'):
        data_dict['OPT_BOUND_LOWER'] = -1e10
    if not data_dict.has_key('OPT_COMBINE_OBJECTIVE'):
        data_dict['OPT_COMBINE_OBJECTIVE'] = "NO"
    if not data_dict.has_key('VALUE_OBJFUNC_FILENAME'):
      data_dict['VALUE_OBJFUNC_FILENAME'] = 'of_eval.dat'
    if not data_dict.has_key('GRAD_OBJFUNC_FILENAME'):
      data_dict['GRAD_OBJFUNC_FILENAME'] = 'of_grad.dat'
    if not data_dict.has_key('NPOIN'):
	data_dict['NPOIN'] = 0
    if not data_dict.has_key('MIN_MAX_GRAD'):
	data_dict['MIN_MAX_GRAD'] = 0
    if not data_dict.has_key('LEARN_RATE'):
	data_dict['LEARN_RATE'] = 0.01    
    if not data_dict.has_key('NUM_CASES'):
	data_dict['NUM_CASES'] = 0
    else : data_dict['NUM_CASES'] = int(data_dict['NUM_CASES'])
    if not data_dict.has_key('CONFIG_I'): data_dict['CONFIG_I'] = 0
    if not data_dict.has_key('STOCHASTIC_GRAD'):
	data_dict['STOCHASTIC_GRAD'] = "NO"
    if not data_dict.has_key('MULTI_MESH'):
	data_dict['MULTI_MESH'] = "NO"
 
    return data_dict
    
#: def read_config()



# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_config(filename,param_dict):
    """ updates an existing config file """
    
    temp_filename = "temp.cfg"
    shutil.copy(filename,temp_filename)
    output_file = open(filename,"w")
    
    # break pointers
    param_dict = copy.deepcopy(param_dict)
    npoin = int(param_dict.NPOIN)
    if not npoin == 0:
	if os.path.isfile('beta_fiml.dat') : 
	    os.remove('beta_fiml.dat')
	output_beta_fiml = open('beta_fiml.dat',"w")
	output_beta_fiml.seek(0)
    
    for raw_line in open(temp_filename):
        # remove line returns
        line = raw_line.strip('\r\n')
        
        # make sure it has useful data
        if not "=" in line:
            output_file.write(raw_line)
            continue
        
        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        old_value  = line[1].strip()
        
        # skip if parameter unwanted
        if not param_dict.has_key(this_param):
            output_file.write(raw_line)
            continue
        
        # start writing parameter
        new_value = param_dict[this_param] 
        output_file.write(this_param + "= ")
        
        # handle parameter types
        for case in switch(this_param):  
              
            # comma delimited list of floats
            if case("DV_VALUE_NEW") : pass
            if case("DV_VALUE_OLD") : 
		n_lists = len(new_value)
		for i_value in range(n_lists):
		    if (npoin == 0 or i_value == 0):
			output_file.write("%s" % new_value[i_value])
		    if (npoin == 0 and i_value+1 < n_lists):
			output_file.write(", ")
		break
	    
	    if case("DV_VALUE") :
		n_lists = len(new_value)
		for i_value in range(n_lists):
		    if (npoin == 0 or i_value == 0):
			output_file.write("%s" % new_value[i_value])
		    if not npoin == 0 :
			output_beta_fiml.write("%s\n" % new_value[i_value])
		    if (npoin == 0 and i_value+1 < n_lists):
			output_file.write(", ")
		break		
            
            # comma delimited list of strings no paren's
            if case("DV_KIND")            : pass
            if case("TASKS")              : pass
            if case("GRADIENTS")          :            
                if not isinstance(new_value,list):
                    new_value = [ new_value ]
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")               
                break            
            
            # comma delimited list of strings inside paren's
            if case("MARKER_EULER")      : pass
            if case("MARKER_FAR")        : pass
            if case("MARKER_PLOTTING")   : pass
            if case("MARKER_MONITORING") : pass
            if case("MARKER_SYM")        : pass            
            if case("DV_MARKER") : 
                if not isinstance(new_value,list):
                    new_value = [ new_value ]                
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value+1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )") 
                break                
            
            # semicolon delimited lists of comma delimited lists
            if case("DV_PARAM") :

                assert isinstance(new_value['PARAM'],list) , 'incorrect specification of DV_PARAM'
                if not isinstance(new_value['PARAM'][0],list): new_value = [ new_value ]
                
                for i_value in range(len(new_value['PARAM'])):

                    output_file.write("( ")
                    this_param_list = new_value['PARAM'][i_value]
                    this_ffd_list = new_value['FFDTAG'][i_value]
                    n_lists = len(this_param_list)
                    
                    if this_ffd_list != []:
                      output_file.write("%s, " % this_ffd_list)
                      for j_value in range(1,n_lists):
                        output_file.write("%s" % this_param_list[j_value])
                        if j_value+1 < n_lists:
                          output_file.write(", ")
                    else:
                      for j_value in range(n_lists):
                        output_file.write("%s" % this_param_list[j_value])
                        if j_value+1 < n_lists:
                          output_file.write(", ")

                    output_file.write(") ")
                    if i_value+1 < len(new_value['PARAM']):
                        output_file.write("; ")
                break
            
            # int parameters
            if case("NUMBER_PART")            : pass
            if case("ADAPT_CYCLES")           : pass
            if case("TIME_INSTANCES")         : pass
            if case("AVAILABLE_PROC")         : pass
            if case("UNST_ADJOINT_ITER")      : pass
            if case("EXT_ITER")               :
                output_file.write("%i" % new_value)
                break
                        
            if case("DEFINITION_DV") :
                n_dv = len(new_value['KIND'])
                if not n_dv:
                    output_file.write("NONE")
                for i_dv in range(n_dv):
                    this_kind = new_value['KIND'][i_dv]
                    output_file.write("( ")
                    output_file.write("%i , " % get_dvID(this_kind) )
                    output_file.write("%s " % new_value['SCALE'][i_dv])
                    output_file.write("| ")
                    # markers                  
                    n_mark = len(new_value['MARKER'][i_dv])
                    for i_mark in range(n_mark):                       
                        output_file.write("%s " % new_value['MARKER'][i_dv][i_mark])
                        if i_mark+1 < n_mark:
                            output_file.write(", ")
                    #: for each marker
                    if not this_kind in ['AOA','MACH_NUMBER']:
                        output_file.write(" | ")
                        # params
                        if this_kind in ['FFD_SETTING','FFD_ANGLE_OF_ATTACK','FFD_CONTROL_POINT','FFD_NACELLE','FFD_GULL','FFD_TWIST_ANGLE','FFD_TWIST','FFD_TWIST_2D','FFD_ROTATION','FFD_CAMBER','FFD_THICKNESS','FFD_CONTROL_POINT_2D','FFD_CAMBER_2D','FFD_THICKNESS_2D','FIML']:
                            n_param = len(new_value['PARAM'][i_dv])
                            output_file.write("%s , " % new_value['FFDTAG'][i_dv])
                            for i_param in range(1,n_param):
                                output_file.write("%s " % new_value['PARAM'][i_dv][i_param])
                                if i_param+1 < n_param:
                                    output_file.write(", ")
                        else:
                            n_param = len(new_value['PARAM'][i_dv])
                            for i_param in range(n_param):
                                output_file.write("%s " % new_value['PARAM'][i_dv][i_param])
                                if i_param+1 < n_param:
                                    output_file.write(", ")
                    
                        #: for each param                    
                    output_file.write(" )")
                    if i_dv+1 < n_dv:
                        output_file.write("; ")
                #: for each dv
                break
            
            if case("OPT_OBJECTIVE"):
                i_name = 0
                for name,value in new_value.iteritems():
                    if i_name>0: output_file.write("; ")
                    output_file.write( "%s * %s" % (name,value['SCALE']) )
                    i_name += 1
                break
            
            if case("OPT_CONSTRAINT"):
                i_con = 0
                for con_type in ['EQUALITY','INEQUALITY']:
                    this_con = new_value[con_type]
                    for name,value in this_con.iteritems():
                        if i_con>0: output_file.write("; ")
                        output_file.write( "( %s %s %s ) * %s" 
                                          % (name, value['SIGN'], value['VALUE'], value['SCALE']) ) 
                        i_con += 1
                    #: for each constraint
                #: for each constraint type
                if not i_con: output_file.write("NONE")
                break
            
            # default, assume string, integer or unformatted float 
            if case():
                output_file.write('%s' % new_value)
                break                         
                
        #: for case
        
        # remove from param dictionary
        del param_dict[this_param]
        
        # next line
        output_file.write("\n")        
        
    #: for each line
    
    # check that all params were used
    for this_param in param_dict.keys():
        if not this_param in ['JOB_NUMBER']:
            print ( 'Warning: Parameter %s not found in config file and was not written' % (this_param) )
        
    if not npoin == 0:
	output_file.close()
	output_beta_fiml.close()
    os.remove( temp_filename )
    
#: def write_config()


def dump_config(filename,config):
    ''' dumps a raw config file with all options in config 
        and no comments
    '''
    
    # HACK - twl
    if config.has_key('DV_VALUE_NEW'):
        config.DV_VALUE = config.DV_VALUE_NEW
        
    config_file = open(filename,'w')
    # write dummy file
    for key in config.keys():
        config_file.write( '%s= 0 \n' % key )
    config_file.close()
    # dump data
    write_config(filename,config)    

