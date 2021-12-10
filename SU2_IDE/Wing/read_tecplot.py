'''
Tecplot ASCII file reader for python
author: awaldm
date: Feb 2014
'''
import numpy as np
import pandas as pd

def parse_vars(vars_in):
    variables = [item.strip() for item in vars_in]
    variables = [item.replace("'","") for item in variables]
    #print variables
    return variables

# this is from http://pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


def read_tecplot(filename, verbose):

    datafile = open(filename, "r")

    struct=1
    variables = []
    zones = []
    header_lines = 0
    varmode = 0
    data_match = 0
    mode = 0
    # mode 0 = reading header
    # mode 1 = reading VARIABLES
    # mode 2 = reading ZONES
    # mode 3 = data
    dfs = []  # list of dataframes
    els = []
    for line in datafile:
        comm_match = line.find("#")
        title_match = line.find("TITLE")
        var_match = line.find("VARIABLES")
        zone_match = line.find("ZONE T")
        dt_match = line.find("DT=(")
        nodes_match = line.find("Nodes=")
        elems_match = line.find("Elements=")

        if title_match >= 0:
            header_lines += 1
            title_line = line
            if verbose:
                print("Title on line " + str(header_lines) + ": " + str(title_line))
            mode = 0
        elif var_match >= 0:  # match lines starting with variable keyword, add the following var names to list
            header_lines += 1
            # split and get words enclosed by double quotations, return only the odd elements
            l = line.split('"')[1::2]
            l = ['"{}"'.format(word) for word in l]
            if verbose:
                print(l)
            variables = l

            if verbose:
                print("Variables on line " + str(header_lines) + ": " + str(variables))
            mode = 1

        elif zone_match >= 0:
            if mode == 3: # that is, if we have been reading a previous zone
                if struct == 1: #only if we had no nodes or other info, otherwise we'd add double the amount of DFs
                    if verbose:
                        print 'found new zone'
                    dfs.append(pd.DataFrame(data))
                mode = 2
#                continue
            else:
                # if we are here, this means that the VARIABLES keyword is finished. We can safely parse the variables now:
                # this makes only sense if this is the first occurrence of the ZONE keyword
                variables = parse_vars(variables)

            # parse the actual line
            header_lines += 1
            zone = line
            zone = line.split("ZONE T=", 1).pop()
            zones.append(zone)
            varmode = 0
            if verbose:
                print("Zone on line " + str(header_lines) + ": " + str(zone))
            mode = 2  # vars are apparently finished, shift mode to 2
            data = []


        elif nodes_match >= 0:
            # match number of nodes between substrings
            a = "Nodes="
            b = ","
            nodes = int(line.split(a)[-1].split(b)[0])
            if verbose:
                print ("found nodes= " + str(nodes))
            struct = 0 # if Nodes keyword is present, we read unstructured data
            header_lines += 1
            #print str(mode)
            if elems_match >= 0:  # in case of Nodes= and Elements= in the same line
                a = "Elements="
                b = ","
                elements = int(line.split(a)[-1].split(b)[0])
                if verbose:
                    print ("found elements= " + str(elements))


        elif elems_match >= 0:
            a = "Elements="
            b = ","
            elements = int(line.split(a)[-1].split(b)[0])
            if verbose:
                print ("found elements= " + str(elements))

            header_lines += 1

        elif dt_match >= 0:
            header_lines += 1
            if verbose:
                print "found DT, mode is " + str(mode)
            data_match = 1
            mode = 3  # DT is for sure the last line of ZONE
        elif (line.find('STRANDID') >= 0) or (line.find('DATAPACKING') >= 0):
            header_lines += 1
            if verbose:
                print "found STRANDID or DATAPACKING, mode remains at " + str(mode)
        # mode remains unchanged, reading data
        elif mode == 3:
            if struct == 1 or len(data) < nodes: # read data
                data.append(line.split())
                if struct == 0 and len(data) == nodes:
                    if verbose:
                        print ("found last unstructured data line, overall length: " + str(len(data)))
                        print (type(data))
                        print len(els)

                    elems=[]
                    dfs.append(pd.DataFrame(data))

            elif struct == 0 and (len(elems) < elements):
                elems.append(line.split())
                if len(elems) == elements:
                    els.append(elems)

            else:
                continue  # finished reading data, now doing elements

        elif (mode==2 and is_number(line.split()[0])): # this is true if we are below the ZONE keyword and at a data line
            mode = 3
            if verbose:
                print 'found first data line'
            data.append(line.split())


        elif mode == 1:  # lines between VARIABLES and ZONE keywords: these may be additional variables or AUXDATA 
            if verbose:
                print("More variables on line " + str(header_lines) + ": " + str(line))
            if line.startswith( 'DATASETAUX' ): # TODO: match other AUX data
                continue
            else:
                variables.append(line)
            header_lines += 1
        else:
            continue
    else: # for ..else meaning we reached EOF in this context
        if mode == 3:
            if struct == 1 or len(data) < nodes:
                dfs.append(pd.DataFrame(data))

    zones = [item.strip() for item in zones]



    # set variables as column names of each pandas dataframe
    for df in dfs:
        df.columns = variables
        df.apply(pd.to_numeric)

    return zones, dfs, els
