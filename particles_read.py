# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import struct

def NStepStr(ns):
    nsstr = str(ns)
    while len(nsstr) < 4:
        nsstr = '0'+nsstr

    return nsstr


def ReadPartData(ns):

    nstepstr = NStepStr(ns)
    fname = r"./problem/out/particles."+nstepstr+".partdbl"
    h_lines = 0 
    val_dict = {}

	#READ HEADER. 
    with open(fname, "rb") as f:
        for line in f:
            if h_lines < 13:
                if h_lines > 0 and h_lines < 13:
                    val_dict.update({line.split()[1].decode('utf8'):[i.decode('utf8') for i in line.split()[2:]]})
                h_lines += 1

    fdata = open(fname, "rb")
 
    cnt = 0
    while (cnt < h_lines):
        fdata.readline()
        cnt += 1
        
    data = fdata.read()
    fdata.close()

    dt = np.dtype({'names':val_dict['field_names'], 'formats':['('+i+',)<d' for i in val_dict['field_dim']]})

    val = np.frombuffer(data, dtype=dt)

    for i in range(len(val_dict['field_names'])):
        name = val_dict['field_names'][i]
        val_dict.update({name:val[name]})

    return val_dict

