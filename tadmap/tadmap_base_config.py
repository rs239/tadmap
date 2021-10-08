#!/usr/bin/env python

PROJ_DIR="/afs/csail.mit.edu/u/r/rsingh/work/vartad" #"../.."

DEBUG_MODE=True #False

def dbg_print(s, *args, **kwargs):
    if DEBUG_MODE or s[:5] != "Flag ":
        print(s, *args, **kwargs)


import logging

tadmap_loglevel = logging.WARNING #WARNING


def tadmap_debug(*args, **kwargs):
    if tadmap_loglevel <= logging.DEBUG: print("DEBUG: ", *args, **kwargs)

def tadmap_info(*args, **kwargs):
    if tadmap_loglevel <= logging.INFO: print("INFO: ", *args, **kwargs)

def tadmap_warning(*args, **kwargs):
    if tadmap_loglevel <= logging.WARNING: print("WARNING: ", *args, **kwargs)

def tadmap_error(*args, **kwargs):
    if tadmap_loglevel <= logging.ERROR: print("ERROR: ", *args, **kwargs)

    
########## for maintenance ###################
# def noop(*args, **kwargs):
#     pass
#
# logging.info = print
# logging.debug = noop
##############################################
