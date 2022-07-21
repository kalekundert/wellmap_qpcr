#!/usr/bin/env python3

from . import biorad

# In the future, I'll want these functions to make an effort at guessing the 
# correct file format.  But for now I only have one, so it's kinda trivial.

def load_cq(path):
    return biorad.load_cq(path)

def load_trace(path):
    return biorad.load_trace(path)

def load_melt(path):
    return biorad.load_melt(path)

