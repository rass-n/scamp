# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 22:25:46 2016

@author: Rasmus Nordfang

Module for simple climate modeling


"""

#from classes import *
from classes import dim
from classes import ocean
from classes import ice
from classes import world
from classes import atmos

from update import *
from ode_solver import *
from salinity import *
from ice_volume import *
from density import *
from heat import *

from animation import ani