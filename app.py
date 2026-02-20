#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 11:58:49 2025

@author: fgraziani
"""

from dash import Dash, DiskcacheManager
import diskcache
import dash_bootstrap_components as dbc
import os
from config import *
from logger import setup_logger
import logging
import os


#Database version. This version is associated with a database model; 
#compatibility with other versions is not guaranteed.
DB_VERSION="1.0.0"

#Mode server : if the tool is not installed locally
SERVER_MODE = is_server_mode()

#Admin mode : if set to true in the conf file, then it is required to log in
#This mode allow to write into database and to load files
#If PanAbyss is installed locally then this option will be ignored
ADMIN_MODE = is_admin_mode()

BLOCK_ADMIN_FUNCTIONNALITIES = SERVER_MODE and not ADMIN_MODE

success_style = {"color": "green", "marginTop": "10px"}
warning_style = {"color": "orange", "marginTop": "10px"}
error_style = {"color": "red", "marginTop": "10px"}

logger = setup_logger(name="panabyss_logger")

cache = diskcache.Cache("./cache")
background_callback_manager = DiskcacheManager(cache)

app = Dash(__name__,
            suppress_callback_exceptions=True,
            #external_stylesheets=[dbc.themes.BOOTSTRAP],
            background_callback_manager=background_callback_manager)

server = app.server


    

