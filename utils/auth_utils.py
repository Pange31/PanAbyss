#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 15:27:29 2025

@author: fgraziani
"""

from functools import wraps

from config import *



AUTHORIZATION = check_authorization()

print(f"Writing operations enable : {AUTHORIZATION}")

#This function is used to check authorization for every writing or admin operations
def require_authorization(func):
    """Check authorization in BO functions."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        if not AUTHORIZATION:
            print(f"Blocked unauthorized admin callback: {func.__name__}")
            return "🚫 Unauthorized access."
        return func(*args, **kwargs)
    return wrapper