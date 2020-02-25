#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#First lets import the ESO phase 3 module
from starsearch.phase3Archive import ESOquery 

#Connect with your ESO account
user = ESOquery('ESOusername') 

#Our package was built to search spectra only of FEROS, UVES, HARPS and ESPRESSO
print(user.instruments)