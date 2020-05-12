#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:35:29 2020

@author: camacho
"""
from starsearch.phase3Archive import ESOquery

user = ESOquery('jdcamacho')

print('\nSearching star...')
star = user.searchStar('hd41248')
print('\nSearching instruments...')
instruments = user.searchInstruments('hd41248')
print('\nSearching ESPRESSO spectra...')
espresso = user.searchInstrumentSpectra('hd41248', instrument='ESPRESSO')