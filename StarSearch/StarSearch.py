#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.eso import Eso

class ESOquery(object):
    """
    Lets check what instruments observed a given star!!!
    """
    def __init__(self, user):
        self.user = str(user)
        self.eso = Eso()
        self.eso.login(self.user)
        
    def searchStarDates(self, star):
        """
        Lets see when the star was observed
        """
        search = self.eso.query_main(column_filters={'target': star})
        result = np.array(search['Release_Date'])
        return result
    
    def searchStarInstrument(self, star):
        """
        Lets see what instrument observed the star
        """
        searchResult = self.eso.query_main(column_filters={'target': star})
        instruments = np.unique(np.array(searchResult['Instrument']), 
                                return_counts=True)
        instrumentDict = dict()
        for i, j in enumerate(instruments[0]):
            instrumentDict[j] = instruments[1][i]
        return instrumentDict