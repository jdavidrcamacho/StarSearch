#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.esocas import Eso
from sys import maxsize
np.set_printoptions(threshold = maxsize)

class ESOquery():
    """
    ESO query class
    
    Parameters
    ----------
    user: str
        User name used in ESO website
        
    Returns
    -------
    
    """
    def __init__(self, user):
        self.user = str(user) #user name 
        self.eso = Eso()
        self.eso.login(self.user) #login to eso archive
        self.eso.ROW_LIMIT = -1 #unlimited number of search results = -1
        self.surveys = self.eso.list_surveys() #list of available surveys
        self.instruments = ['FEROS', 'UVES', 'HARPS', 'ESPRESSO']
        
        
    def searchInstruments(self, star):
        """
        Checks which instruments where used for the given star and how many
        observations it made
        
        Parameters
        ----------
        star: str
            Name of the star
            
        Returns
        -------
        instrumentDict: dict
            Instruments and number of observations
        """
        searchResult = self.eso.query_main(column_filters={'target': star})
        instruments = np.unique(np.array(searchResult['Instrument']), 
                                return_counts=True)
        instrumentDict = dict()
        for i, j in enumerate(instruments[0]):
            instrumentDict[j] = instruments[1][i]
        return instrumentDict
    
    
    def searchStar(self, star):
        #url = "http://archive.eso.org/wdb/wdb/adp/phase3_spectral/form"
        search = self.eso.query_surveys(instrument = 'HARPS', 
                                        target = star) 
        

        # checkInstruments = self.searchInstrument(star)
        # esoInst = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])
        # for i, j in enumerate(esoInst):
        #     print('\n*** Searching for {0} results ***\n'.format(j))
        #     if j in checkInstruments:
        #         self._searchAndDownload(star, j, downloadPath, date, calib)
        #     else:
        #         print('No {0} data\n'.format(j))
        # print('\n*** Done ***\n')
        # return 0

        return search
