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
        self.eso.ROW_LIMIT = -1 #unlimited number of search results
        
        
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
    
    
    def retriveStarData(self, star, downloadPath = None):
        """
        Lets download the data
        """
        checkInstruments = self.searchStarInstrument(star)
        esoInst = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])
        
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***'.format(j))
            if j in checkInstruments:
                print('Downloading {0} data\n'.format(j))
                table = self.eso.query_main(column_filters = {'instrument': j, 
                                                              'target': star})
                if downloadPath:
                    self.eso.retrieve_data(table['Dataset ID'], 
                                           destination = downloadPath)
                else:
                    self.eso.retrieve_data(table['Dataset ID'])
            else:
                print('No {0} data\n'.format(j))
        return 0