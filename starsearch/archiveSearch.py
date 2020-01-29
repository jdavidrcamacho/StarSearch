#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.eso import Eso

class ESOquery(object):
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
        self.user = str(user)
        self.eso = Eso()
        self.eso.login(self.user)
        self.eso.ROW_LIMIT = -1 #unlimited number of search results
        
        
    def searchDates(self, star):
        """
        Checks the day the data was released 
        
        Parameters
        ----------
        star: str
            Name of the star
            
        Returns
        -------
        result: array
            Array with the date of the data release
        """
        search = self.eso.query_main(column_filters={'target': star})
        result = np.array(search['Release_Date'])
        return result
    
    
    def searchInstrument(self, star):
        """
        Checks which instruments where used
        
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
    
    
    def _searchAndDownload(self, star, instrument, downloadPath, date):
        """
        Download ESO data
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Instrument we are searching the darta
        date: str
            Download only the data past a certain date
            
        Returns
        -------
        """
        print('Downloading {0} data\n'.format(star))
        table = self.eso.query_main(column_filters = {'instrument': instrument, 
                                                      'target': star})
        if downloadPath:
            self.eso.retrieve_data(table['Dataset ID'], 
                                   destination = downloadPath)
        else:
            self.eso.retrieve_data(table['Dataset ID'])
        
    
    def getStar(self, star, downloadPath = None , date = None):
        """
        Download ESO data
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
            
        Returns
        -------
        
        """
        checkInstruments = self.searchStarInstrument(star)
        esoInst = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])
        
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(j, star, downloadPath, date)
#                print('Downloading {0} data\n'.format(j))
#                table = self.eso.query_main(column_filters = {'instrument': j, 
#                                                              'target': star})
#                if downloadPath:
#                    self.eso.retrieve_data(table['Dataset ID'], 
#                                           destination = downloadPath)
#                else:
#                    self.eso.retrieve_data(table['Dataset ID'])
            else:
                print('No {0} data\n'.format(j))
        return 0