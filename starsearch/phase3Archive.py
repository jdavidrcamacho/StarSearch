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
        self.user = user #user name 
        self.eso = Eso()
        self.eso.login(self.user) #login to eso archive
        self.eso.ROW_LIMIT = -1 #unlimited number of search results = -1
        self.surveys = self.eso.list_surveys() #list of available surveys
        self.instruments = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])
        
        
    def searchInstruments(self, star):
        """
        Checks which instruments observed the given star and how many
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
    
    
    def searchStar(self, star, instrument=None):
        #url = "http://archive.eso.org/wdb/wdb/adp/phase3_spectral/form"
        if instrument:
            search = self.eso.query_surveys(instrument = instrument, 
                                            target = star) 
        else:
            search = self.eso.query_surveys(instrument = self.instruments, 
                                            target = star) 
        return search


    def _searchAndDownload(self, star, instrument, downloadPath=None, date=None):
        starARCFILE = np.array(self.searchStar(star, instrument)['ARCFILE'])
        if downloadPath:
            self.eso.retrieve_data(datasets = starARCFILE, 
                                   destination = downloadPath)
            return 0
        else:
            self.eso.retrieve_data(datasets = starARCFILE)
        return 0



    def getHARPSdata(self, star, downloadPath = None , date = None):
        """
        Download HARPS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data younger than a certain date
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstruments(star)
        esoInst = np.array(['HARPS'])
        for _, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0