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
        self.eso.login(self.user) #login in eso 
        self.eso.ROW_LIMIT = -1 #unlimited number of search results
        
        
    def searchReleaseDate(self, star):
        """
        Checks the date the data was released 
        
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
    
    
    def _searchAndDownload(self, star, instrument, downloadPath, date, calib):
        """
        Download ESO spectra of a given star; to be used is getStarData()
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Instrument we are searching the darta
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'none' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        print('Downloading {0} data\n'.format(instrument))
        table = self.eso.query_main(column_filters = {'instrument': instrument, 
                                                      'target': star})
        if downloadPath:
            self.eso.retrieve_data(table['Dataset ID'], 
                                   destination = downloadPath,
                                   with_calib = calib)
        else:
            self.eso.retrieve_data(table['Dataset ID'], with_calib = calib)
        return 0
        
        
    def getALLdata(self, star, downloadPath = None , date = None, calib = 'none'):
        """
        Download ESO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstrument(star)
        esoInst = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0
    
    
    def getFEROSdata(self, star, downloadPath = None , date = None, 
                     calib = 'none'):
        """
        Download FEROS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstrument(star)
        esoInst = np.array(['FEROS'])
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0
    
    
    def getUVESdata(self, star, downloadPath = None , date = None, 
                    calib = 'none'):
        """
        Download UVES spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstrument(star)
        esoInst = np.array(['UVES'])
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0
    
    
    def getHARPSdata(self, star, downloadPath = None , date = None, 
                    calib = 'none'):
        """
        Download HARPS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstrument(star)
        esoInst = np.array(['HARPS'])
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0
    
    
    def getESPRESSOdata(self, star, downloadPath = None , date = None, 
                        calib = 'none'):
        """
        Download ESPRESSO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstrument(star)
        esoInst = np.array(['ESPRESSO'])
        for i, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0