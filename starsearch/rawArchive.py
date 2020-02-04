#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.eso import Eso
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
        self.user = str(user)
        self.eso = Eso()
        self.eso.login(self.user) #login in eso 
        self.eso.ROW_LIMIT = -1 #unlimited number of search results
        self.instruments = np.array(['FEROS', 'UVES', 'HARPS', 'ESPRESSO'])

        
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
            
            
        MJD-OBS is the modified Julian Date (JD - 2400000.5) of the start of the observation.
        """
        search = self.eso.query_main(column_filters={'target': star})
        result = np.array(search['Release_Date'])
        return result
    
    
    def searchObservationDate(self, star):
        """
        Checks the modified Julian Date (JD - 2400000.5) of the start of the 
        observation
        
        Parameters
        ----------
        star: str
            Name of the star
            
        Returns
        -------
        result: array
            Array with the start date of the observations
            
        """
        search = self.eso.query_main(column_filters={'target': star})
        result = np.array(search['MJD-OBS'])
        return result
    
    
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
    
    
    def _searchAndDownload(self, star, instrument, downloadPath, date, calib):
        """
        Download ESO spectra of a given star; to be used is getStarData()
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Instrument we are searching the darta
        date: float
            Download spectra younger than date (in modified Julian Date)
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
            return 0
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
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(self.instruments):
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
        date: float
            Download spectra younger than date (in modified Julian Date)
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstruments(star)
        esoInst = np.array(['FEROS'])
        for _, j in enumerate(esoInst):
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
        checkInstruments = self.searchInstruments(star)
        esoInst = np.array(['UVES'])
        for _, j in enumerate(esoInst):
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
        checkInstruments = self.searchInstruments(star)
        esoInst = np.array(['HARPS'])
        for _, j in enumerate(esoInst):
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
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past a certain date
        calib : str
            Retrieve associated calibration files: 'None' (default), 'raw' for
            raw calibrations, or 'processed' for processed calibrations.
            
        Returns
        -------
        
        """
        checkInstruments = self.searchInstruments(star)
        esoInst = np.array(['ESPRESSO'])
        for _, j in enumerate(esoInst):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, calib)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        return 0
