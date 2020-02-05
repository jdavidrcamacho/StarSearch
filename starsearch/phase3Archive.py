#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.eso import Eso
from astroquery.esocas import Eso
from astropy.time import Time
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
        WARNING: The number of observations might be different from the number
        of spectra available to download in phase 3
        
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
    

    def searchNumberByInstrument(self, star, instrument=None):
        """
        Checks how many spectra is available to downloan for a given instrument
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
            
        Returns
        -------
        instrumentDict: dict
            Instruments and number of observations
        """
        searchStar = self.searchStar(star)
        instrumentDict = dict()
        if instrument:
            number = np.sum(searchStar['Instrument'] == instrument)
            instrumentDict[instrument] = number
            return instrumentDict
        else:
            for _, j in enumerate(self.instruments):
                number = np.sum(searchStar['Instrument'] == j)
                instrumentDict[j] = number
            return instrumentDict
    

    def searchStar(self, star, instrument=None):
        """
        Return phase 3 ESO query for given star and a given instrument

        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
            
        Returns
        -------
        instrumentDict: dict
            Instruments and number of observations
        """
        if instrument:
            search = self.eso.query_surveys(instrument = instrument, 
                                            target = star) 
        else:
            search = self.eso.query_surveys(instrument = self.instruments, 
                                            target = star) 
        return search

    def searchObservationDate(self, star, instrument=None):
        """
        Searches for the modified Julian Date (JD - 2400000.5) of the 
        observation
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
            
        Returns
        -------
        result: array
            Array with the start date of the observations
        """
        if instrument:
            search = self.eso.query_surveys(instrument = instrument, 
                                            target = star) 
        else:
            search = self.eso.query_surveys(instrument = self.instruments, 
                                            target = star) 
        result = Time(search['Date Obs'], format='isot', scale='utc').mjd
        return result


    def searchSNR(self, star, instrument=None):
        """
        Return the signal-to-noise ration on the available spectra from a given
        instrument on a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments

        Returns
        -------
        SNR: array
            Array with the signal-to-noise ratios
        """
        if instrument:
            search = self.eso.query_surveys(instrument = instrument, 
                                            target = star) 
        else:
            search = self.eso.query_surveys(instrument = self.instruments, 
                                            target = star) 
        SNR = np.array(search['SNR (spectra)'])
        return SNR


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
        ------
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