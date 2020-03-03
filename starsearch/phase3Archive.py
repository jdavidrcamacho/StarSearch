#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
from astropy.time import Time
from sys import maxsize, stdout
from datetime import datetime
np.set_printoptions(threshold = maxsize)
from starsearch.core import Eso


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
        super(ESOquery, self).__init__()
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
    
    
    def searchInstrumentSpectra(self, star, instrument=None):
        """
        Checks how many spectra is available to downloand for a given instrument
        
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
        
        
    def searchStar(self, star, instrument=None, date=None, SNR=None):
        """
        Return phase 3 ESO query for given star and a given instrument
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if not date: 
            date = Time('1990-01-23')
        date = Time(date)
        if not SNR: 
            SNR = 30
        if instrument:
            search = self.eso.query_surveys(surveys = instrument, 
                                            target = star)
        else:
            search = self.eso.query_surveys(surveys = list(self.instruments), 
                                            target = star)
        search.remove_rows(Time(search['Date Obs']) < date) #Date criteria
        search.remove_rows(search['SNR (spectra)'] < SNR) #SNR critetia
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
    
    
    def searchByDate(self, star, instrument=None, date=None):
        """
        Searches for spectra of a given star past a certain date of observation
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if not date: 
            date = Time('1990-01-23')
        date = Time(date)
        search = self.searchStar(star, instrument)
        search.remove_rows(Time(search['Date Obs']) < date)
        return search
    
    
    def searchBySNR(self, star, instrument=None, SNR = None):
        """
        Return spectra with a signal-to-noise ratio higher than a certain 
        SNR value
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        SNR: array
            Array with the signal-to-noise ratios
        """
        if not SNR: 
            SNR = 30
        search = self.searchStar(star, instrument)
        search.remove_rows(search['SNR (spectra)'] < SNR)
        return search
    
    
    def _searchAndDownload(self, star, instrument = None, 
                            downloadPath = None, date = None, SNR = None):
        starARCFILE = np.array(self.searchStar(star, instrument, 
                                               date, SNR)['ARCFILE'])
        if downloadPath:
            self.eso.retrieve_data(datasets = starARCFILE, 
                                   destination = downloadPath)
        else:
            self.eso.retrieve_data(datasets = starARCFILE)
            
            
    def _getData(self, star, instrument, downloadPath = None , date = None,
                 SNR = None):
        """
        Downloads spectra from a given instrument of a given star 
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data younger than a certain date
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array([instrument])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getALLdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Download ESO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(self.instruments):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getFEROSdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Download FEROS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data oast than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array(['FEROS'])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getUVESdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Download UVES spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array(['UVES'])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getHARPSdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Download HARPS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array(['HARPS'])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getESPRESSOdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Download ESPRESSO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 30
            
        Returns
        -------
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array(['ESPRESSO'])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def readFILE(self, filename, instrument = None, downloadPath = None, 
                 date = None, SNR = None):
        """
        Load SWEET-Cat stars catalogue. Basically a copy of np.loadtxt() but 
        thinking in using it on SWEET-Cat.
        Each row in the text file must have the same number of columns.
    
        Parameters
        ----------
        filename : str
            File, filename, or generator to read. 
            filename = '~/WEBSITE_online.rdb'
        instrument: str
            Name of the instrument, None for default instruments
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 30
            
        Returns
        -------
        starsInArchive : array
            List of stars found on ESO archive
        starsNotInArchive: array
            List of stars not found on ESO archive
        """
        stars = np.loadtxt(filename, dtype=str, delimiter='\t', 
                           usecols=[0], skiprows=0)
        starsInArchive = np.array([])
        starsNotInArchive = np.array([])
        for _, j in enumerate(stars):
            print('*************')
            print('*', j)
            print('*************')
            try:
                self.searchStar(j, instrument, date, SNR)
                starsInArchive = np.append(starsInArchive, j)
                print('Star found in ESO archive.\n')
            except:
                starsNotInArchive = np.append(starsNotInArchive, j)
                print('Star not found in ESO archive!\n')
        return starsInArchive, starsNotInArchive
    
    
    def getFILEdata(self, filename, downloadPath = None, 
                    date = None, SNR = None):
        """
        Search and downloads spectra of SWEET-Cat stars catalogue.
    
        Parameters
        ----------
        filename : str
            File, filename, or generator to read. 
            filename = '~/WEBSITE_online.rdb'
        downloadPath: str
            Adress where to download data
            If None: downloadPath = '~/'
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 30
            
        Returns
        -------
        """
        if downloadPath is None:
            downloadPath = '~/'
        savePath = downloadPath
        stars = np.loadtxt(filename, dtype=str, delimiter='\t', 
                           usecols=[0], skiprows=0)
        for _, j in enumerate(stars):
            print('*************')
            print('*', j)
            print('*************')
            try:
                if not os.path.exists('{0}/{1}'.format(downloadPath, j)):
                    os.makedirs('{0}/{1}'.format(downloadPath, j))
                downloadPath = '{0}/{1}'.format(downloadPath, j)
                self.getALLdata(j, downloadPath, date, SNR)
            except:
                print('Star not found in ESO archive!\n')
            downloadPath = savePath
    
    
    def _remove_planet(self, name):
        """
        Remove the trailing b, c, d, etc in the stellar name, no sure if it is
        going to be necessary
        
        Parameters
        ----------
        name: str
            Name of the star+planet
            
        Returns
        -------
        name: str
            Name  of the star
        """
        planets = 'bcdefghijB'
        for planet in planets:
            if name.endswith(' %s' % planet):
                return name[:-2]
        # some exoplanets have .01 or .02 in the name 
        if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
            return name[:-3]                
        return name
    
    
    def summaryStar(self, star, instrument=None, date=None, SNR=None, 
                    savetofile = False):
        """
        Return a summary of the spectra available of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
            
        Returns
        -------
        """
        if not date: 
            date = Time('1990-01-23')
        date = Time(date)
        if not SNR: 
            SNR = 30
        if instrument:
            search = self.eso.query_surveys(surveys = instrument, 
                                            target = star)
        else:
            search = self.eso.query_surveys(surveys = list(self.instruments), 
                                            target = star)
        search.remove_rows(Time(search['Date Obs']) < date) #Date criteria
        search.remove_rows(search['SNR (spectra)'] < SNR) #SNR critetia
        #organizing our stuff
        fileName = np.array(search['ARCFILE'])
        spectrograph = np.array(search['Instrument'])
        obsDate = np.array(search['Date Obs'])
        snr = np.array(search['SNR (spectra)'])
        now = datetime.now()
        if savetofile:
            f = open("{0}_{1}.txt".format(star,now.strftime("%Y-%m-%dT%H:%M:%S")),"a")
        else:
            f = stdout 
        #number of spectra
        print('*** SUMMARY of {0} *** \n'.format(star), file = f)
        print('Total number of spectra found: {0}'.format(snr.size), file = f)
        value, count = np.unique(spectrograph, return_counts=True)
        for i in range(value.size):
            print('{0} spectra: {1}'.format(value[i], count[i]), file = f)
        print(file = f)
        #signal-to-noise
        maxSNRpos = np.argmax(snr)
        print('Maximum SNR: {0}'.format(snr[maxSNRpos]), file = f)
        print('Observation date: {0}'.format(obsDate[maxSNRpos]), file = f)
        print('Instrument: {0}'.format(spectrograph[maxSNRpos]), file = f)
        print('File name: {0} \n'.format(fileName[maxSNRpos]), file = f)
        minSNRpos = np.argmin(snr)
        print('Minimum SNR: {0}'.format(snr[minSNRpos]), file = f)
        print('Observation date: {0}'.format(obsDate[minSNRpos]), file = f)
        print('Instrument: {0}'.format(spectrograph[minSNRpos]), file = f)
        print('File name: {0} \n'.format(fileName[minSNRpos]), file = f)
        #spectra found
        print('ARCFILE\tInstrument\tObservationDate\tSNR', file = f)
        for i, j in enumerate(fileName):
            print('{0}\t{1}\t{2}\t{3}'.format(j,spectrograph[i],obsDate[i], snr[i]),
                file = f)
        if savetofile:
            f.close()
    
    
    def summaryList(self, filename, instrument=None, date=None, SNR=None,
                    savetofile = False):
        """
        Return a summary of the spectra available of the stars in a given list
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 30
        savetofile: bool
            Save summary in a .txt file
            Default: False
        Returns
        -------
        """   
        stars = np.loadtxt(filename, dtype=str, delimiter='\t', 
                           usecols=[0], skiprows=0)
        now = datetime.now()
        if savetofile:
            f = open("{0}_{1}.txt".format(filename,
                    now.strftime("%Y-%m-%dT%H:%M:%S")),"a")
        else:
            f = stdout 
        for _, j in enumerate(stars):
            try:
                self.summaryStar(j, instrument=instrument, date=date, SNR=SNR,
                                 savetofile=savetofile)
            except:
                print('{0} not found in archive'.format(j), file = f)
        if savetofile:
            f.close()

