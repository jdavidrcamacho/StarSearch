#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
from sys import maxsize, stdout
np.set_printoptions(threshold = maxsize)

from astropy.time import Time
from datetime import datetime
from starsearch.core import Eso
from starsearch.utils import HMS2deg
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

class ESOquery():
    """
    ESO query class
    
    Parameters
    ----------
    user: str
        User name used for ESO website login
    store_password : bool
        Optional, stores the password securely in your keyring
        Default: store_password = False
        
    Returns
    -------
    """
    def __init__(self, user, store_password=False):
        super(ESOquery, self).__init__()
        self.user = user #user name 
        self.eso = Eso()
        if store_password:
            self.eso.login(self.user, store_password=True) #login to eso
        else:
            self.eso.login(self.user)
        self.eso.ROW_LIMIT = -1 #unlimited number of search results = -1
        self.surveys = self.eso.list_surveys() #list of available surveys
        self.instruments = np.array(['FEROS', 'HARPS', 'ESPRESSO'])
        self.UVES = np.array(['UVES'])
        
        
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
    
    
    def searchInstrumentSpectra(self, star, instrument = None):
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
        
        
    def searchStar(self, star, instrument = None, date = None, SNR = None, 
                   dec = None):
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
            If None: SNR = 10
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if not date: 
            date = Time('1990-01-23')
        date = Time(date)
        if not SNR: 
            SNR = 10
        if not dec:
            dec = 30
        if instrument:
            search = self.eso.query_surveys(surveys=instrument, target=star)
        else:
            search = self.eso.query_surveys(surveys = list(self.instruments), 
                                            target = star)
        try:
            search.remove_rows(Time(search['Date Obs']) < date) #Date criteria
            search.remove_rows(search['SNR (spectra)'] < SNR) #SNR critetia
        except:
            pass
        return search
    
    
    def searchObservationDate(self, star, instrument = None):
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
    
    
    def searchByDate(self, star, instrument = None, date = None):
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
    
    
    def searchBySNR(self, star, instrument = None, SNR = None):
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
            If None: SNR = 10
            
        Returns
        -------
        SNR: array
            Array with the signal-to-noise ratios
        """
        if not SNR: 
            SNR = 10
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
            If None: SNR = 10
            
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
            If None: SNR = 10
            
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
            If None: SNR = 10
            
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
            If None: SNR = 10
            
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
            If None: SNR = 10
            
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
            If None: SNR = 10
            
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
    
    
    def getFILEdata(self, filename, header = 0, downloadPath = None, 
                    date = None, SNR = None):
        """
        Search and downloads spectra of SWEET-Cat stars catalogue.
    
        Parameters
        ----------
        filename : str
            File, filename, or generator to read. 
            Example: filename = '~/WEBSITE_online.rdb'
        header: int
            Number of header lines to skip
            Default: header = 0
        downloadPath: str
            Adress where to download data
            If None: downloadPath = '~/'
        date: str
            Download only the data past than a certain date (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 10
            
        Returns
        -------
        """
        if downloadPath is None:
            downloadPath = '~/'
        savePath = downloadPath
        stars = np.loadtxt(filename, dtype = str, delimiter='\t', 
                           usecols = [0], skiprows = header)
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
    
    
    def summaryStar(self, star, instrument = None, date = None, SNR = None, 
                    saveFile = False, savePath = None, 
                    printFiles = False, fromList = False):
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
            If None: SNR = 10
        saveFile: bool
            Save summary in a .txt file
            Default: saveFile = False
        savePath: str
            Path to save the generated .txt file
        printFiles: bool
            Print each spectra file and respective parameters
            Default: printFiles = False
        fromList: Bool
            If we are checking a list or not
        Returns
        -------
        """
        now = datetime.now()
        if savePath:
            os.chdir(savePath)
        if saveFile:
            f = open("{0}_{1}.txt".format(star,now.strftime("%Y-%m-%dT%H:%M:%S")),"a")
        else:
            f = stdout 
        if fromList:
            f = fromList
        print('*** SUMMARY of {0} ***'.format(star), file = f)
        search = self.searchStar(star, instrument = instrument, date = date,
                                 SNR = SNR)
        #organizing our stuff
#        try:
        fileName = np.array(search['ARCFILE'])
        spectrograph = np.array(search['Instrument'])
        obsDate = np.array(search['Date Obs'])
        snr = np.array(search['SNR (spectra)'])
        #number of spectra
        print('Total number of spectra found: {0}'.format(snr.size), file = f)
        value, count = np.unique(spectrograph, return_counts=True)
        for i in range(value.size):
            specSNR = search[search['Instrument']==value[i]]['SNR (spectra)']
            quadSum = np.sqrt(np.sum(specSNR**2))
            print('{0} spectra: {1}'.format(value[i], count[i]), '|',
                  'SNR Quadratic Sum: {0}'.format(quadSum), file = f)
        print(file = f)
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
        if printFiles:
            print('ARCFILE\tInstrument\tObservationDate\tSNR', file = f)
            for i, j in enumerate(fileName):
                print('{0}\t{1}\t{2}\t{3}'.format(j, spectrograph[i], 
                      obsDate[i], snr[i]), file = f)
        if saveFile:
            f.close()
#    except:
#        print('{0} not found in archive\n'.format(star), file = f)
            
        
    def summaryList(self, filename, header = 0, instrument = None, date = None, 
                    SNR = None, dec = None, saveFile = False, savePath = None, 
                    printFiles = False):
        """
        Return a summary of the spectra available of the stars in a given list
        
        Parameters
        ----------
        filename: str
            File, filename, or generator to read. 
            Example: filename = '~/something.txt'
        header: int
            Number of header lines to skip
            Default: header = 0
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 10
        dec: float
            Search for stars with declination lower that dec
            If None: dec = 30
        saveFile: bool
            Save summary in a .txt file
            Default: saveFile = False
        savePath: str
            Path to save the generated .txt files
        printFiles: bool
            Print each spectra file and respective parameters
            Default: printFiles = False
            
        Returns
        -------
        noSpectra: arr
            Array with the stars with no spectra on ESO archives
        """   
        stars = np.loadtxt(filename, dtype = str, delimiter = '\t', 
                           skiprows = header)
        if not dec:
            dec = 30
        for i, j in enumerate(stars):
            checkLash = j.find('/') #to separate the stars that have two names
            stars[i] = ' '.join(j.split(' ', -1)[:-2]) #remove last 2 columns
            if checkLash != -1:
                newStar = j.split('/')
                stars[i] = newStar[0]
                stars = np.append(stars, ' '.join(newStar[1].split(' ', -1)[:-2]))
        now = datetime.now()
        if saveFile:
            f = open("summary_{0}.txt".format(now.strftime("%Y-%m-%dT%H:%M:%S")), "a")
        else:
            f = stdout
        noSpectra = [] #to add the stars with no spectra on archive
        for _, j in enumerate(stars):
            try:
                self.summaryStar(j, instrument=instrument, date=date, SNR=SNR, 
                                 fromList=f);
            except:
                noSpectra.append(j)
                print('{0} not found in archive\n'.format(j), file = f)
        f.close()
        if saveFile:
            f1 = open("{0}_noSpectra.txt".format(filename[0:-4]), "a")
            f2 = open("{0}_checkStar.txt".format(filename[0:-4]), "a")
            for _, j in enumerate(noSpectra):
                try:
                    jSearch = Simbad.query_object(j)
                    RAandDEC = HMS2deg(jSearch['RA'][0], jSearch['DEC'][0])
                    if float(RAandDEC[1]) > dec:
                        pass
                    else:
                        print('{0}\t{1}degress'.format(j, RAandDEC[1]), 
                              file = f1)
                except:
#                    try:
#                    jSearch = Vizier.query_object(j)
#                    RAandDEC = HMS2deg(jSearch['RA'][0], jSearch['DEC'][0])
#                    if float(RAandDEC[1]) > dec:
#                        pass
#                    else:
#                        print('{0}\t{1} degress'.format(j, RAandDEC[1]), 
#                              file = f)
                    print('{0} not found on SIMBAD'.format(j),
                          file = f2)
            f1.close(); f2.close()
        return noSpectra


