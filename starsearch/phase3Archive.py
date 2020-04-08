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
Simbad.add_votable_fields('flux(V)')

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
        Checks how many spectra is available to downloand for a given 
        instrument
        
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
        
        
    def searchStar(self, star, instrument = None, date = None, SNR = None):
        """
        Return phase 3 ESO query for given star and a given instrument
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 0
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if not date: 
            date = Time('1990-01-23')
        date = Time(date)
        if not SNR: 
            SNR = 0
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
            Date to search for obvervations ('YYYY-MM-DD')
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
        Downloaded spectra
        """
        checkInstruments = self.searchInstruments(star)
        print('check', checkInstruments)
        for _, j in enumerate(np.array([instrument])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getALLdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Downloads ESO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPatch: str
            Adress where to download data
        date: str
            Download only the data past a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 10
            
        Returns
        -------
        Downloaded spectra
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(self.instruments):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getFEROSdata(self, star, downloadPath = None , date = None, 
                     SNR = None):
        """
        Download FEROS spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data oast than a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 10
            
        Returns
        -------
        Downloaded FEROS spectra
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
            Download only the data past than a certain date (YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 10
            
        Returns
        -------
        Downloaded UVES spectra
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
            Download only the data past than a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 10
            
        Returns
        -------
        Downloaded HARPS spectra
        """
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(np.array(['HARPS'])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getESPRESSOdata(self, star, downloadPath = None , date = None, 
                        SNR = None):
        """
        Download ESPRESSO spectra of a given star
        
        Parameters
        ----------
        star: str
            Name of the star
        downloadPath: str
            Adress where to download data
        date: str
            Download only the data past than a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 10
            
        Returns
        -------
        Downloaded ESPRESSO data
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
            Download only the data past than a certain date ('YYYY-MM-DD')
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
        Search and downloads spectra for a file with stars
    
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
            Download only the data past than a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 10
            
        Returns
        -------
        Downloaded spectra
        """
        if downloadPath is None:
            downloadPath = '~/'
        stars = np.loadtxt(filename, dtype = str, delimiter='\t', 
                           usecols = [0], skiprows = header)
        for _, j in enumerate(stars):
            print('*************')
            print('*', j)
            print('*************')
            try:
                if not os.path.exists('{0}/{1}'.format(downloadPath, j)):
                    os.makedirs('{0}/{1}'.format(downloadPath, j))
                downloadPath4Star = '{0}/{1}'.format(downloadPath, j)
                self.getALLdata(j, downloadPath4Star, date, SNR)
            except:
                print('Star not found in ESO archive!\n')
    
    
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
        Print a bunch of information of the available star's spectra
        """
        now = datetime.now()
        if savePath:
            os.chdir(savePath)
        if saveFile:
            f = open("{0}_{1}.txt".format(star,
                     now.strftime("%Y-%m-%dT%H:%M:%S")),"a")
        else:
            f = stdout 
        if fromList:
            f = fromList
        print('*** SUMMARY of {0} ***'.format(star), file = f)
        search = self.searchStar(star, instrument = instrument, date = date,
                                 SNR = SNR)
        #organizing our stuff
        try:
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
        except:
            print('{0} not found in archive\n'.format(star), file = f)
            
            
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
            If None: dec = 180
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
        if savePath:
            os.chdir(savePath)
        stars = np.loadtxt(filename, skiprows = header, usecols=(0),
                           dtype = str, delimiter = '\t')
        if not dec:
            dec = 180
        for i, j in enumerate(stars):
            checkLash = j.find('/') #to separate the stars that have two names
            stars[i] = ' '.join(j.split(' ', -1)[:-2]) #remove last 2 columns
            if checkLash != -1:
                newStar = j.split('/')
                stars[i] = newStar[0]
                stars = np.append(stars, ' '.join(newStar[1].split(' ', 
                                                  -1)[:-2]))
        now = datetime.now()
        if saveFile:
            f = open("summary_{0}.txt".format(now.strftime("%Y-%m-%dT%H:%M:%S")),
                     "a")
        else:
            f = stdout
        noSpectra = [] #to add the stars with no spectra on archive
        for _, j in enumerate(stars):
            try:
                print(j)
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
                    print('{0} not found on SIMBAD'.format(j),
                          file = f2)
            f1.close(); f2.close()
        return noSpectra
    
    
    def searchSWEETCatSpectra(self, filename, table, mag=12, noHigher = True, 
                               savePath = None, download = False):
        """
        Compares the SWEET-Cat spectra with the one available on the ESO
        archive
        
        Parameters
        ----------
        filename: str
            File with the SWEET-Cat spectra info
        table: str
            SWEET-cat table downloaded from the website (table.txt)
        mag: float
            Maximum magnitude of the stars we are searching
            Default: mag = 12
        savePath: str
            Path to save the generated .txt files
        download: bool
            True if we want to download the spectra
            Default: download = False
            
        Returns
        -------
        01_spectra: file
            List of stars of ESO archive and respective signal-to-noise ratio
        02_spectraTBC: file
            List of stars to check manually (TBC = To Be Checked) 
            Either has a different name or no magnitude found on SIMBAD
        03_noSpectraFound: file
            List of stars with either no spectra on ESO archive or with 
            magnitude outside our specifications
        stars: array
            Array with the star from 01_spectra file
        stars2download: array
            Array with the stars to be downloaded
        """
        #name of the stars
        SWEETstars= np.loadtxt(table, usecols=(0), delimiter='\t', dtype=np.str)
        for i, j in enumerate(SWEETstars):
            #to remove empty spaces
            SWEETstars[i] = j.replace(" ", "")
        #magnitudes
        SWEETmags= np.loadtxt(table, usecols=(4), delimiter='\t', dtype=np.str)
        if savePath:
            os.chdir(savePath)
        else:
            savePath = 'spectra/'
        print()
        spectra = np.loadtxt(filename, skiprows=2, usecols=(0), 
                             delimiter='\t', unpack=True, dtype=np.str)
        spectraInst = []
        print()
        rv, sn, sn2 = np.loadtxt(filename, skiprows=2, usecols=(1,2,3), 
                                 delimiter='\t', unpack=True)
        t = datetime.now()
        f1 = open("01_spectra_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f1)
        print('-------\t--\t--\t---\t------\t---', file = f1)
        f2 = open("02_spectraTBC_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f2)
        print('-------\t--\t--\t---\t------\t---', file = f2)
        f3 = open("03_noSpectraFound_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('star\tmag', file = f3)
        print('----\t---', file = f3)
        stars = []
        starSN, starSN2, starQuadSum, starName = [], [], [], []
        for i, j in enumerate(spectra):
            nameSpliting = np.array(j.split('_'))
            spectra[i] = nameSpliting[0] #to fix the stars names
            spectraInst.append(nameSpliting[1]) #instrument
            try:
                #search magnitudes in SIMBAD
                starMag = Simbad.query_object(spectra[i])['FLUX_V'][0] #<- ERROR HERE
                
                if np.float(starMag) is np.nan:
                    try:
                        position= np.where( SWEETstars == spectra[i])
                        starMag = np.round(float(SWEETmags[position]), 2)
                    except ValueError:
                        starMag = '--'
                        print('spectra i', spectra[i])
                        print('{0}\t{1}\t{2}\t{3}\t--\t'.format(spectra[i], 
                              rv[i], sn[i], sn2[i]), starMag, file = f2)
            except:
                try:
                    position= np.where( SWEETstars == spectra[i])
                    starMag = np.round(float(SWEETmags[position]), 2)
                except:
                    starMag = '--'
                    print('{0}\t{1}\t{2}\t{3}\t--\t'.format(spectra[i], 
                          rv[i], sn[i], sn2[i]), starMag, file = f2)
            if type(starMag) is (np.float32 or np.float64) and starMag <= mag:
                try:
                    search = self.searchStar(spectra[i], SNR = 1)
                    spectrograph = np.array(search['Instrument'])
                    #number of spectra
                    value, count = np.unique(spectrograph, return_counts=True)
                    for k in range(value.size):
                        specSNR = search[search['Instrument']==value[k]]['SNR (spectra)']
                        quadSum = np.round(np.sqrt(np.sum(specSNR**2)), 2)
                        stars.append(spectra[i])
                        if value[k] == spectraInst[i]:
                            print('{0}_{1}\t{2}\t{3}\t{4}\t{5}\t'.format(spectra[i], 
                                  value[k], rv[i], sn[i], sn2[i], quadSum),
                                    starMag, file = f1)
                        else:
                            print('{0}_{1}\t--\t--\t--\t{2}\t'.format(spectra[i], 
                                  value[k], quadSum),
                                    starMag, file = f1)
                    starName.append(spectra[i]); starQuadSum.append(quadSum)
                    starSN.append(sn[i]); starSN2.append(sn2[i]); print()
                except TypeError:
                    print('{0}\t{1}'.format(spectra[i], starMag), file=f3)
        
        f1.close(); f2.close(); f3.close()
        stars2download = [] # compare SN, SN2 and quadSum
        for i, j in enumerate(starQuadSum):
            if (j>(starSN[i] and starSN2[i]) and (starSN[i]<300 or starSN2[i]<300)):
                stars2download.append(starName[i])
                if download:
                    self._downloadSWEETCatSpectra(starName[i], savePath)
        #to remove duplicates
        stars, stars2download = np.unique(stars), np.unique(stars2download)
        return stars, stars2download
        
        
    def searchSWEETCatDatabase(self, filename, mag, noHigher = True, 
                               savePath = None, download = False):
        """
        Compares the SWEET-Cat spectra with the one available on the ESO
        archive
        
        Parameters
        ----------
        filename: str
            File with the SWEET-Cat spectra info
        mag: float
            Magnitude of the stars 
        noHigher: bool
            True if we want to search stars with a magnitude lower than mag, 
            False is we want to search stars with magnitude higher than mag.
            Default: noHigher = True
        savePath: str
            Path to save the generated .txt files
        download: bool
            True if we want to download the spectra
            Default: download = False
            
        Returns
        -------
        01_spectra: file
            List of stars of ESO archive and respective signal-to-noise ratio
        02_spectraTBC: file
            List of stars to check manually (TBC = To Be Checked) 
            Either has a different name or no magnitude found on SIMBAD
        03_noSpectra: file
            List of stars with either no spectra on ESO archive or with 
            magnitude outside our specifications
        stars: array
            Array with the star from 01_spectra file
        stars2download: array
            Array with the stars to be downloaded
        """
        if savePath:
            os.chdir(savePath)
        else:
            savePath = 'spectra/'
        spectra = np.loadtxt(filename, skiprows=2, usecols=(0), 
                             delimiter='\t', unpack=True, dtype=np.str)
        rv, sn, sn2 = np.loadtxt(filename, skiprows=2, usecols=(1,2,3), 
                                 delimiter='\t', unpack=True)
        t = datetime.now()
        f1 = open("01_spectra_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f1)
        print('-------\t--\t--\t---\t------\t---', file = f1)
        f2 = open("02_spectraTBC_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f2)
        print('-------\t--\t--\t---\t------\t---', file = f2)
        f3 = open("03_noSpectra_{0}.txt".format(t.strftime("%Y-%m-%dT%H:%M:%S")),
                  "a")
        print('star\tmag', file = f3)
        print('----\t---', file = f3)
        stars = []
        starSN, starSN2, starQuadSum, starName = [], [], [], []
        for i, j in enumerate(spectra):
            spectra[i] = j.split('_')[0] #to fix the stars names
            if noHigher:
                try: 
                    starMag = Simbad.query_object(spectra[i])['FLUX_V'][0]
                    starMag = np.round(starMag, 2)
                    if starMag <= mag:
                        search = self.searchStar(spectra[i], SNR = 1)
                        spectrograph = np.array(search['Instrument'])
                        #number of spectra
                        value, count = np.unique(spectrograph, return_counts=True)
                        for k in range(value.size):
                            specSNR = search[search['Instrument']==value[k]]['SNR (spectra)']
                            quadSum = np.round(np.sqrt(np.sum(specSNR**2)), 2)
                            print('{0}spectra:'.format(spectra[i]))
                            print('{0} spectra: {1}'.format(value[k], count[k]), 
                                  '|', 'SNR Quadratic Sum: {0}'.format(quadSum))
                            stars.append(spectra[i])
                            print('{0}_{1}\t{2}\t{3}\t{4}\t{5}\t'.format(spectra[i], 
                                  value[k], rv[i], sn[i], sn2[i], quadSum),
                                    starMag, file = f1)
                        starName.append(spectra[i]); starQuadSum.append(quadSum)
                        starSN.append(sn[i]); starSN2.append(sn2[i]); print()
                    else:
                        print('{0}\tmag ='.format(spectra[i]), starMag, file=f3)
                except:
                    try:
                        starMag = Simbad.query_object(spectra[i])['FLUX_V'][0]
                        starMag = np.round(starMag, 2)
                    except:
                        starMag = '--'
                    print('{0}\t{1}\t{2}\t{3}\tNo\t'.format(spectra[i], 
                           rv[i],sn[i],sn2[i]), starMag, file = f2)
            else:
                try: 
                    starMag = Simbad.query_object(spectra[i])['FLUX_V'][0]
                    starMag = np.round(starMag, 2)
                    if starMag >= mag:
                        search = self.searchStar(spectra[i], SNR = 1)
                        spectrograph = np.array(search['Instrument'])
                        #number of spectra
                        value, count = np.unique(spectrograph, return_counts=True)
                        for k in range(value.size):
                            specSNR = search[search['Instrument']==value[k]]['SNR (spectra)']
                            quadSum = np.round(np.sqrt(np.sum(specSNR**2)), 2)
                            print('{0} spectra:'.format(spectra[i]))
                            print('{0} spectra: {1}'.format(value[k], count[k]), 
                                  '|', 'SNR Quadratic Sum: {0}'.format(quadSum))
                            stars.append(spectra[i])
                            print('{0}_{1}\t{2}\t{3}\t{4}\t{5}\t'.format(spectra[i], 
                                  value[k], rv[i], sn[i], sn2[i], quadSum),
                                    starMag, file = f1)
                        starName.append(spectra[i]); starQuadSum.append(quadSum)
                        starSN.append(sn[i]); starSN2.append(sn2[i]); print()
                    else:
                        print('{0}\tmag ='.format(spectra[i]), starMag, file=f3)
                except:
                    try:
                        starMag = Simbad.query_object(spectra[i])['FLUX_V'][0]
                        starMag = np.round(starMag, 2)
                    except:
                        starMag = '--'
                    print('{0}\t{1}\t{2}\t{3}\tNo\t'.format(spectra[i], 
                           rv[i],sn[i],sn2[i]), starMag, file = f2)
        f1.close(); f2.close(); f3.close()
        stars2download = [] # compare SN, SN2 and quadSum
        for i, j in enumerate(starQuadSum):
            if j > (starSN[i] and starSN2[i]):
                stars2download.append(starName[i])
        #to remove duplicates
        stars, stars2download = np.unique(stars), np.unique(stars2download)
        if download:
            self._downloadSWEETCatSpectra(stars2download, savePath)
        return stars, stars2download
    
    
    def _downloadSWEETCatSpectra(self, stars, savePath):
        """
        Function to download the spectra found by searchSWEETCatDatabase()
        
        Parameters
        ----------
        stars: array
            Array with the stars to download
        savePath: string
            Address where to save the spectra
            
        Returns
        -------
        A bunch of folders (1 per star) with the spectra
        """
        downloadPath = savePath
        for _, j in enumerate(stars):
            print('*************')
            print('*', j)
            print('*************')
            try:
                if not os.path.exists('{0}/{1}'.format(downloadPath, j)):
                    os.makedirs('{0}/{1}'.format(downloadPath, j))
                downloadPath4Star = '{0}/{1}'.format(downloadPath, j)
                self.getALLdata(j, downloadPath4Star)
            except:
                print('Star not found in ESO archive!\n')
        return 0
