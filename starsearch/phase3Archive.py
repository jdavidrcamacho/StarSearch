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
        #user name
        self.user = user
        #login to ESO
        self.eso = Eso()
        self.eso.login(self.user, store_password=store_password)
        #unlimited number of search results = -1
        self.eso.ROW_LIMIT = -1
        #list of available surveys
        self.surveys = self.eso.list_surveys()
        #default instruments for our package
        self.instruments = np.array(['FEROS', 'HARPS', 'ESPRESSO'])
        #In the future we might include UVES
        self.UVES = np.array(['UVES'])
        
        
    def searchStar(self, star, instrument = None, date = None, SNR = None):
        """
        Return phase 3 ESO query for selected star. 
        Includes other options such as instrument, date, and signal-to-noise.
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses our default instruments
        date: str
            Date to search for obvervation ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 1
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if date is None: 
            date = Time('1990-01-23')
        else:
            date = Time(date)
        if SNR is None: 
            SNR = 1

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
    
    
    def searchInstruments(self, star, phase3 = True):
        """
        Eithers checks which number of reduced spectra is available per 
        intrument (phase3=True) or checks which number of raw observations are
        in the ESO ARCHIVE (phase3=False).
        Most likelly the results will be different.
        
        Parameters
        ----------
        star: str
            Name of the star
        phase3: bool
            True searches instruments in the phase3 archive, False searches 
            instruments in the raw archive
            Default: phase3 = True
        Returns
        -------
        instrumentDict: dict
            Instruments and number of observations
        """
        if phase3 is True:
            searchResult = self.eso.query_surveys(target = star)
        else:
            searchResult = self.eso.query_main(column_filters={'target': star})
        instruments = np.unique(np.array(searchResult['Instrument']), 
                                return_counts=True)
        instrumentDict = dict()
        for i, j in enumerate(instruments[0]):
            instrumentDict[j] = instruments[1][i]
        return instrumentDict
    
    
    def searchInstrumentSpectra(self, star, instrument = None):
        """
        Checks how many reduced spectra is available to downloand for a given 
        instrument
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses our default instruments
            
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
        
        
    def searchObservationDate(self, star, instrument = None):
        """
        Searches for the modified Julian Date (JD - 2400000.5) of the 
        observation
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses our default instruments
            
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
        Searches for spectra of the selected star past the given date of 
        observation
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses our default instruments
        date: str
            Date to search for obvervations ('YYYY-MM-DD')
            If None: date = '1990-01-23'
            
        Returns
        -------
        search: table
            Result of the query on ESO arquive
        """
        if date is None: 
            date = Time('1990-01-23')
        else:
            date = Time(date)
        search = self.searchStar(star, instrument)
        search.remove_rows(Time(search['Date Obs']) < date)
        return search
    
    
    def searchBySNR(self, star, instrument = None, SNR = None):
        """
        Searches for spectra of the selected star with a signal-to-noise ratio
        higher than a given SNR value
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses our default instruments
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 1
            
        Returns
        -------
        SNR: array
            Array with the signal-to-noise ratios
        """
        if SNR is None: 
            SNR = 1
        search = self.searchStar(star, instrument)
        search.remove_rows(search['SNR (spectra)'] < SNR)
        return search
    
    
    def _searchAndDownload(self, star, instrument, downloadPath, date, SNR):
        """
        Auxiliary function to search and downlaod spectra
        
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
        starARCFILE = np.array(self.searchStar(star, instrument, 
                                               date, SNR)['ARCFILE'])
        if downloadPath:
            self.eso.retrieve_data(datasets = starARCFILE, 
                                   destination = downloadPath)
        else:
            self.eso.retrieve_data(datasets = starARCFILE)
            
            
    def _getData(self, star, instrument, downloadPath, date = None, SNR = None):
        """
        Auxiliary function to download spectra given an instrument and star 
        
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
        for _, j in enumerate(np.array([instrument])):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                self._searchAndDownload(star, j, downloadPath, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def getALLdata(self, star, downloadPath = None , date = None, SNR = None):
        """
        Downloads HARPS, FEROS, and ESPRESSO reduced spectra of a selected star
        
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
        if downloadPath is None:
            downloadPath = '~/'
        checkInstruments = self.searchInstruments(star)
        for _, j in enumerate(self.instruments):
            print('\n*** Searching for {0} results ***\n'.format(j))
            if j in checkInstruments:
                downloadPathInst = '{0}/{1}'.format(downloadPath, j)
                if not os.path.exists(downloadPathInst):
                    os.makedirs(downloadPathInst)
                self._searchAndDownload(star, j, downloadPathInst, date, SNR)
            else:
                print('No {0} data\n'.format(j))
        print('\n*** Done ***\n')
        
        
    def GetInstrumentData(self, star, instrument, downloadPath = None , 
                          date = None, SNR = None):
        """
        Download [INSTRUMENT] spectra of a selected star
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument; Must be one of "ESPRESSO, UVES, HARPS, FEROS"; 
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
        Downloaded [INSTRUMENT] spectra
        """
        if downloadPath is None:
            downloadPath = '~/'
        checkInstruments = self.searchInstruments(star)
        print('\n*** Searching for {0} results ***\n'.format(instrument))
        if instrument.upper() in checkInstruments:
            self._searchAndDownload(star, instrument, downloadPath, date, SNR)
        else:
            print('No {0} data for {1}\n'.format(instrument, star))
        print('\n*** Done ***\n')
        
        
    def getFILEdata(self, filename, header = 0, column = 0,
                    downloadPath = None, date = None, SNR = None):
        """
        Search and downloads spectra from FEROS, HARPS, amd ESPRESSO from a
        given file with stars
        WARNING
    
        Parameters
        ----------
        filename : str
            File name 
        header: int
            Number of header lines to skip
            Default: header = 0
        column: int
            Number of the column to read
            Default: column = 0
        downloadPath: str
            Adress where to download data
            If None: downloadPath = '~/'
        date: str
            Download only the data past than a certain date ('YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio.
            If None: SNR = 1
            
        Returns
        -------
        Downloaded spectra
        """
        if downloadPath is None:
            downloadPath = '~/'
        stars = np.loadtxt(filename, dtype = str, delimiter='\t', 
                           usecols = [column], skiprows = header)
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
                    saveFile = False, savePath = '', 
                    printFiles = False, fromList = False):
        """
        Return a summary of the reduced spectra available of the selected star
        
        Parameters
        ----------
        star: str
            Name of the star
        instrument: str
            Name of the instrument
            If None: Uses default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 1
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
        Prints a bunch of information of the available star's spectra
        """
        now = datetime.now()
        if saveFile:
            f = open(os.path.join(savePath, "{0}_{1}.txt".format(star,
                     now.strftime("%Y-%m-%dT%H:%M:%S"))),"a")
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
            print('Total number of spectra found: {0:8.1f}'.format(snr.size), file = f)
            value, count = np.unique(spectrograph, return_counts=True)
            for i in range(value.size):
                specSNR = search[search['Instrument']==value[i]]['SNR (spectra)']
                quadSum = np.sqrt(np.sum(specSNR**2))
                print('\n{0} spectra: {1}'.format(value[i], count[i]), '|',
                      'SNR Quadradratic Sum: {0:8.1f}'.format(quadSum), file = f)
                comparison = search[search['Instrument']==value[i]]['SNR (spectra)']
                maxSNRpos = np.argmax(comparison)
                print('Maximum SNR: {0:8.1f}'.format(comparison[maxSNRpos]), file=f)
                minSNRpos = np.argmin(comparison)
                print('Minimum SNR: {0:8.1f}'.format(comparison[minSNRpos]), file=f)
            #spectra found
            if printFiles:
                print('ARCFILE\tInstrument\tObservationDate\tSNR', file = f)
                for i, j in enumerate(fileName):
                    print('{0}\t{1}\t{2}\t{3}'.format(j, spectrograph[i], 
                          obsDate[i], snr[i]), file = f)
            if saveFile:
                f.close()
        ## This expection is not correct.
        except:
            print('{0} not found in archive\n'.format(star), file = f)
            
            
    def summaryList(self, starList, instrument = None, date = None, SNR = None, 
                    saveFile = False, savePath = '', dec = None,
                    printFiles = False):
        """
        Return a summary of the spectra available of the stars in a given list
        or numpy array
        
        Parameters
        ----------
        starList: array
            List or numpy array of stars to look at
        header: int
            Number of header lines to skip
            Default: header = 0
        column: int
            Number of the column to read
            Default: column = 0
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 1
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
        stars = np.array(starList) # in case the input is a list
        if not dec:
            dec = 180
        now = datetime.now()
        if saveFile:
            storage_name = "summary_{0}.txt".format(now.strftime("%Y-%m-%dT%H:%M:%S"))
            f = open(os.path.join(savePath, storage_name), "a")
        else:
            f = stdout
        noSpectra = [] #to add the stars with no spectra on archive
        for j in stars:
            try:
                self.summaryStar(j, instrument=instrument, date=date, SNR=SNR, 
                                 fromList=f);
            except:
                noSpectra.append(j)
                print('{0} not found in archive\n'.format(j), file = f)
        if saveFile:
            with open(os.path.join(savePath,
                            "{0}_noSpectra.txt".format(starList[0:-4])), "a"):
                for j in noSpectra:
                    try:
                        jSearch = Simbad.query_object(j)
                        RAandDEC = HMS2deg(jSearch['RA'][0], jSearch['DEC'][0])
                        if float(RAandDEC[1]) > dec:
                            pass
                        else:
                            print('{0}\t{1}degress'.format(j, RAandDEC[1]))
                    except:
                        with open(os.path.join(savePath,
                            "{0}_checkStar.txt".format(starList[0:-4])), "a"):
                            print('{0} not found on SIMBAD'.format(j))
        return noSpectra
        
    
    def summaryFile(self, filename, header = 0, column = 0, instrument = None, 
                    date = None, SNR = None, dec = None, saveFile = False, 
                    savePath = None, printFiles = False):
        """
        Return a summary of the spectra available of the stars in a given file
        
        Parameters
        ----------
        filename: str
            File, filename, or generator to read. 
            Example: filename = '~/something.txt'
        header: int
            Number of header lines to skip
            Default: header = 0
        column: int
            Number of the column to read
            Default: column = 0
        instrument: str
            Name of the instrument, None for default instruments
        date: str
            Date to search for obvervation (FORMAT: 'YYYY-MM-DD')
            If None: date = '1990-01-23'
        SNR: float
            Signal to noise ratio. 
            If None: SNR = 1
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
        stars = np.loadtxt(filename, skiprows = header, usecols=(column),
                           dtype = str, delimiter = '\t')
        noSpectra = self.summaryList(np.array(stars), instrument = instrument, 
                                     date = date, SNR = SNR, saveFile = saveFile, 
                                     savePath = savePath, dec = dec, 
                                     printFiles = printFiles)
        return noSpectra
    
    
############################ SWEET-cat functions ##############################
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
    
    
    def searchSWEETCatSpectra(self, filename, table, mag=12, dec = 30, 
                              noHigher = True, savePath = None, download = False):
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
        dec: float
            Search for stars with declination lower that dec
            Default: dec = 30
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
        stars: arra    y
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
        #declinations
        SWEETmags= np.loadtxt(table, usecols=(3), delimiter='\t', dtype=np.str)

        if savePath is None:
            os.mkdir("spectra")  #create this folder to store the data
            savePath = 'spectra/'

        print()
        spectra = np.loadtxt(filename, skiprows=2, usecols=(0), 
                             delimiter='\t', unpack=True, dtype=np.str)
        spectraInst = []    
        print()
        rv, sn, sn2 = np.loadtxt(filename, skiprows=2, usecols=(1,2,3), 
                                 delimiter='\t', unpack=True)
        t = datetime.now()

        build_name = lambda name, t: "{name}_{1}.txt".format(name, t.strftime("%Y-%m-%dT%H:%M:%S"))

        spectra_name  = os.path.join(savePath, build_name("01_starsWithSpectra",t))

        spectra_anom  = os.path.join(savePath, build_name("02_spectraAnomalies",t))
        no_spectra  = os.path.join(savePath, build_name("03_noSpectraFound",t))
        to_download  = os.path.join(savePath, build_name("04_starsToDownload",t))

        stars = []
        starSN, starSN2, starQuadSum, starName = [], [], [], []
        stars2download = [] #to compare SN, SN2 and quadSum

        with open(spectra_name, mode = 'a') as f1, open(spectra_anom, mode = 'a') as f2, open(no_spectra, mode = 'a') as f3, open(to_download, mode = 'a') as f4: 
            print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f1)
            print('-------\t--\t--\t---\t------\t---', file = f1)
            print('spectra\tRV\tSN\tSN2\tquadSN\tmag', file = f2)
            print('-------\t--\t--\t---\t------\t---', file = f2)
            print('star\tmag', file = f3)
            print('----\t---', file = f3)
            print('spectra\tinstr\tquadSN', file = f4)
            print('-------\t-----\t------', file = f4)

            for i, j in enumerate(spectra):
                nameSpliting = np.array(j.split('_'))
                spectra[i] = nameSpliting[0] #to fix the stars names
                spectraInst.append(nameSpliting[1]) #instrument
                #search star on SIMBAD
                starSearch = Simbad.query_object(spectra[i])
                try:
                    RAandDEC = HMS2deg(starSearch['RA'][0], starSearch['DEC'][0])
                    starDEC = RAandDEC[1]
                except:
                    starMag = '--'
                    print('{0}\t{1}\t{2}\t{3}\t--\t'.format(spectra[i], 
                        rv[i], sn[i], sn2[i]), starMag, file = f2)
                    pass
                if float(starDEC) > dec:
                    pass #star with declination higher than what we want
                else:
                    #we now look at the magnitudes on Simbad
                    try:
                        starMag = starSearch['FLUX_V'][0]
                    #or on SWEET-cat
                    except:
                        try:
                            position= np.where( SWEETstars == spectra[i])
                            starMag = np.round(float(SWEETmags[position]), 2)
                        except ValueError:
                            #Not found on SIMBAD and on SWETT-cat
                            starMag = '--'
                            print('{0}\t{1}\t{2}\t{3}\t--\t'.format(spectra[i], 
                                  rv[i], sn[i], sn2[i]), starMag, file = f2)
                            starMag = np.nan
                    #having magnitudes now we compare spectra
                    if np.isfinite(starMag) and starMag <= mag:
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
                                    print('{0}\t{1}\t{2}'.format(spectra[i], value[k], quadSum), file=f4)
                                if (quadSum>(sn[i] and sn2[i]) and (sn[i]<300 or sn2[i]<300)):
                                        stars2download.append(spectra[i])
                                        print('{0}\t{1}\t{2}'.format(spectra[i], value[k], quadSum), file=f4)
                            starName.append(spectra[i]); starQuadSum.append(quadSum)
                            starSN.append(sn[i]); starSN2.append(sn2[i]); print()
                        except TypeError:
                            print('{0}\t{1}'.format(spectra[i], starMag), file=f3)
        if download:
            for i, j in enumerate(stars2download):
                    self._downloadSWEETCatSpectra(j, savePath)
        #to remove duplicates
        stars, stars2download = np.unique(stars), np.unique(stars2download)
        return stars, stars2download
    
    
    def downloadSWEETCatSpectra(self, starsFile, savePath = None):
        """
        Function to download the spectra 04_starsToDownload file generated by
        searchSWEETCatSpectra() function
        
        Parameters
        ----------
        starsFile: str
            Path to 04_starsToDownload file
        savePath: string
            Address where to save the spectra
            
        Returns
        -------
        A bunch of folders (1 per star) with the spectra
        """

        if savePath is None:
            os.mkdir('spectra')
            savePath = 'spectra/'
        downloadPath = savePath
        stars = np.loadtxt(starsFile, skiprows = 2, usecols = (0), 
                           delimiter = '\t', dtype = np.str)
        inst = np.loadtxt(starsFile, skiprows = 2, 
                          usecols = (1), delimiter = '\t', dtype = np.str)
        for i, j in enumerate(stars):
            print('*************')
            print('*', j)
            print('*************')
            try:
                downloadPath4Star = os.path.join(savePath, '{0}/{1}'.format(downloadPath, j))
                if not os.path.exists(downloadPath4Star):
                    os.makedirs(downloadPath4Star)
                self._getData(j, inst[i], downloadPath=downloadPath4Star)
            except:
                print('Star not found in ESO archive!\n')
        return 0
    
    
    def _downloadSWEETCatSpectra(self, stars, savePath):
        """    
        Function to download the spectra found by searchSWEETCatDatabase(), to
        be used if download = True
        
        Parameters
        ----------
        stars: str
            Star's name we want to download
        savePath: string
            Address where to save the spectra
            
        Returns
        -------
        A bunch of folders (1 per star) with the spectra
        """
        downloadPath = savePath
        np.loadtxt
        print('*************')
        print('*', stars)
        print('*************')
        try:
            if not os.path.exists('{0}/{1}'.format(downloadPath, stars)):
                os.makedirs('{0}/{1}'.format(downloadPath, stars))
            downloadPath4Star = '{0}/{1}'.format(downloadPath, stars)
            self.getALLdata(stars, downloadPath4Star)
        except:
            print('Star not found in ESO archive!\n')
        return 0
    
    
### END
