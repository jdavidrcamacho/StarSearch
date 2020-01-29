#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.eso import Eso
from sys import maxsize
np.set_printoptions(threshold = maxsize)

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
        self.user = str(user) #user name 
        self.eso = Eso()
        self.eso.login(self.user) #login to eso archive
        self.eso.ROW_LIMIT = -1 #unlimited number of search results
        self.surveys = self.eso.list_surveys() #list of surveys
        
    def searchSurveys(self, star):
        return 0