#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
#from astroquery.esocas import Eso
from astroquery.eso import Eso


class ESOquery(object):
    """
    Lets check what instruments observed a given star!!!
    """
    
    def __init__(self, user):
        self.user = str(user)
        self.eso = Eso()
        self.eso.login(self.user)

    def searchStarDates(self, star):
        """
        Lets see when the star was observed
        """
        searchResult = self.eso.query_main(column_filters={'target': star})
        return searchResult['Release_Date']
    
    