#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def split_str(str):
    """ 
    Splits strings, taken from: https://stackoverflow.com/a/57109577
    
    Parameters
    ----------
    str
        String to separate (duh!)
        
    Returns
    -------
    res_arr: list
        List of new strings
    """
    ref_dict = {
        '\x07':'a',
        '\x08':'b',
        '\x0C':'f',
        '\n':'n',
        '\r':'r',
        '\t':'t',
        '\x0b':'v',
    }
    res_arr = []
    temp = ''
    for i in str :
        if not i == '\\':
            if i in ref_dict:
                if not temp == "":
                    res_arr.append(temp)
                res_arr.append(ref_dict[i])
                temp = ''
            else:    
                temp += i
        else:
            if not temp == '':
                res_arr.append(temp)
            temp = ''
    res_arr.append(temp)
    return res_arr


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
    
    