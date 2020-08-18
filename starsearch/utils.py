def split_str(str):
    """ 
    Splits strings
    Taken from: https://stackoverflow.com/a/57109577
    
    Parameters
    ----------
    str
        String to separate (duh!)
        
    Returns
    -------
    res_arr: list
        List of new strings
    """
    ref_dict = {'\x07':'a', '\x08':'b', '\x0C':'f', '\n':'n', '\r':'r', 
                '\t':'t', '\x0b':'v'}
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
    #Some exoplanets have .01 or .02 in the name 
    if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
        return name[:-3]
    return name


def HMS2deg(ra='', dec=''):
    """
    To convert hours, minutes and seconds into decimal degrees
    Taken from: http://www.bdnyc.org/2012/10/decimal-deg-to-hms/
    
    Parameters
    ----------
    ra: str
        Right Ascension in hours, minutes, seconds
    dec: str
        Declination in hours, minutes, seconds
    
    Returns
    -------
    RA: str
        Right Ascencion in degrees
    DEC: str
        Declination in degrees
    """
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC
