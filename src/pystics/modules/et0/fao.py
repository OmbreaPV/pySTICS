import numpy as np

def deg2rad(degrees):
    """
    Convert angular degrees to radians

    :param degrees: Value in degrees to be converted.
    :return: Value in radians
    :rtype: float
    """
    return degrees * (np.pi / 180.0)

# Internal constants
# Latitude
_MINLAT_RADIANS = deg2rad(-90.0)
_MAXLAT_RADIANS = deg2rad(90.0)

# Solar declination
_MINSOLDEC_RADIANS = deg2rad(-23.5)
_MAXSOLDEC_RADIANS = deg2rad(23.5)

def _check_doy(doy):
    """
    Check day of the year is valid.
    """
    if not 1 <= doy <= 366:
        raise ValueError(
            'Day of the year (doy) must be in range 1-366: {0!r}'.format(doy))
    
def _check_latitude_rad(latitude):
    if not _MINLAT_RADIANS <= latitude <= _MAXLAT_RADIANS:
        raise ValueError(
            'latitude outside valid range {0!r} to {1!r} rad: {2!r}'
            .format(_MINLAT_RADIANS, _MAXLAT_RADIANS, latitude))
    

def _check_sol_dec_rad(sd):
    """
    Solar declination can vary between -23.5 and +23.5 degrees.

    See http://mypages.iit.edu/~maslanka/SolarGeo.pdf
    """
    if not _MINSOLDEC_RADIANS <= sd <= _MAXSOLDEC_RADIANS:
        raise ValueError(
            'solar declination outside valid range {0!r} to {1!r} rad: {2!r}'
            .format(_MINSOLDEC_RADIANS, _MAXSOLDEC_RADIANS, sd))



def sol_dec(day_of_year):
    """
    Calculate solar declination from day of the year.

    Based on FAO equation 24 in Allen et al (1998).

    :param day_of_year: Day of year integer between 1 and 365 or 366).
    :return: solar declination [radians]
    :rtype: float
    """
    _check_doy(day_of_year)
    return 0.409 * np.sin(((2.0 * np.pi / 365.0) * day_of_year - 1.39))


def sunset_hour_angle(latitude, sol_dec):
    """
    Calculate sunset hour angle (*Ws*) from latitude and solar
    declination.

    Based on FAO equation 25 in Allen et al (1998).

    :param latitude: Latitude [radians]. Note: *latitude* should be negative
        if it in the southern hemisphere, positive if in the northern
        hemisphere.
    :param sol_dec: Solar declination [radians]. Can be calculated using
        ``sol_dec()``.
    :return: Sunset hour angle [radians].
    :rtype: float
    """
    _check_latitude_rad(latitude)
    _check_sol_dec_rad(sol_dec)

    cos_sha = -np.tan(latitude) * np.tan(sol_dec)
    # If tmp is >= 1 there is no sunset, i.e. 24 hours of daylight
    # If tmp is <= 1 there is no sunrise, i.e. 24 hours of darkness
    # See http://www.itacanet.org/the-sun-as-a-source-of-energy/
    # part-3-calculating-solar-angles/
    # Domain of acos is -1 <= x <= 1 radians (this is not mentioned in FAO-56!)
    return np.acos(min(max(cos_sha, -1.0), 1.0))


def inv_rel_dist_earth_sun(day_of_year):
    """
    Calculate the inverse relative distance between earth and sun from
    day of the year.

    Based on FAO equation 23 in Allen et al (1998).

    :param day_of_year: Day of the year [1 to 366]
    :return: Inverse relative distance between earth and the sun
    :rtype: float
    """
    _check_doy(day_of_year)
    return 1 + (0.033 * np.cos((2.0 * np.pi / 365.0) * day_of_year))