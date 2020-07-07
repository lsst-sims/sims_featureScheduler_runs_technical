import numpy as np
import matplotlib.pylab as plt
from astropy.coordinates import get_sun, SkyCoord, EarthLocation, AltAz, get_moon
from astropy.time import Time
from lsst.sims.utils import Site, raDec2Hpid, m5_flat_sed
from astropy import units as u
from lsst.sims.skybrightness_pre import SkyModelPre
import healpy as hp
from lsst.sims.seeingModel import SeeingModel
from lsst.sims.downtimeModel import ScheduledDowntimeData
from lsst.sims.featureScheduler.surveys import generate_dd_surveys



def place_obs(nights, space=2):
    """Take some input nights and flag nights as observe or not
    """
    # say 1 for an observation that can float 1.7 days, 2 for one that's good for 12 hours?
    result = np.zeros(nights.size, dtype=int) 
    # Place observations before and after each large gap
    gaps = nights*0
    gaps[1:] = nights[1:] - nights[0:-1]
    gaps[0] = space
    after_gap = np.where(gaps >= space)[0]
    result[after_gap] = 1
    pre_gap = after_gap[1:]-1
    result[pre_gap] = 2
    last_obs_indx = 0

    for i, night in enumerate(nights):
        if result[i] > 0:
            last_obs_indx = i
        elif (night-nights[last_obs_indx]) >= space:
            result[i] = 1
            last_obs_indx = i
    return result

def generate_ddf(ddf_name, nyears=10, space=2):
    previous_ddf = generate_dd_surveys()
    survey_names = np.array([survey.survey_name for survey in previous_ddf])
    survey_indx = np.where(survey_names == ddf_name)[0].max()
    ddf_ra = previous_ddf[survey_indx].ra * u.rad
    ddf_dec = previous_ddf[survey_indx].dec* u.rad

    site = Site('LSST')
    location = EarthLocation(lat=site.latitude, lon=site.longitude, height=site.height)

    mjd = np.arange(59853.5, 59853.5+365.25*nyears, 20./60/24.)
    times = Time(mjd, format='mjd', location=location)

    airmass_limit = 2.5  # demand airmass lower than this
    twilight_limit = -18.  # Sun below this altitude in degrees
    dist_to_moon_limit = 30.  # minimum distance to keep from moon degrees
    zenith_limit = 10.  # Need to be this far away from zenith to start (20 min = 5 deg)
    g_m5_limit = 23.5  # mags

    season_gap = 20.  # days. Count any gap longer than this as it's own season
    season_length_limit = 80  # Days. Demand at least this many days in a season

    # How long to keep attempting a DDF
    expire_dict = {1: 36./24., 2: 0.5}

    sun_coords = get_sun(times)
    moon_coords = get_moon(times)

    sched_downtime_data = ScheduledDowntimeData(Time(mjd[0], format='mjd'))
    observatory_up = np.ones(mjd.size, dtype=bool)
    for dt in sched_downtime_data():
        indx = np.where((mjd >= dt['start'].mjd) & (mjd <= dt['end'].mjd))[0]
        observatory_up[indx] = False

    lst = times.sidereal_time('mean')

    sun_altaz = sun_coords.transform_to(AltAz(location=location))

    # generate a night label for each timestep
    sun_rise = np.where((sun_altaz.alt[0:-1] < 0) & (sun_altaz.alt[1:] > 0))[0]
    night = np.zeros(mjd.size, dtype=int)
    night[sun_rise] = 1
    night = np.cumsum(night) + 1  # 1-index for night

    sun_down = np.where(sun_altaz.alt < twilight_limit*u.deg)[0]

    ddf_coord = SkyCoord(ra=ddf_ra, dec=ddf_dec)
    ddf_altaz = ddf_coord.transform_to(AltAz(location=location, obstime=times))
    ddf_airmass = 1./np.cos(np.radians(90.-ddf_altaz.az.deg))
    zenith = AltAz(alt=90.*u.deg, az=0.*u.deg)
    ddf_zenth_dist = zenith.separation(ddf_altaz)

    nside = 32
    ddf_indx = raDec2Hpid(nside,ddf_coord.ra.deg, ddf_coord.dec.deg)
    sm = SkyModelPre()

    g_sb = mjd*0+np.nan

    indices = np.where((sun_altaz.alt < twilight_limit*u.deg) & (ddf_airmass > airmass_limit))[0]
    # In theory, one could reach into the sky brightness model and do a much faster interpolation
    # There might be an airmass limit on the sky brightness.
    for indx in sun_down:
        g_sb[indx] = sm.returnMags(mjd[indx], indx=[ddf_indx], filters='g', badval=np.nan)['g']

    dist_to_moon = ddf_coord.separation(moon_coords)
    seeing_model = SeeingModel()
    ddf_approx_fwhmEff = seeing_model(0.7, ddf_airmass)
    # I think this should pluck out the g-filter. Really should be labled
    ddf_approx_fwhmEff = ddf_approx_fwhmEff['fwhmEff'][1].ravel()

    ddf_m5 = m5_flat_sed('g', g_sb, ddf_approx_fwhmEff, 30., ddf_airmass, nexp=1.)

    # demand sun down past twilight, ddf is up, and observatory is open, and not too close to the moon
    good = np.where((ddf_airmass < airmass_limit) & (sun_altaz.alt < twilight_limit*u.deg) &
                    (ddf_airmass > 0) & (observatory_up == True) & (dist_to_moon > dist_to_moon_limit*u.deg) &
                    (ddf_zenth_dist > zenith_limit*u.deg) &
                    (ddf_m5 > g_m5_limit))

    potential_nights = np.unique(night[good])
    night_gap = potential_nights[1:] - potential_nights[0:-1]
    big_gap = np.where(night_gap > season_gap)[0]+1
    season = potential_nights*0
    season[big_gap] = 1
    season = np.cumsum(season)

    u_seasons = np.unique(season)
    season_lengths = []
    for se in u_seasons:
        in_se = np.where(season == se)
        season_lengths.append(np.max(potential_nights[in_se]) - np.min(potential_nights[in_se]))
    season_lengths = np.array(season_lengths)

    good_seasons = u_seasons[np.where(season_lengths > season_length_limit)[0]]
    gn = np.isin(season, good_seasons)
    potential_nights = potential_nights[gn]
    season = season[gn]

    obs_attempts = []
    for sea in np.unique(season):
        night_indx = np.where(season == sea)
        obs_attempts.append(place_obs(potential_nights[night_indx], space=space))
    obs_attempts = np.concatenate(obs_attempts)

    mjd_observe = []
    m5_approx = []
    for indx in np.where(obs_attempts > 0)[0]:
        in_night_indx = np.where(night == potential_nights[indx])[0]
        best_depth_indx = np.min(np.where(ddf_m5[in_night_indx] == np.nanmax(ddf_m5[in_night_indx]))[0])
        mjd_start = mjd[in_night_indx[best_depth_indx]]
        m5_approx.append(ddf_m5[in_night_indx[best_depth_indx]])
        mjd_end = mjd_start + expire_dict[obs_attempts[indx]]
        mjd_observe.append((mjd_start, mjd_end))

    result = np.zeros(len(mjd_observe), dtype=[('mjd_start', '<f8'),
                      ('mjd_end', '<f8'), ('label', '<U10')])
    mjd_observe = np.array(mjd_observe)
    result['mjd_start'] = mjd_observe[:, 0]
    result['mjd_end'] = mjd_observe[:, 1]
    result['label'] = ddf_name

    return result #, ddf_ra, ddf_dec #, previous_ddf[survey_indx].observations, m5_approx


if __name__ =="__main__":

    roll = False
    if roll:
        ddf_names = ['DD:XMM-LSS', 'DD:ECDFS', 'DD:COSMOS', 'DD:EDFS']
        filename = 'test_sched.npz'
    else:
        ddf_names = ['DD:ELAISS1', 'DD:XMM-LSS', 'DD:ECDFS', 'DD:COSMOS', 'DD:EDFS']
        filename = 'test_sched_noroll.npz'

    results = []
    for ddf_name in ddf_names:
        print('generating %s' % ddf_name)
        mjd = generate_ddf(ddf_name, nyears=10.0)
        results.append(mjd)

    if roll:
        ddf_names = ['DD:ELAISS1']
        for ddf_name in ddf_names:
            print('generating %s' % ddf_name)
            mjd = generate_ddf(ddf_name, nyears=1.0, space=0.8)
            results.append((mjd, ddf_name))

    final = np.concatenate(results)

    np.savez(filename, results=final)

