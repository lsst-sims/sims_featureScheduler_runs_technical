import numpy as np
from lsst.sims.featureScheduler.utils import empty_observation
from lsst.sims.featureScheduler.surveys import BaseSurvey
import copy
from lsst.sims.utils import _approx_RaDec2AltAz


def basic_sequence(ra, dec, survey_name='', sequence='urgizy', nvis=[8, 20, 10, 20, 26, 20],
                   exptime=30., u_exptime=30., nexp=1):
    """Generate the list of observations that should happen in a ddf sequence
    """
    observations = []
    for num, filtername in zip(nvis, sequence):
        for j in range(num):
            obs = empty_observation()
            obs['filter'] = filtername
            if filtername == 'u':
                obs['exptime'] = u_exptime
            else:
                obs['exptime'] = exptime
            obs['RA'] = np.radians(ra)
            obs['dec'] = np.radians(dec)
            obs['nexp'] = nexp
            obs['note'] = survey_name
            observations.append(obs)
    return np.array(observations)


def ddf_info(exptime=30., u_exptime=30., nexp=1):
    """
    Basic info for each DDF
    """
    ddf_dict = {}
    ha_limits_dict = {}

    # ELAIS S1
    RA = 9.45
    dec = -44.
    ddf_dict['DD:ELAISS1'] = basic_sequence(RA, dec, 'DD:ELAISS1',
                                            exptime=exptime,
                                            u_exptime=u_exptime, nexp=nexp)
    ha_limits_dict['DD:ELAISS1'] = ([0., 1.5], [21.5, 24.])

    # XMM-LSS
    RA = 35.708333
    dec = -4-45/60.
    ddf_dict['DD:XMM-LSS'] = basic_sequence(RA, dec, 'DD:XMM-LSS',
                                            exptime=exptime,
                                            u_exptime=u_exptime, nexp=nexp)
    ha_limits_dict['DD:XMM-LSS'] = ([0., 1.5], [21.5, 24.])

    # Extended Chandra Deep Field South
    RA = 53.125
    dec = -28.-6/60.
    ddf_dict['DD:ECDFS'] = basic_sequence(RA, dec, 'DD:ECDFS',
                                          exptime=exptime,
                                          u_exptime=u_exptime, nexp=nexp)
    ha_limits_dict['DD:ECDFS'] = ([0.5, 3.0], [20., 22.5])

    # COSMOS
    RA = 150.1
    dec = 2.+10./60.+55/3600.
    ddf_dict['DD:COSMOS'] = basic_sequence(RA, dec, 'DD:COSMOS', exptime=exptime, u_exptime=u_exptime, nexp=nexp)
    ha_limits_dict['DD:COSMOS'] = ([0., 2.5], [21.5, 24.])

    # Euclid Fields
    # I can use the sequence kwarg to do two positions per sequence
    filters = 'urgizy'
    nviss = [8, 20, 10, 20, 26, 20]
    survey_name = 'DD:EDFS'
    # Note the sequences need to be in radians since they are using observation objects directly
    RAs = np.radians([58.97, 63.6])
    decs = np.radians([-49.28, -47.60])
    sequence = []

    for filtername, nvis in zip(filters, nviss):
        for ra, dec in zip(RAs, decs):
            for num in range(nvis):
                obs = empty_observation()
                obs['filter'] = filtername
                if filtername == 'u':
                    obs['exptime'] = u_exptime
                else:
                    obs['exptime'] = exptime
                obs['RA'] = ra
                obs['dec'] = dec
                obs['nexp'] = nexp
                obs['note'] = survey_name
                sequence.append(obs)
    ddf_dict[survey_name] = np.array(sequence)
    ha_limits_dict[survey_name] = ([0., 1.5], [22.5, 24.])

    return ddf_dict, ha_limits_dict


def read_times(filename='schedule_59002.txt', end_time=36.):
    end_time = end_time/24.  # to Days
    names = ['mjd_start', 'label']
    types = [float, '|U10']
    data = np.genfromtxt(filename, dtype=list(zip(names, types)), skip_header=1)

    info, ha_dict = ddf_info()
    ddf_names = list(info.keys())
    # convert the names
    orig_names = np.unique(data['label'])
    for name in orig_names:
        indx = np.where(data['label'] == name)
        for new_name in ddf_names:
            if name.lower() in new_name.lower():
                data['label'][indx] = new_name

    order = np.argsort(data['mjd_start'])
    data = data[order]
    mjd_end = np.empty(data.size, dtype=[('mjd_end', float)])
    mjd_end = data['mjd_start'] + end_time

    names = ['mjd_start', 'mjd_end', 'label']
    types = [float, float, '|U10']
    # I can never figure out how to stack numpy arrays
    result = np.empty(data.size, dtype=list(zip(names, types)))
    result['mjd_start'] = data['mjd_start']
    result['label'] = data['label']
    result['mjd_end'] = mjd_end
    return result


class Scheduled_ddfs(BaseSurvey):
    def __init__(self, times_array, sequence_dict, ha_dict, basis_functions=[],
                 detailers=[], ignore_obs=None, reward_value=100, flush_time=40.):

        super(Scheduled_ddfs, self).__init__(basis_functions=basis_functions,
                                             detailers=detailers, ignore_obs=ignore_obs)

        self.times_array = times_array
        self.sequence_dict = sequence_dict
        self.flush_time = flush_time/60./24.
        self.observation_complete = np.zeros(self.times_array.size, dtype=bool)
        self.reward_value = reward_value
        self.ha_dict = ha_dict

    def _check_feasibility(self, conditions):
        result = False
        # Check if there is a sequence that wants to go
        # XXX--can probably use searchsorted here for faster results
        in_mjd = np.where((conditions.mjd > self.times_array['mjd_start']) &
                          (conditions.mjd < self.times_array['mjd_end']) &
                          (self.observation_complete == False))[0]
        for indx in in_mjd:
            ra = self.sequence_dict[self.times_array[indx]['label']]['RA'][0]
            # HA
            target_HA = (conditions.lmst - ra*12/np.pi) % 24
            in_HA_range = False
            for ha_range in self.ha_dict[self.times_array[indx]['label']]:
                lres = np.min(ha_range) <= target_HA < np.max(ha_range)
                in_HA_range = in_HA_range | lres

            if in_HA_range:
                result = True
                self.observations = copy.deepcopy(self.sequence_dict[self.times_array[indx]['label']])
                self.indx = indx
                return result

        return result

    def calc_reward_function(self, conditions):
        result = -np.inf
        if self._check_feasibility(conditions):
            result = self.reward_value
        return result

    def generate_observations_rough(self, conditions):
        result = []
        if self._check_feasibility(conditions):
            result = self.observations

            # Set the flush_by
            result['flush_by_mjd'] = conditions.mjd + self.flush_time

            # remove filters that are not mounted
            mask = np.isin(result['filter'], conditions.mounted_filters)
            result = result[mask]
            # Put current loaded filter first
            ind1 = np.where(result['filter'] == conditions.current_filter)[0]
            ind2 = np.where(result['filter'] != conditions.current_filter)[0]
            result = result[ind1.tolist() + (ind2.tolist())]

            # convert to list of array. Arglebargle, don't understand why I need a reshape there
            final_result = [row.reshape(1,) for row in result]
            result = final_result

            # Assume if we've gotten here, the observations are on the way to being executed,
            # Mark them as complete. Could get fancier, but let's see if simple works
            # XXX--might need to have a way to mark complete if we restart?
            self.observation_complete[self.indx] = True

        return result
