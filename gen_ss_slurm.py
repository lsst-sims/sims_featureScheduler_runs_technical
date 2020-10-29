import glob 
import os


if __name__ == "__main__":
    db_files = glob.glob('*10yrs.db')
    db_files.sort()
    orbit_path = '/gscratch/astro/lynnej/orbits/'

    orbit_files = {'granvik_5k': os.path.join(orbit_path, 'granvik_5k/granvik_5k.txt'),
                   'granvik_pha_5k': os.path.join(orbit_path, 'granvik_pha_5k/granvik_pha_5k.txt'),
                   'l7_5k': os.path.join(orbit_path, 'l7_5k/l7_5k.txt'),
                   'mba_5k': os.path.join(orbit_path, 'mba_5k/mba_5k.txt'),
                   'oort': os.path.join(orbit_path, 'oort/oort.txt'),
                   'sdo_5k': os.path.join(orbit_path, 'sdo_5k/sdo_5k.txt'),
                   'trojan_5k': os.path.join(orbit_path, 'trojan_5k/trojan_5k.txt'),
                   }

    pops = list(orbit_files.keys())

    defaultRanges = {'PHA': [16, 28, 0.2, 22],
                     'NEO': [16, 28, 0.2, 22],
                     'MBA': [16, 26, 0.2, 20],
                     'Trojan': [14, 22, 0.2, 18],
                     'TNO': [6, 12, 0.2, 8],
                     'SDO': [4, 12, 0.2, 7],
                     'Oort': [4, 20, 0.5, 5]}

    pop_labels = {'granvik_5k': 'NEO',
                  'granvik_pha_5k': 'PHA',
                  'l7_5k': 'TNO',
                  'mba_5k': 'MBA',
                  'oort': 'Oort',
                  'sdo_5k': 'SDO',
                  'trojan_5k': 'Trojan'}

    # Not bothering to re-type it
    hmarks = {}
    ranges = {}
    for pop in pops:
        hmarks[pop] = defaultRanges[pop_labels[pop]][-1]
        ranges[pop] = tuple(defaultRanges[pop_labels[pop]][0:3])

    inner_outer = {'PHA': 'inner', 'NEO': 'inner',
                   'MBA': 'inner', 'Trojan': 'inner',
                   'TNO': 'outer', 'SDO': 'outer',
                   'Oort': 'outer'}

    runs = [file.replace('.db', '') for file in db_files]
    runs = [run for run in runs if 'tracking' not in run]

    with open('ss1_script.sh', 'w') as f:
        for pop in pops:
            for filename in db_files:
                print('makeLSSTobs.py --opsimDb %s --orbitFile %s' % (filename, orbit_files[pop]), file=f)

    with open('ss2_script.sh', 'w') as f:
        for pop in pops:
            for run in runs:
                print('run_moving_calc.py --characterization %s --obsFile ' % inner_outer[pop_labels[pop]] +
                      '%s__%s_obs.txt ' % (run, pop) +
                      '--opsimDb %s.db ' % run + '--orbitFile %s ' % orbit_files[pop] +
                      '--outDir %s ' % run +
                      '--opsimRun %s ' % run +
                      '--hMin %0.2f --hMax %0.2f --hStep %0.2f ' % ranges[pop] +
                      '--metadata %s ' % pop_labels[pop],
                      file=f)

    with open('ss3_script.sh', 'w') as f:
        for run in runs:
            for pop in pops:
                print('run_moving_fractions.py --workDir %s ' % run +
                      '--metadata %s ' % pop_labels[pop] +
                      '--hMark %i' % hmarks[pop],
                      file=f)
