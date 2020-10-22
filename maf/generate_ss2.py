import glob

dbfiles = glob.glob('*.db')
dbfiles.sort()

pops = ['granvik_5k', 'l7_5k']

orbit_files = {'granvik_5k':'/gscratch/astro/lynnej/orbits/granvik_5k/granvik_5k.txt',
               'l7_5k': '/gscratch/astro/lynnej/orbits/l7_5k/l7_5k.txt'}

runs = [file.replace('.db', '') for file in dbfiles]

with open('ss2_script.sh', 'w') as f:
    for pop in pops:
        for run in runs:
            print('run_moving_calc.py --characterization inner --obsFile '+
                  '%s__%s_obs.txt ' % (run, pop) + 
                  '--opsimDb %s.db ' % run + '--orbitFile '+orbit_files[pop] + 
                  '--outDir %s_ss ' % run +
                  '--runName %s' % run,
                  file=f)