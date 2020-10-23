import glob

dbfiles = glob.glob('*.db')
dbfiles.sort()

pops = ['granvik_5k', 'l7_5k']

orbit_files = {'granvik_5k': '/gscratch/astro/lynnej/orbits/granvik_5k/granvik_5k.txt',
               'l7_5k': '/gscratch/astro/lynnej/orbits/l7_5k/l7_5k.txt'}

ranges = {'granvik_5k': (16, 18, 0.2),
          'l7_5k': (4, 12, 0.2)}

metadatas = {'granvik_5k': 'NEO',
             'l7_5k': 'TNO'}

hmarks = {'granvik_5k': 22,
          'l7_5k': 8}

runs = [file.replace('.db', '') for file in dbfiles]

runs = [run for run in runs if 'tracking' not in run]

with open('ss2_script.sh', 'w') as f:
    for pop in pops:
        for run in runs:
            print('run_moving_calc.py --characterization inner --obsFile '+
                  '%s__%s_obs.txt ' % (run, pop) + 
                  '--opsimDb %s.db ' % run + '--orbitFile '+orbit_files[pop] + 
                  ' --outDir %s_%s ' % (run, pop) +
                  '--opsimRun %s ' % run +
                  '--hMin %0.2f --hMax %0.2f --hStep %0.2f' % ranges[pop],
                  file=f)

with open('ss3_script.sh', 'w') as f:
    for run in runs:
        for pop in pops:
            print('run_moving_fractions.py --workDir %s_%s ' % (run, pop) + 
                  '--metadata %s ' % metadatas[pop] + 
                  '--hMark %i' % hmarks[pop],
                  file=f)
