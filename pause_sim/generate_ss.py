import glob

dbfiles = glob.glob('*.db')
dbfiles.sort()

pops = ['granvik_5k', 'l7_5k']

orbit_files = {'granvik_5k': '/gscratch/astro/lynnej/orbits/granvik_5k/granvik_5k.txt',
               'l7_5k': '/gscratch/astro/lynnej/orbits/l7_5k/l7_5k.txt'}

ranges = {'granvik_5k': (16, 22.2, 0.2),
          'l7_5k': (4, 12, 0.2)}

metadatas = {'granvik_5k': 'NEO',
             'l7_5k': 'TNO'}

hmarks = {'granvik_5k': 22,
          'l7_5k': 8}

runs = [file.replace('.db', '') for file in dbfiles]

runs = [run for run in runs if 'tracking' not in run]

with open('ss_script.sh', 'w') as f:
    for pop in pops:
        for run, filename in zip(runs, dbfiles):
            s1 = ('makeLSSTobs.py --opsimDb %s --orbitFile %s' % (filename, orbit_files[pop]))
            s2 = ('run_moving_calc.py --characterization inner --obsFile ' +
                  '%s__%s_obs.txt ' % (run, pop) +
                  '--opsimDb %s.db ' % run + '--orbitFile '+orbit_files[pop] +
                  ' --outDir %s_%s ' % (run, pop) +
                  '--opsimRun %s ' % run +
                  '--hMin %0.2f --hMax %0.2f --hStep %0.2f ' % ranges[pop] +
                  '--metadata %s ' % metadatas[pop])
            s3 = ('run_moving_fractions.py --workDir %s_%s ' % (run, pop) +
                  '--metadata %s ' % metadatas[pop] +
                  '--hMark %i' % hmarks[pop])
            print(s1 + ' ; ' + s2 + ' ; ' + s3, file=f)
