

With branch tickets/OPSIM-771:

python baseline.py --survey_length 30 --verbose
progress = 100.05%Skipped 0 observations
Flushed 0 observations from queue for being stale
Completed 21053 observations
ran in 8 min = 0.1 hours
Writing results to  baseline_nexp2_v1.7_0yrs.db

python baseline.py --survey_length 90 --verbose
progress = 100.06%Skipped 0 observations
Flushed 0 observations from queue for being stale
Completed 55120 observations
ran in 31 min = 0.5 hours
Writing results to  baseline_nexp2_v1.7_0yrs.db

with master:

python baseline.py --survey_length 30 --verbose
progress = 100.05%Skipped 0 observations
Flushed 0 observations from queue for being stale
Completed 21048 observations
ran in 6 min = 0.1 hours
Writing results to  baseline_nexp2_v1.7_0yrs.db

python baseline.py --survey_length 90 --verbose
progress = 100.06%Skipped 0 observations
Flushed 0 observations from queue for being stale
Completed 55125 observations
ran in 17 min = 0.3 hours
Writing results to  baseline_nexp2_v1.7_0yrs.db


-----------------------------------------

Running on inga.astro.washington.edu, for full 10 years

on master (old skybrightness_pre files)
python baseline.py --verbose
progress = 100.01%Skipped 0 observations
Flushed 47 observations from queue for being stale
Completed 2044491 observations
ran in 470 min = 7.8 hours
Writing results to  baseline_nexp2_v1.7_10yrs.db


