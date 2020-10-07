For SNe, it would be nice not to have very long gaps with no g-band observations

So, let's add some basis functions that will supress taking g in consecutive nights, and boost taking go if there has been a long gap.

typical spot gets 86 g observations. 208 in i and 207 in r

so 8 and 20 per season. 

if season is 2.5 hours. --> 38 days. that seems short. 
so g cadence of 4.75 days
r,i of 1.9 days

so then I would say supress g if it's been observed in the last 3 days. Enhance after 5 days. Should we do a ramp up and down to prevent step function badness?


Thinking more--it is going to be tough to make this behave where it goes after areas that haven't been observed in g for a long time, and also not slam to high airmass.  

Might need to make a survey object that looks to see if a blob of observations filling a gap should be done. Could actially have it look ahead on the night, then broadcast when it will want to observe so other things respect it!

-------------
Notes on the survey object

* we need to include the moon mask (even if there are lots of HEALpix around the moon, we can't get them)
* There's also the issue of the zenith mask. That needs to be dealt with

So we need something a bit more sophisticated. Maybe the feature of last observed in g, the moon mask. Then grow a blob. If the blob is big enough, figure out when it would be best observed in the night (near meridian, small hour angle, moon down if possible), and schedule it.