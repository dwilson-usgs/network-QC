# network-QC
Python codes for seismic network QC

These codes are the Python version of Matlab code found here:
https://github.com/dwilson-usgs/temp-network-QC

There have also been improvements and additions to the Matlab version.

---------------------------------------------------------
* NetQCamps.py
```
usage: NetQCamps.py [-h] -nets NETS [-stas STAS] [-chans CHANS]
                    [-minlat MINLAT] [-maxlat MAXLAT] [-minlon MINLON]
                    [-maxlon MAXLON] [-n NET] [-s STA] [-r RAD] [-t1 TIME]
                    [-t2 TIME2]
```
Check station amplitudes relative to other local stations

optional arguments:
```
   -h, --help      show this help message and exit
   -nets NETS      Networks for analysis
   -stas STAS      Stations for analysis
   -chans CHANS    Channels for analysis
   -minlat MINLAT  Min Latitude for QC
   -maxlat MAXLAT  Max Latitude for QC
   -minlon MINLON  Min Longitude for QC
   -maxlon MAXLON  Max Longitude for QC
   -n NET          Option to specify net individual station. If your use this flag,
                        you must also use the -s flag.
   -s STA          Option to specify individual station
   -r RAD          Option to specify radius (in degrees) from individual station when using
                        the -n and -s flags (default is 3 degrees)
   -t1 TIME        Start Time for analysis (default is 24hrs ago)
   -t2 TIME2       End Time for analysis (default is now)
```
Examples:

    To check amplitudes of newly installed aftershock stations:
    python3 NetQCamps.py -nets GS -t1 2016-10-01 -t2 2016-10-02 -minlat 36 -maxlat 37 -minlon -97.5 -maxlon -96.2

    To check amplitudes of stations in the vicinity of station N4-G62A over the last 24 hours:
    python3 NetQCamps.py -nets N4,US,NE -n N4 -s G62A

---------------------------------------------------------

**Disclaimer:**

>This software is preliminary or provisional and is subject to revision. It is 
being provided to meet the need for timely best science. The software has not 
received final approval by the U.S. Geological Survey (USGS). No warranty, 
expressed or implied, is made by the USGS or the U.S. Government as to the 
functionality of the software and related material nor shall the fact of release 
constitute any such warranty. The software is provided on the condition that 
neither the USGS nor the U.S. Government shall be held liable for any damages 
resulting from the authorized or unauthorized use of the software.

---------------------------------------------------------
