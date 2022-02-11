#!/usr/bin/env python

from obspy.clients.fdsn import Client
#from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as ml
from obspy.core import UTCDateTime
import datetime as datetime
#from scipy.interpolate import griddata, interp1d
from scipy.interpolate import griddata
from obspy.geodetics import gps2dist_azimuth
import matplotlib
from scipy import signal
#from obspy.taup import TauPyModel
#from scipy import signal as scisig
#from obspy.signal import spectral_estimation as spec
from matplotlib import cm
#import pickle
#import gc
#from numba import jit
from matplotlib.transforms import blended_transform_factory
import matplotlib.patches as mpatches
#from scipy.optimize import leastsq
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import cartopy
from cartopy import geodesic
import shapely
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import argparse
import copy
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

def process_stream(st):
    st.detrend()
    st.remove_response(output="ACC")
    st.detrend()
    st.filter('bandpass', freqmin=fmin, freqmax=fmax)
    winl=int(np.floor(len(st.data)/24))
    nn=[]
    for n in range(23):
        nn.append(np.std(st.data[n*winl:((n+2)*winl)-1]))
    nn=np.asarray(nn)
    return np.median(nn)

def get_loc_list(sta):
    chanlist=[]
    for chan in sta:
        if chan.location_code not in chanlist:
            chanlist.append(chan.location_code)
    return chanlist

def get_meds(zn,zla,zlo,nsta=5):
    # for each noise measure, compute ratio to median of nsta nearest stations
    meds=[]
    if nsta<len(zn):
        tmpnsta=nsta
    else:
        tmpnsta=len(zn)-1
        
    for n in range(len(zn)):
        tmpdists=[]
        for m in range(len(zn)):
            if m==n:
                tmpdists.append(0.0)
            else:
                epi_dist, az, baz = gps2dist_azimuth(zla[n],zlo[n],zla[m],zlo[m])
                tmpdists.append(epi_dist/1000.)
        ii=np.argsort(tmpdists)
        tmp2=[]
        for m in range(tmpnsta+1):
            tmp2.append(zn[ii[m]])
        meds.append(zn[n]/np.median(tmp2))
            
    return meds

def find_lag(y1, y2):
    
    y1.data=y1.data/np.mean(abs(y1.data))
    y2.data=y2.data/np.mean(abs(y2.data))
    if y1.stats.sampling_rate > y2.stats.sampling_rate:
        y1.resample(y2.stats.sampling_rate)
    elif y1.stats.sampling_rate < y2.stats.sampling_rate:
        y2.resample(y1.stats.sampling_rate)
    sr=y1.stats.sampling_rate
    # trick to weight higher amplitudes
    maxsc=np.max(np.abs(y1.data))/np.log(2)
    y1.data=y1.data*(np.exp(np.abs(y1.data)/maxsc)-.0)
    maxsc=np.max(np.abs(y2.data))/np.log(2)
    y2.data=y2.data*(np.exp(np.abs(y2.data)/maxsc)-.0)
    y1.data=y1.data/np.mean(abs(y1.data))
    y2.data=y2.data/np.mean(abs(y2.data))
    
    n = len(y1.data)
    #print(y1.stats.sampling_rate,y2.stats.sampling_rate,sr)
    #corr = signal.correlate(y2.data, y1.data, mode='same') / np.sqrt(signal.correlate(y1.data, y1.data, mode='same')[int(n/2)] * signal.correlate(y2.data, y2.data, mode='same')[int(n/2)])
    corr = signal.correlate(y2.data, y1.data, mode='same') / np.sqrt(np.max(signal.correlate(y1.data, y1.data, mode='same')) * np.max(signal.correlate(y2.data, y2.data, mode='same')))
    if 0:
                            plt.figure(1)
                            plt.plot(y1.data,label=rstat.code)
                            plt.plot(y2.data,label=stat.code)
                            plt.legend()
                            plt.show()
    delay_arr = np.linspace(-0.5*n/sr, 0.5*n/sr, n)
    delay = delay_arr[np.argmax(corr)]
    return delay, np.max(corr)

parser = argparse.ArgumentParser(description='Check station amplitudes relative to other local stations')

parser.add_argument("-nets", action="store", dest="Nets",
                    required=True, help="Networks for analysis (required)")
parser.add_argument("-stas", action="store", dest="Stas",
                    default="*", help="Stations for analysis (default is *)")
parser.add_argument("-chans", action="store", dest="Chans",
                    default="BHZ", help="Channels for analysis (default is BHZ)")
parser.add_argument("-locs", action="store", dest="Locs",
                    default="*", help="location codes for analysis (default is *)")

parser.add_argument("-minlat", action="store", dest="MinLat",
                    default=-90, help="Min Latitude for QC")
parser.add_argument("-maxlat", action="store", dest="MaxLat",
                    default=90, help="Max Latitude for QC")
parser.add_argument("-minlon", action="store", dest="MinLon",
                    default=-180, help="Min Longitude for QC")
parser.add_argument("-maxlon", action="store", dest="MaxLon",
                    default=180, help="Max Longitude for QC")


parser.add_argument("-n", action="store", dest="Net",
                    default=999, help="Option to specify an individual reference network. If your use this flag, you must also use the -s flag.")
parser.add_argument("-s", action="store", dest="Sta",
                    default=999, help="Option to specify individual reference station. If you don't use this option, the first station in nets,stas will be used.")
parser.add_argument("-r", action="store", dest="Rad",
                    default=3, help="Option to specify radius (in degrees) from individual station when using the -n and -s flags (default is 3 degrees)")

parser.add_argument("-t1", action="store", dest="Time",
                    default=999,  help="Start Time for analysis (default is one week ago)")
parser.add_argument("-t2", action="store", dest="Time2",
                    default=999,  help="End Time for analysis (default is now)")

args = parser.parse_args()


net1=args.Net
sta1=args.Sta

rad=args.Rad
nets=args.Nets
stas=args.Stas
chans=args.Chans
myloc=args.Locs

minlat=args.MinLat
minlon=args.MinLon
maxlat=args.MaxLat
maxlon=args.MaxLon
tt = args.Time
tt2 = args.Time2

if tt==999:
    starttime = UTCDateTime()-3600*24*7
else:
    starttime=UTCDateTime(tt)

if tt2==999:
    endtime = UTCDateTime()
else:
    endtime=UTCDateTime(tt2)

plotrad=False
plotbox=False
#try:
#    xmlf = read_inventory(fil)
#except:
#    sys.exit("Couldn't read in file 1.")
    

model = TauPyModel(model="iasp91")
client = Client("IRIS")
#client = Client("https://service.scedc.caltech.edu")

###########
# to check timing, grab events at a good range for P-wave arrivals,
# use taup to get predicted arrivals, grab data, filter and
# cross correlate to see if predicted arrivals difference matches.
###########

# time for station analysis
#starttime = UTCDateTime("2021-01-21 00:00:00")
#endtime = UTCDateTime("2021-01-21 23:59:00")
#minlat=43.3
#minlon=-75
#maxlat=47
#maxlon=-68.3
#chans=['LHZ','LHN','LHE']
#lat=37.6
#lon=-84.6
#dep = 6
#rad=12
#Vp=5
#Vs=Vp/1.8

#stas= "BO04,MT01,BO03,FAR1,MT10,AC04,GO04,GO04,ROC1,AFO1,LCO"
#nets="US,N4"
#stas = "*"
#nets="IW"

#chans="LH*"
#chans="HHZ"

plotsect=0
debug = True

#fmin=1/8.
#fmax=1/4.
#gwidth=300

rsta=False
if net1 != 999:
    plotrad=True
    if sta1 != 999:
        rsta=True
        inv1=client.get_stations(network=net1,station=sta1,starttime=starttime, endtime=endtime,channel=chans, level="response")
        inventory = client.get_stations(network=nets,station=stas,starttime=starttime, endtime=endtime,
                                    maxradius=np.float(rad),latitude=inv1[0][0].latitude,longitude=inv1[0][0].longitude,channel=chans, location=myloc, level="response")
        rstat=inv1[0][0]
        rcnet=inv1[0].code
    else:
        sys.exit("You must specify both -n and -s for the individual station option.")
else:
    inventory = client.get_stations(network=nets,station=stas,starttime=starttime, endtime=endtime,
                                    minlatitude=minlat,minlongitude=minlon, maxlatitude=maxlat,maxlongitude=maxlon,channel=chans, location=myloc, level="response")
    rstat=inventory[0][0]
    rcnet=inventory[0].code

print(inventory)

slatsz =[]
slonsz =[]
slats1 =[]
slons1 =[]
slats2 =[]
slons2 =[]

znoise=[]
zstats=[]
h1noise=[]
h1stats=[]
h2noise=[]
h2stats=[]

r1=25
r2=90
mag=5.2
flow=.35
fhigh=1.5
pretim=60
winl=120

# Fetch a catalog relative to your reference station
try:
    cat = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=mag, latitude=rstat.latitude,
                                     longitude=rstat.longitude, minradius=r1, maxradius=r2)
except:
    sys.exit("No earthquakes found")
    
print(cat.__str__(print_all=True))
rdts=[] # this will be len(cat) * number of stations successfully processed
rcors=[] # this will be len(cat) * number of stations successfully processed
rtims=[] # this will be len(cat) * number of stations successfully processed
rmags=[] # this will be len(cat) * number of stations successfully processed
rstanums=[]  # this will be len(cat) * number of stations successfully processed

refrcors=[] #this will be len(cat)
othrcors=[] #this will be len(cat)
evals2=[] #this will be len(cat)
ecors=[]
evalstd=[] #this will be len(cat)


mags=[] #this will be len(cat)
tims=[] #this will be len(cat)
tims2=[] #this will be len(cat) of successful events
minetims=[] #this will be len(cat)
rlabs=[] #this will be len of stations

for i in range(len(cat)):
    minetims.append(UTCDateTime.utcnow())
print(len(minetims))
evnum=0
for evt in cat:
    try:
        evnum+=1
        #print("Event #%i"%evnum)
        #print(cat)

        mag1=evt.magnitudes[0].mag
        
        lat=evt.origins[0].latitude
        lon=evt.origins[0].longitude
        tim=evt.origins[0].time
        edeps=evt.origins[0].depth / 1000
        epi_dist, eaz, baz = gps2dist_azimuth(lat,lon, rstat.latitude, rstat.longitude)
        epi_dist = epi_dist / 1000
        rdistance_in_degree=epi_dist/111.1949
        rarrivals = model.get_travel_times(source_depth_in_km=edeps,
                  distance_in_degree=rdistance_in_degree,
                  phase_list=["P"])
        tr1=rarrivals[0].time-pretim
        if evnum==1:
            minetims[evnum-1]=tim+tr1
        elif tim+tr1 < minetims[evnum-1]:
            minetims[evnum-1]=tim+tr1
        #tr2=v2*distance_in_degree
        rst = client.get_waveforms(rcnet, rstat.code, rstat[0].location_code, rstat[0].code, tim+tr1, tim+tr1+winl, attach_response=True)
        rst.merge(fill_value=0)
        rst.trim(starttime=tim+tr1, endtime=tim+tr1+winl,pad=True,fill_value=0)
        # remove response and filter
        rst.detrend(type='linear')
        rst.taper(max_percentage=0.3)
        rst.remove_response(output="DISP")
        rst.detrend(type='simple')
        rst.taper(max_percentage=0.3)
        rst.filter(type='bandpass',freqmin=flow,freqmax=fhigh)
        
        
        mags.append(mag1)
        tt=tim+tr1
        tims.append(tt.datetime)
        stanums=[]
        event_dat=[]
        ddegs=[]
        stanum=0 # the reference trace will be stanum 0
        # loop through and collect all the data for this event               
        for cnet in inventory:
            for stat in cnet:
                if (rstat.code != stat.code):
                    locs=get_loc_list(stat)
                    #print(stat.code,chan)
                    for loc in locs:
                        stanum+=1
                        try:
                            
                            epi_dist, eaz, baz = gps2dist_azimuth(lat,lon, stat.latitude, stat.longitude)
                            epi_dist = epi_dist / 1000
                            distance_in_degree=epi_dist/111.1949
                            arrivals = model.get_travel_times(source_depth_in_km=edeps,
                                      distance_in_degree=distance_in_degree,
                                      phase_list=["P"])
                            tr1=arrivals[0].time-pretim
                            if tim+tr1 < minetims[evnum-1]:
                                minetims[evnum-1]=tim+tr1
                            #tr2=v2*distance_in_degree
                            st = client.get_waveforms(cnet.code, stat.code, loc, chans, tim+tr1, tim+tr1+winl, attach_response=True)
                            st.merge(fill_value=0)
                            st.trim(starttime=tim+tr1, endtime=tim+tr1+winl,pad=True,fill_value=0)
                            # remove response and filter
                            st.detrend(type='linear')
                            st.taper(max_percentage=0.3)
                            st.remove_response(output="DISP")
                            st.detrend(type='simple')
                            st.taper(max_percentage=0.3)
                            st.filter(type='bandpass',freqmin=flow,freqmax=fhigh) # should use multiband filter for best correlation.
                            event_dat.append(st.copy())
                            stanums.append(stanum)
                            ddegs.append(distance_in_degree)
                        except:
                            print("Could not fetch %s-%s-%s-%s" % (cnet.code, stat.code, loc, chans))
        # loop through, correlating every station with the reference
        dts=[]
        cors=[]
        cors2=[]
        nn=-1
        for st in event_dat:
            nn+=1
           
            if 1:
                rstc=copy.deepcopy(rst)
                stc=copy.deepcopy(st)
                dt,cor = find_lag(rstc[0],stc[0])
                cors.append(cor)
                if cor >=.7:
                    dts.append(-dt)
                    rdts.append(-dt)
                    rcors.append(cor)
                    rtims.append(tt.datetime)
                    rmags.append(mag1)
                    rstanums.append(stanums[nn])
                    cors2.append(cor)

                    #mags.append(mag1)
                    #tt=tim+tr1
                    #tims.append(tt.datetime)
                print('%s, Delay is %3.2f seconds, %s relative to %s, corr=%3.2f'%(tim,-dt,rstat.code, st[0].stats.station,cor))
                if 0:
                    #maxsc=np.max([np.max(np.abs(rst[0].data)),np.max(np.abs(st[0].data))])/1.5
                    plt.figure(1)
                    #plt.plot(rst[0].data*np.exp(np.abs(rst[0].data)/maxsc),label=rstat.code)
                    #plt.plot(st[0].data*np.exp(np.abs(st[0].data)/maxsc),label=stat.code)
                    plt.plot(rst[0].data,label=rstat.code)
                    plt.plot(st[0].data,label=st[0].stats.station)
                    plt.title('%s, Delay is %3.2f seconds, %s relative to %s, corr=%3.2f'%(tim,-dt,rstat.code, st[0].stats.station,cor))
                    plt.legend()
                    plt.show()
            else:
                print("Could not fetch %s-%s" % (cnet.code, stat.code))
                
        refrcors.append(np.mean(cors))
        if len(dts)>=1:
            evals2.append(np.mean(dts))
            evalstd.append(np.std(dts))
            tims2.append(tt.datetime)
            ecors.append(np.mean(cors2))
        # loop through all other stations to compute average correlation
        
        cors=[]
        for st1 in event_dat:
            for st2 in event_dat:
                if st1[0].stats.station != st2[0].stats.station:
                    stc1=copy.deepcopy(st1)
                    stc2=copy.deepcopy(st2)
                    dt,cor = find_lag(stc1[0],stc2[0])
                    cors.append(cor)
                
        othrcors.append(np.median(cors))

        
        if len(dts) and plotsect:  # plot event section
            tr_scale=(np.max(ddegs)-np.min(ddegs))/10
            plt.figure(2)
            n=-1
            for st in event_dat:
                n+=1
                # trick to weight higher amplitudes
                maxsc=np.max(np.abs(st[0].data))/np.log(2)
                st[0].data=st[0].data*(np.exp(np.abs(st[0].data)/maxsc)-.0)
                st[0].data=tr_scale*(st[0].data/np.std(st[0].data*3))+ ddegs[n]
                dgs=np.arange(0,len(st[0].data))*st[0].stats.delta -pretim
                plt.plot(dgs,st[0].data,label=st[0].stats.station,alpha=.7)
            maxsc=np.max(np.abs(rst[0].data))/np.log(2)
            rst[0].data=rst[0].data*(np.exp(np.abs(rst[0].data)/maxsc)-.0)
            stdata=tr_scale*(rst[0].data/np.std(rst[0].data*3))+ rdistance_in_degree
            dgs=np.arange(0,len(stdata))*rst[0].stats.delta -pretim
            plt.plot(dgs,stdata,'k',label=rstat.code)
            plt.xlabel('time relative to predicted arrival (s)')
            plt.ylabel('distance to earthquake (degrees)')
            plt.legend()
            plt.show()
    except:
        print("Couldn't fetch reference station data")
    


stanum=0 # the reference trace will be stanum 0
# loop through and plot station results               
for cnet in inventory:
    for stat in cnet:
        if (rstat.code != stat.code):
            locs=get_loc_list(stat)
            #print(stat.code,chan)
            for loc in locs:
                stim=[]
                dts=[]
                stanum+=1
                for n in range(len(rtims)):
                    if rstanums[n]==stanum:
                        dts.append(rdts[n])
                        stim.append(rtims[n])
                if len(dts)>=1:
                                
                    plt.figure(6)
                    plt.plot(stim,dts,'d',label=rcnet+"-"+rstat.code+"-"+rstat[0].location_code+" vs. "+cnet.code+"-"+stat.code+"-"+loc)
                    plt.xlabel('date')
                    plt.ylabel('time delay (s)')
                    #plt.plot(st[0].data,label=stat.code)
                    plt.grid()
                    plt.gcf().autofmt_xdate()
                    plt.legend()
        
etims2=np.array(tims2)
evals2=np.array(evals2)
evalstd=np.array(evalstd)

plt.figure(7)
plt.plot(etims2,evals2,'kd',label=rcnet+"-"+rstat.code+"-"+rstat[0].location_code+" mean")
plt.plot(etims2,evals2+evalstd,'k+',label="mean +/- 1 std")
plt.plot(etims2,evals2-evalstd,'k+')

plt.xlabel('date')
plt.ylabel('mean error (s)')
#plt.plot(st[0].data,label=stat.code)
plt.grid()
plt.gcf().autofmt_xdate()
plt.legend()

plt.figure(8)
plt.plot(tims,np.asarray(refrcors)/np.asarray(othrcors),'k.',label=rcnet+"-"+rstat.code+"-"+rstat[0].location_code+" mean correlation relative to network mean")
#plt.plot(tims,othrcors,'kd',label="rest of network mean correlation")


plt.xlabel('date')
plt.ylabel('correlation value')
#plt.plot(st[0].data,label=stat.code)
plt.grid()
plt.gcf().autofmt_xdate()
plt.legend()

plt.figure(9)
plt.plot(rmags,rcors,'k.')
#plt.plot(tims,othrcors,'kd',label="rest of network mean correlation")


plt.xlabel('magnitude')
plt.ylabel('correlation value')
#plt.plot(st[0].data,label=stat.code)
plt.grid()

plt.figure(10)
plt.plot(rcors,rdts,'k.', label="all measurements")
plt.plot([np.min(rcors),np.max(rcors)],[np.mean(rdts)+3*np.std(rdts),np.mean(rdts)+3*np.std(rdts)],'k--', label="95% measurement confidence")
plt.plot([np.min(rcors),np.max(rcors)],[np.mean(rdts)-3*np.std(rdts),np.mean(rdts)-3*np.std(rdts)],'k--')
plt.plot(ecors,evals2,'go', label="event means")
plt.plot([np.min(rcors),np.max(rcors)],[np.mean(rdts)+3*np.std(rdts)/np.sqrt(5),np.mean(rdts)+3*np.std(rdts)/np.sqrt(5)],'g--', label="95% event means confidence")
plt.plot([np.min(rcors),np.max(rcors)],[np.mean(rdts)-3*np.std(rdts)/np.sqrt(5),np.mean(rdts)-3*np.std(rdts)/np.sqrt(5)],'g--')
#plt.plot(tims,othrcors,'kd',label="rest of network mean correlation")

plt.title('%i measurements, %i events, time mean = %3.2f, std = %3.2f'%(len(rdts),len(ecors),np.mean(rdts),np.std(rdts)))
plt.ylabel('delay time')
plt.xlabel('correlation value')
plt.legend()
#plt.plot(st[0].data,label=stat.code)
plt.grid()

plt.show()

                

