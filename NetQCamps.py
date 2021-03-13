#!/usr/bin/env python

from obspy.clients.fdsn import Client
#from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as ml
from obspy.core import UTCDateTime
#from scipy.interpolate import griddata, interp1d
from scipy.interpolate import griddata
from obspy.geodetics import gps2dist_azimuth
import matplotlib
from obspy.taup import TauPyModel
#from scipy import signal as scisig
#from obspy.signal import spectral_estimation as spec
from matplotlib import cm
#import pickle
#import gc
#from numba import jit
from matplotlib.transforms import blended_transform_factory
#from scipy.optimize import leastsq
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import cartopy

import argparse

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

parser = argparse.ArgumentParser(description='Check station amplitudes relative to other local stations')

parser.add_argument("-nets", action="store", dest="Nets",
                    required=True, help="Networks for analysis")
parser.add_argument("-stas", action="store", dest="Stas",
                    default="*", help="Stations for analysis")
parser.add_argument("-chans", action="store", dest="Chans",
                    default="LH*", help="Channels for analysis")

parser.add_argument("-minlat", action="store", dest="MinLat",
                    default=-90, help="Min Latitude for QC")
parser.add_argument("-maxlat", action="store", dest="MaxLat",
                    default=90, help="Max Latitude for QC")
parser.add_argument("-minlon", action="store", dest="MinLon",
                    default=-180, help="Min Longitude for QC")
parser.add_argument("-maxlon", action="store", dest="MaxLon",
                    default=180, help="Max Longitude for QC")



parser.add_argument("-n", action="store", dest="Net",
                    default=999, help="Option to specify net individual station. If your use this flag, you must also use the -s flag.")
parser.add_argument("-s", action="store", dest="Sta",
                    default=999, help="Option to specify individual station")
parser.add_argument("-r", action="store", dest="Rad",
                    default=3, help="Option to specify radius from individual station")

parser.add_argument("-t1", action="store", dest="Time",
                    default=999,  help="Start Time for analysis (default is 24hrs ago)")
parser.add_argument("-t2", action="store", dest="Time2",
                    default=999,  help="End Time for analysis (default is now)")

args = parser.parse_args()


net1=args.Net
sta1=args.Sta
rad=args.Rad
nets=args.Nets
stas=args.Stas
chans=args.Chans
minlat=args.MinLat
minlon=args.MinLon
maxlat=args.MaxLat
maxlon=args.MaxLon
tt = args.Time
tt2 = args.Time2

if tt==999:
    starttime = UTCDateTime()-3600*24
else:
    starttime=UTCDateTime(tt)

if tt2==999:
    endtime = UTCDateTime()
else:
    endtime=UTCDateTime(tt2)


#try:
#    xmlf = read_inventory(fil)
#except:
#    sys.exit("Couldn't read in file 1.")
    

model = TauPyModel(model="iasp91")
client = Client("IRIS")
#client = Client("https://service.scedc.caltech.edu")


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
dep = 6
#rad=12
Vp=5
Vs=Vp/1.8

#stas= "BO04,MT01,BO03,FAR1,MT10,AC04,GO04,GO04,ROC1,AFO1,LCO"
#nets="US,N4"
#stas = "*"
#nets="IW"

#chans="LH*"
#chans="HHZ"

debug = True

fmin=1/8.
fmax=1/4.
#gwidth=300


if net1 != 999:
    if sta1 != 999:
        inv1=client.get_stations(network=net1,station=sta1,starttime=starttime, endtime=endtime)
        inventory = client.get_stations(network=nets,station=stas,starttime=starttime, endtime=endtime,
                                    maxradius=rad,latitude=inv1[0][0].latitude,longitude=inv1[0][0].longitude,channel=chans, level="response")
    else:
        sys.exit("You must specify both -n and -s for the individual station option.")
else:
    inventory = client.get_stations(network=nets,station=stas,starttime=starttime, endtime=endtime,
                                    minlatitude=minlat,minlongitude=minlon, maxlatitude=maxlat,maxlongitude=maxlon,channel=chans, level="response")

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

for cnet in inventory:
    for stat in cnet:
                       
        locs=get_loc_list(stat)
        #print(stat.code,chan)
        for loc in locs:
            try:
                st = client.get_waveforms(cnet.code, stat.code, loc, chans, starttime, endtime, attach_response=True)
                #print(st)
                #st.plot()
                st.merge(fill_value=0)
                #if debug == True:
                #    st.plot()
                #print(st)
                print("Fetched %i traces, %i samples, at %3.1f sample rate"%(len(st),st[0].stats.npts,st[0].stats.sampling_rate))
                #st.decimate(int(np.floor(st[0].stats.sampling_rate / 50)))
                #print(st)
                st.trim(starttime=starttime, endtime=endtime,pad=True,fill_value=0)
                #if debug == True:
                #    st.plot()
                #print(st)
                print("Trimmed to %i samples, at %3.1f sample rate"%(st[0].stats.npts,st[0].stats.sampling_rate))
                for tr in st:
                    n2=process_stream(tr)
                    mychan=tr.stats.channel
                    if "Z" in mychan[2]:
                            znoise.append(n2)
                            zstats.append(stat.code)
                            slatsz.append(stat.latitude)
                            slonsz.append(stat.longitude)
                    elif "N" in mychan[2] or "1" in mychan[2]:
                            h1noise.append(n2)
                            h1stats.append(stat.code)
                            slats1.append(stat.latitude)
                            slons1.append(stat.longitude)
                    elif "E" in mychan[2] or "2" in mychan[2]:
                            h2noise.append(n2)
                            h2stats.append(stat.code)
                            slats2.append(stat.latitude)
                            slons2.append(stat.longitude)
                    
                            
            except:
                print("Could not fetch %s-%s" % (cnet.code, stat.code))

#with open('Jul15SV0Noise.pickle', 'wb') as f3:
#    pickle.dump([zstats, znoise, h1noise, h1stats, h2noise, h2stats, slats, slons], f3)
#f3.close()
#slats = np.asarray(slats)
#slons = np.asarray(slons)

znoise = np.asarray(znoise)
h1noise = np.asarray(h1noise)
h2noise = np.asarray(h2noise)

medz=np.median(znoise,axis=0)
medh1=np.median(h1noise,axis=0)
medh2=np.median(h2noise,axis=0)

if 1:
    zmeds=get_meds(znoise,slatsz,slonsz)
    h1meds=get_meds(h1noise,slats1,slons1)
    h2meds=get_meds(h2noise,slats2,slons2)
    print(len(zmeds),len(h1meds),len(h2meds))
if 1:
    fig=plt.figure(num=4, figsize=(8,8))
    plt.plot( slonsz,zmeds,'kd',markersize=8,linewidth=4,label='vertical')
    plt.plot( slons1,h1meds,'g+',markersize=8,linewidth=4,label='h1')
    plt.plot( slons2,h2meds,'bx',markersize=8,linewidth=4,label='h2')                
    plt.xlabel('Longitude')
    plt.ylabel('Microseismic amp relative to local average')
    ax1 = fig.axes[0]
    transform = blended_transform_factory(ax1.transData, ax1.transAxes)
    for n in range(len(slonsz)):
        ax1.text(slonsz[n], 1.0, zstats[n], rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10)
    plt.grid()
    plt.legend()
    #plt.show()

if 1:
    fig=plt.figure(num=5, figsize=(8,8))
    
    nn=0
    #transform = blended_transform_factory(ax.transData, ax.transAxes)
    #transform = blended_transform_factory(ax.transAxes,ax.transData)
    labels=[]
    for n in range(len(slonsz)):
        sig=1.5
        #if ((zmeds[n] < np.mean(zmeds)-sig*np.std(zmeds)) or (zmeds[n] > np.mean(zmeds)+sig*np.std(zmeds)) or
        #        (h1meds[n] < np.mean(h1meds)-sig*np.std(h1meds)) or (h1meds[n] > np.mean(h1meds)+sig*np.std(h1meds)) or
        #        (h2meds[n] < np.mean(h2meds)-sig*np.std(h2meds)) or (h2meds[n] > np.mean(h2meds)+sig*np.std(h2meds))):
        if ((zmeds[n] < .55) or (zmeds[n] > 1/.55) or
                (h1meds[n] < .55) or (h1meds[n] > 1/.55) or
                (h2meds[n] < .55) or (h2meds[n] > 1/.55)):    
            nn=nn+1
            if nn==1:
                plt.semilogx( np.max((.04,zmeds[n])),-nn,'kd',markersize=8,linewidth=4,label='vertical')
                plt.semilogx( np.max((.04,h1meds[n])),-nn,'g+',markersize=8,linewidth=4,label='h1')
                plt.semilogx( np.max((.04,h2meds[n])),-nn,'bx',markersize=8,linewidth=4,label='h2')
            else:
                plt.semilogx( np.max((.04,zmeds[n])),-nn,'kd',markersize=8,linewidth=4)
                plt.semilogx( np.max((.04,h1meds[n])),-nn,'g+',markersize=8,linewidth=4)
                plt.semilogx( np.max((.04,h2meds[n])),-nn,'bx',markersize=8,linewidth=4)
            ax = fig.axes[0]
            ax.text(10,-nn, zstats[n], rotation=0,
                va="center", ha="right" ,zorder=10)
            #labels.append(zstats[n])
    if nn>0:
        xx=ax.get_xlim()
        ax.set_xlim(left=xx[0],right=np.max((xx[1],12.)))
        ax.set_yticklabels([''])
        #plt.xlabel('Longitude')
        plt.xlabel('Microseismic amp relative to local average')
        
        plt.grid()
        plt.legend(loc='upper left')
        #plt.show()

if 1:
    x=np.asarray(slonsz)
    y=np.asarray(slatsz)
    #z=np.asarray(20*np.log10(znoise))
    z=np.asarray(znoise/zmeds)
    #extent=[np.min(y)-.2, np.min(x)-.2, np.max(y)+.2, np.max(x)+.2]
    extent=[np.min(x)-.2, np.max(x)+.2, np.min(y)-.2, np.max(y)+.2]
    central_lon = np.median(x)
    central_lat = np.median(y)
    nbins=300
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]

    zi = griddata( (x,y), z, (xi, yi), method='linear')
    

#with open('DetectGrid%s%s.pickle'%(titl,starttime.strftime('%Y%m%d%H')),'wb') as f:
#    pickle.dump([xi,yi,zi],f)
#f.close()
    
if 1:
    
    c1=min(z)
    c2=max(z)
    if c2>c1:
        plt.figure(7, figsize=(6,6))
        #c1=-145
        #c2=-111
        #c1=1.0
        #c2=2.9
        #ax = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
        ax = plt.axes(projection=ccrs.Mercator(central_lon))
        ax.set_extent(extent)

        ax.add_feature(cartopy.feature.OCEAN)
        ax.add_feature(cartopy.feature.LAND, edgecolor='black')
        ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
        ax.add_feature(cartopy.feature.STATES, edgecolor='black')

        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        #print(xi.shape,c1,c2)
        plt.contourf(xi, yi, zi.reshape(xi.shape), np.arange(c1, c2, (c2-c1)/10), cmap=plt.cm.plasma, transform=ccrs.PlateCarree() )
        gridlines=ax.gridlines(draw_labels=True, color='gray', alpha=.8, linestyle=':')
        gridlines.xlabels_top=False
        gridlines.ylabels_right=False
        gridlines.xlocator = mticker.FixedLocator(np.arange(np.floor(np.min(x)-.2), np.ceil(np.max(x)+.2),2))
        gridlines.ylocator = mticker.FixedLocator(np.arange(np.floor(np.min(y)-.2), np.ceil(np.max(y)+.2),2))

        # Add color bar
        plt.clim(c1,c2)
        
        cbar=plt.colorbar()
        cbar.set_label('Locally weighted noise (5 station median)')

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        #plt.title("%s - %s "% (letter.upper(), content.upper()))
        #for sta in Sdict:
        #    plt.plot(Sdict[sta]['lon'],Sdict[sta]['lat'], 'kd', markersize=4.5, transform=ccrs.PlateCarree())
        plt.plot(slonsz,slatsz, 'kd', markersize=4.5, transform=ccrs.PlateCarree())
    plt.show()
