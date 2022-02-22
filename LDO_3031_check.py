#!/usr/bin/env python

from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
from obspy import read
import datetime as datetime
import sys
import argparse
import copy
import os

def get_loc_list(sta):
    chanlist=[]
    for chan in sta:
        if chan.location_code not in chanlist:
            chanlist.append(chan.location_code)
    return chanlist

def elv_cor(stel,stlat):
    # more accurate elevation correction
    P_b = 101325 # reference pressure (Pa)
    T = 293 -25*(1-np.cos(np.abs(stlat)*np.pi/180))# reference temperature (K), 20 C = 293 K
    #L_L = temperature lapse rate (K/m) in ISA
    hh = stel # height at which pressure is calculated (m)
    h_b = 0 # height of reference level b (meters; e.g., hb = 11 000 m)
    R= 8.3144598 # = universal gas constant: 8.3144598 J/(mol·K)
    g0=9.798 + (1-np.cos(np.abs(stlat)*np.pi/180))*(9.863-9.798) # = gravitational acceleration: 9.80665 m/s2, 9.798 - 9.863
    M = .0289644 # = molar mass of Earth's air: 0.0289644 kg/mol
    P= P_b * np.exp( (-g0 * M * (hh - h_b)) / (R * T))
    return P_b - P

parser = argparse.ArgumentParser(description='Fetch LDO channels and compare 30-LDO vs. 31-LDO')

parser.add_argument("-nets", action="store", dest="Nets",
                    default="IU", help="Networks for analysis (default is IU)")
parser.add_argument("-stas", action="store", dest="Stas",
                    default="*", help="Stations for analysis (default is *)")
"""
parser.add_argument("-chans", action="store", dest="Chans",
                    default="LDO", help="Channels for analysis (default is LDO)")
parser.add_argument("-locs", action="store", dest="Locs",
                    default="30,31", help="location codes for analysis (default is *)")

parser.add_argument("-n", action="store", dest="Net",
                    default=999, help="Option to specify an individual reference network. If your use this flag, you must also use the -s flag.")
parser.add_argument("-s", action="store", dest="Sta",
                    default=999, help="Option to specify individual reference station. If you don't use this option, the first station in nets,stas will be used.")
parser.add_argument("-r", action="store", dest="Rad",
                    default=20, help="Option to specify radius (in degrees) from individual station when using the -n and -s flags (default is 3 degrees)")
"""
parser.add_argument("-t1", action="store", dest="Time",
                    default=999,  help="Start Time for analysis (default is two days before end time)")
parser.add_argument("-t2", action="store", dest="Time2",
                    default=999,  help="End Time for analysis (default is now)")

args = parser.parse_args()

plotme=True
writefile=True
frommsd=False

#net1=args.Net
#sta1=args.Sta

#rad=args.Rad
nets=args.Nets
stas=args.Stas

#chans=args.Chans
#Myloc=args.Locs
chans="LDO"
Myloc="30,31"

tt = args.Time
tt2 = args.Time2


if tt2==999:
    endtime = UTCDateTime()
    #endtime = starttime+3600*24
else:
    endtime=UTCDateTime(tt2)

if tt==999:
    starttime = endtime-3600*24*2
    #starttime = UTCDateTime()-3600*24*2
    #starttime = UTCDateTime("2020-01-16 00:00:00")
else:
    starttime=UTCDateTime(tt)

client = Client("IRIS")
try:
    inventory = client.get_stations(network=nets,station=stas,starttime=starttime, endtime=endtime,
                            channel=chans, location=Myloc, level="response")
except:
    print("No channels found for %s-%s-%s-%s %s"%(nets,stas,chans,Myloc,starttime.strftime("%Y-%m-%d %H:%M")))
    sys.exit()

try:
    fls = os.popen('ls -ld /msd/IU_ANMO')
    if 'cannot' not in fls:
        frommsd=True
    #print(fls)
except:
    frommsd=False
    
#print(inventory)

sncls=[]
elevs=[]
lats=[]
cormeans=[]
errs=[]
if writefile:
    f=open("LDO_out.csv","w+")
ecor=11.897 # Pa/m
Pa2atm=1/101325
plt.figure(3)
plt.plot([101325*.97, 101325*1.03],[101325*.97, 101325*1.03],'b--',label="1 to 1 line")
print("Channel,  elevation, coef0, coef1, mean raw counts, Corrected Pressure (Pa), Pressure metric")
for cnet in inventory:
    for stat in cnet:
        if 1:
            elv = stat.elevation
            
            #locs=get_loc_list(stat)
            #print(stat.code,chan)
            val30=[]
            val31=[]
            for loc in ["30","31"]:
                
                try:
                    
                    inv1=inventory.select(cnet.code, stat.code, loc, chans)
                    mychan=inv1[0][0][0].response
                    #print(inv1[0][0][0])
                    #print("/msd/%s_%s/%s/%s_%s*"%(cnet.code,stat.code,starttime.strftime("%Y/%j"),loc,inv1[0][0][0].code))
                    #cfs=mychan.instrument_polynomial.coefficients
                    #print(cfs)
                    try:
                        st = client.get_waveforms(cnet.code, stat.code, loc, chans, starttime, endtime, attach_response=True)
                    except:
                        
                        if frommsd:
                            try:
                                st=read("/msd/%s_%s/%s/%s_%s*"%(cnet.code,stat.code,starttime.strftime("%Y/%j"),loc,inv1[0][0][0].code))
                                print("Couldn't fetch %s-%s-%s-%s, loaded from /msd instead"%(cnet.code, stat.code, loc, chans))
                            except:
                                print("Couldn't fetch %s-%s-%s-%s, or load from /msd "%(cnet.code, stat.code, loc, chans))
                                break
                        
                    #st.merge(fill_value=0)
                    #st.trim(starttime=starttime, endtime=endtime,pad=True,fill_value=0)
                    # remove response and filter
                    if mychan.instrument_polynomial:
                        cfs=mychan.instrument_polynomial.coefficients
                    else:
                        
                        print("%s-%s-%s-%s not a polynomial, sensitivity is %9.2f counts per %s"%(cnet.code, stat.code, loc, chans, mychan.instrument_sensitivity.value, mychan.instrument_sensitivity.input_units))
                        if "kPa" in mychan.instrument_sensitivity.input_units:
                            cfs=[0.,1000./mychan.instrument_sensitivity.value]
                        else:
                            cfs=[0.,1./mychan.instrument_sensitivity.value]

                        
                    if len(st)>1:
                        print("warning, %i gaps for %s-%s-%s-%s"%(len(st),cnet.code, stat.code, loc, chans))
                    trstd=np.std(st[0].data)
                    trmean=np.mean(st[0].data)
                    cormean=trmean*cfs[1] + cfs[0] + elv_cor(elv,stat.latitude)
                    err=(cormean-cfs[0])/(101325-cfs[0])
                    err=cormean*Pa2atm
                    elevs.append(elv)
                    lats.append(stat.latitude)
                    errs.append(err)
                    cormeans.append(cormean)
                    sncls.append("%s-%s-%s-%s"%(cnet.code,stat.code,loc,st[0].stats.channel))
                    #print("%s-%s-%s-%s, %5.1f, %5.3e, %5.3e, %5.3e, %5.3e, %5.3e,%5.2f"%(cnet.code,stat.code,loc,st[0].stats.channel,elv,cfs[0],cfs[1],trmean,cormean,trstd,err))
                    print("%s-%s-%s-%s, %5.1f, %5.3e, %5.3e, %5.3e, %5.3e, %5.2f"%(cnet.code,stat.code,loc,st[0].stats.channel,elv,cfs[0],cfs[1],trmean,cormean,err))
                    #if writefile:
                    #    f.write("%s-%s,%s,%s, %5.1f, %5.3e, %5.3e, %5.3e, %5.3e, %5.2f \n"%(cnet.code,stat.code,loc,st[0].stats.channel,elv,cfs[0],cfs[1],trmean,cormean,err))
                    if "30" in loc:
                        val30.append(cormean)
                    else:
                        val31.append(cormean)
                   
                except:
                    print("Couldn't process %s-%s-%s-%s"%(cnet.code, stat.code, loc, chans))
            if len(val30) > 0 and len(val31) > 0:
                if writefile:
                        f.write("%s-%s, %5.1f, %5.3f \n"%(cnet.code,stat.code,elv,val31[0]/val30[0]))
                if val30[0]/val31[0] > 1.005 or val31[0]/val30[0] > 1.005:
                    plt.plot(val30[0],val31[0],'o',label=stat.code)
                else:
                    plt.plot(val30[0],val31[0],'ko')
if writefile:
    f.close()
if plotme:
    plt.figure(1)
    plt.plot(cormeans,elevs,'ko')
    plt.plot([101325, 101325],[0,np.max(elevs)],'b--',label="1 atm")
    for nn in range(len(errs)):
        if cormeans[nn] < 97000 or cormeans[nn]>104000:
            plt.plot(cormeans[nn],elevs[nn], 'o',label=sncls[nn])
    plt.ylabel('elevation (m)')
    plt.xlabel('Corrected Pressure (Pa)')
    plt.legend()
    plt.title('%s-%s-%s-%s %s-%s'%(nets,stas,chans,Myloc,starttime.strftime("%Y %m/%d"),endtime.strftime("%m/%d")))
    #plt.show()
    plt.figure(2)
    plt.plot(lats,cormeans,'ko')
    plt.plot([-90,90],[101325, 101325],'b--',label="1 atm")
    plt.plot([-60,-60],[101325*.99, 101325],'g:',label="predicted lows")
    plt.plot([0,0],[101325*.99, 101325],'g:')
    plt.plot([60,60],[101325*.99, 101325],'g:')
    plt.plot([-30,-30],[101325, 101325*1.01],'r:',label="predicted highs")
    plt.plot([30,30],[101325, 101325*1.01],'r:')
    for nn in range(len(errs)):
        if cormeans[nn] < 97000 or cormeans[nn]>104000:
            plt.plot(lats[nn],cormeans[nn], 'o',label=sncls[nn])
    plt.xlabel('latitude')
    plt.ylabel('Corrected Pressure (Pa)')
    plt.legend()
    plt.title('%s-%s-%s-%s %s-%s'%(nets,stas,chans,Myloc,starttime.strftime("%Y %m/%d"),endtime.strftime("%m/%d")))
    plt.figure(3)
    #plt.plot(cormeans,elevs,'ko')
    
    
    plt.ylabel('Corrected 31-LDO Pressure (Pa)')
    plt.xlabel('Corrected 30-LDO Pressure (Pa)')
    plt.legend()
    plt.title('%s-%s-%s-%s %s-%s'%(nets,stas,chans,Myloc,starttime.strftime("%Y %m/%d"),endtime.strftime("%m/%d")))
    
    plt.show()


