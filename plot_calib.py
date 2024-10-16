import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math

def variance(data, ddof=0):
    n = len(data)
    mean = sum(data) / n
    return sum((x - mean) ** 2 for x in data) / (n - ddof)

def scalingrelation(datax,datay,bins,option='median'):
    y=[]
    x=[]
    err=[]
    dim=len(bins)-1
    for i in range(dim):
        aa=datay[(datax > bins[i]) & (datax < bins[i+1])]
        if (len(aa) ==0):
            continue
        if (option =='mean'):
            meanv=np.average(aa)
            sd = np.sqrt(variance(aa))#/(len(aa)-1)) #np.std(aa)#
            x.append(bins[i]/2+bins[i+1]/2)
            y.append(meanv)
            err.append([sd,sd])
        if (option =='median'):
            x.append(bins[i]/2+bins[i+1]/2)
            y.append(np.median(aa))
            err.append(np.abs(np.percentile(aa,[30,70])-np.median(aa)))
    x=np.array(x)
    y=np.array(y)
    err=np.array(err).T
    return x,y,err

#HIMF

def jones2018_himf(logM, phi=0.0045, log_M0=9.94, alpha=-1.25, e=2.718281828):
    schechter = math.log(10) * phi*((10**(logM-log_M0))**(alpha+1))*(e**(-1* 10**(logM-log_M0)))
    return schechter

def obs_himf(redshift,iplot):
    if (redshift=='0' or redshift=='0.0' or redshift==0):
        x,y=readfile('/home/lzxie/Work/obs_data/zwaan2005_himf_z_0')
        ax[i].plot(x,np.log10(y),color='k',marker='^',alpha=0.5,linestyle='none',label='Zwaan et al. 2005')
        x,y=readfile('/home/lzxie/Work/obs_data/haynes2011_himf_z_0')
        ax[i].plot(x,np.log10(y),color='k',marker='s',alpha=0.5,linestyle='none',label='Haynes et al. 2011')
        

def readfile(filename,option='nonerror'):
    xx=[]
    yy=[]
    if (option == 'err'):
        err1=[]
        err2=[]
    with open(filename,'r') as f:
        while True:
            line=f.readline()
            if not line:
                break
            if line.strip().startswith("#"):
                continue
            if (option == 'err'):
                x,y,e1,e2=[float(i) for i in line.split()]
            else:
                x,y=[float(i) for i in line.split()]
            xx.append(x)
            yy.append(y)
            if (option == 'err'):
                err1.append(e1)
                err2.append(e2)
    if (option == 'err'):
        xx,yy,err1,err2=np.array(xx),np.array(yy),np.array(err1),np.array(err2)
        return xx,yy,err1,err2
    else:
        xx,yy=np.array(xx),np.array(yy)
        return xx,yy
    
def xfraction(sfr,mstar,mmin,mmax,binsize,sfr_limit):
    xx=[]
    yy=[]
    m=mmin
    if (max(mstar) < mmax):
        mmax = max(mstar)
    while(m < mmax):
        tmp=sfr[(mstar > m) & (mstar < m+binsize)]
        if(len(tmp)==0):
            continue
#        quench = tmp[tmp<sfr_limit]
#        fq=len(quench)*1.0/len(tmp)
        fq=np.sum(tmp < sfr_limit)*1.0/len(tmp)
        yy.append(fq)
        xx.append(m+binsize/2)
        m=m+binsize
    xx=np.array(xx)
    yy=np.array(yy)
    return xx,yy    

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math

def scalingrelation(datax,datay,bins,option='median'):
    y=[]
    x=[]
    err=[]
    dim=len(bins)-1
    for i in range(dim):
        aa=datay[(datax > bins[i]) & (datax < bins[i+1])]
        if (len(aa) ==0):
            continue
        if (option =='mean'):
            meanv=np.average(aa)
            sd = np.sqrt(variance(aa))#/(len(aa)-1)) #np.std(aa)#
            x.append(bins[i]/2+bins[i+1]/2)
            y.append(meanv)
            err.append([sd,sd])
        if (option =='median'):
            x.append(bins[i]/2+bins[i+1]/2)
            y.append(np.median(aa))
            err.append(np.abs(np.percentile(aa,[30,70])-np.median(aa)))
    x=np.array(x)
    y=np.array(y)
    err=np.array(err).T
    return x,y,err

#HIMF

def jones2018_himf(logM, phi=0.0045, log_M0=9.94, alpha=-1.25, e=2.718281828):
    schechter = math.log(10) * phi*((10**(logM-log_M0))**(alpha+1))*(e**(-1* 10**(logM-log_M0)))
    return schechter

def obs_himf(redshift,iplot):
    if (redshift=='0' or redshift=='0.0' or redshift==0):
        x,y=readfile('/home/lzxie/Work/obs_data/zwaan2005_himf_z_0')
        ax[i].plot(x,np.log10(y),color='k',marker='^',alpha=0.5,linestyle='none',label='Zwaan et al. 2005')
        x,y=readfile('/home/lzxie/Work/obs_data/haynes2011_himf_z_0')
        ax[i].plot(x,np.log10(y),color='k',marker='s',alpha=0.5,linestyle='none',label='Haynes et al. 2011')
        

def readfile(filename,option='nonerror'):
    xx=[]
    yy=[]
    if (option == 'err'):
        err1=[]
        err2=[]
    with open(filename,'r') as f:
        while True:
            line=f.readline()
            if not line:
                break
            if line.strip().startswith("#"):
                continue
            if (option == 'err'):
                x,y,e1,e2=[float(i) for i in line.split()]
            else:
                x,y=[float(i) for i in line.split()]
            xx.append(x)
            yy.append(y)
            if (option == 'err'):
                err1.append(e1)
                err2.append(e2)
    if (option == 'err'):
        xx,yy,err1,err2=np.array(xx),np.array(yy),np.array(err1),np.array(err2)
        return xx,yy,err1,err2
    else:
        xx,yy=np.array(xx),np.array(yy)
        return xx,yy
    
def xfraction(sfr,mstar,mmin,mmax,binsize,sfr_limit):
    xx=[]
    yy=[]
    m=mmin
    if (max(mstar) < mmax):
        mmax = max(mstar)
    while(m < mmax):
        tmp=sfr[(mstar > m) & (mstar < m+binsize)]
        if(len(tmp)==0):
            continue
#        quench = tmp[tmp<sfr_limit]
#        fq=len(quench)*1.0/len(tmp)
        fq=np.sum(tmp < sfr_limit)*1.0/len(tmp)
        yy.append(fq)
        xx.append(m+binsize/2)
        m=m+binsize
    xx=np.array(xx)
    yy=np.array(yy)
    return xx,yy    


def MassFunction(data,volume,xmin,xmax,binsize):
    m=xmin
    x,y=[],[]
    while (m<xmax):
        x.append(m+binsize/2)
        y.append(np.sum((data>= m) & (data < m+binsize))/binsize/volume) 
        m=m+binsize
    x=np.array(x)
    y=np.array(y)
    return x,y


li_m=[8.32609 ,     8.45942  ,    8.59275  ,    8.72609  ,    8.85942    ,  8.99275  ,    9.12609    ,  9.25942   ,   9.39275 ,     9.52609  ,    9.65942 ,    9.79275,
      9.92609  ,    10.0594    ,  10.1928    ,  10.3261    ,  10.4594   ,   10.5928   ,   10.7261    ,  10.8594  ,    10.9928  ,    11.1261   ,   11.2594    ,  11.3928,
      11.5261  ,    11.6594    ,  11.7928   ,   11.9261    ,  12.0594   ,   12.1928]
li_smf=[0.0172218 ,   0.0164026 ,   0.0155730  ,  0.0147218 ,   0.0138358  ,  0.0129000  ,  0.0118980  ,  0.0108136  , 0.00963336  , 0.00835177  , 0.00708076 ,  0.00706432,
   0.00696476 ,  0.00675671 ,  0.00641269 ,  0.00590743 ,  0.00522585  , 0.00437514 ,  0.00339864 ,  0.00238457,   0.00145687 , 0.000781879 , 0.000370797 , 0.000150015,
  4.89036e-05 , 1.18864e-05 , 1.93838e-06 , 1.83754e-07 , 8.33251e-09 , 1.38666e-10]

bardry_m=[ 6.25000 , 6.75000,7.10000 , 7.30000 ,  7.50000, 7.70000,7.90000,8.10000,8.30000, 8.50000, 8.70000 ,8.90000
           ,9.10000, 9.30000, 9.50000,9.70000,9.90000,10.1000 ,10.3000,10.5000 ,10.7000,10.9000 ,11.1000,11.3000
           ,11.5000,11.7000, 11.9000]
bardry_smf=np.array([0.0311000 ,0.0181000,0.0179000,0.0431000,0.031600, 0.034800,0.027300,0.02830,0.023500,0.019200, 0.018000,
            0.0143000,0.0102000 ,0.0095900, 0.007420, 0.006210,0.0057100,0.005510,0.00548000,0.005120, 0.0035500, 
            0.00241000 ,  0.00127000 , 0.000338000,4.20000e-05 , 2.10000e-05 , 4.20000e-05])
bardry_err=np.array([0.0216000, 0.0066000,0.0057000,0.0087000, 0.00900000, 0.00840000 , 0.00420000 ,  0.00280000 ,  0.00300000,
            0.00120000 ,  0.00260000,   0.00170000, 0.000600000,  0.000550000  ,0.000410000 , 0.000370000 , 0.0003500,
            0.00034000,  0.000340000 , 0.000330000 , 0.0002700, 0.000230000 , 0.000160000 , 8.50000e-05,3.00000e-05,
            2.10000e-05 , 3.00000e-05])

xu2022_m = np.linspace(8.2, 11.6, 18)
xu2022_mf = np.array([-1.36,-1.35,-1.47, -1.62, -1.76, -1.9,-1.99, -2.08, -2.08, -2.14, -2.23, -2.30, -2.43, -2.60, -2.85, -3.20, -3.80, -4.49])
xu2022_mf_err_b = np.array([0.22, 0.15, 0.16, 0.10, 0.15, 0.15, 0.14, 0.11, 0.08, 0.08, 0.05, 0.05, 0.04, 0.03, 0.03, 0.05, 0.08, 0.12])
xu2022_mf_err_t = np.array([0.20, 0.11, 0.12, 0.11, 0.11, 0.14, 0.08, 0.07, 0.06, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.04, 0.05])

# Graham & Scott 2015
Mbulge_g=np.array([   47,   14,    9,   19,   10,     8,  9.1, 1.2, 5.4,  8.4, 14.7, 
                   2.57, 2.29, 1.32, 2.95, 1.26, 2.88, 2.57, 2.14, 1.32, 1.74, 24.2, 
                   16.8, 18.4, 26.6, 66.7, 18.4, 60.8, 15.3, 9.64, 42.1, 15.2, 53.3, 
                   102, 11.7, 32.7, 63.5, 10.3, 27.2, 36.2, 32.1, 164])  * 0.1   # ; 10^10 Msun
Mbh_g   =np.array([1.5e4,8.1e3,1.5e3,1.7e3, 75.0, 2.6e3, 42.6, 2.0, 16., 12.2,  5.1,  
                   5.0,  2.5,  0.8, 12.6,  1.0,  1.6,  2.5,  5.0,  1.3,  1.6, 71.0, 
                   210., 54.0, 45.0,269.0,217.0, 167., 124., 39.0, 100.,  6.0, 39.8, 
                   331, 39.8, 166., 155., 15.5, 135., 195., 490., 1000]) * 1e-5 #; 10^10 Msun
Mbulge_g = np.log10(Mbulge_g)+10
Mbh_g = np.log10(Mbh_g)+10
#Plot,alog10(Mbulge)+10,alog10(Mbh)+10, psym=5,yrange=[4,12],xrange=[8,12]
# Scott 2013
Mbh_s   = np.array([ 39, 11,0.45, 25, 24,0.044, 1.40, 0.73, 9.0, 58, 0.10, 8.3, 0.39, 0.42, 
           0.084, 0.66, 0.73,  15, 4.7, 0.083, 0.14, 0.15, 0.4, 0.12, 1.7, 0.024, 
           8.8, 0.14, 2.0, 0.073,0.77, 4.0,0.17,0.34, 2.4,0.058, 3.1, 1.3, 2.0, 97, 
           8.1, 1.8,0.65, 0.39,  5, 3.3,  4.5,0.075,0.68, 1.2, 0.13, 4.7, 0.59, 6.4,
           0.79, 3.9, 47, 1.8, 0.06, 0.016, 210, 0.014, 7.4, 1.6, 6.8, 2.6, 11, 37,
           5.9,0.31, 0.1,3.7,0.55, 13,0.11]) * 1e-2
Mbulge_s= np.array([ 69, 37, 1.4, 55, 27,  2.4, 0.46,  1.0,  19, 23, 0.61, 4.6,   11,  1.9,   
           4.5,  1.4, 0.66, 4.7,  26,   2.0, 0.39, 0.35, 0.3,  3.5, 6.7,  0.88, 1.9, 
           0.93,1.24,  0.86, 2.0, 5.4, 1.2, 4.9, 2.0, 0.66, 5.1, 2.6, 3.2,100, 1.4,
           0.88, 1.3, 0.56, 29, 6.1, 0.65,  3.3, 2.0, 6.9,  1.4, 7.7, 0.90, 3.9, 1.8, 
           8.4, 27, 6.0, 0.43,   1.0, 122,  0.30,  29,  11,  20, 2.8, 24, 78, 96, 3.6, 
           2.6, 55, 1.4, 64, 1.2])
#oplot,alog10(Mbulge)+10,alog10(Mbh)+10, psym=6,color=2000
Mbh_s = np.log10(Mbh_s)+10
Mbulge_s = np.log10(Mbulge_s)+10

#SMF
data_dir='/home/lzxie/obsdata/SMF_sf_q.csv'
zlist=['0.5','1.0','1.5','2.0','2.5','3.0']
labels=[]
for z in zlist:
    labels.append('sf_z'+z+'_x')
    labels.append('sf_z'+z+'_y')
for z in zlist:
    labels.append('qg_z'+z+'_x')
    labels.append('qg_z'+z+'_y')
obs_smf = pd.read_csv(data_dir,header=2,names=labels)
data_dir='/home/lzxie/obsdata/SMF_all.csv'
zlist=['0.0','0.5','1.0','1.5','2.0','2.5','3.0']
labels_all=[]
for z in zlist:
    labels_all.append('all_z'+z+'_x')
    labels_all.append('all_z'+z+'_y')
obs_smf_all = pd.read_csv(data_dir,header=2,names=labels_all)

def HI_H2_xGASS(gdata='HI'):
    if (gdata =='H2_HI'):
        data = pd.read_csv('/home/lzxie/obsdata/H2_HI_mstar_xGASS.csv')
    elif (gdata =='HI'):
        data = pd.read_csv('/home/lzxie/obsdata/HI_mstar_xGASS.csv')
    elif (gdata =='H2'):
        data = pd.read_csv('/home/lzxie/obsdata/H2_mstar_xGASS.csv')
    return data

#quench fraction
def Fq_sdss(option='all'):
    if (option =='all'):
        obsxx=np.array([ 9.54, 9.67, 9.80, 9.93, 10.06, 10.19, 10.32, 10.45, 10.58, 10.71, 10.84, 10.97, 11.10, 11.23, 11.36, 11.49, 11.61, 11.74, 11.87]) #11.48, 11.60, 11.73, 11.87, 
        obsyy=np.array([0.1726, 0.1980, 0.2284, 0.2589, 0.3198, 0.3909, 0.4670, 0.5431, 0.6091, 0.6599, 0.7056, 0.7563, 0.8223, 0.8680, 0.9289, 0.9645, 0.9442, 0.9949, 1.000]) # , 0.9036, 0.8071, 0.7157, 0.6599]
        obsyy_t=np.array([0.1726, 0.1980, 0.2284, 0.2589, 0.3198, 0.3909, 0.4670, 0.5431, 0.6091, 0.6599, 0.7056, 0.7563, 0.8223, 0.8480, 0.9089, 0.9817, 0.9870, 1.000, 1.000]) #, 0.9036, 0.8071, 0.7157, 0.6599]
        obsyy_b=np.array([0.1726, 0.1980, 0.2284, 0.2589, 0.3198, 0.3909, 0.4670, 0.5431, 0.6091, 0.6599, 0.7056, 0.7563, 0.8223, 0.8880, 0.9509, 0.90517, 0.8017, 0.716, 0.664]) #, 0.9036, 0.8071, 0.7157, 0.6599]
    elif (option == 'central'):
        obsxx=np.array([ 9.54, 9.67, 9.80, 9.93, 10.06, 10.19, 10.32, 10.45, 10.58, 10.71, 10.84, 10.97, 11.10, 11.23, 11.36, 11.49, 11.61, 11.74, 11.87]) #11.48, 11.60, 11.73, 11.87, 
        obsyy=np.array([0.06588,0.09014,0.1398,0.1640,0.2239,0.3093,0.4000,0.5000,0.5604,0.6305,0.6800,0.7400,0.8000,0.8653,0.9267,0.9565,0.9469,0.9925,0.9948]) #,0.8881,0.7694,0.6660,0.6036]
        obsyy_t=np.array([0.06588,0.09014,0.1398,0.1640,0.2239,0.3093,0.4000,0.5000,0.5604,0.6305,0.6800,0.7400,0.8000,0.8870,0.9563,0.842,0.9923,1.000,1.000]) #,0.8881,0.7694,0.6660,0.6036]
        obsyy_b=np.array([0.06588,0.09014,0.1398,0.1640,0.2239,0.3093,0.4000,0.5000,0.5604,0.6305,0.6800,0.7400,0.7960,0.8436,0.9010,0.8994,0.7852,0.6927,0.6318])
    elif (option == 'satellite'):
        obsxx=np.array([9.535,9.658,9.790,9.931,10.05,10.19,10.31,10.44,10.58,10.71,10.83,10.95,11.07,11.22,11.35,11.48]) #,11.22,11.35,11.47]
        obsyy=np.array([0.3250, 0.3657, 0.3703, 0.4153, 0.4713, 0.5116, 0.5778, 0.6079, 0.6479, 0.6933, 0.7139, 0.7495, 0.7902, 0.8046, 0.8041, 0.8954]) #, 0.7536, 0.6919, 0.6509]
        obsyy_t=np.array([0.3250, 0.3657, 0.3703, 0.4153, 0.4713, 0.5116, 0.5778, 0.6079, 0.6479, 0.6933, 0.7139, 0.7495, 0.7902, 0.8046, 0.8041, 0.8954]) #, 0.7536, 0.6919, 0.6509]
        obsyy_b=np.array([0.3250, 0.3657, 0.3703, 0.4153, 0.4713, 0.5116, 0.5778, 0.6079, 0.6479, 0.6933, 0.7139, 0.7495, 0.7902, 0.8046, 0.8041, 0.8954]) #, 0.7536, 0.6919, 0.6509]
        
    obsyy_b = obsyy -obsyy_b
    obsyy_t = obsyy_t-obsyy
    return obsxx,obsyy, obsyy_b, obsyy_t


ytitle = ['log (MHI/M*)','log (MH2/M*)','log (MH2/MHI)']
gtype =['all','central','satellite']
colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][:6]
sfr_limit = -10.6

datadir = '/home/cossim/Jiutian/M300/GAEA-SAM/'
fname = 'Gal_testpara_'

zlist = ['0.0','0.5','1.0']#,'2.0']#,'3.0']
snaplist = [127,98, 83,]# 66]#, 56]

filelist = [18,19,20,21]
filenr = 467
Nfile = 2
firstfile=36
volume = (300.0/0.677)**3 /1000 * Nfile
ifile = 6

iz = 0
nrow = 6
nline = 3
fig,ax = plt.subplots(nrow,nline,figsize=(nline*6+1,nrow*6))
for z,snap in zip(zlist,snaplist):
    ddir = datadir+'%03d/' % (snap)
    

    ddir = datadir+'%03d/' % (snap)
    for ifile in filelist:
        for filenr in range(Nfile):
            filename = ddir + fname +'%d_%03d_%03d.hdf5' % (ifile, snap, filenr+firstfile)
            with h5py.File(filename,'r') as f:
                gals = f['Galaxies']
                gals = gals[(gals['StellarMass'] > 1e-3)]
                if (filenr ==0):
                    allgal = gals
                else:
                    allgal = np.concatenate((allgal, gals),axis = None)

        mstar = np.log10(allgal['StellarMass']/0.677)+10
    
        x,y = MassFunction(mstar,volume,7,12,0.2)
        y = np.log10(y)
        ax[0,iz].plot(x,y,linewidth=3,alpha=0.6, label='M300 '+str(ifile),color = colors[ifile-filelist[0]])
        
        ####################   
        mhi=np.log10((allgal['ColdGas']-allgal['H_mol'])*0.74/0.677)+10
        x_hi,y_hi = MassFunction(mhi,volume,6,12,0.2)
        ax[1,iz].plot(x_hi,np.log10(y_hi),linewidth=2,linestyle='solid', label='M300 '+str(ifile),color = colors[ifile-filelist[0]])
        #####################
        if (z=='0.0'):
            mh2 = np.log10(allgal['H_mol']*0.74/0.677)+10
            ggtype = allgal['Type']
        
            sel1 = np.argwhere((mstar > 10) & (mh2 - mstar > -1.6)).squeeze()
            sel2 = np.argwhere((mstar < 10) & (mh2 - mstar > -1.8)).squeeze()
            select = np.concatenate((sel1,sel2),axis = None)
            mhi = mhi[select]
            mh2 = mh2[select]
            mstar = mstar[select]
            ggtype = ggtype[select]

            for igas,yt in enumerate(ytitle):
                if (yt == 'log (MHI/M*)'):
                    ss = mhi -mstar
                elif (yt =='log (MH2/M*)'):
                    ss = mh2 -mstar
                elif (yt == 'log (MH2/MHI)'):
                    ss = mh2-mhi

                for j,gt in enumerate(gtype):
                    if (gt == 'all'):
                        x,y,err = scalingrelation(mstar,ss,np.linspace(8,12),option='median')

                    elif (gt =='central'):
                        x,y,err = scalingrelation(mstar[(ggtype ==0)], ss[(ggtype ==0)], np.linspace(8,12,10))

                    elif ( gt =='satellite'):
                        x,y,err = scalingrelation(mstar[(ggtype !=0)], ss[(ggtype !=0)], np.linspace(8,12,10))
                        
                    ax[igas+2,j].plot(x,y, label = 'model, median'+str(ifile),color = colors[ifile-filelist[0]],linewidth=2)
                    


            mstar = np.log10(allgal['StellarMass']/0.677)+10
            sfr = np.log10(allgal['Sfr']) - mstar
            sfr[(allgal['Sfr'] == 0 )] = -12
            Gtype = allgal['Type']

            for i,option in enumerate(['all','central','satellite']):
                if (option =='all'):
                    x,y = xfraction(sfr,mstar,7,12,0.5,sfr_limit)
                elif (option =='central'):
                    x,y = xfraction(sfr[(Gtype == 0)], mstar[(Gtype == 0)],7,12,0.5,sfr_limit)
                elif (option == 'satellite'):
                    x,y = xfraction(sfr[(Gtype != 0)], mstar[(Gtype != 0)],7,12,0.5,sfr_limit)

                ax[5,i].plot(x,y,linewidth = 2, label='M300 ' + str(ifile),color = colors[ifile-filelist[0]])



    obsxname = 'all_z'+z+'_x'
    obsyname = 'all_z'+z+'_y'
    xobs = obs_smf_all[obsxname].values
    yobs = obs_smf_all[obsyname].values
    xobs = xobs[np.isfinite(yobs)]
    yobs = yobs[np.isfinite(yobs)]
    ax[0,iz].plot(xobs,yobs,alpha=0.6,linewidth = 4,color='k',label='Obs')
    
    if (z == '0.0'):
        plt.scatter(li_m,np.log10(li_smf),label='Li et al. 2009 SDSS',linewidth=4,alpha=0.6,color='k')
        ax[0,iz].fill_between(xu2022_m, xu2022_mf-xu2022_mf_err_b, xu2022_mf + xu2022_mf_err_t, color='cyan', alpha=0.3, label = 'Xu & Jing 2022')

    if (z=='2.0'):
        fontx = np.array([9.88300     , 10.0830    , 10.2830     , 10.4830    , 10.6830    ,  10.8830     , 11.0830     , 11.2830     , 11.4830])
        fonty = np.array([0.00707946  , 0.00199526 , 0.000891251 , 0.00104713 , 0.000758578,  0.000338844 , 0.000354813 , 0.000109648 , 2.75423e-05])
        ax[0,iz].scatter(fontx,np.log10(fonty),label='Font 2006', linewidth=4,alpha=0.6)
    #plt.ylim(-5,-1)
    ax[0,iz].set_xlabel('$log(M_{\star}/M_{\odot})$',fontsize=12)
    ax[0,iz].set_ylabel('$\log \Phi /(\log M_{\odot} Mpc^3)$',fontsize=12)
    ax[0,iz].set_title('z=~'+z)
    ax[0,iz].legend()
    #plt.savefig('figure/SMF_z'+z+'.pdf')
    
    #HIMF
    obs_himf(zlist[i],i)
    if (z=='0.0'):
        xx=np.linspace(6,12,20)
        #ax[0].plot(xx,np.log10(jones2018_himf(xx)),label='Jones et al. 2018',color='k',linestyle=':',linewidth=3)
        ax[1,iz].plot(xx+math.log10(0.7)*2,np.log10(jones2018_himf(xx)/(0.7*0.7*0.7)),
                 label='Jones et al. 2018',color='k',linestyle=':',linewidth=3)
    ax[1,iz].set_xlim(8,12)
    ax[1,iz].set_ylim(-6,0)
    ax[1,iz].set_xlabel('$M_{HI}/M_{\odot}$',fontsize=12)
    ax[1,iz].set_ylabel('$\Phi /Mpc^3/dM_{HI}$',fontsize=12)
    ax[1,iz].legend()
    
    iz=iz+1
    print(iz)

    
                #ax[i,j].fill_between(x,y-err[0,:],y+err[1,:],color = colors[icolor],alpha= 0.5)
                
for igas,yt in enumerate(ytitle):
    
    if (yt == 'log (MHI/M*)'):
        data = HI_H2_xGASS(gdata='HI')
    elif (yt =='log (MH2/M*)'):
        data = HI_H2_xGASS(gdata = 'H2')
    elif (yt == 'log (MH2/MHI)'):
        data = HI_H2_xGASS(gdata = 'H2_HI')
            
    obsx = data['Mstar']
    for j,gt in enumerate(gtype):
        if (gt == 'all'):
            obsy = data['mean_all']
        elif (gt =='central'):
            obsy = data['mean_central']
        elif ( gt =='satellite'):
            obsy = data['mean_satellite']
        ax[igas+2,j].plot(obsx,obsy,label = 'xGASS', color= 'k', linestyle='dotted',linewidth=4)
        ax[igas+2,j].set_ylabel(yt)
        ax[igas+2,j].set_xlabel('log $M*/M_{\odot}$')
    
    ax[igas+2,j].set_title(gt)


for i,option in enumerate(['all','central','satellite']):
    obsxx,obsyy,obse1,obse2 = Fq_sdss(option)
    ax[5,i].errorbar(obsxx,obsyy,[obse1,obse2], label='SDSS',marker = 'o')
    ax[5,i].set_xlabel('$M_{\star}/M_{\odot}$',fontsize=16)
    if (i==0):
        ax[5,i].set_ylabel('$f_q$',fontsize=16)
    ax[5,i].set_title(option)
    ax[5,i].set_ylim(0.,1.0)


fig.tight_layout()
plt.savefig('figure/300.pdf')