import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import FormatStrFormatter
import mpl_toolkits.axisartist as axisartist


import matplotlib
print(matplotlib.__version__)

#x,y,z = np.loadtxt("borex_chi3.dat",unpack=True)

#plt.show()


#N = 1000
#_y = np.log10(y)
#_x = np.log10(x)
#chi2_min = np.min(z)
#
#z = z - chi2_min
#
#xi = np.linspace(_x.min(),_x.max(),N)
#yi = np.linspace(_y.min(),_y.max(),N)
#
#zi = scipy.interpolate.griddata((_x,_y),z,(xi[None,:],yi[:,None]))
#
#fig,ax = plt.subplots()
#ax.set_xlabel(r"$tan^{2}\theta_{12}$",size=20)
#ax.set_ylabel(r"$\Delta m^{2}_{21}$ eV$^{2}$",size=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.contour(10**xi,10**yi,zi,levels=[2.61,5.99,10.0],colors=['r','b','k'])
#plt.show()


class data_set:
    def __init__(self,_data,_level):
        self.data = _data
        self.level = _level

        self.prepare_data()


    def prepare_data(self):
        x = self.data[0]
        y = self.data[1]
        z = self.data[2]

        N = 1000
        _x = np.log10(x)
        _y = np.log10(y)

        chi2_min = np.min(z)
        print(chi2_min)

        z = z - chi2_min

        xi = np.linspace(_x.min(),_x.max(),N)
        yi = np.linspace(_y.min(),_y.max(),N)

        zi = interp.griddata((_x,_y),z,(xi[None,:],yi[:,None]))

        cs = plt.contour(10**xi,10**yi,zi,[self.level])
        p = cs.collections[0].get_paths()[0]

        v = p.vertices
        self.X = v[:,0]
        self.Y = v[:,1]

        return 0


data = np.loadtxt("borex_chi2.dat",unpack=True)


#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

cs_data_3s = data_set(data,11.83)
cs_data_95 = data_set(data,5.99)
cs_data_1s = data_set(data,2.61)

fig = plt.figure(figsize=(5,4))
#ax = fig.add_subplot(axes_class=axisartist.Axes)

ax = fig.add_subplot()

plt.text(0.1,0.8,"Borexino",transform=ax.transAxes,font="Times New Roman",fontsize=20)
plt.text(0.1,0.7,"No decay",transform=ax.transAxes,font="Times New Roman",fontsize=15)


ax.set_ylabel(r"$\Delta m^{2}_{21}$ [eV$^{2}$]",fontsize=20)
ax.set_xlabel(r"tan$^{2}\theta_{12}$",fontsize=20)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-2,1)
ax.set_ylim(1e-6,1e-4)
ax.minorticks_on()

ax.plot(cs_data_3s.X,cs_data_3s.Y,color='k',linewidth=1.0,label='$3 \sigma$ C.L.')
ax.plot(cs_data_95.X,cs_data_95.Y,color='r',linewidth=1.0,label='95 % C.L.')
ax.plot(cs_data_1s.X,cs_data_1s.Y,color='b',linewidth=1.0,label=r'$1 \sigma$ C.L.')
ax.tick_params(axis='x', which='minor',direction='in')
ax.tick_params(axis='y', which='minor',direction='in')
ax.tick_params(axis='x', which='major',direction='in')
ax.tick_params(axis='y', which='major',direction='in')

plt.legend(loc=3,edgecolor=None)
#plt.show()
plt.savefig("Borex_decall.png",dpi=300,bbox_inches='tight')

