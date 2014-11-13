from common import *

#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)


list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()

nkpath = 500

kpath = Path(acell,nkpath)

HW = []
HWr= []


path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/modules/'
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

for q in kpath.list_k:
    list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q)

    HW.append(list_hw)
    HWr.append(list_hw_renormalized)

HW = N.array(HW)
HWr= N.array(HWr)

#================================================================================

fig = plt.figure(figsize=(14,10))

ax  = fig.add_subplot(111)

First = True
for hw , hwr in zip(HW.transpose(),HWr.transpose()):

    if First:
            label1 = 'Dubay and Kresse (2003)'
            label2 = 'Renormalized'
            First  = False
    else:
            label1 = '__nolabel__'
            label2 = '__nolabel__'

    ax.plot(kpath.list_x,1000*hwr,'b--',lw=4,label=label2)
    ax.plot(kpath.list_x,1000*hw,'r-',lw=4,label=label1)


ax.set_title(r'phonon dispersion of Graphene')
ax.set_ylabel(r'$\hbar\omega{\bf q}$ (meV)')

ax.set_xticks(kpath.list_xticks)
ax.set_xticklabels(kpath.list_labels)

ax.set_xlim(kpath.list_x[0],kpath.list_x[-1])


L = ax.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax.grid(True,linestyle='-',color='black',axis='x',alpha=1.)


emin = -10.
emax = 220.

ax2  = ax.twinx()

ax2.grid(True,linestyle='-',color='grey',axis='y',alpha=0.5)
ax.set_ylim([emin,emax])

yticks = N.arange(0.,1800.,100.)

ax2.set_yticks(yticks)

ax2.set_ylim([eV_to_cm1*emin/1000,eV_to_cm1*emax/1000])

ax2.set_ylabel(r'$\hbar\omega{\bf q}$ (cm$^{-1}$)')


fig.subplots_adjust(    left    =       0.15,
                        bottom  =       0.10,
                        right   =       0.85,
                        top     =       0.92,
                        wspace  =       0.50,
                        hspace  =       0.50)
plt.show()

