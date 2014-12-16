import os
import numpy as N
from Scientific.IO.NetCDF import NetCDFFile as Dataset

#================================================================================
# Define some job bundles, which contain all the necessary info to plot them
#================================================================================

jobs_bundle = []

dic_job = {}
mu     =-0.400
Gamma  = 0.050
nq = (32,128,5)
nk = (8,256,2)
dic_job['nq']                  =  nq
dic_job['nk']                  =  nk
dic_job['mu']                  =  mu
dic_job['Gamma']               =  Gamma
dic_job['data_dir']            = 'Gamma=%i_meV_mu=%i_meV_nq=%i_%i_%i_nk=%i_%i_%i/'%(1000*Gamma,1000*mu,nq[0],nq[1],nq[2],nk[0],nk[1],nk[2])
dic_job['filename_minus_tmpl'] = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i'%(nq[0],nq[1],nq[2],nk[0],nk[1],nk[2])+'_iq_minus=%i_nu=%i.nc'
dic_job['filename_plus_tmpl']  = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i'%(nq[0],nq[1],nq[2],nk[0],nk[1],nk[2])+'_iq_plus=%i_nu=%i.nc'
dic_job['combined_filename']   = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i_mu=%i_meV_Gamma=%i_meV'%(nq[0],nq[1],nq[2],nk[0],nk[1],nk[2],1000*mu,1000*Gamma)+'.nc'
jobs_bundle.append(dic_job)





error_file = open('error_file_mu=%i.txt'%(1000*mu),'w')
log_file   = open('convert_mu=%i.log'%(1000*mu),'w')

print >> error_file, " The following files had problems in them"
error_file.flush()


list_nu = [2,3,4,5]

First = True

for dic_job in jobs_bundle:

        data_dir = dic_job['data_dir']
        Gamma    = dic_job['Gamma']

        filename_minus_tmpl = data_dir+dic_job['filename_minus_tmpl']
        filename_plus_tmpl  = data_dir+dic_job['filename_plus_tmpl']

        list_iq = []
        for filename in os.listdir(data_dir):
                if 'iq_plus' in filename:
                        iq = int(filename.strip().split('=')[3][:-3])
                        list_iq.append(iq)

        number_of_q = N.array(list_iq).max()

        list_iq = N.arange(1,number_of_q +1)

        Nph = len(list_iq)


        # Extract basic properties
        iq = 1
        nu = 2

        filename = filename_plus_tmpl%(iq,nu)

        ncfile   = Dataset(filename,'r')
        list_hw  = ncfile.variables['list_hw'][:]

        Area  = ncfile.Area
        acell = ncfile.acell
        mu    = ncfile.mu
        Gamma = ncfile.Gamma_width
        beta  = ncfile.beta
        knmax_coarse = ncfile.nmax_coarse
        knmax_fine   = ncfile.nmax_fine
        knmax_block  = ncfile.n_blocks_coarse_to_fine


        ncfile.close()


        combined_filename = dic_job['combined_filename']


        #--------------------------------------------
        # Write to netcdf file 
        #--------------------------------------------
        Cncfile   = Dataset(combined_filename,'w')

        # --- set various attributes, identifying the parameters of the computation ----
        setattr(Cncfile,'mu',mu) 
        setattr(Cncfile,'beta',beta) 
        setattr(Cncfile,'acell',acell) 
        setattr(Cncfile,'Area',Area) 
        setattr(Cncfile,'Gamma_width',Gamma) 


        setattr(Cncfile,'k_nmax_coarse',knmax_coarse)
        setattr(Cncfile,'k_nmax_fine',knmax_fine)
        setattr(Cncfile,'k_nmax_block',knmax_block)

        setattr(Cncfile,'q_nmax_coarse',nq[0])
        setattr(Cncfile,'q_nmax_fine',nq[1])
        setattr(Cncfile,'q_nmax_block',nq[2])



        # --- Create dimensions ----
        Cncfile.createDimension("number_of_frequencies",list_hw.shape[0])
        Cncfile.createDimension("xy",2)
        Cncfile.createDimension("L",2)
        Cncfile.createDimension("uv",2)
        Cncfile.createDimension("phonon_alpha_kappa",6)
        Cncfile.createDimension("nu",len(list_nu))
        Cncfile.createDimension("number_of_q_points",Nph)


        # --- Write data ----
        HW     = Cncfile.createVariable("list_hw",'d',('number_of_frequencies',))
        HW[:]  = N.real(list_hw)

        NU     = Cncfile.createVariable("phonon_nu",'d',('nu',))
        NU[:]  = list_nu


        Qplus  = Cncfile.createVariable( "list_q_plus",'d',('number_of_q_points','xy'))
        Qminus = Cncfile.createVariable("list_q_minus",'d',('number_of_q_points','xy'))


        HWPH   = Cncfile.createVariable("phonon_frequencies",'d',('nu','number_of_q_points'))

        REPHp  = Cncfile.createVariable("Re_E_phonon_plus",'d',('nu','number_of_q_points','phonon_alpha_kappa',))
        IEPHp  = Cncfile.createVariable("Im_E_phonon_plus",'d',('nu','number_of_q_points','phonon_alpha_kappa',))

        REPHm  = Cncfile.createVariable("Re_E_phonon_minus",'d',('nu','number_of_q_points','phonon_alpha_kappa',))
        IEPHm  = Cncfile.createVariable("Im_E_phonon_minus",'d',('nu','number_of_q_points','phonon_alpha_kappa',))


        RHp    = Cncfile.createVariable("Re_H_plus",'d',('nu','number_of_q_points','xy','L','uv','number_of_frequencies'))
        IHp    = Cncfile.createVariable("Im_H_plus",'d',('nu','number_of_q_points','xy','L','uv','number_of_frequencies'))

        RHm    = Cncfile.createVariable("Re_H_minus",'d',('nu','number_of_q_points','xy','L','uv','number_of_frequencies'))
        IHm    = Cncfile.createVariable("Im_H_minus",'d',('nu','number_of_q_points','xy','L','uv','number_of_frequencies'))

        #REPH[:] = N.real(CS.E_ph)
        #IEPH[:] = N.imag(CS.E_ph)
#
#

        First = True
        for inu, nu in enumerate(list_nu):
                print >> log_file, 'nu = %i'%nu
        log_file.flush()

                for iq in list_iq:
                    print >> log_file, '\t iq = %i of %i'%(iq, len(list_iq))
            log_file.flush()

                        # careful! iq starts at 1; indices should start at zero


                        #=============== q plus ======================================= 
                        filename = filename_plus_tmpl%(iq,nu)

            try:
                ncfile   = Dataset(filename,'r')


                HWPH[N.int(inu),N.int(iq-1)]  = ncfile.phonon_frequency[0]

                RHp[N.int(inu),N.int(iq-1),:] = ncfile.variables['Re_H'][:]
                IHp[N.int(inu),N.int(iq-1),:] = ncfile.variables['Im_H'][:]

                REPHp[N.int(inu),N.int(iq-1),:] = ncfile.variables['Re_E_phonon'][:]
                IEPHp[N.int(inu),N.int(iq-1),:] = ncfile.variables['Im_E_phonon'][:]

                if First:
                    Qplus[N.int(iq-1),:]  = ncfile.variables['q_phonon'][:]


                            ncfile.close()

            except:
                print >> error_file , filename 
                error_file.flush()



                        #=============== q minus ======================================= 
                        filename = filename_minus_tmpl%(iq,nu)
            try:
                ncfile   = Dataset(filename,'r')

                RHm[N.int(inu),N.int(iq-1),:] = ncfile.variables['Re_H'][:]
                IHm[N.int(inu),N.int(iq-1),:] = ncfile.variables['Im_H'][:]

                REPHm[N.int(inu),N.int(iq-1),:] = ncfile.variables['Re_E_phonon'][:]
                IEPHm[N.int(inu),N.int(iq-1),:] = ncfile.variables['Im_E_phonon'][:]



                if First:
                    Qminus[N.int(iq-1),:]  = ncfile.variables['q_phonon'][:]



                            ncfile.close()
            except:
                print >> error_file , filename 
                error_file.flush()




                First = False


        Cncfile.close()

print >> error_file, "END of error file"
error_file.close()
log_file.close()
#================================================================================
