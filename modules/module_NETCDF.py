#================================================================================
#
#       Module NETCDF
#       ===============
#
#       This module implements the necessary mechanics to write the data to
#       netCDF files!
#
#================================================================================
from module_Constants import *

from Scientific.IO.NetCDF import NetCDFFile as Dataset

def write_splmake(splmake_tuple,filename,mu,beta,delta):
    """
    Write the spline parameters to a netcdf file
    """
    #--------------------------------------------
    # Write to netcdf file 
    #--------------------------------------------
    ncfile   = Dataset(filename,'w')

    # --- set various attributes, identifying the parameters of the computation ----
    setattr(ncfile,'mu',mu) 
    setattr(ncfile,'beta',beta) 
    setattr(ncfile,'delta',delta) 

    setattr(ncfile,'d3',splmake_tuple[2]) 

    # --- Create dimensions ----
    d1 = splmake_tuple[0].shape[0]
    d2 = splmake_tuple[1].shape[0]
    ncfile.createDimension("d1",d1)
    ncfile.createDimension("d2",d2)
 
    D1 = ncfile.createVariable("array1",'d',('d1',))
    D2 = ncfile.createVariable("array2",'d',('d2',))

    D1[:]    = N.real(splmake_tuple[0])
    D2[:]    = N.real(splmake_tuple[1])

    ncfile.close()

    return

def read_splmake(filename):
    """
    Read the spline parameters from a netcdf file
    """
    #--------------------------------------------
    # read to netcdf file 
    #--------------------------------------------
    ncfile   = Dataset(filename,'r')

    file_mu     = ncfile.mu
    file_beta   = ncfile.beta
    file_delta  = ncfile.delta

    d3 = ncfile.d3

    D1 = ncfile.variables['array1'][:]
    D2 = ncfile.variables['array2'][:]

    splmake_tuple = (D1,D2,d3)

    return splmake_tuple, file_mu, file_beta, file_delta

def write_to_file(CS, nmax_coarse, nmax_fine, nblocks, hw_ph,filename):

    """
    Writes the results of a computation to a netcdf file.
    Takes a Compute_Loop_Function object as input; it is assumed that 
    this object has already computed what we wish to write!
    """

    #--------------------------------------------
    # Write to netcdf file 
    #--------------------------------------------
    ncfile   = Dataset(filename,'w')

    # --- set various attributes, identifying the parameters of the computation ----
    setattr(ncfile,'mu',CS.mu) 
    setattr(ncfile,'beta',CS.beta) 
    setattr(ncfile,'acell',acell) 
    setattr(ncfile,'Area',Area) 
    setattr(ncfile,'nmax_coarse',nmax_coarse) 
    setattr(ncfile,'nmax_fine',nmax_fine) 
    setattr(ncfile,'n_blocks_coarse_to_fine',nblocks) 
    setattr(ncfile,'Green_Gamma_width',CS.Green_Gamma_width) 
    setattr(ncfile,'kernel_Gamma_width',CS.kernel_Gamma_width) 
    setattr(ncfile,'phonon_frequency',hw_ph) 

    # --- Create dimensions ----
    ncfile.createDimension("xy",2)
    ncfile.createDimension("L_AB",2)
    ncfile.createDimension("phonon_alpha_kappa",6)

    # --- Write data ----
    Q      = ncfile.createVariable("q_phonon",'d',('xy',))
    REPH   = ncfile.createVariable("Re_E_phonon",'d',('phonon_alpha_kappa',))
    IEPH   = ncfile.createVariable("Im_E_phonon",'d',('phonon_alpha_kappa',))

    Re_R = ncfile.createVariable("Re_R",'d',('xy','L_AB'))
    Im_R = ncfile.createVariable("Im_R",'d',('xy','L_AB'))
    Re_I = ncfile.createVariable("Re_I",'d',('xy','L_AB'))
    Im_I = ncfile.createVariable("Im_I",'d',('xy','L_AB'))

    Q[:]    = CS.q
    REPH[:] = N.real(CS.E_ph)
    IEPH[:] = N.imag(CS.E_ph)

    Re_R[:,:] = N.real(CS.Rq)
    Im_R[:,:] = N.imag(CS.Rq)
    Re_I[:,:] = N.real(CS.Iq)
    Im_I[:,:] = N.imag(CS.Iq)

    ncfile.close()


