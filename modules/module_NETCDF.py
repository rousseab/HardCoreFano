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

def write_to_file(CS,nmax_coarse, nmax_fine, nblocks ,hw_ph,filename):
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
    setattr(ncfile,'Gamma_width',CS.Gamma) 
    setattr(ncfile,'phonon_frequency',hw_ph) 


    # --- Create dimensions ----
    ncfile.createDimension("number_of_frequencies",CS.list_hw.shape[0])
    ncfile.createDimension("xy",2)
    ncfile.createDimension("L",2)
    ncfile.createDimension("uv",2)
    ncfile.createDimension("phonon_alpha_kappa",6)


    # --- Write data ----
    Q      = ncfile.createVariable("q_phonon",'d',('xy',))
    REPH   = ncfile.createVariable("Re_E_phonon",'d',('phonon_alpha_kappa',))
    IEPH   = ncfile.createVariable("Im_E_phonon",'d',('phonon_alpha_kappa',))
    HW     = ncfile.createVariable("list_hw",'d',('number_of_frequencies',))

    RH     = ncfile.createVariable("Re_H",'d',('xy','L','uv','number_of_frequencies'))
    IH     = ncfile.createVariable("Im_H",'d',('xy','L','uv','number_of_frequencies'))


    Q[:]    = CS.q
    REPH[:] = N.real(CS.E_ph)
    IEPH[:] = N.imag(CS.E_ph)
    HW[:]   = N.real(CS.list_hw)

    RH[:,:,:,:] = N.real(CS.Hq)
    IH[:,:,:,:] = N.imag(CS.Hq)



    ncfile.close()
