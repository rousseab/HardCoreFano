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

def write_splmake(      re_SR_spline, im_SR_spline, re_SI_spline, im_SI_spline,
                        re_TR_spline, im_TR_spline, re_TI_spline, im_TI_spline,
                        filename, spline_order, list_hw_ext, mu, 
                        beta, kernel_Gamma_width, Green_Gamma_width)
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
    setattr(ncfile,'kernel_Gamma_width',kernel_Gamma_width) 
    setattr(ncfile,'Green_Gamma_width',Green_Gamma_width) 

    setattr(ncfile,'spline_order',spline_order) 

    # --- Create dimensions ----
    dim_spline = len(re_SR_spline)
    dim_hw_ext = len(list_hw_ext)
    dim_eta    = 2

    ncfile.createDimension("dim_eta",dim_eta)
    ncfile.createDimension("dim_hw_ext",dim_hw_ext)
    ncfile.createDimension("dim_spline",len(list_hw_ext))

    # Create variables
    HW = ncfile.createVariable("list_hw_ext",'d',('dim_hw_ext',))
    HW[:] = list_hw_ext

    RE_SR = ncfile.createVariable("re_SR_spline",'d',('dim_spline',))
    RE_SR[:] =  re_SR_spline 

    IM_SR = ncfile.createVariable("im_SR_spline",'d',('dim_spline',))
    IM_SR[:] =  im_SR_spline 

    RE_SI = ncfile.createVariable("re_SI_spline",'d',('dim_spline',))
    RE_SI[:] =  re_SI_spline 

    IM_SI = ncfile.createVariable("im_SI_spline",'d',('dim_spline',))
    IM_SI[:] =  im_SI_spline 

    RE_TR = ncfile.createVariable("re_TR_spline",'d',('dim_eta','dim_hw_ext','dim_spline',))
    RE_TR[:] =  re_TR_spline 

    IM_TR = ncfile.createVariable("im_TR_spline",'d',('dim_eta','dim_hw_ext','dim_spline',))
    IM_TR[:] =  im_TR_spline 

    RE_TI = ncfile.createVariable("re_TI_spline",'d',('dim_eta','dim_hw_ext','dim_spline',))
    RE_TI[:] =  re_TI_spline 

    IM_TI = ncfile.createVariable("im_TI_spline",'d',('dim_eta','dim_hw_ext','dim_spline',))
    IM_TI[:] =  im_TI_spline 

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

    file_kernel_Gamma_width = ncfile.kernel_Gamma_width
    file_Green_Gamma_width  = ncfile.Green_Gamma_width

    spline_order = ncfile.spline_order[0]

    list_keys = ['Re_fKR', 'Im_fKR', 'Re_dfKR', 'Im_dfKR', 'Re_fKI', 'Im_fKI', 'Re_dfKI', 'Im_dfKI']

    splmake_tuple_dict = {}

    for key in list_keys:
        D1 = ncfile.variables['%s_array1'%key][:]
        D2 = ncfile.variables['%s_array2'%key][:]

        splmake_tuple = (D1,D2,spline_order)

        splmake_tuple_dict[key] = splmake_tuple 

    return splmake_tuple_dict, file_mu, file_beta, file_kernel_Gamma_width, file_Green_Gamma_width  

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

    return

