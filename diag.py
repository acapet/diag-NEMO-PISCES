"""
    The script will list all instances of *ptrc*.nc files and issue additional *diag*nc files with same filename structure and attributes.
    
    Args:
        dir (str): directoy of NMO-PISCES outputs (default : './')
        diaglist (list of str): identifiers of required diagnostics. Each should correspond to an entry in the diagnsotic dictionnary. 
""" 
   # The main arguments is a directory.
# The script will list all instances of *ptrc*.nc files and issue additional *diag*nc files with same filename structure and attributes.


import argparse
import os
from glob import glob
import xarray as xr
import DiagFunctions_NEMOPISCES as diag

parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, default='./')
parser.add_argument('--diaglist', nargs='+', default=['TPPI','pocF200','TrophicEfficiencyI','zchlmax','nitracline','ratioLargeM'])
args = parser.parse_args()

indir =args.dir

# List of diagnostics to be computed. Kept as hard-coded until further requirements are defined. 
dlist=args.diaglist

print('Selected Diags : ')
diag.diaglist(dlist)


flist_p = glob(indir + '*ptrc*nc')
flist_d = [ f.replace('ptrc','diad') for f in flist_p ] 
flist_g = [ f.replace('ptrc','gridT') for f in flist_p ] 
flist_o = [ f.replace('ptrc','diag') for f in flist_p ] 


for i, (fp, fd, fg, fo) in enumerate(zip(flist_p, flist_d, flist_g, flist_o) ):
    x_p=xr.load_dataset(fp)
    x_d=xr.load_dataset(fd)
    x_g=xr.load_dataset(fg)

    #FIXME TPP is now provided as time instant.
    # For now I just overwrite and assume time_centered coordinates instead.
    x_d=x_d.assign(time_counter=x_g['time_counter'].data )
    x_d=x_d.assign(time_counter_bounds=(('time_counter', 'axis_nbounds'), x_g['time_counter_bounds'].data) )

    # Ensure we got all we may need
    x_a=xr.merge([x_p,x_g,x_d])

    # Need to define a cell height variable to use xgcm 
    x_a['h']=(x_a['deptht_bounds'][:,1:]-x_a['deptht_bounds'][:,:-1]).squeeze()
    x_a['h'].attrs={'units'     : 'm',
              'long_name' : 'cells height',
              'valid_min' : -1e20,
              'valid_max' : 1e20,
              'cell_methods' : 'time: mean',
              'coordinates': 'lon lat deptht'}

    bibi=diag.add2D(x_a,dlist, verbose=False)
    bibi[dlist].to_netcdf(fo)

    print(fo + ' completed')
