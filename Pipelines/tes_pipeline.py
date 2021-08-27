from os import scandir
from sys import path
import warnings
import numpy as np
import pandas as pd
import xarray as xr
warnings.filterwarnings("ignore", category=FutureWarning, message=r'.+resample')
warnings.filterwarnings("ignore", category=RuntimeWarning, message='Mean of')
file_prefix = 'EXD/TES_data/'

attrs = {'latitude': {'LongName':'''Latitude coordinate in geodetic \
                                    cartesian coordinate system''',
                      'Units':'degrees'},
         'longitude': {'LongName':'''Longitude coordinate in geodetic \
                                     cartesian coordinate system''',
                      'Units':'degrees'},
         'date': {'LongName':'Datetime object containing month and year',
                  'Units':'month'},
         'surface_alt': {'LongName':'Altitude at planet surface',
                         'Units':'meters'},
         'species': {'LongName':'''Average CO2 volume mixing ratio for \
                                   range of altitudes''',
                     'Units':'volume mixing ratio relative to dry air'}}

def int_to_datetime(xarr):
    '''
    Extraxt month and year from date as int
    '''
    date = xarr['YYYYMMDD'].values[0].astype(str)
    date_fmt = np.vectorize(lambda x: f'{x[:4]}-{x[4:6]}')
    return date_fmt(date).astype(np.datetime64).reshape(1)

def reformat_array(xarr, initial_index=0, mean=False, smooth=True):
    '''
    Generate a new xarray with subset of data that can be combined.
    '''

    coords = {'latitude': (('geo_index'), xarr['Latitude'].values,
                           attrs['latitude']),
              'longitude': (('geo_index'), xarr['Longitude'].values,
                           attrs['longitude']),
              'date': (('date'), int_to_datetime(xarr), attrs['date'])}

    if smooth:
        species = smooth_column(xarr, mean)
    else:
        species = xarr.Species.where(xarr.Species > 0)
        if mean:
            species = species.mean(axis=1).values.reshape(-1,1)
        else:
            species = species.bfill('Grid_Pressure').values.T[0].reshape(-1,1)
    variables = {'species':(('geo_index', 'date'), species,
                            attrs['species']),
                 'surface_alt':(('geo_index'), xarr['SurfaceAltitude'].values,
                                attrs['surface_alt'])}

    new_xarr = xr.Dataset(data_vars=variables, coords=coords)
    new_xarr = new_xarr.assign_coords({'geo_index':new_xarr.geo_index + initial_index})
    xarr.close()
    return new_xarr

def smooth_column(xarr, mean):
    '''
    Smooth the retreived species ratio according to the equation:
    x_smooth = A*x_retreived + (I-A)*x_apriori, where A is the averaging kernel, a square matrix
    '''
    A = xarr.AveragingKernel.values
    I = np.identity(14)
    x_retreived = xarr.Species.values.reshape(-1, 14, 1)
    x_apriori = xarr.ConstraintVector.values.reshape(-1, 14, 1)
    x_smooth = np.matmul(A, x_retreived) + np.matmul(I-A, x_apriori)
    x_s = pd.DataFrame(x_smooth.squeeze())
    if mean:
        return x_s.where(x_s > 0).mean(axis=1).values.reshape(-1,1)
    else:
        return x_s.where(x_s > 0).bfill(axis=1)[0].values.reshape(-1,1)

def coarsen(xarr, degrees):
    '''
    Aggregate longitude and latitude into n by n degree blocks
    '''
    rounded_latitude = xarr.latitude//degrees * degrees
    rounded_longitude = xarr.longitude//degrees * degrees
    xarr = xarr.assign_coords({'latitude':(('geo_index'),rounded_latitude),
                               'longitude':(('geo_index'),rounded_longitude)})

    xarr = xarr.set_index(geo_index=['longitude','latitude'])
    xarr = xarr.groupby('geo_index') \
                .mean().reset_index('geo_index') \
                .rename({'geo_index_level_0':'longitude',
                         'geo_index_level_1':'latitude'})

    return xarr

def get_data(agg_year = True, degrees=None, save=False, mean=False, smooth=True):
    '''
    Load data from list of files and combine into single xarray
    '''
    ds_list = []
    offset = 0
    for file in scandir(file_prefix):
        filename = file.name
        if filename.endswith('.nc') and not 'all_data' in filename:
            ds = reformat_array(xr.open_dataset(file_prefix + filename),
                                offset, mean, smooth)
            offset += ds.geo_index.size
            ds_list.append(ds)

    full_ds = xr.combine_by_coords(ds_list)

    if agg_year:
        full_ds = full_ds.resample(date='AS-MAR').mean()
        if isinstance(full_ds, xr.core.dataarray.DataArray):
            full_ds = full_ds.to_dataset()
        full_ds = full_ds.assign_coords({'geo_index':(('geo_index'),
                                                range(full_ds.geo_index.size))})

    if degrees:
        full_ds = coarsen(full_ds, degrees)

    if save:
        if degrees:
            path = f'tes_data/all_data_by{degrees}'
        else:
            path = f'tes_data/all_data'
        if not smooth:
            path += '_rough'
        if mean:
            path += '_mean.nc'
        else:
            path += '.nc'
        full_ds.to_netcdf(path)
    return full_ds

def save_years(xarr, degrees=None, mean=False, smooth=True):
    kwargs = {'squeeze':False, 'restore_coord_dims':True}
    for label, dataset in [*xarr.groupby('date', **kwargs)]:
        if degrees:
            path = f'tes_data/all_data{str(label)[:4]}_by{degrees}'
        else:
            path = f'tes_data/all_data{str(label)[:4]}'
        if mean:
            path += '_mean.nc'
        else:
            path += '.nc'
        dataset.to_netcdf(path)

def load_years(years='2005', degrees=0, mean=False, smooth=True):
    if not isinstance(years, list):
        years = [years]
    if degrees > 0:
        filenames = [f'all_data{year}_by{degrees}' for year in years]
    else:
        filenames = [f'all_data{year}' for year in years]
    if not smooth:
        path += '_rough'
    if mean:
        filenames = [path + '_mean.nc' for path in filenames]
    else:
        filenames = [path + '.nc' for path in filenames]

    found_files = []
    for file in scandir('tes_data'):
        if file.name in filenames:
            found_files.append(file.name)
    if found_files == []:
        raise ValueError('Cannot find files for specified years and degrees')
    for file in filenames:
        if file not in found_files:
            print(f'Warning: cannot find file for {file[8:-3]}')

    ds_list = []
    for filename in filenames:
        ds_list.append(xr.load_dataset('tes_data/' + filename))

    return xr.combine_nested(ds_list, concat_dim='date')
