import os
import sys
import warnings
from re import search, findall
import xarray as xr
import numpy as np
import pandas as pd
import dask
from dask.diagnostics import ProgressBar
from pyproj import CRS, Transformer
from time import perf_counter
warnings.filterwarnings("ignore", category=RuntimeWarning, message='Mean of')

#create a transform to convert sinusoidal coordinates to longitude and latitude
crs = CRS.from_string('esri:54008') #sinusoidal
crs_geo = CRS.from_string('epsg:4326') #conventional lat/lon
transform = Transformer.from_crs(crs, crs_geo)

original_prefix = 'EXD/veg_data/'

def parse_metadata(attrs):

    #find year associated with data
    date_search = r'(?s)(?:BEGINNINGDATE.+?VALUE\s+=\s\"(\S+)-.+?END_OBJECT)'
    date = search(date_search, attrs['CoreMetadata.0']).groups()[0]

    #find upper left and lower right coordinates for grid
    up_left, low_right = findall(r'(?:Mtrs=(\S+))', attrs['StructMetadata.0'])
    up_left = up_left.strip('()').split(',')
    low_right = low_right.strip('()').split(',')
    upper_left = tuple([float(x) for x in up_left])
    lower_right = tuple([float(x) for x in low_right])

    north_search = r'(?s)(?:NORTHBOUND.+?VALUE\s+=\s(\S+).+?END_OBJECT)'
    south_search = r'(?s)(?:NORTHBOUND.+?VALUE\s+=\s(\S+).+?END_OBJECT)'
    north = search(north_search, attrs['ArchiveMetadata.0']).groups()[0]
    south = search(south_search, attrs['ArchiveMetadata.0']).groups()[0]
    if round(float(north)) <= -39 or round(float(south)) >= 51:
        in_bounds = False
    else:
        in_bounds = True

    return date, upper_left, lower_right, in_bounds


def reformat_array(xarr):
    # xarr = xarr.rename({'YDim:MOD44B_250m_GRID':'y',
                        # 'XDim:MOD44B_250m_GRID':'x'})

    date, upper_left, lower_right, in_bounds = parse_metadata(xarr.attrs)
    if not in_bounds:
        xarr.close()
        return None

    grid_shape = xarr.Percent_Tree_Cover.shape

    coords = {'longitude': (('longitude'), range(grid_shape[0])),
              'latitude': (('latitude'), range(grid_shape[1])),
              'date': (('date'), np.asarray([date[:7]], dtype='M'))}

    tree_cover = xarr.Percent_Tree_Cover
    tree_cover = tree_cover.where(tree_cover < 200)
    non_tree_veg = xarr.Percent_NonTree_Vegetation
    non_tree_veg = non_tree_veg.where(non_tree_veg < 200)
    non_veg = xarr.Percent_NonVegetated
    non_veg = non_veg.where(non_veg != 200, 100).where(non_veg < 200)

    variables = {'tree_cover': (('latitude', 'longitude'), tree_cover),
                 'non_tree_veg': (('latitude', 'longitude'), non_tree_veg),
                 'non_veg': (('latitude', 'longitude'), non_veg)}
    xarr.close()
    new_xarr = xr.Dataset(data_vars=variables, coords=coords)
    new_xarr = new_xarr.coarsen(latitude=4, longitude=4,
                                boundary='pad', keep_attrs=True).mean()

    grid_shape = new_xarr.tree_cover.shape
    longitude = np.linspace(upper_left[0], lower_right[0], grid_shape[1],
                          endpoint=False)
    latitude = np.linspace(upper_left[1], lower_right[1], grid_shape[0],
                          endpoint=False)
    new_xarr = transform_coords(new_xarr, longitude, latitude)
    return new_xarr

def transform_coords(xarr, lon, lat):
    lat, lon = [arr.flatten() for arr in np.meshgrid(lat, lon)]
    longitude, latitude = transform.transform(lon, lat)
    xarr = xarr.stack(geo_index=['longitude','latitude']) \
                .reset_index('geo_index') \
                .assign_coords({'longitude':(('geo_index'),longitude),
                               'latitude':(('geo_index'),latitude)}) \
                .dropna('geo_index', how='all')
    return xarr

def save_latitude_files(xarr, folder, degrees):

    if not xarr:
        return False

    bins = pd.interval_range(start=-40, end=50, freq=degrees, closed='left')
    kwargs = {'squeeze':False, 'restore_coord_dims':True}
    for label, dataset in xarr.groupby_bins('latitude', bins, **kwargs):

        dataset = dataset.assign_coords({'latitude': (('latitude'),
                                                        [label.left])})
        bounds = [str(int(bound)) for bound in [label.left, label.right]]
        path = folder + bounds[0] + '_' + bounds[1] + '.nc'

        if os.path.exists(path):
            ds = xr.load_dataset(path)
            geo_range = range(ds.geo_index.size, ds.geo_index.size +
                                                    dataset.geo_index.size)
            geo_index = (('geo_index'),geo_range)
            dataset = dataset.assign_coords({'geo_index':geo_index})
            xr.combine_by_coords([ds, dataset]).to_netcdf(path)
        else:
            geo_coords = (('geo_index'),range(dataset.geo_index.size))
            dataset.assign_coords({'geo_index':geo_coords}).to_netcdf(path)

    return True

def latitude_files(year, filenames, degrees):
    current_prefix = file_prefix + f'{year}/'
    tag_prefix = current_prefix + 'tagged/'
    for prefix in [current_prefix, tag_prefix]:
        if not os.path.exists(prefix):
            os.mkdir(prefix)
    checked, written, valid, total = 0, 0, 0, len(filenames)
    print(f'Grouping and saving {year} files by latitude...')
    t1 = perf_counter()
    for filename in filenames:
        name = tag_prefix + filename.strip(original_prefix) + '.txt'
        if not os.path.exists(name):
            xarr = reformat_array(xr.open_dataset(filename))
            if xarr:
                saved = save_latitude_files(xarr, current_prefix, degrees)
                if saved:
                    written += 1
                    valid += 1
                with open(name, 'wt') as f:
                    f.write('eat your veggies')
        else:
            valid += 1
        checked += 1
        if not checked%10:
            t2 = perf_counter()
            time = f'{(t2-t1)//60}min {round((t2-t1)%60, 4)}s'
            print(f'''Checked {checked} and wrote {written} of {total} files \
                      in {time}\n{valid}/{checked} files were in range''')
    t2 = perf_counter()
    time = f'{(t2-t1)//60}min {round((t2-t1)%60, 4)}s'
    print(f'''Grouping and saving {year} files by latitude complete.\n\
              Checked {checked} and wrote {written} of {total} files\
              in {time}'\n{valid}/{checked} files were in range''')

def scan_files(years):
    file_list = {}
    for i, year in enumerate(years):
        years[i] = 'MOD44B.A' + year
        file_list[year] = []

    for file in os.scandir(original_prefix):
        filename = file.name
        for year in years:
            if filename.startswith(year):
                file_list[year[-4:]].append(original_prefix + filename)

    return file_list

def aggregate(xarr, degrees):
    bins = pd.interval_range(start=-180, end=180, freq=degrees, closed='left')
    labels = [bin.left for bin in bins]
    kwargs = {'squeeze':False, 'restore_coord_dims':True}
    xarr = xarr.groupby_bins('longitude', bins, **kwargs, labels=labels) \
                .mean('geo_index', keep_attrs=True) \
                .rename({'longitude_bins':'longitude'})
    if xarr.longitude.dtype == np.dtype('O'):
        lon = [bin.left for bin in xarr.longitude.values]
        xarr = xarr.assign_coords({'longitude':(('longitude'), lon)})
    return xarr

def aggregate_data(years='2005', degrees=10):

    global file_prefix
    file_prefix = original_prefix
    for subdir in ['latitudes/', f'{degrees}degrees/']:
        file_prefix += subdir
        if not os.path.exists(file_prefix):
            os.mkdir(file_prefix)

    if not isinstance(years, list):
        years = [years]

    file_list = scan_files(years)

    for year, filenames in file_list.items():
        print(f'Loading {year} data...')
        latitude_files(year, filenames, degrees)
        print(f'Opening {year} latitude files...')
        year_ds = xr.open_mfdataset(f'{file_prefix}{year}/*.nc',
                                    preprocess=lambda x: aggregate(x, degrees),
                                    parallel=True, coords='minimal')
        print(f'Loaded, grouped by longitude, and combined {year} data')
        y = year_ds.to_netcdf(f'veg_data/all_data{year}by{degrees}.nc',
                            engine='h5netcdf', compute=False)
        y.compute()
        print(f'Saved {year} data to file')
        year_ds.close()

def load_years(years='2005', degrees=0):
    if not isinstance(years, list):
        years = [years]
    if degrees > 0:
        filenames = [f'all_data{year}by{degrees}.nc' for year in years]
    else:
        filenames = [f'all_data{year}.nc' for year in years]

    found_files = []
    for file in os.scandir('veg_data'):
        if file.name in filenames:
            found_files.append(file.name)
    if found_files == []:
        raise ValueError('Cannot find files for specified years and degrees')
    for file in filenames:
        if file not in found_files:
            print(f'Warning: cannot find file for {file[8:12]}')

    ds_list = []
    for filename in filenames:
        ds_list.append(xr.load_dataset('veg_data/'+ filename))

    return xr.combine_nested(ds_list, concat_dim='date')

def main(argv):
    if len(str(argv[-1])) < 4:
        degrees = int(argv.pop())
    else:
        degrees = 0
    years = [str(arg) for arg in argv]

    pbar = ProgressBar(minimum=1.0)
    pbar.register()

    aggregate_data(years=years, degrees=degrees)

if __name__ == '__main__':
    main(sys.argv[1:])
