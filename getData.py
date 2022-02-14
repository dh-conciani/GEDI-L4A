from google.colab import drive
drive.mount('/content/drive')

!pip install geopandas

!pip install contextily

import geopandas as gpd
grsm_poly = gpd.read_file('/content/drive/MyDrive/GEDI_L4A/gooo.geojson') 
grsm_poly.crs
## Arquivo no drive ## https://drive.google.com/file/d/1exiy_3VpjKPTw2pRnd25-WPk1Jv-tiVC/view?usp=sharing

grsm_poly.crs = "epsg:4674"
grsm_poly = grsm_poly.to_crs('epsg:31983')

## create regular grid from json
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

xmin, ymin, xmax, ymax = grsm_poly.to_crs(epsg=3857).total_bounds
print(ymin, ymin, xmax, ymax)

## set extent from each grid (in meter?)
length = 500000
wide = 500000

cols = list(np.arange(xmin, xmax + wide, wide))
rows = list(np.arange(ymin, ymax + length, length))

polygons = []
for x in cols[:-1]:
    for y in rows[:-1]:
        polygons.append(Polygon([(x,y), (x+wide, y), (x+wide, y+length), (x, y+length)]))

grid = gpd.GeoDataFrame({'geometry':polygons})
grid.crs = 'epsg:31983'

## plot
import contextily as ctx
ax=grid.plot(figsize=(20, 10), alpha=0.3, edgecolor='red')
ctx.add_basemap(ax)

import h5py
import numpy as np
import requests
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import orient
import pandas as pd
hf = h5py.File('/content/drive/MyDrive/GEDI_L4A/full_orbits2/GEDI04_A_2021075083514_O12789_04_T11129_02_002_01_V002.h5', 'r')
hf.keys()
## Arquivo ## https://drive.google.com/file/d/1U8Ey1C2eSD8bRYHmdNCIJKesBao4Q1fO/view?usp=sharing

beam0110 = hf.get('BEAM0110')
beam0110.keys()

lat_l = []
lon_l = []
beam_n = []
for var in list(hf.keys()):
    if var.startswith('BEAM'):
        beam = hf.get(var)
        lat = beam.get('lat_lowestmode')[:]
        lon = beam.get('lon_lowestmode')[:]
        lat_l.extend(lat.tolist()) # latitude
        lon_l.extend(lon.tolist()) # longitude
        n = lat.shape[0] # number of shots in the beam group
        beam_n.extend(np.repeat(str(var), n).tolist())
geo_arr = list(zip(beam_n,lat_l,lon_l))
l4adf = pd.DataFrame(geo_arr, columns=["beam", "lat_lowestmode", "lon_lowestmode"])
print(l4adf)

## compute the length of the grid
numbers = range(0, len(grid))

## create the index 
grid_index = []
for number in numbers:
  grid_index.append(number)

## criar geo ponto
l4agdf = gpd.GeoDataFrame(l4adf, geometry=gpd.points_from_xy(l4adf.lon_lowestmode, l4adf.lat_lowestmode))
## definir crs de entrada  
l4agdf.crs = 'EPSG:4326'
## reprojetar
l4agdf = l4agdf.to_crs(epsg=3857)

## process for each grid cell 
for i_grid in grid_index:
  print('grid ', i_grid+1, ' of ', len(grid))
  print('i_grid ==', i_grid)
    
  ## subset for the grid [i_grid]
  l4agdf_gsrm = l4agdf[l4agdf['geometry'].within(grid.geometry[i_grid])]  
  print('recorte')
  print(l4agdf_gsrm)  
  print('=============')

  ## quando não houver orbita para o tile, pula processo
  if len(l4agdf_gsrm) == 0:
	  print("vazio, pular")
  else:
    print('com dado -->')
    import h5py
    import numpy as np
    from glob import glob
    from os import path

    indir = '/content/drive/MyDrive/GEDI_L4A/full_orbits2'
    outdir = '/content/drive/MyDrive/GEDI_L4A/subsets2'

    for infile in glob(path.join(indir, '*.h5')):
        outfile = path.join(outdir, 'grid' + str(i_grid) + '_' + path.basename(infile))
        hf_in = h5py.File(infile, 'r')
        hf_out = h5py.File(outfile, 'w')
        
        # copy ANCILLARY and METADATA groups
        var1 = ["/ANCILLARY", "/METADATA"]
        for v in var1:
            hf_in.copy(hf_in[v],hf_out)
        
        # loop through BEAMXXXX groups
        for v in list(hf_in.keys()):
            if v.startswith('BEAM'):
                beam = hf_in[v]
                # find the shots that overlays the area of interest (GRSM)
                lat = beam['lat_lowestmode'][:]
                lon = beam['lon_lowestmode'][:]
                i = np.arange(0, len(lat), 1) # index

                geo_arr = list(zip(lat,lon, i))
                l4adf_i = pd.DataFrame(geo_arr, columns=["lat_lowestmode", "lon_lowestmode", "i"])
                l4agdf_i = gpd.GeoDataFrame(l4adf_i, geometry=gpd.points_from_xy(l4adf_i.lon_lowestmode, l4adf_i.lat_lowestmode))
                l4agdf_i.crs = 'EPSG:4326'
                ## reprojetar
                l4agdf_i = l4agdf_i.to_crs(epsg=3857)

                l4agdf_gsrm_i = l4agdf_i[l4agdf_i['geometry'].within(grid.geometry[i_grid])]  
                indices = l4agdf_gsrm_i.i

                # copy BEAMS to the output file
                for key, value in beam.items():
                    if isinstance(value, h5py.Group):
                        for key2, value2 in value.items():
                            group_path = value2.parent.name
                            group_id = hf_out.require_group(group_path)
                            dataset_path = group_path + '/' + key2
                            hf_out.create_dataset(dataset_path, data=value2[:][indices])
                            for attr in value2.attrs.keys():
                                hf_out[dataset_path].attrs[attr] = value2.attrs[attr]
                    else:
                        group_path = value.parent.name
                        group_id = hf_out.require_group(group_path)
                        dataset_path = group_path + '/' + key
                        hf_out.create_dataset(dataset_path, data=value[:][indices])
                        for attr in value.attrs.keys():
                            hf_out[dataset_path].attrs[attr] = value.attrs[attr]

        hf_in.close()
        hf_out.close()

    ## step 4 - create above-ground biomass plot 
    lat_l = []
    lon_l = []
    agbd = []
    agbd_se=[]

    ## set directory to read files 
    outdir = '/content/drive/MyDrive/GEDI_L4A/subsets2'
    ## for each file in folder:
    for subfile in glob(path.join(outdir, '*.h5')):
        hf_in = h5py.File(subfile, 'r')
        for v in list(hf_in.keys()):
            if v.startswith('BEAM'):
                beam = hf_in[v]
                lat_l.extend(beam['lat_lowestmode'][:].tolist()) 
                lon_l.extend(beam['lon_lowestmode'][:].tolist()) 
                agbd.extend(beam['agbd'][:].tolist())
                agbd_se.extend(beam['agbd_se'][:].tolist())  
        hf_in.close()
    geo_arr = list(zip(agbd,agbd_se,lat_l,lon_l))
    df = pd.DataFrame(geo_arr, columns=["agbd", "agbd_se", "lat_lowestmode", "lon_lowestmode"])
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon_lowestmode, df.lat_lowestmode))
    gdf


    ## plot
    #grsm_df = pd.DataFrame([[-9999,-9999,0,-9999,-9999,grsm_epsg4326.geometry.item()]], columns=["agbd", "agbd_se", "l4_quality_flag","lat_lowestmode", "lon_lowestmode", "geometry"])
    #gdf = pd.concat([gdf, grsm_df])
    #gdf.crs="EPSG:4326"
    #gdf_epsg3857 = gdf.to_crs(epsg=3857)
    #ax4=gdf_epsg3857[-1:].plot(color='white', edgecolor='red', alpha=0.3, linewidth=5, figsize=(22, 7))
    #gdf_epsg3857[gdf_epsg3857['agbd'] != -9999][:-1].plot(ax=ax4, column='agbd', alpha=0.1, linewidth=0, legend=True)
    #ctx.add_basemap(ax4)

    ## plot2
    #ax4=gdf_epsg3857[-1:].plot(color='white', edgecolor='red', alpha=0.3, linewidth=5, figsize=(22, 7))
    #gdf_epsg3857[gdf_epsg3857['agbd_se'] != -9999][:-1].plot(ax=ax4, column='agbd_se', alpha=0.1, linewidth=0, legend=True)
    #ctx.add_basemap(ax4)

    ## step 4c
    from glob import glob
    from os import path
    import pandas as pd
    import h5py

    outdir = '/content/drive/MyDrive/GEDI_L4A/subsets2'
    subset_df = pd.DataFrame()
    for subfile in glob(path.join(outdir, '*.h5')):
        hf_in = h5py.File(subfile, 'r')
        for v in list(hf_in.keys()):
            if v.startswith('BEAM'):
                col_names = []
                col_val = []
                beam = hf_in[v]
                # copy BEAMS 
                for key, value in beam.items():
                    # looping through subgroups
                    if isinstance(value, h5py.Group):
                        for key2, value2 in value.items():
                            if (key2 != "shot_number"):
                                # xvar variables have 2D
                                if (key2.startswith('xvar')):
                                    for r in range(4):
                                        col_names.append(key2 + '_' + str(r+1))
                                        col_val.append(value2[:, r].tolist())
                                else:
                                    col_names.append(key2)
                                    col_val.append(value2[:].tolist())
                    
                    #looping through base group
                    else:
                        # xvar variables have 2D
                        if (key.startswith('xvar')):
                            for r in range(4):
                                col_names.append(key + '_' + str(r+1))
                                col_val.append(value[:, r].tolist())
                        else:
                            col_names.append(key)
                            col_val.append(value[:].tolist())
                
                # create a pandas dataframe        
                beam_df = pd.DataFrame(map(list, zip(*col_val)), columns=col_names) 
                # Inserting BEAM names
                beam_df.insert(0, 'BEAM', np.repeat(str(v), len(beam_df.index)).tolist())
                # Appending to the subset_df dataframe
                subset_df = subset_df.append(beam_df)
        hf_in.close()

    # Setting 'shot_number' as dataframe index. shot_number column is unique
    subset_df = subset_df.set_index('shot_number') 
    subset_df.head()

    ## save dataframe
    subset_df.to_csv('/content/drive/MyDrive/GEDI_L4A/exports/GRID' + str(i_grid) + '.csv') # Export to CSV

    ## save as shapefile 
    subset_gdf = gpd.GeoDataFrame(subset_df, geometry=gpd.points_from_xy(subset_df.lon_lowestmode, subset_df.lat_lowestmode))
    subset_gdf.crs = "EPSG:4326"

    # convert object types columns to strings. object types are not supported
    for c in subset_gdf.columns:
        if subset_gdf[c].dtype == 'object':
            subset_gdf[c] = subset_gdf[c].astype(str)

    # Export to GeoJSON
    subset_gdf.to_file('/content/drive/MyDrive/GEDI_L4A/exports/GRID' + str(i_grid) + '.geojson', driver='GeoJSON')
    # Export to ESRI Shapefile
    subset_gdf.to_file('/content/drive/MyDrive/GEDI_L4A/exports/GRID' + str(i_grid) + '.shp')

    ## clean files from subsets
    import os
    path = "/content/drive/MyDrive/GEDI_L4A/subsets2"
    imgs = os.listdir(path)
    for img in imgs:
        os.remove(f'{path}/{img}') 

    ## double-check
        imgs = os.listdir(path)
    for img in imgs:
        os.remove(f'{path}/{img}') 
