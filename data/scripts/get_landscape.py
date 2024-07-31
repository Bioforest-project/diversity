# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
amazonia = snakemake.input
nc = snakemake.output

# test
# amazonia = "results/data/amazonia"
amazonia = "results/data/guaviare"
# amazonia = "results/data/capricho"
# nc = "results/data/tmf/tmf.nc"
nc = "results/data/forest/deforestation_0.05.nc"
res = 0.05        

# libs
import geopandas as gp
import ee
import xarray as xr
import numpy as np
import xesmf as xe

# grid
area = gp.read_file(amazonia).dissolve()
lon_min = np.floor(area.bounds.minx[0])
lon_max = np.ceil(area.bounds.maxx[0])
lat_min = np.floor(area.bounds.miny[0])
lat_max = np.ceil(area.bounds.maxy[0])
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(lat_min, lat_max, res), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(lon_min, lon_max, res), {"units": "degrees_east"}),
    }
)

# tmf
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
leg = ee.Geometry.Rectangle(lon_min, lat_min, lon_max, lat_max)
ic = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
ds = xr.open_mfdataset(
        [ic],
        engine='ee',
        projection=ic.first().select(0).projection(),
        geometry=leg
).sel(time=2).drop('time')
ds2 = ds.to_array("year", name="tmf").to_dataset().assign_coords(
    year=np.arange(1990,2023)
)
# ds2 = ds2.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
# ds2 = ds2.rio.write_crs("epsg:4362")
# ds2 = ds2.rio.clip(area.geometry.values, area.crs)
# ds2.to_netcdf(nc)

# deforest
deforest = ds2[["lon", "lat", "year"]]
deforest["deforestation"] = (ds2.tmf == 3).astype(int)
res_tmf = ((ds2['lat'].max() - ds2['lat'].min())/(ds2['lat'].count()-1.)).values
factor = int(np.floor(res/res_tmf))
deforest_res = deforest.coarsen(lon=factor, lat=factor, boundary="trim").mean().chunk({"lon": 50, "lat": 50, "year": 5})
regridder = xe.Regridder(deforest_res, ds_out, "bilinear")
deforest_r = regridder(deforest_res, keep_attrs = True)
deforest_r.to_netcdf(nc)


# deforest = deforest.chunk({"lon": 100, "lat": 100, "year": 1})
  
# val 1. Undisturbed Tropical moist forest (TMF)  
# val 2. Degraded TMF  
# val 3. Deforested land  
# val 4. Forest regrowth  
# val 5. Permanent or seasonal Water  
# val 6. Other land cover  

# plot
# from matplotlib import pyplot as plt
# area.plot()
# ds2.sel(year=1990)["tmf"].plot()
# plt.show()
