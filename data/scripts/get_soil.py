# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
sites_tab = snakemake.input[0]
site = snakemake.params.site
out_file = snakemake.output[0]

# test
sites_tab = "config/sites.tsv"
site = "paracou"
out_file = "soil/paracou_soil.tsv"

# libs
import pandas as pd
import ee
import xarray as xr

sites = pd.read_table(sites_tab)
sites = sites[sites["site"]==site]
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
leg = ee.Geometry.Rectangle(sites.lon[0], sites.lat[0], sites.lon[0], sites.lat[0])

def get_var(var):
    im = ee.Image("projects/soilgrids-isric/" + var)
    ic = ee.ImageCollection.fromImages([im])                       
    ds = xr.open_mfdataset([ic], engine='ee', geometry=leg)
    return ds.to_dataframe()

all_tabs = [get_var("bdod_mean"),
            get_var("cec_mean"),
            get_var("cfvo_mean"),
            get_var("clay_mean"),
            get_var("sand_mean"),
            get_var("silt_mean"),
            get_var("nitrogen_mean"),
            get_var("phh2o_mean"),
            get_var("soc_mean"),
            get_var("ocd_mean"),
            get_var("ocs_mean")]

tab = pd.concat(all_tabs, axis=1)
tab.insert(0, "site", site)
tab.to_csv(out_file, sep="\t", index=False)
