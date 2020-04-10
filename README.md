# satellite-gage
LRL remote sensing river gage development

This work was supported by the United States Agency for International Development, Southern Africa Regional Mission, Fixed Amount Award 72067419FA00001. This work reflects the work of the authors and does not necessarily reflect the views of USAID or the United States Government.

All material is public; material is provided as-is without any warranty.  We request attribution for both the developer and the funding source.  Please notify the organization if you find this material particularly useful.  Thank you.

For more information, please visit: duq.edu/limpopo

The codes in this repository were written in 2020 for R.  There are several portions of the program that were written to perform various tasks then combined into the single code, master.r.  For reference, the following codes were written to perform the following tasks:
crop.r      crops a GeoTIFF file (we used Planet Labs satellite images) to a specific area of analysis
peakfind.r  analyzes normalized difference water index (NDWI) data to determine the threshold of the edge of water
meta_ext.r  extracts the metadata from Planet Labs in the accompanying XML file, specifically for the top of atmosphere reflectance coefficient
