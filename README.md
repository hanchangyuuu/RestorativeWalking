# Restorative Walking Analysis

Run the Jupyter notebooks in order, then run the R script from the directory containing the model CSV.

| File | What it does |
| --- | --- |
| `1_1DataMBQ.ipynb` | Filters MOASIS walking episodes and scores restoration. Writes `active_move_respond_CH.geojson`. |
| `1_2DataNDVI.ipynb` | Samples routes every 10 m and estimates episode-level NDVI from monthly rasters. Writes `tripleg_ndvi_exposure_step10m_b200.csv`. |
| `1_3DataALL.ipynb` | Joins episodes, NDVI, and baseline participant data. Exports the model CSV. |
| `Version5_A3.R` | Fits mixed-effects models for greenness and restoration, including moderation by age and social context. Saves an HTML report and figures. |

For the **200 m** analysis, set the NDVI input in `1_3DataALL.ipynb` to `tripleg_ndvi_exposure_step10m_b200.csv` and its output to `temp_for_model_200.csv`. Then run `Rscript 'Version5_A3.R'`.
