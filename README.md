# cellutils
Utility functions for High throughput microscopy machine learning pipelines

## Installation

TODO:include `requirements.txt` file

```
git clone git@github.com:halliganbs/cellutils.git

cd cellutils

pip install .
```

## Utils

- `make_well` - Converts Object level measurements to well level means
- `split_metadata` - Splits filename into metadata columns
- `zpad` - Pads zeros into well ids ie. A1 -> A01
- `zprime` - Calculates ZPrime of a plate
- `get_data_cols` - Sorts out measurement columns
- `char_range` - Generates a letters in an iterable range
- `add_controls` - Adds control rows to platemap file
- `synergy_convert` - Converts HPD300 despense platemap/viability score into synergy format for processing
- `find_number_components` - Calculates optimal number of components needed to meet provided variance threshold (defaults to 90%)

## Distance

- `euclidean` - Calculates euclidean distance for each row
- `mahalanobis` - Calculates mahalobis distance for the plate

## Plotting

- `plot_plate` - Plots well level plate
- `plot_reg_facet` - Plots facet grid of regression plots
- `synergy_plot` - Plots synergy score plot with option Bliss Synergy

## Score

- `score_plate` - Scores a database file using a provided model and column names