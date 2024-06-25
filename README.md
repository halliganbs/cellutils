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

## Distance

- `euclidean` - Calculates euclidean distance for each row
- `mahalanobis` - Calculates mahalobis distance for the plate

## Plotting

- `plot_plate` - Plots well level plate
