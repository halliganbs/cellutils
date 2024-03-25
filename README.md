# cellutils
Utility functions for High throughput microscopy machine learning pipelines

## Utils

- `make_well` - Converts Object level measurements to well level means
- `split_metadata` - Splits filename into metadata columns
- `zpad` - Pads zeros into well ids ie. A1 -> A01
- `zprime` - Calculates ZPrime of a plate

## Distance

- `euclidean` - Calculates euclidean distance for each row
- 'mahalanobis` - Calculates mahalobis distance for the plate

## Plotting

- `plot_plate` - Plots well level plate
