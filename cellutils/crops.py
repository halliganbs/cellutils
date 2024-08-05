import os
import numpy as np
import pandas as pd

from skimage import measure, io
from skimage.morphology import label

# Function to find the centroid of the mask for a given ObjectNumber
def find_centroid(mask_path, object_number):
    mask_image = io.imread(mask_path, as_gray=True)
    labeled_mask = label(mask_image)
    props = measure.regionprops(labeled_mask)
    
    for prop in props:
        if prop.label == object_number:
            centroid_y, centroid_x = prop.centroid
            return int(centroid_x), int(centroid_y)
    return None, None

# Function to find the diameter of the mask for a given ObjectNumber
def find_diameter(mask_path, object_number):
    mask_image = io.imread(mask_path, as_gray=True)
    labeled_mask = label(mask_image)
    props = measure.regionprops(labeled_mask)
    
    for prop in props:
        if prop.label == object_number:
            return prop.equivalent_diameter
    return None
 
# Filter cells based on their diameter (remove false masks)
def filter_cells(df, min_diameter=100, wellid='Image_Metadata_WellID', field='Image_Metadata_Field', obj_num_col='ObjectNumber'):
    filtered_rows = []
    
    for index, row in df.iterrows():
        well_id = row[wellid]
        field_id = row[field]
        object_number = row[obj_num_col]
        
        # Format the FieldID with 'F' prefix
        field_id_formatted = f"F{field_id}"
        
        # Get the mask path, adjusting for A03 in the path
        mask_path = f"Y:/CV8000/Sophia/BIS009_20240701_102850/masks/PECCU_{well_id}_T0001{field_id_formatted}L01A03Z01C03.tif" #change to your mask path
        
        diameter = find_diameter(mask_path, object_number)
        if diameter and diameter > min_diameter:
            filtered_rows.append(row)
    
    return pd.DataFrame(filtered_rows)