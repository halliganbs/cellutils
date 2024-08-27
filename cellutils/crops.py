import os
import numpy as np
import pandas as pd

from skimage import measure, io
from skimage.morphology import label


# Function to create a cropped image with the cell
def create_cropped_image(base_path, well_id, field_id_formatted, object_number, padding=5):
    # Load the images
    dna_path = f"{base_path}{well_id}_T0001{field_id_formatted}L01A01Z01C01.tif"
    cmo_path = f"{base_path}{well_id}_T0001{field_id_formatted}L01A02Z01C02.tif"
    claudin2_path = f"{base_path}{well_id}_T0001{field_id_formatted}L01A03Z01C04.tif"
    
    try:
        dna_image = io.imread(dna_path)
        cmo_image = io.imread(cmo_path)
        claudin2_image = io.imread(claudin2_path)
        
        # Load the mask image to find the bounding box of the object
        mask_path = f"Y:/CV8000/Sophia/BIS009_20240701_102850/masks/PECCU_{well_id}_T0001{field_id_formatted}L01A03Z01C03.tif" #change to your mask path
        mask_image = io.imread(mask_path, as_gray=True)
        labeled_mask = label(mask_image)
        props = measure.regionprops(labeled_mask)
        
        for prop in props:
            if prop.label == object_number:
                minr, minc, maxr, maxc = prop.bbox
                break
        
        # Calculate the crop boundaries with padding
        minr = max(minr - padding, 0)
        minc = max(minc - padding, 0)
        maxr = min(maxr + padding, mask_image.shape[0])
        maxc = min(maxc + padding, mask_image.shape[1])
        
        # Crop the images
        cropped_dna = dna_image[minr:maxr, minc:maxc]
        cropped_cmo = cmo_image[minr:maxr, minc:maxc]
        cropped_claudin2 = claudin2_image[minr:maxr, minc:maxc]
        
        return cropped_dna, cropped_cmo, cropped_claudin2
    except Exception as e:
        print(f"Error creating cropped image: {e}")
        return None, None, None

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
def filter_cells(df, maskpath, min_diameter=100, wellid='Image_Metadata_WellID', field='Image_Metadata_Field', obj_num_col='ObjectNumber',
                 mask_names = {
        'plate_type':"PECCU",
        'timepoint':"T0001",
        'location':'L01',
        'wid':'A03',
        'zstack':'Z01',
        'channel':'C03',
    },
                 ):
    filtered_rows = []
    
    for index, row in df.iterrows():
        well_id = row[wellid]
        field_id = row[field]
        object_number = row[obj_num_col]
        
        # Format the FieldID with 'F' prefix
        field_id_formatted = f"F{field_id}"
        
        # Get the mask path, adjusting for A03 in the path
        mask_name = f"{mask_names['plate_type']}_{well_id}_{mask_names['timepoint']}{field_id_formatted}{mask_names['location']}{mask_names['wid']}{mask_names['zstack']}{mask_names['channel']}.tif"
        mask_path = os.path.join(maskpath, mask_name)
        diameter = find_diameter(mask_path, object_number)
        if diameter and diameter > min_diameter:
            filtered_rows.append(row)
    
    return pd.DataFrame(filtered_rows)

# Randomly sample 5 objects from each cluster
def sample_objects(df, cluster_col, n_samples=5):
    sampled_df = df.groupby(cluster_col).apply(lambda x: x.sample(n=n_samples, random_state=1)).reset_index(drop=True)
    return sampled_df