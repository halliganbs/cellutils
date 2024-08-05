import os
import sqlite3
import numpy as np
import pandas as pd

from cellutils.crops import *

# Fetch data from table with cluster labels
QUERY = """
SELECT a.leiden_cluster, a.ObjectNumber, b.Image_Metadata_WellID, b.Image_Metadata_Field, b.ImageNumber
FROM MyExpt_Per_Object_UMAP a
JOIN MyExpt_Per_Image b ON a.ImageNumber = b.ImageNumber
""" #change UMAP table name and column names accodingly


def get_crops(
    db_path, output_dir='output_single_cell_crops',
    base_path="Y:/CV8000/Sophia/BIS009_20240701_102850/",
    mask_path="Y:/CV8000/Sophia/BIS009_20240701_102850/masks/",
    mask_names = {
        'plate_type':"PECCU",
        'timepoint':"T0001",
        'location':'L01',
        'wid':'A03',
        'zstack':'Z01',
        'channel':'C03',
    },
    query=QUERY,):
    
    if not os.path.exists(output_dir):
        os.path.mkdir(output_dir, exist_ok=True)
    conn = sqlite3.connect(database=db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    filtered_df = filter_cells(df, min_diameter=100)
    sampled_df = sample_objects(filtered_df, 'leiden_cluster')
    # Create and save single-cell crops for sampled objects
    for index, row in sampled_df.iterrows():
        well_id = row['Image_Metadata_WellID']
        field_id = row['Image_Metadata_Field']
        object_number = row['ObjectNumber']
        
        # Format the FieldID with 'F' prefix
        field_id_formatted = f"F{field_id}"
        
        # Get the mask path, adjusting for A03 in the path
        mask_path = f"{mask_names['plate_type']}_{mask_names['well_id']}_{mask_names['timepoint']}{field_id_formatted}{mask_names['location']}{mask_names['wid']}{mask_names['zstack']}{mask_names['channel']}.tif"
        # Find the centroid of the cell in the mask
        centroid_x, centroid_y = find_centroid(mask_path, object_number)
        if centroid_x is None or centroid_y is None:
            continue  # Skip this object if centroid calculation failed
        
        # Create the cropped images
        cropped_dna, cropped_cmo, cropped_claudin2 = create_cropped_image(base_path+mask_names['plate_type'], well_id, field_id_formatted, object_number)
        if cropped_dna is not None and cropped_cmo is not None and cropped_claudin2 is not None:
            # Save the cropped images as a multi-channel TIFF stack
            output_path = os.path.join(output_dir, f"cluster_{row['leiden_cluster']}_object_{index}.tif") #change cluster label accordingly
            io.imsave(output_path, np.stack([cropped_dna, cropped_cmo, cropped_claudin2], axis=0))
            print(f'Saved: {output_path}')