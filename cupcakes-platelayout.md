## Summary
- takes a csv formatted with SAMPLE ID, ROW, COLUMN and outputs into a new csv formatted like a grid
- prevent copying and pasting of sample names for sample tracking and plating
- visualize samples in either 384-well format or 96-well format
### USAGE
- if running in same directory as files

```
python3 make384map.py 'OPS_010_Echo_Worksheet.csv' 'output.csv'
```
- if running in different directory as files

```
python3 make384map.py '~/worksheets/input.csv' '~/grids/output.csv'
```

## Installation
- pandas is a software library written for the Python programming language for data manipulation and analysis.
- numpy is a library for the Python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays.
- The value of sys.argv is a Python list of command line arguments that were passed to the Python script.

```
import pandas as pd
import numpy as np
import sys
df = pd.read_csv(sys.argv[1])
```

## Convert numbers to int type so we can sort
```
df["COLUMN"] = df["COLUMN"].astype(int)
```
## sort_values tells python to first sort by the ROW column then the COLUMN column
```
sort_df = df.sort_values(["ROW", "COLUMN"])
```

## get all the letters that denote the rows
```
row_elements = df["ROW"].unique()
```
## create empty dataframe with rows labeled by the letters and columns 1-24
```
new_df = pd.DataFrame(index=row_elements, columns=np.arange(1,25))
```
## iterate over the letters
- for el in row_elements:
- gets all the SAMPLE ID values that match the letter
```
IDs = df.loc[sort_df['ROW'] == el]['SAMPLE ID'].values
    col_ids = df.loc[sort_df['ROW'] == el]['COLUMN'].values
    ```
- sorted the dataframe earlier, so we know enumerating over IDs will be in the correct order
- assigns row, column, coordinate in new dataframe to matching SAMPLE ID

```
for ID, col_id in zip(IDs, col_ids):
        new_df.loc[el, col_id] = ID    
```

## final sort on the index (letters) to get everything in the right order
```
new_df = new_df.sort_index()
```

## writing the final output to a csv
```
new_df.to_csv(sys.argv[2])
```
