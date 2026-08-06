# Examine Parquet and NPZ Artifacts

Start one Python session and replace the file placeholders:

```text
python3
>>> import pyarrow.parquet as pq
>>> parquet = pq.ParquetFile("FILE.parquet")
>>> table = parquet.read()
```

## Parquet

A Parquet file is a column-oriented, compressed table format designed for efficient storage and analytical queries. Indices are zero-based.

### Schema, dimensions, and row groups

```text
>>> parquet.schema_arrow
>>> parquet.metadata.num_rows
>>> parquet.metadata.num_columns
>>> parquet.metadata.num_row_groups
```

These return, respectively, the schema, number of rows, number of columns, and
number of row groups.

### First and last rows

```text
>>> table.slice(0, 5)
>>> table.slice(max(0, table.num_rows - 5), 5)
```

### Extract a column, row, or cell

```text
>>> table["COLUMN"]
>>> table.slice(123, 1)
>>> table["COLUMN"][123]
```

The three expressions return a column, row 123, and the value in `COLUMN` at
row 123.

## NPZ

NPZ is a container of named NumPy arrays. Indices are zero-based. Use `ARRAY` to select one array.

```text
python3
>>> import numpy as np
>>> npz = np.load("FILE.npz", allow_pickle=False)
```

### Array names, dimensions, and dtypes

```text
>>> npz.files
>>> npz["ARRAY"].shape
>>> npz["ARRAY"].dtype
```

### First and last rows/elements

```text
>>> np.atleast_1d(npz["ARRAY"])[:5]
>>> np.atleast_1d(npz["ARRAY"])[-5:]
```

### Extract an array row, column, or cell

These examples assume `ARRAY` is 2-D.

```text
>>> npz["ARRAY"][123]
>>> npz["ARRAY"][:, 4]
>>> npz["ARRAY"][123, 4]
```

A 1-D array has one index only; for example, `npz["ARRAY"][123]`.
