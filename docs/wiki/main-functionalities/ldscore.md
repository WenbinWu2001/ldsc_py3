



### Effective CM Coordinates and Metadata Export

With `--ld-wind-cm`, `ldscore` uses CM coordinates interpolated from the provided genetic map when one is supplied; otherwise, it uses the `.bim` CM values, which must be informative. When `--export-ref-metadata` is requested, the exported sidecar records these effective CM coordinates—interpolated map values when a map is provided, or the original `.bim` values otherwise—without modifying the input `.bim` file.