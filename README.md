# On technical biases in iST

## Data

- CosMx data on human tonsil are available from [Zenodo](https://zenodo.org/records/16782547)
- CosMx data on human NSCLS are available from [Bruker](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/)
- Xenium data on human CRC are available through [OSTA.data](https://bioconductor.org/packages/OSTA.data)

*Data are retrieved programmatically in the corresponding script(s).
An exception is the NSCLS dataset, which requires logging in for access,
and is expected to have been downloaded to the `data/cosmx-lung` folder.

## Code

- **cosmx_tonsil_fov.R** and **xenium_colon_fov.R**:
Per-cell target probe counts per area vs. tissue-wide xy-coordinates
and, for CosMx, average counts per area across local coordinates.
- **cosmx_lung_gc.R** and **cosmx_lung_iso.R**:
Median target probe detection vs. mean GC percentage across the probe set and
number of distinct isoforms in the probe set, respectively, stratified by sample.
- **cosmx_lung_z.R**: Transcript counts by z-slice vs. FOV and visualized spatially, respectively.