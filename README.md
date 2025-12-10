# PCI Checking and Audit Tool

A Python desktop application for detecting and analyzing Physical Cell Identifier (PCI) conflicts in LTE/4G cellular networks. This tool performs geospatial analysis to identify cells with identical PCI values within specified buffer zones, a critical task for network optimization and interference management.

## Overview

PCI conflicts occur when multiple cell sites use the same Physical Cell Identifier within overlapping coverage areas, leading to interference, handover failures, and degraded network performance. This tool automates the detection and analysis of such conflicts using buffer-based spatial analysis.

## Key Features

### Automated PCI Conflict Detection
- **Spatial Join Analysis**: Uses buffer zones to identify overlapping cell coverage
- **Unique Identifier System**: Creates composite IDs (buffer + EARFCN + PCI) to detect conflicts
- **Conflict Flagging**: Automatically marks cells with "PCI Not OK" when duplicates are found
- **Configurable Buffer Radius**: Adjustable buffer distance (default 2 km) to match network planning requirements

### Distance Metrics Calculation

#### Standard Metrics (Always Calculated)
- **Q - Distance_Site_to_Buffer**: Distance between cell site and buffer center point
- **R - Distance_same_uniq_m**: Nearest neighbor distance within conflict groups (only for problematic cells)

#### Optional Group Metrics (Toggle-able)
- **S - Nearest_same_uniq_SiteID**: Identifies the nearest conflicting site
- **T - Nearest_pair_m_per_uniq**: Minimum pairwise distance within each conflict cluster
- **U - Cluster_diameter_m**: Maximum distance between any two cells in a conflict group
- **V - Pair_ok_vs_buffer**: Boolean check if nearest pair distance ≤ 2× buffer radius

### Multiple Output Formats
- **XLSX**: Excel-compatible spreadsheet with all attributes and calculated metrics
- **GeoPackage (.gpkg)**: Geospatial format for GIS visualization and further analysis
- **Timestamped Naming**: Automatic filename generation with date/time stamps

### Performance Features
- **Efficient Haversine Calculations**: Vectorized distance computations using NumPy
- **Optimized Pairwise Analysis**: Matrix-based approach for large conflict groups
- **Optional Metrics**: Disable S/T/U/V calculations for faster processing on large datasets
- **Low Memory Operations**: Handles large CSV files efficiently

## Technical Details

### Spatial Processing
- **Coordinate System**: WGS84 (EPSG:4326) for input/output
- **Projection**: EPSG:3857 (Web Mercator) for accurate metric buffering
- **Buffer Strategy**: Creates circular buffer zones around reference points
- **Join Predicate**: Uses spatial intersection to find overlapping cells

### Distance Algorithms
- **Haversine Formula**: Accurate great-circle distances on Earth's surface
- **Pairwise Distance Matrix**: Efficient N×N calculations for group analysis
- **Self-Distance Filtering**: Excludes diagonal (self-comparisons) and near-zero duplicates
- **Nearest Neighbor**: Per-row minimum distance calculation with index tracking

### Data Requirements

#### Buffer CSV (Reference Points)
Required columns:
- `Longitude`: Decimal degrees
- `Latitude`: Decimal degrees
- `buffer`: Unique identifier for each buffer point

#### PCI Template CSV (Cell Sites)
Required columns:
- `Long`: Cell site longitude (decimal degrees)
- `Lat`: Cell site latitude (decimal degrees)
- `earfcn`: E-UTRA Absolute Radio Frequency Channel Number
- `PCI`: Physical Cell Identifier (0-503)

Optional but recommended:
- `Site ID` / `Site_ID` / `SiteID` / `Site Name`: For nearest neighbor identification

### Conflict Detection Logic

```
uniq_id = buffer_id + "_" + earfcn + "_" + PCI

If uniq_count > 1:
    → PCI conflict detected
    → Calculate nearest neighbor distances
    → Flag as "PCI Not OK"
```

## Use Cases

### Network Planning
- Pre-deployment PCI planning to avoid conflicts
- Validation of proposed site configurations
- Coverage overlap analysis

### Network Optimization
- Post-deployment PCI audit and cleanup
- Interference troubleshooting
- Handover optimization

### RF Engineering
- Drive test data validation
- Cell reselection analysis
- Co-site interference detection

### Quality Assurance
- Regular network health checks
- Compliance verification
- Performance baseline establishment

## Workflow

1. **Prepare Input Data**
   - Export buffer reference points (planned coverage centers)
   - Export cell site configuration data with PCI assignments

2. **Configure Tool**
   - Load Buffer CSV and PCI Template CSV
   - Set buffer radius based on cell coverage (typically 1-3 km)
   - Optionally enable detailed group metrics (S/T/U/V)

3. **Run Analysis**
   - Tool creates buffer zones around reference points
   - Performs spatial join to find cells within buffers
   - Calculates distance metrics and identifies conflicts

4. **Review Results**
   - Excel file contains all cells with conflict flags
   - GeoPackage allows visualization in QGIS/ArcGIS
   - Distance metrics help prioritize conflict resolution

## Output Structure

### XLSX Columns
- **Original Attributes**: All columns from both input CSVs
- **uniq_id**: Composite identifier (buffer_earfcn_PCI)
- **uniq_count**: Number of cells sharing the same uniq_id
- **PCI_check**: "PCI Not OK" flag for conflicts
- **Q - Distance_Site_to_Buffer**: Buffer-to-site distance (meters)
- **R - Distance_same_uniq_m**: Nearest neighbor in conflict group (meters)
- **S - Nearest_same_uniq_SiteID**: Site ID of nearest conflicting cell (optional)
- **T - Nearest_pair_m_per_uniq**: Minimum pair distance in group (optional)
- **U - Cluster_diameter_m**: Maximum span of conflict cluster (optional)
- **V - Pair_ok_vs_buffer**: Distance validation flag (optional)

### GeoPackage
- Polygon geometries representing buffer zones
- All XLSX attributes included
- CRS: EPSG:4326
- Layer name: "intersection"

## Performance Optimization

### For Large Datasets
- Disable S/T/U/V metrics (use checkbox) for 2-3× faster processing
- Only problematic cells are analyzed for nearest neighbors
- Vectorized operations minimize Python loops

### Memory Efficiency
- Streams CSV data in manageable chunks
- Uses NumPy arrays for distance calculations
- Minimal pandas DataFrame copies

## Requirements

- Python 3.8+
- pandas, numpy
- geopandas
- shapely
- tkinter (usually included with Python)
- openpyxl (for Excel export)

## Installation

```bash
pip install pandas numpy geopandas shapely openpyxl
```

## Usage

1. Launch the application: `python PCI_GUI.py`
2. Browse and select Buffer CSV file
3. Browse and select PCI Template CSV file
4. Set buffer radius in kilometers (default: 2 km)
5. Optionally enable group metrics calculation
6. Click "Run" to process
7. Click "Open Output Folder" to view results

## Output Location

Results are saved in an `output` subdirectory within the PCI Template CSV folder:
```
/path/to/PCI_Template.csv
/path/to/output/
    ├── intersect_result_with_check_1_20250210_143022.xlsx
    └── intersect_result_20250210_143022.gpkg
```

## Best Practices

### Buffer Radius Selection
- **Urban areas**: 1-1.5 km (smaller cells, higher density)
- **Suburban areas**: 2-3 km (medium cell sizes)
- **Rural areas**: 3-5 km (larger cell coverage)

### Data Quality
- Ensure accurate coordinates in both input files
- Verify EARFCN and PCI values are correct
- Remove duplicate rows before processing

### Conflict Resolution
- Prioritize conflicts with smallest Distance_same_uniq_m
- Check cells with Cluster_diameter_m < buffer radius
- Focus on cells with Pair_ok_vs_buffer = False

## Troubleshooting

### Common Issues
- **Missing columns**: Verify CSV headers match required names exactly
- **No conflicts detected**: Increase buffer radius or check coordinate accuracy
- **Slow processing**: Disable S/T/U/V metrics or reduce dataset size
- **Memory errors**: Process data in smaller geographic chunks

### Validation Checks
- Output row count should equal or exceed PCI template rows that fall within buffers
- uniq_count = 1 means no conflict for that cell
- Distance_same_uniq_m should be NaN for non-conflicting cells

## Author

© 2025 - baim.muhammad@gmail.com

## License

[Specify your license here]

---

## Related Tools

This tool is part of a network planning suite that includes:
- **GCELL Pie-Wedge Generator**: Creates sector polygons for coverage visualization
- **PCI Checking and Audit Tool**: Detects PCI conflicts (this tool)

Both tools integrate seamlessly into QGIS/ArcGIS workflows for comprehensive network analysis.
