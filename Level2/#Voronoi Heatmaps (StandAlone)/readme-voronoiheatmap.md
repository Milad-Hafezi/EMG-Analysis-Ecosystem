# Voronoi Heatmap

## Description
The Voronoi Heatmap creates visual representations of muscle excitability across the cortical surface by coloring Voronoi cells based on stimulation threshold values. This component transforms threshold current data into color-coded maps that reveal the spatial organization of muscle representation in the sensorimotor cortex, providing intuitive visualization of functional topography.

## Inputs
- Voronoi map (.svg file) processed by the Voronoi Map Generation scripts
- Verdict file (.xlsx) containing threshold current values for each site-muscle combination
- Configuration settings for visualization parameters

## Outputs
- Set of color-coded visualizations in multiple formats (fig, eps, svg):
  - **Heatmaps**: Individual muscle maps with Voronoi cells colored by threshold current
    - Darker blue indicating lower currents (more excitable sites)
    - White indicating non-responsive sites
    - Gray indicating unfair test sites
  - **Boolean Maps**: Visualizations of overlapping representations between muscle pairs
  - **Green Plots**: Composite maps showing number of muscles activated per site
  - **Centroid Analysis**: Statistical visualization of muscle representation distribution

## Visualization Types
1. **Muscle Heatmaps**: 
   - One map per muscle channel
   - Blue color gradient based on threshold current
   - Threshold sequence: 5μA, 10μA, 20μA, 40μA, 80μA, 160μA, 320μA

2. **Boolean Maps**: 
   - Shows sites where two muscles can be co-activated
   - Color intensity based on maximum threshold between muscles
   - Useful for identifying functional relationships

3. **Green Plots**:
   - Shows number of muscles activated at each site
   - Darker green indicates more muscles
   - Numeric labels showing exact count
   - Separate maps for different threshold levels

4. **Centroid Analysis**:
   - Statistical representation of muscle distribution
   - Standard and weighted centroids
   - Elliptical visualization of spatial variance
   - Color-coded by muscle

## Usage
1. Run the Voronoi Map Generation scripts to create the spatial framework
2. Run the Voronoi Heatmap component
3. Select the Voronoi map and verdict file when prompted
4. Configure visualization parameters if needed
5. Heatmaps and statistical visualizations will be generated in their respective directories

## Statistical Outputs
- **Centroid_logs**: Contains spatial statistics for each muscle representation
  - Standard and weighted centroids (x,y coordinates)
  - Spatial variance measures
  - Available in both .mat and .xlsx formats

## Notes
- The color coding is consistent across maps for easy comparison
- Statistical analysis provides quantitative measures of spatial organization
- Multiple output formats support both analysis and publication needs
- Special handling of merged sites ensures accurate representation
- All visualizations maintain anatomical orientation and scale information