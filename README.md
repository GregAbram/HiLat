# HiLat project data

## Big Picture

This repo contains scripts and job files to convert raw HiLat netcdf files to files amenable for visualization, as well as a few custom modules to use in visualziation.

To use: Clone this directory to \$SCRATCH direcory and cd into it.   From there, create a symlink where the netcdf data resides, named hf_output.link.  E.g.:

    ln -s /scratch/08830/dcomeau/e3sm_arctic_data/hf_output hf_output.link
    
Note that the target directory must have appropriate permissions

Once that is done, you can use the instructions below to create subdirs of the current directory containing VTK-compatible versions of teh data.

Once *that* is done, create a 'match' directory, cd into it and run:

    bash ../match.sh
       
This will create links to the actual data files but which are aligned temporally.  Then, still in tat directory,

     python ../mkpvds.py

This will create pvd files aggregating all the timestep data.   These will then be used in Paraview.

##Atmosphere


The 3 hourly atmosphere data data in a set of files containing 30 days of data: 

>20220526.WCYCL1950.arcticx4v1pg2_oARRM60to10.HFoutput.cori-knl.eam.h4.0001-MM-DD-00000.nc

Where MM and DD are the month and day of the period contained in the file.  These files contain 2D data 'tables' with rows for timesteps and columns for individual values.  Two 1D tables give the lat and lon to locate each data point.

'time' variable is days since 0001-01-01 in floats; time*24 is hours

With a python interpreter that supports netCDF4, this info comes from the following code:

    from netCDF4 import Dataset
    from glob import glob
    import numpy as np
    
    dset = '20220526.WCYCL1950.arcticx4v1pg2_oARRM60to10.HFoutput.cori-knl.eam.h4.0001-MM-DD-00000.nc'
    ds = Dataset(dset)
    time = np.array(ds.variables['time']) * 24
    print(time[0], time[-1])
    
 Note: on stampede2, **module purge; module load intel/18 impi python3 pnetcdf phdf5**

Variables in the atmosphere data are:

+ PRECC Convective precipitation rate (liq + ice) m/s
+ PRECT Total (convective and large-scale) precipitation rate (liq + ice) m/s
+ PS Surface pressure Pa
+ PSL Sea level pressure Pa
+ Q300 Specific Humidity at 300 mbar pressure surface kg/kg
+ Q500 Specific Humidity at 500 mbar pressure surface kg/kg
+ Q850 Specific Humidity at 850 mbar pressure surface kg/kg
+ QREFHT Reference height humidity kg/kg
+ T300 Temperature at 300 mbar pressure surface K
+ T500 Temperature at 500 mbar pressure surface K
+ T850 Temperature at 850 mbar pressure surface K
+ TMQ Total (vertically integrated) precipitable water kg/m2
+ TREFHT Reference height temperature K
+ U300 Zonal wind at 300 mbar pressure surface m/sa 
+ U500 Zonal wind at 500 mbar pressure surface m/s
+ U850 Zonal wind at 850 mbar pressure surface m/s
+ UBOT Lowest model level zonal wind m/s
+ V300 Meridional wind at 300 mbar pressure surface m/s
+ V500 Meridional wind at 500 mbar pressure surface m/s
+ V850 Meridional wind at 850 mbar pressure surface m/s
+ VBOT Lowest model level meridional wind m/s
+ Z300 Geopotential Z at 300 mbar pressure surface m
+ Z500 Geopotential Z at 500 mbar pressure surface m
+ Z850 Geopotential Z at 850 mbar pressure surface m
 
Code to get this info from an atmosphere data file... preferred this to ncdump -h

    ds = Dataset(dsets[0])
    for v in ds.variables.keys():
      var = ds.variables[v]
      dims = var.get_dims()
      if len(dims) == 2 and dims[0].name == 'time' and dims[1].name == 'ncol':
        print('+', var.name, var.long_name, var.units)
   

##atm.py and atm.job

The 'grid' is in arcticx4v1pg2_scrip.nc.   This is a netcdf 'table' file with a 2D variables *grid\_corner\_lat* and *grid\_corner\_lon*.   The i'th row of these contain the lat(lon) of the corners of the i'th quad.  The unique lat/lon points from these corner variables from the grid file correspond approximately to the data points in the data files.  This correspondence is made exact by rounding the lat/lon points to 1 part in 1000000.

**atm.py** creates surface data (as a VTU object) from the above files.  It first pulls the grid file, creating a spherical quad mesh consisting of the unique lat/lon points of the mesh file converted to spherical xyz points and a map associating each (lat/lon) pair with its indices in the point list.  Note that this python script uses MPI to run in parallel on subsets of the input files.

Then, for each data file, the map is used to associate data values with the vertices of the VTU.  The results are written to separate files for each timestep - that is, each input data file containing 240 timesteps (30 days, 3 hourly) results in 240 output files.

**atm.job** is an sbatch wrapper that calls vtkpython to run atm.py.  Parameters are appropriate parameters: all the above variables, a regex pointing at the input data files, the input grid file, and the output template.   It **ibrun**s vtkpython with 3 processes to run atm.py 3-way parallel.


##wind.py and wind.job

**wind.py** takes the meridional and zonal components of the wind at 4 levels to create 4 wind vectors on the XYZ sphere.   To do so it creates U and V tangent vectors and multiplies these by the input meridional and zonal components.   Note that this takes no parameters; it assumes the input are VTU files contained in an 'atm' subdirectory.  This has to match the output template of the atm step.  The output files are placed in a 'wind' subdirectory.   Again, MPI enabled.

wind.py receives no options.  The 'vectors' array in the Python code contain the output wind variables and the zonal and meridional components that comprise each.  For each point on the mesh, meridional and zonal vectors are computed and scaled by the meridional and zonal components at the point to get the wind direction in 3D coordinates.   Given the above set of available data variables, the current vectors array is:

    vectors = {
      'W300': ['U300', 'V300'],
      'W500': ['U500', 'V500'],
      'W850': ['U850', 'V850'],
      'WBOT': ['UBOT', 'VBOT']
    }

so the wind vectors are available at 4 different pressure surfaces.

**wind.job** is an sbatch wrapper for vtkpython to run wind.py 3-way parallel.

##advec.py

**advec.py** computes *pathlines* of the time-varying wind data. To quote Wikipedia, "Pathlines are the trajectories that individual fluid particles follow. These can be thought of as "recording" the path of a fluid element in the flow over a certain period. The direction the path takes will be determined by the streamlines of the fluid at each moment in time.".  

In this case, advec.py determines a set of initial seed points in the 'sample' procedure.    This selects a given set of points from the vector field vertices using a random sampling of points above some given Z value, acting as a minimum latitude, but in 3-space.

advec.py traces pathlines from these intial points iteratively in a doubly-nested loop, with an outer loop through the vector field time series and an inner loop within each time-series interval.  With each iteration, a new point is added to each pathline until some maximum number of points has been added, at which time they are 'retired' to a list of terminated lines.  (The initial cound is randomized so that all the pathlines are not retired simultaneously.  Each iteration also removes the earliest point of each terminated pathline so that they shorten and terminate synchronously with the extending active path lines.

Periodically (by default matching the time steps of the input vector field) the last N points of each pathline (both active and terminated) are produced as output points.   If there are more than N points in the active pathlines, they are truncated to the last N points.

When a pathline is terminatrd (moved to the terminated list) it is re-initialized at the original corresponding seed point and begins anew.

Since we cannot interpolate points on the spherical surface itself, the vector fields are transformed to lat/lon space.   With each iteration, the endpoints of the active pathlines are also transformed to lat/lon space, where they interpolate the vector field to produce step vectors in 3D space.   A new set of endpoints are then created and added to the active pathlines in 3D space.

The arguments to advec.py are:

-  -N substeps        substeps in each timestep interval(10)
-  -P points          max number of points in any streamline (100)
-  -T term            max number of samples from any one seed (200)
-  -S nseeds          appx number of seeds (500)
-  -O template        output file template ('tt_%05d.vtu')
-  -A variable        wind vector (WBOT)
-  -X n               (testing seed n)
-   vtu               file containing vector field time series

###A note on visualizing advec pathlines
advec.py adds a timestamp to each point along the pathlines.   By passing the pathlines into a Python Calculator using:

> age = np.max(inputs[0].PointData['timestamp']) - inputs[0].PointData['timestamp']

you convert the timestamp on each point to its age in the current timestep.   If you put a continuous colormap on that variable, enable opacity on surfaces, and reverse the default opacity (to get the least old - that is, newest) points opaque, you get pathlines that fade off into oblivion.   Tou can adjust the falloff by lowering the top end of the color range, effectively moving the ransparent point backward toward the newest end of the pathlines.

###Another note on visualizing advec pathlines`
The atmosphere data and its resulting pathlines are on the unit sphere.   You can use the Calculator to scale the coordinates to match the ocean and seaice data, which are in real-world meters.   In the resulting coordinates, the velocity units one the path lines are correct.

##Ocean/SeaIce Grid

Ocean and seaice data are processed to VTK data files (unstructured grids, .vtu extension) by the **grid.py** python code.

Although these datasets contain hexagonal cells, we instead produce a triangle set using the dual of that cell set: we form triangles by linking the centers of the three hexagonal cells that share internal vertices of the hex grid.   Thus the triangle-set vertices are the hex-set cell centers, found in the grid file 'xCell', 'yCell' and 'zCell' variables.   The triangles correspond to the 'cellsOnVertex' variable of the grid file when it contains 3 valid indices.  There is additional logic to aquire a mapping from hex-set vertices to hex-set centers, so that we can average hex-set vertex dependent data to get hex-set cell centers, and thus triangle-set vertices.

  
Code to get the contents of Ocean/SeaIce netcdf files::

	import sys
	from netCDF4 import Dataset

    ds = Dataset(sys.argv[1])
    for v in ds.variables.keys():
      va = ds.variables[v]
      dims = [str(d.name) for d in va.get_dims()]
      print v, ': (', ','.join(dims), ')'
      if 'long_name' in dir(va):
        ln = va.long_name
      else:
        ln = "None"
      if 'units' in dir(va):
        u = va.units
      else:
        u = "None"
      print v, ': (', ','.join(dims), ')', ln, ud

###Grid file

The grid file is: ocean.ARRM60to10.180715.nc containing 

    SSHGradientMeridional : ( Time,nCells ) Meridional gradient of SSH reconstructed at cell centers m m^{-1}
    SSHGradientZonal : ( Time,nCells ) Zonal gradient of SSH reconstructed at cell centers m m^{-1}
    angleEdge : ( nEdges ) Angle the edge normal makes with local eastward direction. radians
    areaCell : ( nCells ) Area of each cell in the primary grid. m^2
    areaTriangle : ( nVertices ) Area of each cell (triangle) in the dual grid. m^2
    atmosphericPressure : ( Time,nCells ) Pressure at the sea surface due to the atmosphere. Pa
    bottomDepth : ( nCells ) Depth of the bottom of the ocean. Given as a positive distance from sea level. m
    boundaryLayerDepth : ( Time,nCells ) CVMix/KPP: diagnosed depth of the ocean surface boundary layer m
    cellsOnCell : ( nCells,maxEdges ) List of cells that neighbor each cell. unitless
    cellsOnEdge : ( nEdges,TWO ) List of cells that straddle each edge. unitless
    cellsOnVertex : ( nVertices,vertexDegree ) List of cells that share a vertex. unitless
    dcEdge : ( nEdges ) Length of each edge, computed as the distance between cellsOnEdge. m
    dvEdge : ( nEdges ) Length of each edge, computed as the distance between verticesOnEdge. m
    edgesOnCell : ( nCells,maxEdges ) List of edges that border each cell. unitless
    edgesOnEdge : ( nEdges,maxEdges2 ) List of edges that border each of the cells that straddle each edge. unitless
    edgesOnVertex : ( nVertices,vertexDegree ) List of edges that share a vertex as an endpoint. unitless
    fCell : ( nCells ) Coriolis parameter at cell centers. s^{-1}
    fEdge : ( nEdges ) Coriolis parameter at edges. s^{-1}
    fVertex : ( nVertices ) Coriolis parameter at vertices. s^{-1}
    filteredSSHGradientMeridional : ( Time,nCells ) Time filtered meridional gradient of SSH m m^{-1}
    filteredSSHGradientZonal : ( Time,nCells ) Time filtered zonal gradient of SSH m m^{-1}
    forcingGroupNames : ( Time,nForcingGroupsMax,StrLen ) None None
    forcingGroupRestartTimes : ( Time,nForcingGroupsMax,StrLen ) None None
    indexToCellID : ( nCells ) List of global cell IDs. unitless
    indexToEdgeID : ( nEdges ) List of global edge IDs. unitless
    indexToVertexID : ( nVertices ) List of global vertex IDs. unitless
    kiteAreasOnVertex : ( nVertices,vertexDegree ) Area of the portions of each dual cell that are part of each cellsOnVertex. m^2
    latCell : ( nCells ) Latitude location of cell centers in radians. radians
    latEdge : ( nEdges ) Latitude location of edge midpoints in radians. radians
    latVertex : ( nVertices ) Latitude location of vertices in radians. radians
    layerThickness : ( Time,nCells,nVertLevels ) layer thickness m
    lonCell : ( nCells ) Longitude location of cell centers in radians. radians
    lonEdge : ( nEdges ) Longitude location of edge midpoints in radians. radians
    lonVertex : ( nVertices ) Longitude location of vertices in radians. radians
    maxLevelCell : ( nCells ) Index to the last active ocean cell in each column. unitless
    meshDensity : ( nCells ) Value of density function used to generate a particular mesh at cell centers. unitless
    nEdgesOnCell : ( nCells ) Number of edges that border each cell. unitless
    nEdgesOnEdge : ( nEdges ) Number of edges that surround each of the cells that straddle each edge. These edges are used to reconstruct the tangential velocities. unitless
    normalBarotropicVelocity : ( Time,nEdges ) barotropic velocity, used in split-explicit time-stepping m s^{-1}
    normalVelocity : ( Time,nEdges,nVertLevels ) horizonal velocity, normal component to an edge m s^{-1}
    refBottomDepth : ( nVertLevels ) Reference depth of ocean for each vertical level. Used in 'z-level' type runs. m
    restingThickness : ( nCells,nVertLevels ) Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations. m
    salinity : ( Time,nCells,nVertLevels ) salinity grams salt per kilogram seawater
    salinitySurfaceValue : ( Time,nCells ) salinity extrapolated to ocean surface PSU
    seaIcePressure : ( Time,nCells ) Pressure at the sea surface due to sea ice. Pa
    simulationStartTime : ( StrLen ) start time of first simulation, with format 'YYYY-MM-DD_HH:MM:SS' unitless
    surfaceVelocityMeridional : ( Time,nCells ) Meridional surface velocity reconstructed at cell centers m s^{-1}
    surfaceVelocityZonal : ( Time,nCells ) Zonal surface velocity reconstructed at cell centers m s^{-1}
    temperature : ( Time,nCells,nVertLevels ) potential temperature degrees Celsius
    temperatureSurfaceValue : ( Time,nCells ) potential temperature extrapolated to ocean surface degrees Celsius
    vertCoordMovementWeights : ( nVertLevels ) Weights used for distribution of sea surface heigh purturbations through multiple vertical levels. unitless
    vertNonLocalFluxTemp : ( Time,nCells,nVertLevelsP1 ) CVMix/KPP: nonlocal boundary layer mixing term for temperature nondimensional
    verticesOnCell : ( nCells,maxEdges ) List of vertices that border each cell. unitless
    verticesOnEdge : ( nEdges,TWO ) List of vertices that straddle each edge. unitless
    weightsOnEdge : ( nEdges,maxEdges2 ) Reconstruction weights associated with each of the edgesOnEdge. unitless
    xCell : ( nCells ) X Coordinate in cartesian space of cell centers. unitless
    xEdge : ( nEdges ) X Coordinate in cartesian space of edge midpoints. unitless
    xVertex : ( nVertices ) X Coordinate in cartesian space of vertices. unitless
    yCell : ( nCells ) Y Coordinate in cartesian space of cell centers. unitless
    yEdge : ( nEdges ) Y Coordinate in cartesian space of edge midpoints. unitless
    yVertex : ( nVertices ) Y Coordinate in cartesian space of vertices. unitless
    zCell : ( nCells ) Z Coordinate in cartesian space of cell centers. unitless
    zEdge : ( nEdges ) Z Coordinate in cartesian space of edge midpoints. unitless
    zVertex : ( nVertices ) Z Coordinate in cartesian space of vertices. unitless
    
with:

    Time = UNLIMITED ; // (1 currently)
    nCells = 619264 ;
    nEdges = 1882616 ;
    nVertices = 1262761 ;
    maxEdges = 8 ;
    TWO = 2 ;
    vertexDegree = 3 ;
    maxEdges2 = 16 ;
    nForcingGroupsMax = 4 ;
    StrLen = 64 ;
    nVertLevels = 80 ;
    nVertLevelsP1 = 81 ;

We note the presence of *layerThickness* and *maxLevelCell* which may be used to 'stack' multi-layer variables.


###Variables
The ocean data is incomplete.   One data file (20220526.WCYCL1950.arcticx4v1pg2_oARRM60to10.HFoutput.cori-knl.mpaso.hist.0001-01-01_00000.nc
) is avalable and contains:

    ssh : ( Time,nCells ) sea surface height m
    kineticEnergyCell : ( Time,nCells,nVertLevels ) kinetic energy of horizontal velocity on cells m^2 s^{-2}
    dThreshMLD : ( Time,nCells ) mixed layer depth based on density threshold m
    temperature : ( Time,nCells,nVertLevels ) potential temperature degrees Celsius
    salinity : ( Time,nCells,nVertLevels ) salinity grams salt per kilogram seawater
    surfaceVelocityZonal : ( Time,nCells ) Zonal surface velocity reconstructed at cell centers m s^{-1}
    surfaceVelocityMeridional : ( Time,nCells ) Meridional surface velocity reconstructed at cell centers m s^{-1}

with:

    nCells = 619264 ;
    Time = UNLIMITED ; // (30 currently)
    nVertLevels = 80 ;

Note the multi-layer salinity and temperature, and that there is no explicit time.   The number of timesteps and file names would indicate daily.

###grid.py

grid.py extracts a triangle grid from a 'meshfile', then for each datafile matching a regex, it creates a set of timestep files by pulling data for each timestep in the datafile, installing it on the grid, and writing it to the output as a vtu file.  Data files are expected to end in .....yyyy-mm-dd.nc, and a dimension from the file should contain the number of timesteps in the file.  The month, year and time interval - the number pwd
of hours bewteen timesteps - are used to create a time-offset in hours for each timestep and this value is passed through a template to create output file names.

Arguments to grid.py are:

-  -meshfile filename   .nc file containing grid information (required)
-  -datafile regex      regex for datafiles (required)
-  -gvars gridvars      comma separated array of data variables (all variables)
-  -vars datavars       comma separated array of data variables (all variables)
-  -o template          template for output file name - one %s for timestep, .vtu extension ('timestep_%08d.vtu')
-  -tint varname        interval between timesteps (24)
-  -tdim varname        name of dimension containing the number of timesteps ('Time')
-  -uvw                 add UVW tangent/radial vectors to output (no)
-  -l layer             for multi-layer data, whch layer (0)


###stack.py

stack.py is similar to grid.py, but it assumes that the mesh file has two special variables:  *layerThickness* and *bottomDepth*.    Variables that are dimensioned (time, ncells, nVertLevels) contain data for each layer.   stack.py uses this data to create a space-filling grid using triangular prisms to join triangles from othe top of the layer to the bottom.    

Generating this space-filling grid is a little tricky.   The data values are *layer-centered*; that is, if there are *n* layerThicknesses, then there will be *n* data values.   However, there are *n+1* interfaces - the top, the *n-1* interfaces between layers and the bottom.   We therefore duplicate the top layer's data to the surface and the bottom layer's data to the bottom and generate inernal interface data as a weighted sum of the above and below layer's data, with weights determined by the thicknesses of the above and below layers.

One other issue makes the code a little complicated.   The layerThicknesses are given in topmost layer thickness first and bottommost last, and are meant to be accumulated upward from the bottom.   

Arguments to stack.py are:

-  -meshfile filename   .nc file containing grid information (required)
-  -datafile regex      regex for datafiles (required)
-  -vars datavars       comma separated array of data variables (all variables)
-  -o template          template for output file name - one %s for timestep, .vtu extension ('timestep_%08d.vtu')
-  -tint varname        interval between timesteps (24)
-  -tdim varname        name of dimension containing the number of timesteps ('Time')

##Ice Data

Unlike the ocean data, each ice data file contains its own 

The variables are:

    xtime : ( Time,StrLen ) None None
    iceAreaCell : ( Time,nCells ) None None
    iceVolumeCell : ( Time,nCells ) None None
    uVelocityGeo : ( Time,nVertices ) True eastwards ice velocity m/s
    vVelocityGeo : ( Time,nVertices ) True northwards ice velocity m/s
    shear : ( Time,nCells ) None None
    divergence : ( Time,nCells ) None None
    congelation : ( Time,nCells ) None None
    frazilFormation : ( Time,nCells ) None None
    latCell : ( nCells ) Latitude location of cell centers in radians. radians
    lonCell : ( nCells ) Longitude location of cell centers in radians. radians
    xCell : ( nCells ) X Coordinate in cartesian space of cell centers. unitless
    yCell : ( nCells ) Y Coordinate in cartesian space of cell centers. unitless
    zCell : ( nCells ) Z Coordinate in cartesian space of cell centers. unitless
    indexToCellID : ( nCells ) List of global cell IDs. unitless
    latEdge : ( nEdges ) Latitude location of edge midpoints in radians. radians
    lonEdge : ( nEdges ) Longitude location of edge midpoints in radians. radians
    xEdge : ( nEdges ) X Coordinate in cartesian space of edge midpoints. unitless
    yEdge : ( nEdges ) Y Coordinate in cartesian space of edge midpoints. unitless
    zEdge : ( nEdges ) Z Coordinate in cartesian space of edge midpoints. unitless
    indexToEdgeID : ( nEdges ) List of global edge IDs. unitless
    latVertex : ( nVertices ) Latitude location of vertices in radians. radians
    lonVertex : ( nVertices ) Longitude location of vertices in radians. radians
    xVertex : ( nVertices ) X Coordinate in cartesian space of vertices. unitless
    yVertex : ( nVertices ) Y Coordinate in cartesian space of vertices. unitless
    zVertex : ( nVertices ) Z Coordinate in cartesian space of vertices. unitless
    indexToVertexID : ( nVertices ) List of global vertex IDs. unitless
    cellsOnEdge : ( nEdges,TWO ) List of cells that straddle each edge. unitless
    nEdgesOnCell : ( nCells ) Number of edges that border each cell. unitless
    nEdgesOnEdge : ( nEdges ) Number of edges that surround each of the cells that straddle each edge. These edges are used to reconstruct the tangential velocities. unitless
    edgesOnCell : ( nCells,maxEdges ) List of edges that border each cell. unitless
    edgesOnEdge : ( nEdges,maxEdges2 ) List of edges that border each of the cells that straddle each edge. unitless
    dvEdge : ( nEdges ) Length of each edge, computed as the distance between verticesOnEdge. m
    dcEdge : ( nEdges ) Length of each edge, computed as the distance between cellsOnEdge. m
    areaCell : ( nCells ) Area of each cell in the primary grid. m^2
    areaTriangle : ( nVertices ) Area of each cell (triangle) in the dual grid. m^2
    cellsOnCell : ( nCells,maxEdges ) List of cells that neighbor each cell. unitless
    verticesOnCell : ( nCells,maxEdges ) List of vertices that border each cell. unitless
    verticesOnEdge : ( nEdges,TWO ) List of vertices that straddle each edge. unitless
    edgesOnVertex : ( nVertices,vertexDegree ) List of edges that share a vertex as an endpoint. unitless
    cellsOnVertex : ( nVertices,vertexDegree ) List of cells that share a vertex. unitless
    kiteAreasOnVertex : ( nVertices,vertexDegree ) Area of the portions of each dual cell that are part of each cellsOnVertex. m^2

with:

    StrLen = 64 ;
    Time = UNLIMITED ; // (31 currently)
    nCells = 619264 ;
    nVertices = 1262761 ;
    nEdges = 1882616 ;
    TWO = 2 ;
    maxEdges = 8 ;
    maxEdges2 = 16 ;
    vertexDegree = 3 ;

We note that the grid information is available in the ice data files; thus no explicit grid file is necessary.
    
    
# Notes

## Empty data cells and invalid values
   
The volumetric data contains invalid values that are in fact used to mark the corresponding cells as empty.   This is the case for the bottom-most cells in much of the data, where the water is relatively shallow (relative to the deepest parts of the ocean).    These can be removed by passing the data through the Threshold filter, specifying the 'salinity' component, and setting a lower threshold of 0.   You'll need to atttend to the data range which is by default tied to the erroneous values in the original unthresholded data.  Salinity is *never* negative.

## Pathline Age

Pathlines at any timestep can be assigned a 'age' attribute by subtracting its 'timestamp' attribute from the maximum timestamp in the timestep.  This can be done by the following equation
in a Python Calculator filter:

    np.max(inputs[0].PointData['timestamp']) - inputs[0].PointData['timestamp']

Note that the *lowest* value corresponds to the lowest age, or the most recent bit of the pathlines.   The apparent length of the pathlines can be adjusted by tweaking the colormap range is.   Tapering the opacity is very effective, but the default opacity map maps the highest value to opaque.   This needs to be reversed.

# Matching Timesteps

We have atmosphere timesteps are 3 hourly while the remainder are 24 hourly.   For convenience in Paraview we need to make matching .vtu files.   Right now, I create a 'match' directory with matching timesteps at a 24 hour rate.

## match.sh 

This drives off the 'seaice' data, linking the seaice, ocean, wind and pathline datasets at the seaice time delta.    If the corresponding ocean, wind and pathlines are not available, they are not linked.   Run this in the match directory.  The links are in:

- ocean_nnnnnnnn.vtu
- seaice_nnnnnnnn.vtu
- wind_nnnnnnnnnn.vtu
- wbot_nnnnnnnnn.vtu
- w300_nnnnnnnnn.vtu
- w500_nnnnnnnnn.vtu
- w850_nnnnnnnnn.vtu

## mkpvds.py

Given a comma separated set of variable names (eg: ocean,seaice,wbot...), a starting hour, a stopping hour, a delta hour and a vtu-name template, generate vtus.   For example, in 

    python mkpvds.py ocean,w300,w500,w850,wbot,wind 24 696 24 %s-month.vtu
    
will create one month daily pvd files.

# Depth Scaling

The volumetric data - ocean3D - is very thin: the ocean is shallow compared to the radius of the Earth.   The DScale filter scales the depth by a factor given on its properties page, defaulting to 100.

# Corner Cuts

A useful visualization technique for the ocean data (eg. ocean3D) is to use a wedge of the Earth, with slices of the ocean3D data visible on the sides of the wedge; blank disks filling the remainder of the sides of the wedge, and the surface data clipped to the wedge.    Four custom Python filters are included to do this.


First, however, these filters ***do not work*** when the data is left in meter cordinates.  We must first normalize the data to the unit sphere using a Python programmable filter containing:

    import numpy as np
    r = np.linalg.norm(inputs[0].Points, axis=1)
    rmax = np.max(r)
    output.Points = output.Points / rmax
    
Do not forget to check the Copy Arrays box!

My apologies for the complexity here... this doesn't fit well in the Paraview paradigm,

## Step 1: Selecting the Wedge Coordinates

The wedge is selected using the included **Picker** filter.   Add it to the *normalized surface* data - eg. **ocean3D.pvd**, and Accept.    Now, with Picker highlighted, pick three points on the surface, beginning with the apex, then two other points in counter-clockwise fashion: first select the small '+' in the dotted box on the toolbar above the RenderView window.   Then select the 'select points on' button (or use the 'd' hot button) three times to select the three points.   **Make sure they are on the front face of the surface!**.   When you have done so, poke the 'Reload Python Module' button on the Picker property page.   This should store your pick points in your home directory as **picks.csv** (you can change this location on the property page).  Now hide the Picker output in the Pipeline browser.

## Step 2: Generating the Wedge Data

Now add the Clipper filter to the normalized surface data and Accept.   It will load the pick points from the csv file and clip out the wedge of the surface data.  Now add the Disks filter to the normalized surface data and Accept.   It will generate the faces of the wedge, offset very slightly inward so as not to interfere with the slices of the volumetric data.

Now attach the Slicer filter to the normalized depth-scaled **volumetric** data and Accept.   This will slice the volumetric data along the wedge sides.   Now you should be able to color the Slicer and Clipper outputs with the same variable and see them.

## Changing the Wedge

After you have done the above steps, you can change the wedge geometry.   First, hide the Clipper, Slicer and Disks outputs and make the Picker visible.  Now remove the old pick points by clocking the trash can on the toolbar, then choosing three new points, then hitting 'Reload Python Module' button on the Clipper, Slicer and Disks filter and making them visible again.   And hiding the Picker as well.

This gyration is necessary because the Picker produces the picks.csv file as a side effect - changes to the pick points do not trigger re-running the others.   We reload the Python code to force an update.


## Note:

You can apply the Clipper filter to *any* normalized surface data.  In our case the ocean data does not cover land areas, so will show holes in the wedge outer surface where there is land.   You can apply the Clipper to a sphere of radius *just less than one* to fill these holes.   These filters only work on vtkUnstructuredGrid data types; if you use a Sphere source to generate the sphere, you'll need to use the Append Datasets filter to convert it from PolyData to UnstructuredGrid.

## TruncatePathLine

A python programmable filter to truncate pathlines by age.  age=0 is the head of the pathline
