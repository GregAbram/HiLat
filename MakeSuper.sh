echo "from paraview.util.vtkAlgorithm import *" >> HiLatFilters.py

for i in `grep -l paraview.util.vtkAlgorithm *py` ; do
  grep -v paraview.util.vtkAlgorithm $i >> HiLatFilters.py
done
  
