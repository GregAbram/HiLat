echo "from paraview.util.vtkAlgorithm import *" > HiLatFilters.py

for i in `grep -l paraview.util.vtkAlgorithm *py` ; do
  if test $i != 'HiLatFilters.py' ; then
    grep -v paraview.util.vtkAlgorithm $i >> HiLatFilters.py
  fi
done
  
