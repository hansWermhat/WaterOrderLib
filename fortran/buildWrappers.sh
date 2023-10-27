rm *.so

python -m numpy.f2py --quiet -c -m waterlib waterlib.f90 
python -m numpy.f2py --quiet -c -m sortlib sortlib.f90
python -m numpy.f2py --quiet -c -m imagelib imagelib.f90

