
Build oceantracker conda environment
________________________________________

To ensure python/ module compatibility, recommendation is to build a conda virtual environment with the given environment.yml file

Due to version dependencies of modules outside oceantracker, it curently only works with python 3.10 and numpy < 1.24, and not yet python 3.11 

1. Install 
    
    Install/update anaconda or miniconda and git if not already available


1. Change dir to where oceantracker files will be stored , eg. within

    ``cd  mycodedir``

2. Ensure git is installed, then 

    git clone https://github.com/oceantracker/oceantracker.git

3. Change dir to oceantracker folder, eg.

    ``cd ./oceantracker``


4. Either: Make conda  environment from given file

This method is strongly recommended as it ensures  version compatibility amongst python modules. This conda link is useful. 

https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

First open a conda prompt window,  or in liunx a ordinary command window. 

Note: In Windows may need to run conda prompt window as administrator to install packages.

   From root dir of oceantracker package run 
     
     ``conda env create -f installing/environment.yml``
    
   Activate new environment

      ``conda activate oceantracker``

    
5. Or : Manually build Conda environment

        ``conda create -n oceantracker python=3.10`` 

    Note: Must use python version 3.10 (not 3.11 yet) and NumPy versions 1.21–1.23

    Activate new environment

        ``conda activate oceantracker``
   
   Then install these packages


      ``conda install -c conda-forge numba``  

      ``conda install -c conda-forge netcdf4``
        
      ``conda install -c anaconda scipy``

      ``conda install -c conda-forge pyproj``

      ``conda install -c anaconda pyyaml``

      ``conda install -c anaconda psutil``

      ``conda install -c conda-forge matplotlib``
      
      ``conda install -c conda-forge fiona``

Notes: 
   - numba should also install compatible version of numpy, if not use conda install -c anaconda numpy 
   - the fiona module, which reads geojson and shape files, can be very slow finding a conda repository, alternative ways at
https://anaconda.org/conda-forge/fiona, or use faster to use pip install fiona
   - oceantracker does work with python 3.11, however netcdf4 does not install properly and fails to find a dll when importing in windows, but using pip install netcdf4 does allow netcdf4 to be imported. Hopefully netcdf at conda will soon be compatible with faster python 3.11  

7. Make oceantracker package findable
   
   From root dir of oceantracker 

   ``pip install --no-deps -e .`` 

8. To add ability to make animation movies if needed

   ``conda install -c conda-forge ffmpeg``

9. To work with iPython/Jupyter notebooks

   ``conda install -c anaconda ipykernel``