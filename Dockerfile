FROM continuumio/miniconda3

# Make a build environment
#RUN conda create --yes --name build-env python=3.8

# Enable conda-forge for ipopt
RUN conda config --add channels conda-forge

# Install dependencies
RUN conda install -y -q conda-build

COPY . .

#RUN conda run conda-build conda

#RUN python setup.py build
#
#RUN python setup.py install
