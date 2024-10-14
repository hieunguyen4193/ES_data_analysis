FROM rocker/rstudio:3.6.3
RUN apt-get update && apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt install libudunits2-dev libgdal-dev libgeos-dev libproj-dev libfontconfig1-dev -y
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk40 \
    libgit2-dev \
  && apt-get clean all && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt install libfontconfig1-dev 
RUN apt-get update && apt-get install libmagick++-6.q16-dev -y

RUN mkdir -p /home/hieunguyen
RUN mkdir -p /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_1.0-6.tar.gz -P /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-0.tar.gz -P /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/Archive/fastICA/fastICA_1.2-1.tar.gz -P /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/Archive/FNN/FNN_1.0.tar.gz -P /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/densityClust_0.3.3.tar.gz -P /home/hieunguyen/download_src
RUN wget https://cran.r-project.org/src/contrib/Archive/XML/XML_3.98-1.9.tar.gz -P /home/hieunguyen/download_src

COPY ./install_R_packages.R /home/hieunguyen/download_src
RUN Rscript /home/hieunguyen/download_src/install_R_packages.R
