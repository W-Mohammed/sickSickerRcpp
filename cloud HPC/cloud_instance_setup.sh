#!/bin/bash
# "shebang" or "hashbang"

# Ensure the script is run as root
if [[ $EUID -ne 0 ]]; then
echo "This script must be run as root" 
exit 1
fi

# Update and upgrade the system
apt update && apt upgrade -y

# Install software-properties-common to manage repositories
apt install software-properties-common -y

# Update and upgrade the system again after installing the add-apt-repository
apt update && apt upgrade -y

# Add the R repository appropriate for the Ubuntu release $(lsb_release -cs)
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" -y

# Import the R repository signing key
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# Update and upgrade the system again after adding new repository
apt update && apt upgrade -y

# Install R base development package
apt install r-base-dev -y

# Install R packages
R -e "install.packages(c('future', 'purrr', 'furrr', 'ggplot2', 'Rcpp', 'RcppArmadillo', 'here', 'microbenchmark', 'bench', 'testthat'))"

# From the terminal, cd to folder containing file, and run:
# sudo bash ./setup_linode.sh