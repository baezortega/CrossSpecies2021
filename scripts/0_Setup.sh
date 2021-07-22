# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 0: Setting up project data and directories

# Adrian Baez-Ortega, 2020-21


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with 
# the path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    bash scripts/0_Setup.sh


# URL for data files (Zenodo, DOI:XXXX)
DATA_URL="http://zenodo.org/record/XXXXXXXX/files/CrossSpecies2021_Data.tar.gz"  ### UPDATE


# Create required directories
mkdir -p data/original/RefGenomes data/processed output


# Download and uncompress data files ('data/original')
echo -e "\nDownloading project data:"
curl $DATA_URL -o data.tar.gz  -L -k -#

echo -e "\n\nDecompressing project data..."
tar xf data.tar.gz
rm data.tar.gz


# Download reference genomes
echo -e "\n\nDownloading reference genomes:"
while read SPECIES URL; do

    echo -e "\n$SPECIES"
    curl $URL -o data/original/RefGenomes/${SPECIES}_genome.fa.gz -L -k -#

done < data/original/CrossSpecies_RefGenomes.txt


echo -e "\nDone\n"
