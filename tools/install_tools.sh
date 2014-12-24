# Picard
wget https://github.com/broadinstitute/picard/releases/download/1.126/picard-tools-1.126.zip
unzip -jf picard-tools-1.126.zip
rm picard-tools-1.126.zip

# Beagle
wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
mv beagle.r1398.jar beagle.jar

# Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
extract Trimmomatic-0.32.zip
mv Trimmomatic-0.32/trimmomatic-0.32.jar trimmomatic-0.32.jar
mkdir ../ancillary
mv Trimmomatic-0.32/* ../ancillary/
rm Trimmomatic-0.32.zip
rm -d Trimmomatic-0.32
mv Trimmomatic-0.32.jar Trimmomatic.jar