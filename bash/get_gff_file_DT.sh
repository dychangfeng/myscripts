cut -f1 Deinococcus_thermus_ftp_path.txt |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.gff.gz|' >genomic_directory
