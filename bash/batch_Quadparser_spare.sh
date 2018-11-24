for i in $(ls $(pwd)/genomes/); do
	python ~/scripts/Quadparser_12_spare_tire.py -f $(pwd)/genomes/$i > $i.spare_G4.bed 
done
