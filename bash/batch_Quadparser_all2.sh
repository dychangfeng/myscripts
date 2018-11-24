for i in $(ls $(pwd)/genomes/); do
	python ~/my_scripts/python_scripts/Quadparser3_12.py -f $(pwd)/genomes/$i > all_G4/$i.all_G4.bed 
done
