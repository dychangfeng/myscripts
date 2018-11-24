for i in $(ls $(pwd)/genomes/); do
	python ~/my_scripts/python_scripts/Quadparser1_7_all_Gseq.py -f $(pwd)/genomes/$i > all_G4_7/$i.all_G4.bed 
done
