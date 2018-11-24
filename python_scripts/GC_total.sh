for i in $(ls $(pwd)/genomes/); do
	python ~/scripts/GC_total_length.py -f $(pwd)/genomes/$i  
done
