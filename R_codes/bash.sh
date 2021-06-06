for z in {1..2}
do
	for y in {1..2}
	do
		for x in {1..3}
		do
	        	sbatch job.sub $z $y $x
		done
	done
done

