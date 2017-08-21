#for file in $(ls | grep ".dat"); do
#	echo
#	echo "const int $(basename $file .dat | sed 's/-/_/g')[] = {"
#	cat $file | sed 's/ \|$/,/g'
#	echo "};"
#done
 
for file in $(ls | grep ".dat"); do
	echo $(basename $file .dat | sed 's/-/_/g'),
done
