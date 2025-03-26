#!/bin/bash

# getting subfolders in temp-parenx
FOLDERS=$(find ./data -mindepth 2 -type d)

# run skeletonization and voronoi for each subfolder
for folder in $FOLDERS
do
	if [[ $folder == *temp-parenx* ]] ; then
		echo "Simplification for $folder started";
        
		echo "Running skeletonize.py..."
		/usr/bin/time -v skeletonize.py $folder"/roads_osm.gpkg" $folder"/skeletonize.gpkg" 2> ${folder}/skeletonize_mem.log
		echo "Skeletonize peak memory: $(grep "Maximum resident set size" ${folder}/skeletonize_mem.log | awk '{print $6}') kB"

		echo "Running voronoi.py..."
		/usr/bin/time -v voronoi.py $folder"/roads_osm.gpkg" $folder"/voronoi.gpkg" 2> ${folder}/voronoi_mem.log
		echo "Voronoi peak memory: $(grep "Maximum resident set size" ${folder}/voronoi_mem.log | awk '{print $6}') kB"
	fi
done

echo "Done."