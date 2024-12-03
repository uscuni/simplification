#!/bin/bash

# getting subfolders in temp-parenx
FOLDERS=$(find ./data -type d -mindepth 2)

# run skeletonization and voronoi for each subfolder
for folder in $FOLDERS
do
	if [[ $folder == *temp-parenx* ]] ; then
		echo "Simplification for $folder started";
		skeletonize.py $folder"/roads_osm.gpkg" $folder"/skeletonize.gpkg"
		voronoi.py $folder"/roads_osm.gpkg" $folder"/voronoi.gpkg"
	fi
done

echo "Done."