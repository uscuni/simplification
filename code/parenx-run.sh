#!/bin/bash

# getting subfolders in temp-parenx
FOLDERS=$(find ./temp-parenx -type d -mindepth 1) 

# run skeletonization and voronoi for each subfolder
for folder in $FOLDERS
do
	echo "Simplification for $folder started";
	skeletonize.py $folder"roads_osm.gpkg" $folder"skeletonize.gpkg"
	voronoi.py $folder"roads_osm.gpkg" $folder"voronoi.gpkg"
done

echo "Done."