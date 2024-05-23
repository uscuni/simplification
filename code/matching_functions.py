import math
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial.distance import directed_hausdorff
from shapely.geometry import MultiLineString
from shapely.ops import linemerge, substring

os.environ["USE_PYGEOS"] = "0"


def _get_angle(linestring1, linestring2):
    """
    Function for getting the smallest angle between two lines.
    Does not take the direction of lines into account:
    i.e. if the angle is larger than 90,
    it is instead expressed as 180 minus the original angle.

    Arguments:
        linestring1 (Shapely LineString): first line
        linestring2 (Shapely LineString: second line

    Returns
    -------
    angle_deg (float): angle expressed in degrees
    """

    arr1 = np.array(linestring1.coords)
    arr1 = arr1[1] - arr1[0]

    arr2 = np.array(linestring2.coords)
    arr2 = arr2[1] - arr2[0]

    angle = np.math.atan2(np.linalg.det([arr1, arr2]), np.dot(arr1, arr2))
    angle_deg = abs(np.degrees(angle))

    if angle_deg > 90:
        angle_deg = 180 - angle_deg

    return angle_deg


def _get_hausdorff_dist(osm_edge, ref_edge):
    """
    Computes the Hausdorff distance (max distance) between two LineStrings.

    Arguments:

        osm_edge (LineString): The first geometry for Hausdorff distance calc.
        ref_edge (LineString): Second geometry for in distance calc.

    Returns:
        hausdorff_dist (float): The Hausdorff distance
    """

    osm_coords = list(osm_edge.coords)
    ref_coords = list(ref_edge.coords)

    hausdorff_dist = max(
        directed_hausdorff(osm_coords, ref_coords)[0],
        directed_hausdorff(ref_coords, osm_coords)[0],
    )

    return hausdorff_dist


def _get_segments(linestring, seg_length):
    """
    Convert a Shapely LineString into segments of a speficied length.
    If a line segment ends up being shorter than
    the specified distance, it is merged with
    the segment before it.

    Arguments:
        linestring (Shapely LineString): Line to be cut into segments
        seg_length (numerical): The length of the segments

    Returns:
        lines (MultiLineString):
        A multilinestring consisting of the line segments.
    """

    org_length = linestring.length

    no_segments = math.ceil(org_length / seg_length)

    start = 0
    end = seg_length
    lines = []

    for _ in range(no_segments):
        assert start != end

        line = substring(linestring, start, end)

        lines.append(line)

        start += seg_length
        end += seg_length

    # If the last segment is too short, merge it with the one before
    # Check that more than one line exist
    # (to avoid cases where the line is too short to create multiple segments)
    if len(lines) > 1:
        for i, line in enumerate(lines):
            if line.length < seg_length / 3:
                new_l = linemerge((lines[i - 1], line))

                lines[i - 1] = new_l

                del lines[i]

    lines = MultiLineString(lines)

    return lines


def _merge_multiline(line_geom):
    """
    Convert a Shapely MultiLinestring into a Linestring

    Arguments:
        line_geom (Shapely LineString or MultiLineString): geometry to be merged

    Returns:
        line_geom (Shapely LineString): original geometry as LineString
    """

    if line_geom.geom_type == "MultiLineString":
        line_geom = linemerge(line_geom)

    assert line_geom.geom_type == "LineString"

    return line_geom


def create_segment_gdf(org_gdf, segment_length):
    """
    Takes a geodataframe with linestrings and converts it into shorter segments.

    Arguments:
        gdf (geodataframe): gdf of linestrings to be converted to shorter segments
        segment_length (numerical): The length of the segments

    Returns:
        segments_gdf (geodataframe): new gdf with segments and new unique ids (seg_id)
    """
    gdf = org_gdf.copy()
    gdf["geometry"] = gdf["geometry"].apply(lambda x: _merge_multiline(x))
    assert gdf.geometry.geom_type.unique()[0] == "LineString"

    gdf["geometry"] = gdf["geometry"].apply(lambda x: _get_segments(x, segment_length))
    segments_gdf = gdf.explode(index_parts=False, ignore_index=True)

    segments_gdf.dropna(subset=["geometry"], inplace=True)

    ids = []
    for i in range(1000, 1000 + len(segments_gdf)):
        ids.append(i)

    segments_gdf["seg_id"] = ids
    assert len(segments_gdf["seg_id"].unique()) == len(segments_gdf)

    return segments_gdf


def _find_matches_from_group(group_id, groups, id_col):
    """
    Helper function for overlay_buffer().
    Returns all values in a column for a group in a grouped dataframe.

    Arguments:
        group_id (str/int): identifier of specific group
        groups (grouped datafame):
        id_col (str): name of column to be returned

    Returns:
        matches (list): column values as list
    """

    group = groups.get_group(group_id)

    matches = list(group[id_col])

    return matches


def overlay_buffer(osm_data, reference_data, dist, ref_id_col, osm_id_col):
    """
    Initial buffer matching function.
    Matches each row in a dataset with reference data (linestrings)
    to all the osm features that are within the specified buffer distance.
    The resulting dataframe contains a column with the unique reference id,
    a column 'matches_id' with a list of ids of OSM features within
    the buffered distance, and a column 'count' with the number of matches
    for each row. For better performance, the raw data should be
    segmentized into linestrings of a uniform length.

    Arguments:
        osm_data (gdf):
        reference_data (gdf):
        dist (numeric): max distance (meters) between potential matches
        ref_id_col (str): name of column with unique edge id in reference data

    Returns:
        reference_buff (df): df with the buffered matches for each reference segment
    """

    assert osm_data.crs == reference_data.crs, "Data not in the same crs!"

    reference_buff = reference_data[[ref_id_col, "geometry"]].copy(deep=True)
    reference_buff.geometry = reference_buff.geometry.buffer(distance=dist)

    # Overlay buffered geometries and osm segments
    joined = gpd.overlay(
        reference_buff, osm_data, how="intersection", keep_geom_type=False
    )

    # Group by id - find all matches for each ref segment
    grouped = joined.groupby(ref_id_col)

    reference_buff["matches_id"] = None

    group_ids = grouped.groups.keys()

    reference_buff["matches_id"] = reference_buff.apply(
        lambda x: _find_matches_from_group(x[ref_id_col], grouped, osm_id_col)
        if x[ref_id_col] in group_ids
        else 0,
        axis=1,
    )

    # Count matches
    reference_buff["count"] = reference_buff["matches_id"].apply(
        lambda x: len(x) if x.isinstance(list) else 0
    )

    # Remove rows with no matches
    reference_buff = reference_buff[reference_buff["count"] >= 1]

    reference_buff.drop("geometry", axis=1, inplace=True)

    return reference_buff


# Function for finding the best out of potential/possible matches
def _find_best_match(
    buffer_matches,
    ref_index,
    osm_edges,
    reference_edge,
    angular_threshold,
    hausdorff_threshold,
):
    """
    Finds the best match out of potential matches identifed with a buffer method.
    Computes angle and hausdorff and find best match within threshold, if any exist.

    Arguments:
        buffer_matches(dataframe): Outcome of buffer intersection step
        ref_index: the index of the reference_edge,
            locating it in the original dataset with reference segments
        osm_edges (geodataframe): osm_edges to be matched to the reference_edge
        reference_edge (linestring):
            edge currently being matched to corresponding edge in osm_edges
        angular_threshold (numerical):
            Threshold for max angle between lines considered a match (in degrees)
        hausdorff_threshold:
            Threshold for max Hausdorff distance
            between lines considered a match (in meters)

    Returns:
        best_osm_index:
            The index of the osm_edge
            identified as the best match.
            None if no match is found.
    """

    potential_matches_ix = osm_edges.loc[
        osm_edges.seg_id.isin(buffer_matches.loc[ref_index, "matches_id"])
    ].index

    osm_edges["angle"] = osm_edges.loc[potential_matches_ix].apply(
        lambda x: _get_angle(x.geometry, reference_edge), axis=1
    )

    # Avoid computing Hausdorff distance to edges that fail the angle threshold
    osm_edges["hausdorff_dist"] = osm_edges.loc[potential_matches_ix].apply(
        lambda x: _get_hausdorff_dist(osm_edge=x.geometry, ref_edge=reference_edge)
        if x.angle <= angular_threshold
        else None,
        axis=1,
    )

    # Find matches within thresholds out of all matches for this referehce geometry
    potential_matches_subset = osm_edges[
        (osm_edges.angle <= angular_threshold)
        & (osm_edges.hausdorff_dist <= hausdorff_threshold)
    ].copy()

    if len(potential_matches_subset) == 0:
        best_osm_ix = None

    elif len(potential_matches_subset) == 1:
        best_osm_ix = potential_matches_subset.index.values[0]

    else:
        # Get match(es) with smallest Hausdorff distance and angular tolerance
        potential_matches_subset["hausdorff_dist"] = pd.to_numeric(
            potential_matches_subset["hausdorff_dist"]
        )
        potential_matches_subset["angle"] = pd.to_numeric(
            potential_matches_subset["angle"]
        )

        best_matches_index = potential_matches_subset[
            ["hausdorff_dist", "angle"]
        ].idxmin()
        best_matches = potential_matches_subset.loc[best_matches_index]

        best_matches = best_matches[~best_matches.index.duplicated(keep="first")]
        # !! Duplicates may appear if the same edge
        # is the one with min dist and min angle

        if len(best_matches) == 1:
            best_osm_ix = best_matches.index.values[0]

        elif len(best_matches) > 1:
            # Take the one with the smallest hausdorff distance
            best_match_index = best_matches["hausdorff_dist"].idxmin()
            best_match = potential_matches_subset.loc[
                best_match_index
            ]  # .copy(deep=True)
            best_match = best_match[~best_match.index.duplicated(keep="first")]

            best_osm_ix = best_match.name

    return best_osm_ix


def find_matches_from_buffer(
    buffer_matches,
    osm_edges,
    reference_data,
    angular_threshold=20,
    hausdorff_threshold=12,
):
    """
    Finds the best/correct matches
    in two datasets with linestrings,
    from an initial matching based on a buffered intersection.
    Arguments:
        buffer_matches (dataframe): Outcome of buffer intersection step
        reference_data (geodataframe): reference data to be matched to osm data
        osm_edges (geodataframe): osm data to be matched to reference data
        angular_threshold (numerical):
            Threshold for max angle between lines considered a match (in degrees)
        hausdorff_threshold:
            Threshold for max Hausdorff distance between lines
            considered a match (in meters)
    Returns:
        matched_data (geodataframe):
            Reference data with additional columns
            specifying the index and ids of matched osm edges
    """

    # Get edges matched with buffer
    matched_data = reference_data.loc[buffer_matches.index].copy(deep=True)

    # Find best match within thresholds of angles and distance
    matched_data["matches_ix"] = matched_data.apply(
        lambda x: _find_best_match(
            buffer_matches,
            ref_index=x.name,
            osm_edges=osm_edges,
            reference_edge=x["geometry"],
            angular_threshold=angular_threshold,
            hausdorff_threshold=hausdorff_threshold,
        ),
        axis=1,
    )

    # Drop rows where no match was found
    matched_data.dropna(subset=["matches_ix"], inplace=True)

    matched_data["matches_ix"] = matched_data["matches_ix"].astype(int)

    ixs = matched_data["matches_ix"].values
    ids = osm_edges["seg_id"].loc[ixs].values
    matched_data["matches_id"] = ids

    print(f"{len(matched_data)} reference segments were matched to OSM edges")

    print(
        f"{len(reference_data)-len(matched_data)} reference segments were not matched"
    )

    return matched_data


def summarize_feature_matches(segments, segment_matches, seg_id_col, edge_id_col, osm):
    """
    Determine whether a feature have been matched to a feature in the other dataset
    based on how big a share of the feature's segments have been matched.
    If less than half of the length of a feature has been mapped,
    the match is considered incomplete.
    If half or more is matched, it is considered fully matched.

    Arguments:
        segments (gdf): gdf with data segments
        segment_matches (df): segment matches (result from find_matches_from_buffer())
        seg_id_col (str): column name with unique id of segments
        edge_id_col (str): column name with unique id of edges
        osm (boolean): whether segments are osm (True) or reference (False)

    Returns:
        matched_ids (list): ids of features that have been consistently matched
        undecided_ids (list): ids of features that have not been completely matched
    """

    # Create column to do join on if segments are from OSM (reference data)
    if osm:
        segment_matches[seg_id_col] = segment_matches["matches_id"]
    # Create dataframe with new and old ids and information on matches
    merged = segments.merge(
        segment_matches[[seg_id_col, "matches_ix", "matches_id"]],
        how="left",
        on=seg_id_col,
    )

    # org_ids = list( merged[edge_id_col].unique() ) # Is edge id correct?

    org_ids = list(merged.loc[merged.matches_id.notna()][edge_id_col].unique())

    matched_ids = []
    undecided_ids = []

    for i in org_ids:
        feature = merged.loc[merged[edge_id_col] == i].copy()

        # First check if all segments belonging to the id have been matched
        if len(feature.loc[feature.matches_id.isna()]) == 0:
            # If yes, mark as matched
            matched_ids.append(i)

        elif len(feature.loc[feature.matches_id.isna()]) == len(feature):
            continue

        else:
            # Set matched as either True or False and compute length
            feature["length"] = feature.geometry.length
            feature["matched"] = None

            not_matched_ix = feature.loc[feature.matches_id.isna()].index
            matched_ix = feature.loc[feature.matches_id.notna()].index
            feature.loc[not_matched_ix, "matched"] = False
            feature.loc[matched_ix, "matched"] = True

            summed = feature.groupby("matched").agg({"length": "sum"})

            majority_value = summed["length"].idxmax()

            if majority_value:
                matched_ids.append(i)

            # compute majority value len.
            # if less than some threshold and False,
            # change to True
            else:
                undecided_ids.append(i)

    return matched_ids, undecided_ids


def update_osm(osm_segments, osm_data, final_matches, attr, edge_id_col, seg_id_col):
    """
    Update osm_dataset based on the attributes
    of the reference segments that each OSM feature's segments
    have been matched to.

    Arguments:
        osm_segments (geodataframe): the osm_segments used in the matching process
        osm_data (geodataframe): original osm data to be updated
        final_matches (geodataframe): the result of the matching process
        attr (str): name of column in final_matches data
            with attribute to be transfered to osm data
        edge_id_col(str): name of column in osm_data
            with unique id of all edges/features
        seg_id_col(str): name of column in osm_data
            with unique id of all segments

    Returns:
        updated_osm (geodataframe):
            osm data with an additional column with attribute from reference data
    """

    ids_attr_dict = _summarize_attribute_matches(
        osm_segments, final_matches, edge_id_col, seg_id_col, attr
    )

    attr_df = pd.DataFrame.from_dict(ids_attr_dict, orient="index")
    attr_df.reset_index(inplace=True)
    attr_df.rename(columns={"index": edge_id_col, 0: attr}, inplace=True)
    attr_df[edge_id_col] = attr_df[edge_id_col].astype(int)

    osm_data["id_merge_col"] = osm_data[edge_id_col].astype(int)
    updated_osm = osm_data.merge(
        attr_df,
        left_on="id_merge_col",
        right_on=edge_id_col,
        how="inner",
        suffixes=("", "_matched"),
    )

    updated_osm.drop(["id_merge_col"], axis=1)

    return updated_osm


def _summarize_attribute_matches(
    osm_segments, segment_matches, edge_id_col, seg_id_col, attr
):
    """
    Create a dictionary with the original feature ids
    and the attribute they have been matched to

    Arguments:
        osm_segments (geodataframe): osm_segments used in the analysis
        final_matches: reference_data with information about corresponding osm segments
        attr (str): name of column in final_matches data
            with attribute to be transfered to osm data

    Returns (dict): a dictionary with the original osmid
        as keys and the matched value as values
    """

    # Create dataframe with new and old ids and information on matches

    segment_matches[seg_id_col] = segment_matches["matches_id"]
    osm_merged = osm_segments.merge(
        segment_matches[[seg_id_col, attr, "matches_id"]],
        how="left",
        on=seg_id_col,
        suffixes=("", "_matched"),
    )

    if attr in osm_segments.columns:
        attr = attr + "_matched"

    org_ids = list(osm_merged.loc[osm_merged.matches_id.notna()][edge_id_col].unique())

    matched_attributes_dict = {}

    for i in org_ids:
        feature = osm_merged.loc[osm_merged[edge_id_col] == i].copy(deep=True)
        feature[attr] = feature[attr].fillna("none")

        matched_values = feature[attr].unique()

        if len(matched_values) == 1:
            matched_attributes_dict[str(i)] = matched_values[0]

        else:
            feature["length"] = feature.geometry.length
            summed = feature.groupby(attr).agg({"length": "sum"})

            majority_value = summed["length"].idxmax()

            matched_attributes_dict[str(i)] = majority_value

    matched_attributes_dict = {
        key: val for key, val in matched_attributes_dict.items() if val != "none"
    }

    return matched_attributes_dict
