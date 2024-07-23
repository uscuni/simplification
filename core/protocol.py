import pathlib

import geopandas
import momepy
import pandas
import shapely

from .utils import read_manual, read_parquet_roads
from .viz.protocol import plot_case

__all__ = [
    "protocol_cases",
    "generate_case",
    "process_case",
]


# testing hack... MUST BE EDITABLE INSTALL < probably need to rethink...
curr_path = pathlib.Path(__file__).resolve()
if curr_path.parts[-1] == "protocol.py":
    top_dir = curr_path.parents[1]
else:
    top_dir = curr_path / "simplification"
protocol_images_dir = top_dir / "protocol_images"


protocol_cases = {
    "01": {
        "types": ["original", "manual"],
        "city": "Douala",
        "coordinates": (9.8198052, 4.0024864),
        "buffer": 200,
        "title": "Case {} - Parallel Roads - {}",
    },
    "02": {
        "types": ["original", "manual"],
        "city": "Douala",
        "coordinates": (9.8213187, 4.0098569),
        "buffer": 100,
        "title": "Case {} - Roundabout - {}",
    },
    "03": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.0683693, 36.3143588),
        "buffer": 100,
        "title": "Case {} - Diverging Roads - {}",
    },
    "04": {
        "types": ["original", "manual"],
        "city": "Auckland",
        "coordinates": (174.7512042, -36.8897222),
        "buffer": 100,
        "title": "Case {} - T-Junction - {}",
    },
    "05": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.0659947, 36.3105479),
        "buffer": 100,
        "title": "Case {} - Simple Intersection - {}",
    },
    "06": {
        "types": ["original", "manual"],
        "city": "Douala",
        "coordinates": (9.7410638, 4.0968701),
        "buffer": 100,
        "title": "Case {} - Cross-shaped Intersection - {}",
    },
    "07": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.1427047, 36.2365752),
        "buffer": 100,
        "title": "Case {} - Intersection - {}",
    },
    "08": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.1398569, 36.2401153),
        "buffer": 200,
        "title": "Case {} - Side Roads - {}",
    },
    "09": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (36.9760168, 36.1273910),
        "buffer": 200,
        "title": "Case {} - Cul-de-sac - {}",
    },
    "10": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.1680681, 36.1939477),
        "buffer": 200,
        "title": "Case {} - Ovalabout - {}",
    },
    "11": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.222222, 36.194167),
        "buffer": 500,
        "title": "Case {} - Cloverleaf Interchange - {}",
    },
    "12": {
        "types": ["original", "manual"],
        "city": "Auckland",
        "coordinates": (174.8400000, -36.9186111),
        "buffer": 300,
        "title": "Case {} - Multi-level Carriageway - {}",
    },
    "13": {
        "types": ["original", "manual"],
        "city": "Liège",
        "coordinates": (5.6155107, 50.6764454),
        "buffer": 100,
        "title": "Case {} - Special Case Roundabouts - {}",
    },
    "14": {
        "types": ["original", "manual"],
        "city": "Bucaramanga",
        "coordinates": (-73.1613889, 7.0644444),
        "buffer": 100,
        "title": "Case {} - Parallel Roads Connected with a Linking Road - {}",
    },
    "15": {
        "types": ["original", "manual"],
        "city": "Bucaramanga",
        "coordinates": (-73.1275000, 7.1252778),
        "buffer": 100,
        "title": "Case {} - Outliers - {}",
    },
    "16": {
        "types": ["original", "manual"],
        "city": "Douala",
        "coordinates": (9.6611711, 4.0880919),
        "buffer": 200,
        "title": "Case {} - Parallel Roads Leading to Different Levels - {}",
    },
    "17": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (37.1797839, 36.2086177),
        "buffer": 100,
        "title": "Case {} - Roundabout with Roads on Different Levels - {}",
    },
    "18": {
        "types": ["original", "manual"],
        "city": "Aleppo",
        "coordinates": (36.9912088, 36.0967749),
        "buffer": 500,
        "title": "Case {} - Partial Cloverleaf Interchange - {}",
    },
    "19": {
        "types": ["original", "manual"],
        "city": "Auckland",
        "coordinates": (174.7616667, -36.7991667),
        "buffer": 300,
        "title": "Case {} - Complicated Freeway Intersection - {}",
    },
}


def generate_case(
    lines: geopandas.GeoDataFrame,
    coordinates: tuple[float],
    buffer: int | float,
    remove_false_nodes: bool = False,
) -> tuple[geopandas.GeoDataFrame]:
    """Generate the roads & vertices of a single Case-Type example
    of our simplification protocol clipped to the bounds produced
    from the buffering the central input coordinates.

    Parameters
    ----------
    lines : geopandas.GeoDataFrame
        Streets within FUA.
    coordinates : tuple[float]
        Coordinates of the central location of the case type in
        ``'EPSG:4326'`` and ``(X, Y)`` format.
    buffer : int | float
        Desired buffer – in meters – around ``coordinates``.
    remove_false_nodes : bool = False
        Remove interstitial nodes from the data in ``input_file``.

    Returns
    -------
    clipped_lines,  clipped_vertices : tuple[geopandas.GeoDataFrame]
        Line and vertex data from ``input_file`` clipped to the bounds
        produced from the buffering the central input coordinates.
    """

    if remove_false_nodes:
        lines = momepy.remove_false_nodes(lines)

    # Extract vertices
    _endpoint = lambda x: shapely.get_point(  # noqa: E731
        lines.geometry.explode(ignore_index=True), x
    )
    _endpoints = pandas.concat([_endpoint(0), _endpoint(-1)]).drop_duplicates()
    vertices = geopandas.GeoDataFrame(geometry=_endpoints, crs=lines.crs)

    # Create a point from the coordinates
    point = [shapely.Point(coordinates)]
    center = geopandas.GeoDataFrame(geometry=point, crs="EPSG:4326").to_crs(lines.crs)

    # Create a square buffer around the central case coordinates
    buffer = center.buffer(buffer, cap_style=3)

    # Clip the data to ``buffer`` around ``coordinates``
    clipped_lines = geopandas.clip(lines, buffer)
    clipped_vertices = geopandas.clip(vertices, buffer)

    return clipped_lines, clipped_vertices


def process_case(
    protocol_case: str,
    protocol_type: str,
    city: int | str,
    coordinates: tuple[float],
    buffer: int | float,
    title: str,
    **savefig_kwargs: dict,
) -> None:
    """
    1. Generate the roads & vertices of a single Case-Type example
        of our simplification protocol clipped to the bounds produced
        from the buffering the central input coordinates.
    2. Plot the results from (1.) then save a ``.png``.

    Parameters
    ----------
    protocol_case : str
        The simplification protocol case to generate and plot.
        See the ``core.protocol.protocol_cases`` dictionary for full details.
    protocol_type : str
        The type of ``protocol_case`` to generate and plot.
        For supported types (for each case) see the ``type`` value of that case
        within ``core.protocol.protocol_cases``. For example,
        ``core.protocol.protocol_cases['01']['type']``.
    city : int | str
        Either the name or FUA code of city.
        See ``core.utils.fua_city`` or ``core.utils.city_fua``.
    coordinates : tuple[float]
        Coordinates of the central location of the case type in
        ``'EPSG:4326'`` and ``(X, Y)`` format.
    buffer : int | float
        Desired buffer – in meters – around ``coordinates``.
    title : str
        Unformatted title in the style: ``Case {} - <DESCRIPTION> - {}'``.
    **savefig_kwargs : dict
        Keyword arguments passed into ``matplotlib.pyplot.savefig()``.
    """

    # ensure valid protocol type
    valid_protocol_types = ["original", "manual"]
    if protocol_type not in valid_protocol_types:
        raise ValueError(
            f"Invalid protocol type passed in: {protocol_type=}. "
            f"Must be in {valid_protocol_types}."
        )

    # collate params
    if protocol_type == "original":
        in_data = read_parquet_roads(city)
        remove_false_nodes = True
    else:
        in_data = read_manual(city, read_parquet_roads(city).crs)
        remove_false_nodes = False

    # generate the type-case roads & intersections
    lines, verts = generate_case(
        in_data, coordinates, buffer, remove_false_nodes=remove_false_nodes
    )
    # plot the type-case
    title = title.format(protocol_case, protocol_type.capitalize())
    image_fpath = protocol_images_dir / f"{protocol_case}_{protocol_type}.png"
    plot_case(lines, verts, title, image_fpath, **savefig_kwargs)
