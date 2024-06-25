"""Geometry-related helpers"""

import shapely

__all__ = [
    "is_within",
]


def is_within(
    line: shapely.LineString, poly: shapely.Polygon, rtol: float = -1e4
) -> bool:
    """Check if the line is within a polygon with a set relative tolerance.

    Parameters
    ----------
    line : shapely.LineString
        Input line to check relationship.
    poly : shapely.Polygon
        Input polygon to check relationship.
    rtol : float (default -1e6)
        The set relative tolerance.

    Returns
    -------
    bool
        ``True`` if ``line`` is either entirely within ``poly`` or if
        ``line`` is within `poly`` based on a relaxed ``rtol`` relative tolerance.
    """
    if line.within(poly):
        return True

    intersection = line.intersection(poly)
    return abs(intersection.length - line.length) <= rtol
