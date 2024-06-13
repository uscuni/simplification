"""Geometry-related helpers"""

__all__ = [
    "is_within",
]


def is_within(line, poly, rtol=1e6):
    """Check if the line is within a polygon with a set tolerance

    Parameters
    ----------
    line : shapely.LineString
    poly : shapely.Polygon
    rtol : float, optional
        tolerance, by default 1e6

    Returns
    -------
    bool
    """
    if line.within(poly):
        return True
    intersection = line.intersection(poly)
    return abs(intersection.length - line.length) <= rtol
