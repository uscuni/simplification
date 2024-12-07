from core import utils


def read_method_outputs(fua, method, proj_crs=None):
    match method:
        case "cityseer":
            gdf = utils.read_cityseer(fua, proj_crs)
        case "orig":
            gdf = utils.read_original(fua)
        case "manual":
            gdf = utils.read_manual(fua, proj_crs)
        case "osmnx":
            gdf = utils.read_osmnx(fua, proj_crs)
        case "parenx-voronoi":
            gdf = utils.read_parenx(fua, "voronoi", proj_crs)
        case "parenx-skeletonize":
            gdf = utils.read_parenx(fua, "skeletonize", proj_crs)
        case "sgeop":
            gdf = utils.read_sgeop(fua, proj_crs)
    return gdf
