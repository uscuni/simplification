# Thoughts (and more questions) on `cityseer` simplification.

## Preliminary Evaluation (2024-{01,02})

* A `simplify` parameter is only available as a baked-in process within [`cityseer.tools.io.osm_graph_from_poly()`](https://github.com/benchmark-urbanism/cityseer-api/blob/30c9e9bf9078b111a297e3aeab4034465444afee/pysrc/cityseer/tools/io.py#L218), consisting of the following workflow by function call:
    1. `graphs.nx_simple_geoms()` – happens even if `simplify=False`
    2. `graphs.nx_remove_filler_nodes()` – happens even if `simplify=False`
    3.  `graphs.nx_remove_dangling_nodes()`
    4.  `graphs.nx_consolidate_nodes()`
        1. `graphs.nx_remove_filler_nodes`
        1. `graphs.nx_merge_parallel_edges()`
    5.   `graphs.nx_split_opposing_geoms()`
         1.  `graphs.nx_merge_parallel_edges()`
    6.  `graphs.nx_consolidate_nodes()`
         1. `graphs.nx_remove_filler_nodes()` 
    7.  `graphs.nx_remove_filler_nodes()`
    8.  `graphs.nx_iron_edges()`
         1.  `graphs.nx_merge_parallel_edges()`
    9.  `graphs.nx_split_opposing_geoms()`
         1.  `graphs.nx_merge_parallel_edges()`
    10.  `graphs.nx_consolidate_nodes()`
         1. `graphs.nx_remove_filler_nodes()` 
    11.  `graphs.nx_remove_filler_nodes()`
    12.  `graphs.nx_iron_edges()`
         1. `graphs.nx_merge_parallel_edges()`
* The `cityseer.tools.io.osm_graph_from_poly()` function will not meet our needs (`Linestring` in; `Linestring` out -> LILO?), but the steps above can be broken out.
* Moreover, within the prep func [`io.nx_from_generic_geopandas()`](https://github.com/benchmark-urbanism/cityseer-api/blob/30c9e9bf9078b111a297e3aeab4034465444afee/pysrc/cityseer/tools/io.py#L1010C5-L1010C30), some hardcoded simplification is already happening with [`graphs.nx_merge_parallel_edges()`](https://github.com/benchmark-urbanism/cityseer-api/blob/30c9e9bf9078b111a297e3aeab4034465444afee/pysrc/cityseer/tools/io.py#L1057), which is maybe a bit troublesome out of the box. Or maybe not?
* In short, there are many logical paths with a lot of parameters within `cityseer` to achieve network simplification, which we already kind of knew.

## Unpacking `io.osm_graph_from_poly()` (2024-03-10)
* `./code/cityseer.ipynb` explores and unpacks the simplification routine inside `io.osm_graph_from_poly()`.
* It is clearly sophisticated, but seems to have redundancies, which may of course be required depending on the area of interest.
* There are 7 unique functions being called within `io.osm_graph_from_poly()` that are directly releated to simplification. 
  * Called 1 time:
    * `graphs.nx_simple_geoms()`
    * `graphs.nx_remove_dangling_nodes()`
  * Called 2 times:
    * `graphs.nx_split_opposing_geoms()`
    * `graphs.nx_iron_edges()`
  * Called 3 times:
    * `graphs.nx_consolidate_nodes()` – 2 sets of params
  * Called 5 times:
    * `graphs.nx_merge_parallel_edges()`
  * Called 6 times:
    * `graphs.nx_remove_filler_nodes()`
