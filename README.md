# simplification

[![Continuous Integration](https://github.com/uscuni/simplification/actions/workflows/testing.yml/badge.svg)](https://github.com/uscuni/simplification/actions/workflows/testing.yml)

State-of-the-Art simplification of street network geometry with Python

## Contents

### `./core/`

The curated & tested code base for the project & publication, including:

* The novel algorithm(s) developed and demonstrated; and
* Utility and plotting functionality to support the publication and algorithm(s)

```
- core
    |-- __init__.py
    |-- geometry.py
    |-- stats.py
    |-- utils.py
    |-- algorithms
    |   |-- __init__.py
    |   |-- common.py
    |   |-- triangles.py
    |-- viz
    |   |-- __init__.py
    |   |-- context.py
    |   |-- h3_hex.py
    |-- tests
    |   |-- __init__.py
    |   |-- conftest.py
    |   |-- test_stats.py
    |   |-- test_utils.py
    |   |-- test_viz.py
    |   |-- baseline_images
    |   |   | -- *controls from image comparison*
    |   |-- data
    |   |   | -- test_data_liege.gpkg
```

### `./code/`

Parameterized notebooks.

* `_addartifacts_8989.ipynb` – ...
* `cityseer_*.ipynb` – explorations of the [`cityseer`](https://github.com/benchmark-urbanism/cityseer-api) package
* `_clip_networks.ipynb` – ...
* `clustering.ipynb` – deadend
* `_convert.ipynb` – ...
* `evaluate_h3cells.ipynb` – ...
* `momepy.ipynb` – exploration of the [`momepy`](https://github.com/pysal/momepy) package
* `osmnx.ipynb` – exploration of the [`osmnx.simplification`](https://github.com/gboeing/osmnx/blob/main/osmnx/simplification.py) module
* `parenx.ipynb` & `parenx-run.sh` - exploration of the [`parenx`](https://github.com/anisotropi4/parenx) package (skeletonization & line voronoi diagrams)
* `simplification_protocol.ipynb` – ...
* `triangles.ipynb` – ...
* `usecases.ipynb` – interesting example case locations; add more as desired

### `./data/`

Curated data in `parquet` format for 6 example urban areas
* 1133 – Aleppo, Syria, Middle East / Asia
* 869 – Auckland, New Zealand, Oceania / Asis
* 809 – Douala, Cameroon, Africa
* 1656 – Liège, Belgium, Europe
* 4617 – Bucaramanga, Colombia, S. America
* 4881 – Salt Lake City, Utah, USA, N. America
* 8989 – Wuhan, China, Far East / Asia

Each FUA directory contains (or will contain) the following items housing bespoke data:
* `manual/`
* `no_degree_2/`
* `parenex/`
* `polygons/`
* `roads_osm.parquet`

### `./envs/`

Install and activate the `conda` (or mamba) environment, which creates an environment name `simplification`:

```
conda env create -f environment.yml
conda activate simplification
```

### `./notes/`

Observations & hightlights from each package.

* `osmnx.md`
* `cityseer.md`

### `./resources/`

Additional resources and previous related research.

### `./usecases/`

Demonstration visualizations on specific types of urban form.

* 809 (Douala)
  * `cityseer` (from `notebooks/cityseer_overview_gaboardi.ipynb`)
    * examples as `douala_{1-5}.png`
* 869 (Auckland)
  * `averagedegree` (from `notebooks/evaluate_h3cells.ipynb`)
  * `totallength` (from `notebooks/evaluate_h3cells.ipynb`)
* 1656 (Liège)
  * `cityseer`
    * `parallel_edges_1_midline_False.mp4`
    * `parallel_edges_1_midline_True.mp4`
    * `parallel_edges_2_midline_False.mp4`
    * `parallel_edges_2_midline_True.mp4`
    * `parallel_edges_3_midline_False.mp4`
    * `parallel_edges_3_midline_True.mp4`
  * `osmnx`
    * `highway.mp4`
    * `intersection.mp4`
    * `parkinglot.mp4`
    * `roandabout.mp4`
  * `points.json` - use case locations
