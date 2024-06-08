# simplification

[![Continuous Integration](https://github.com/uscuni/simplification/actions/workflows/testing.yml/badge.svg)](https://github.com/uscuni/simplification/actions/workflows/testing.yml)
[![codecov](https://codecov.io/gh/uscuni/simplification/branch/main/graph/badge.svg)](https://codecov.io/gh/uscuni/simplification)

State-of-the-Art simplification of street network geometry with Python

## Contents

### `./code/`

Parameterized notebooks and code.

* `cityseer_*.ipynb` – explorations of the [`cityseer`](https://github.com/benchmark-urbanism/cityseer-api) package
* `clustering.ipynb` – deadend
* `momepy.ipynb` – exploration of the [`momepy`](https://github.com/pysal/momepy) package
* `osmnx.ipynb` – exploration of the [`osmnx.simplification`](https://github.com/gboeing/osmnx/blob/main/osmnx/simplification.py) module
* `parenx.ipynb` & `parenx-run.sh` - exploration of the [`parenx`](https://github.com/anisotropi4/parenx) package (skeletonization & line voronoi diagrams)
* `usecases.ipynb` – interesting example case locations; add more as desired
* `utils.py` – reusable functions

### `./data/`

Curated data in `parquet` format for 6 example urban areas
* 1133 – Aleppo, Syria, Middle East / Asia
* 869 – Auckland, New Zealand, Oceania / Asis
* 809 – Douala, Cameroon, Africa
* 1656 – Liège, Belgium, Europe
* 4617 – Bucaramanga, Colombia, S. America
* 4881 – Salt Lake City, Utah, USA, N. America
* 8989 – Wuhan, China, Far East / Asia

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
  * `cityseer` (from `code/cityseer_overview_gaboardi.ipynb`)
    * examples as `douala_{1-5}.png`
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
