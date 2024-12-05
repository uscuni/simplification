# simplification

[![Continuous Integration](https://github.com/uscuni/simplification/actions/workflows/testing.yml/badge.svg)](https://github.com/uscuni/simplification/actions/workflows/testing.yml)

State-of-the-Art simplification of street network geometry with Python

## Contents

### `./core/`

The curated & tested code base for the project & publication, including:

* The novel algorithm(s) developed and demonstrated; and
* Utility and plotting functionality to support the publication and algorithm(s)

### `./notebooks/`

Parameterized notebooks.

* `/preprocessing/`: workflows used to preprocess [raw data](https://github.com/martinfleis/urban-block-artifacts) (clipping and removing degree 2 nodes)
* `/methods/`: exploration of different simplification methods, including the new method proposed here, `neatnet`
* `/evaluation/`: comparative evaluation of each simplification method
* `/usecases/`: collection of use cases
* `/typology/`: exploration of different types of face artifacts, foundational to the `neatnet` algorithm
* `/archive/`: archived notebooks

### `./data/`

Curated data in `parquet` format for 6 example urban areas

| FUA  | City                                   |
| ---  | ---                                    |
| 1133 | Aleppo, Syria, Middle East / Asia      |
| 869  | Auckland, New Zealand, Oceania / Asia  |
| 809  | Douala, Cameroon, Africa               |
| 1656 | Liège, Belgium, Europe                 |
| 4617 | Bucaramanga, Colombia, S. America      |
| 4881 | Salt Lake City, Utah, USA, N. America  |
| 8989 | Wuhan, China, Far East / Asia          |


Each FUA directory contains (or will contain) the following items housing bespoke data:
* `manual/`
* `no_degree_2/`
* `original/`
* `parenx/`
* `polygons/`

### `./notes/`

Observations & hightlights from each package.

* `osmnx.md`
* `cityseer.md`

### `./protocol_images/`

Demonstration plots generated through `simplification_protocol.ipynb`.

***This directory only exists once that notebook has been run locally.***

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

### Reproducible environment

This project uses reproducible multi-platform environments using [Pixi](https://pixi.sh). To create an environment able to run all the code in the repository, clone the repository and install the locked environment using:

```sh
pixi install
```

#### Development

If you would like to run the tests, use the `tests` Pixi environment:

```sh
pixi run -e tests pytest
```

Installing pre-commit hook to the env:

```sh
pixi run pre-commit install
```