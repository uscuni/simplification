# simplification

State-of-the-Art simplification of street network geometry with Python

## Contents

### `./code/`

Parameterized notebooks and code.
* `momepy.ipynb` – exploration of the `momepy` package
* `osmnx.ipynb` – exploration of the `osmnx.simplification` module
* `cityseer.ipynb` – exploration of the `cityseer` package

### `./data/`

Curated data in `parquet` format for 6 example urban areas
* 1133 – Aleppo, Syria, Middle East / Asia
* 869 – Auckland, New Zealand, Oceania / Asis
* 809 – Douala, Cameroon, Africa
* 1656 – Liège, Belgium, Europe
* 4617 – Bucaramanga, Colombia, S. America
* 4881 – Salt Lake City, Utah, USA, N. America

### `./envs/`

Install and activate either of the `conda` (or mamba) environments as specified in `./envs/`:
* `py311_uscuni_simplification.yml`
* `py312_uscuni_simplification.yml`

```
conda env create -f <ENV>.yml
conda activate <ENV>
```

### `./notes/`

Observations & hightlights from each package.

* `osmnx.md`
* `cityseer.md`

### `./resources/`

Additional resources and previous related research.

### `./usecases/`

Demonstrations on specific types of urban form.

* 809 (Douala)
  * `cityseer`
    * examples as `douala_{1-5}.png`
* 1656 (Liège) – `osmnx`
  * highway
  * intersection
  * parkinglot
  * roandabout
