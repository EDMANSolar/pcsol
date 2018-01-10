pcsol
====

R code for the paper: "Clear sky solar irradiance models: a review of seventy models", Fernando Antonanzas, Ruben Urraca, Jesus Polo, Oscar Perpiñán, Rodrigo Escobar (under review in Renewable & Sustainable Energy Reviews).

Installation
----

You may either clone the repository with `git` or download a ZIP file with its content.

- Clone with git:
`git clone git://github.com/EDMANSolar/pcsol`

- Download a [ZIP file](https://github.com/EDMANSolar/pcsol/archive/master.zip). 

Usage
----

Once cloned or downloaded, the repository provides a main function
named `clearSky`. This function is an interface to the set of seventy
models. These models can be inspected in the file `clearSky.R` located
in the `R_code` folder.

Use the next code to load the required package,
[`solaR`](https://oscarperpinan.github.io/solar/), configure the
working directory, and load the functions, and the list of models: 

```R
setwd('NAME_OF_YOUR_FOLDER') ## replace the text

source('R_code/clearSky.R')

source('R_code/csMother.R')
```

The function `clearSky` has three arguments:
- `meteo`: a time series of meteorological measurements, including the
  variables required by the corresponding model. This time series must
  be a `zoo` object (see [zoo
  package](https://cran.r-project.org/web/packages/zoo/)).

- `loc`: coordinates of the location where the model is to be
  evaluated. It must be a `list` or `data.frame` with elements `lon`
  (longitude) and `lat` (latitude).

- `model`: name of the model to be evaluated. It must be included in
  the set of seventy models implemented in this code. This list can be
  obtained with `print(csModels)`.
  
The repository includes the two datasets used in the paper:

```R
## cabauw data.frame
load('data_example/cabauw.RData')
## carpentras data.frame
load('data_example/carpentras.RData')
```

The coordinates of these stations are available in the `stations.csv`
file located in the `data_example` folder:

```R
BSRNcc <- read.csv('data_example/stations.csv')
```

For example, the next code evaluates the ASHRAE1972 model in Cabauw:

```R
cabauwASHRAE <- clearSky(cabauw, BSRNcc[1, ], "ASHRAE1972")
```

The result is a `zoo` object with three components, `G0`, `D0`, and
`Bn`. Next code display a comparison with the measurements:

```R
plot(cabauwASHRAE$G0, cabauw$G0, xlab = 'model', ylab = 'measurements')
```

The file `evalCode.R` in the `R_code` folder evaluates the whole set
of models in the two stations, and creates target diagrams to show the
model statistics using the [`tdr`
package](https://github.com/oscarperpinan/tdr).
