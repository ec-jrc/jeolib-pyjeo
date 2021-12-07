# pyjeo

## Description

pyjeo is a library for image processing for geospatial data implemented in
JRC Ispra.

## License

pyjeo is released under an
[GPLv3](http://www.gnu.org/licenses/gpl-3.0.html) license (see
[COPYING](COPYING))

## Reference

Please refer to pyjeo using this [reference](https://doi.org/10.3390/ijgi8100461).

## Dependencies

* [miallib](https://github.com/ec-jrc/jeolib-miallib)
* [jiplib](https://github.com/ec-jrc/jeolib-jiplib)

## Install

### From source
Once the dependencies miallib and jiplib are installed, enter the 
directory of the repository, and run:

```
sudo python setup.py install
```

To install the package without `sudo` right, you can install it with the
`--user` flag.

```
python setup.py install --user
```

### pyjeo in Docker

A [Dockerfile](https://github.com/ec-jrc/jeolib-pyjeo/blob/master/docker/Dockerfile_deb10_pyjeo_public) based on a debian10 image is provided under the docker directory in this repository

Create the pyjeo docker image using:
```
docker build -t deb10_pyjeo_public:0.1.4 -f Dockerfile_deb10_pyjeo_public .
```

Run the pyjeo docker image and execute:
```
docker run --rm deb10_pyjeo_public:0.1.4 python3 -c "import pyjeo as pj; jim = pj.Jim(ncol = 10, nrow = 10, nband = 3); print(jim.properties.nrOfBand())"
```

## Test the installation

From the directory of the repository, run:

```
python -W ignore -m unittest -v tests
```

To run tests only for one module:

```
python -W ignore -m unittest -v tests/test_classify.py
```

## Documentation

An online documentation is available [here](https://jeodpp.jrc.ec.europa.eu/services/processing/pyjeohelp/)

To build the latest documentation from source:

- Install these dependencies:

* python3-sphinx
* sphinx-rtd-theme
* sphinxcontrib-bibtex

Go to directory `doc` and run `make html`.

```
cd doc
make html
```

The documentation is generated in html format in `_build/html` and can be read with your browser (open `index.html`).

A documentation in pdf can be obtained via:

```
cd doc
make latex
```

The documentation is generated in pdf format in `_build/latex`.

## Versions

### Getting the right version

`master` branch is the development branch of `pyjeo`. It contains the newest
features, but it can also contain some API changes against previous versions.
Therefore, it is recommended to use more stable releases: To get them, please
see our [tags](../../tags).

### Development

All development should be done in the development branch (`master`). If
a commit fixes also an issue present in some of the releases, it should be
cherry-picked to the corresponding branch. Later, a patch version will be
released from these cherry-picked fixes.

An example showing how to cherry-pick commit `commit_hash` into branch
`0.5.x` (to fix an issue in release `0.5.0`):

```
git checkout 0.5.x
git cherry-pick commit_hash
git push
git checkout master
```

## See the code coverage

```
python -W ignore -m coverage run --source=pyjeo -m unittest tests
python -m coverage report -m
```
