# pyJEO

## Description

pyJEO is a library for image processing for geospatial data implemented in
JRC Ispra.

## License

pyJEO is released under an
[EUPL](https://joinup.ec.europa.eu/collection/eupl) license (see
[LICENSE.txt](LICENSE.txt))

## Dependencies

* mialib
* jiplib
* numpy

## Install

From the directory of the repository, run:

```
sudo python setup.py install
```

To install without sudo right, you can install with the --user

```
python setup.py install --user
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

## Build documentation

Dependencies for the documentation build:

* python3-sphinx
* sphinx-rtd-theme
* sphinxcontrib-bibtex

Go to directory `doc` and run `make html`.

```
cd doc
make html
```

## See the code coverage

```
python -W ignore -m coverage run --source=pyjeo -m unittest tests
python -m coverage report -m
```
