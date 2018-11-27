# pyJEO

pyJEO is a library for image processing for geospatial data implemented in 
JRC Ispra. 

# License

pyJEO is released under an
[EUPL](https://joinup.ec.europa.eu/collection/eupl) license (see
[LICENSE.txt](LICENSE.txt))

# Dependencies

 * mialib
 * jiplib
 * numpy

# Install

From the directory of the repository, run:
```
sudo python setup.py install
```

# Test the installation

From the directory of the repository, run:
```
python -m unittest -v tests
```

# Build documentation

Go to directory `doc` and run `make html`.
```
cd doc
make html
```
