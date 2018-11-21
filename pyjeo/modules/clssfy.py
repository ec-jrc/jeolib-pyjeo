try:
    import pyjeo as _pj
except:
    try:
        from jeodpp import pyjeo as _pj
    except:
        import jeodpp.pyjeo as _pj


def classify(jim_object, method, model, **kwargs):
    """Supervised classification of a raster dataset.

    The classifier must have been trained via the train() method.
    The classifier can be selected with the key 'method'.

    :param jim_object: a Jim object
    :param method: Classification method (svm or ann)
    :param model: Model filename to save trained classifier
    :param kwargs: See below
    :return: a Jim object

    :keyword arguments:
        :band: Band index (starting from 0). The band order must correspond
            to the band names defined in the model. Leave empty to use all
            bands
        :extent: Data type for output image
        :eo: Special extent options controlling rasterization
        :mask: Only classify within specified mask
        :msknodata: Mask value(s) in mask not to consider for
            classification
        :nodata: Nodata value to put where image is masked as no data
    """
    kwargs.update({'method': method, 'model': model})
    return _pj.Jim(jim_object.classify(kwargs))


def sml(jim_object, **kwargs):
    """Supervised classification of a raster dataset.

    Using symbolic machine learning algorithm SML. For training, one or
    more reference raster datasets with categorical values is expected as
    a JimList. The reference raster dataset is typically at a lower spatial
    resolution than the input raster dataset to be classified. Unlike
    the Jim.classify(), the training is performed not prior to the
    classification, but in the same process as the classification.

    :param jim_object: a Jim object
    :param kwargs: See below
    :return: a Jim object

    :keyword arguments:
        :band: Band index (starting from 0). The band order must correspond
            to the band names defined in the model. Leave empty to use all
            bands
        :class: List of classes to extract from the reference. Leave empty
            to extract two classes only (1 against rest)
        :otype: Data type for output image
    """
    kwargs.update({'method': 'sml'})
    return _pj.Jim(jim_object.classify(kwargs))


class _Classify():
    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def classify(self, method, model, **kwargs):
        """Supervised classification of a raster dataset.

        The classifier must have been trained via the train() method.
        The classifier can be selected with the key 'method'.

        Modifies the instance on which the method was called.

        :param method: Classification method (svm or ann)
        :param model: Model filename to save trained classifier
        :param kwargs: See below

        :keyword arguments:
            :band: Band index (starting from 0). The band order must correspond
                to the band names defined in the model. Leave empty to use all
                bands
            :extent: Data type for output image
            :eo: Special extent options controlling rasterization
            :mask: Only classify within specified mask
            :msknodata: Mask value(s) in mask not to consider for
                classification
            :nodata: Nodata value to put where image is masked as no data
        """
        kwargs.update({'method': method, 'model': model})
        self._jim_object._set(self._jim_object.classify(kwargs))

    def sml(self, **kwargs):
        """Supervised classification of a raster dataset.

        Using symbolic machine learning algorithm SML. For training, one or
        more reference raster datasets with categorical values is expected as
        a JimList. The reference raster dataset is typically at a lower spatial
        resolution than the input raster dataset to be classified. Unlike
        the Jim.classify(), the training is performed not prior to the
        classification, but in the same process as the classification.

        Modifies the instance on which the method was called.

        :param kwargs: See below

        :keyword arguments:
            :band: Band index (starting from 0). The band order must correspond
                to the band names defined in the model. Leave empty to use all
                bands
            :class: List of classes to extract from the reference. Leave empty
                to extract two classes only (1 against rest)
            :otype: Data type for output image
        """
        kwargs.update({'method': 'sml'})
        self._jim_object._set(self._jim_object.classify(kwargs))
