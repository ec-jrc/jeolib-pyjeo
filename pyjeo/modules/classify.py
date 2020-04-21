"""Module for operations connected to classification."""

import pyjeo as _pj


def reclass(jim_object,
            classes: list,
            reclasses: list,
            otype=None):
    """Reclassify a Jim object, replacing all pixels in the set classes.

    Replace all pixels in the set class to the corresponding values in
    reclasses

    :param jim_object: a Jim object
    :param classes: list of source values that need to be replaced
    :param reclasses: list of target values to which the pixels should be
        replaced
    :param otype: Data type
    :return: Jim object with replaced values

    see :py:meth:`~_Classify.reclass` for an example
    """
    kwargs = {'class': classes, 'reclass': reclasses}
    if otype:
        kwargs.update({'otype': otype})
    ret_jim = _pj.Jim(jim_object)
    ret_jim._jipjim.d_reclass(kwargs)
    return ret_jim


def sml(jim_object,
        reflist=None,
        classes: list = None,
        model: str = None,
        **kwargs):
    """Perform supervised classification of a Jim object using SML.

    For training, one or more reference raster datasets with categorical
    values is expected as a JimList. The reference raster dataset is
    typically at a lower spatial resolution than the input raster dataset
    to be classified. Unlike the Jim.classify(), the training is
    performed not prior to the classification, but in the same process as
    the classification.

    :param jim_object: a multi-plane Jim object
    :param reflist: JimList of reference raster datasets containing with
        reference classes
    :param classes: List of classes to extract from the reference
        (leave empty to extract all classes in reference)
    :param model: Model filename for trained classifier
    :return: multi-band Jim object, where each band represents the probability
        for each class

    see :py:meth:`~_Classify.sml` for an example how to use this function
    """
    if model is not None:
        kwargs.update({'method': 'sml', 'model': model})
        return _pj.Jim(jim_object._jipjim.classify(kwargs))
    elif reflist is not None:
        if classes is not None:
            kwargs.update({'class': classes})
        return _pj.Jim(jim_object._jipjim.classifySML(
            reflist, kwargs))
    else:
        raise ValueError('No training found through model or reflist')


class _Classify(_pj.modules.JimModuleBase):
    """Define all classification methods for Jims."""

    def reclass(self,
                classes: list,
                reclasses: list,
                otype=None):
        """Reclassify a Jim object.

        Replace all pixels in the set classes to the corresponding values in
        reclasses.

        :param classes: list of source values that need to be replaced
        :param reclasses: list of target values to which the pixels should be
            replaced
        :param otype: Data type

        Reclass a Jim object, replacing all values in [0,1,2] to
        [250,251,252]::

          jim.classify.reclass(classes=[0,1,2],reclasses[250,251,252])
        """
        kwargs = {'class': classes, 'reclass': reclasses}
        if otype:
            kwargs.update({'otype': otype})
        self._jim_object._jipjim.d_reclass(kwargs)

    def sml(self,
            reflist=None,
            classes: list = None,
            model: str = None,
            **kwargs):
        """Perform supervised classification of a Jim object using SML.

        For training, one or more reference raster datasets with categorical
        values is expected as a JimList. The reference raster dataset is
        typically at a lower spatial resolution than the input raster dataset
        to be classified. Unlike the Jim.classify(), the training is
        performed not prior to the classification, but in the same process as
        the classification.

        Modifies the instance on which the method was called.

        :param reflist: JimList of reference raster datasets containing with
            reference classes
        :param model: Model filename for trained classifier
        :param reflist: JimList of reference raster datasets containing with
            reference classes
        :param classes: list of classes to extract from the reference.
            (leave empty to extract all classes in reference)

        ::

          jim_sml = pj.classify.sml(reference, classes=[0,1,2], jim)

        The result is a multi-band :py:class:`Jim` object where the number of
        bands equals the number of classes and each band represents
        the probability for the respective class. To create a discrete
        classification result, based on the maximum probability for each
        class::

          jim_sml.geometry.band2plane()
          jim_sml.np()[0] = np.argmax(jim_sml.np(), axis=0)
          jim_sml.geometry.cropPlane(0)

        The result contains the indices in the range(0, number of classes).
        Use :py:meth:`~_Classify.reclass` to convert the indices to the actual
        class numbers.
        """
        if model is not None:
            kwargs.update({'method': 'sml', 'model': model})
            self._jim_object._set(self._jim_object._jipjim.classify(kwargs))
        elif reflist is not None:
            if classes is not None:
                kwargs.update({'class': classes})
            self._jim_object._set(self._jim_object._jipjim.classifySML(
                reflist._jipjimlist, kwargs))
        else:
            raise ValueError('No training found through model or reflist')

    def trainSML(self,
                 reference,
                 output: str = None,
                 **kwargs):
        """Train a supervised symbolic machine learning (SML) classifier.

        Train it based on a reference :py:class:`JimList` object

        :param reference: (list of) Jim objects containing reference raster
            dataset(s)
        :param output: output filepath where trained model will be written
            (leave empty to return as a string)
        :param kwargs: See below

        keyword arguments:

        =========  ============================================================
        band       (list of) band index (indices) (starting from 0)
        classes    list of classes to extract from reference. Leave empty to
                   extract two classes only: 1 against rest
        srcnodata  No data value in source (do not consider for training)
        =========  ============================================================

        Example: train a Jim object (jim) with four Sentinel-2 bands using
        a single reference image (e.g., corine land cover)::

          corine = pj.Jim('/path/to/corine.tif')
          jim.classify.trainSML(
              corine, output='/path/to/model.txt', classes=[12, 40])

        Use :py:meth:`~_Classify.sml` to perform the classification


        Typically, the SML algorithm performs a discretization, in order to
        reduce the dynamic range of the input and limit the number of unique
        sequences. For instance, to reduce the dynamic range of the input to
        NBIT bits::

          NBIT=3
          jim += 2 ** (8 - NBIT) - jim.stats.getStats(['min'])['min']
          jim *= (2 ** 8 - 1) / jim.stats.getStats(['max'])['max']
          jim >>= 8-NBIT

        """
        # convert list of Jim to JimList
        if not isinstance(reference, _pj.JimList):
            reference = _pj.JimList([reference])
        kwargs.update({'method': "sml"})
        if 'classes' in kwargs:
            kwargs.update({'class': kwargs.pop('classes')})
        if output:
            kwargs.update({'model': output})
            self._jim_object._jipjim.trainSML(reference._jipjimlist, kwargs)
        else:
            raise ValueError('Error: output for model not set')


class _ClassifyList(_pj.modules.JimListModuleBase):
    """Define all classification methods for JimLists."""

    pass


class _ClassifyVect(_pj.modules.JimVectModuleBase):
    """Define all classification methods for JimVects."""
