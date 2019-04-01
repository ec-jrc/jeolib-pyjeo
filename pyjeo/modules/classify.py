"""Module for operations connected to classification."""

import pyjeo as _pj


def classify(jim_object, method, model, **kwargs):
    """Supervised classification of a raster dataset.

    The classifier must have been trained via the train() method.
    The classifier can be selected with the key 'method'.

    :param jim_object: a Jim object

    :param method: Classification method (:ref:`svm <svm>` or :ref:`ann <ann>`)
    :param model: Model filename for trained classifier
    :param kwargs: See below
    :return: a Jim object with the classification result

    :keyword arguments:
        :band: (list of) band index(indices), starting from 0. The band order must correspond
            to the band names defined in the model. Leave empty to use all bands
        :srcnodata: No data values in input image not to consider for classification
        :dstnodata: No data value to put where input image is no data
        :priors: list of prior probabilities (one for each class)

    see :py:meth:`~_Classify.classify` for an example how to use this function
    """
    kwargs.update({'method': method, 'model': model})
    return _pj.Jim(jim_object._jipjim.classify(kwargs))


def sml(jim_object, **kwargs):
    """Supervised classification of a Jim object, using symbolic machine learning algorithm SML.
    For training, one or more reference raster datasets with categorical values is expected as
    a JimList. The reference raster dataset is typically at a lower spatial
    resolution than the input raster dataset to be classified. Unlike
    the Jim.classify(), the training is performed not prior to the
    classification, but in the same process as the classification.

    Modifies the instance on which the method was called.

    :param jim_object: a Jim object
    :param kwargs: See below
    return: multi-band Jim object, where each band represents the probability for each class

    :keyword arguments:
        :band: Band index (starting from 0). The band order must correspond
            to the band names defined in the model. Leave empty to use all
            bands
        :classes: List of classes to extract from the reference. Leave empty
            to extract two classes only (1 against rest)

    see :py:meth:`~_Classify.sml` for an example how to use this function
    """
    kwargs.update({'method': 'sml'})
    if 'classes' in kwargs:
        kwargs.update({'class': kwargs.pop('classes')})
    return _pj.Jim(jim_object._jipjim.classify(kwargs))


def reclass(jim_object, classes, reclasses, otype=None):
    """Reclassify a Jim object, replacing all pixels in the set classes to the corresponding values in reclasses

    :param jim_object: a Jim object
    :param classes: list of source values that need to be replaced
    :param reclasses: list of target values to which the pixels should be replaced
    return: Jim object with replaced values

    see :py:meth:`~_Classify.reclass` for an example
    """
    kwargs = {'class': classes, 'reclass': reclasses}
    if otype:
        kwargs.update({'otype': otype})
    retJim=_pj.Jim(jim_object)
    retJim._jipjim.d_reclass(kwargs)
    return retJim


class _Classify():

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def trainSML(self, jimlist, output, **kwargs):
        """Train a supervised symbolic machine learning (SML) classifier based on a reference :py:class:`JimList` object

        :param jimlist: list of Jim objects (can be a list of a single Jim object) containing reference raster datasets
        :param output: output filepath where model will be written
        :param kwargs: See below

        keyword arguments:

        =========  =====================================================================================================================
        band       (list of) band index (indices) (starting from 0)
        classes    list of classes to extract from reference. Leave empty to extract two classes only: 1 against rest
        srcnodata  No data value in source (do not consider for training)
        =========  =====================================================================================================================

        Example: train a Jim object (jim) with four Sentinel-2 bands using a reference image (jim_ref).

          jim.classify.trainSML(reflist,output='/path/to/model.txt',classes=sorted(classDict.values()))

        Use :py:meth:`~_Classify.sml` to perform the classification


        Typically, the SML algorithm performs a discretization, in order to reduce the dynamic range of the input and limit the number of unique sequences. For instance, to reduce the dynamic range of the input to NBIT bits::

          NBIT=3
          jim+=2**(8-NBIT)-jim.stats.getStats(['min'])['min']
          jim*=(2**8-1)/jim.stats.getStats(['max'])['max']
          jim>>=8-NBIT

        .. note::
          To train a SVM or ANN classifier, use :py:meth:`~_ClassifyVect.train`

        """
        #convert list of Jim to JimList
        if not isinstance(jimlist, _pj.JimList):
            jimlist=_pj.JimList(jimlist)
        kwargs.update({'method': "sml", 'model': output})
        if 'classes' in kwargs:
            kwargs.update({'class': kwargs.pop('classes')})
        self._jim_object._jipjim.train(jimlist._jipjimlist,kwargs)

    def classify(self, method, model, **kwargs):
        """Supervised classification of a raster dataset.

        The classifier must have been trained via the train() method.
        The classifier can be selected with the key 'method'.

        Modifies the instance on which the method was called.

        :param method: Classification method (:ref:`svm <svm>` or :ref:`ann <ann>`)
        :param model: Model filename for trained classifier
        :param kwargs: See below

        :keyword arguments:
            :band: (list of) band index(indices), starting from 0. The band order must correspond
                to the band names defined in the model. Leave empty to use all bands
            :srcnodata: No data values in input image not to consider for classification
            :dstnodata: No data value to put where input image is no data
            :priors: list of prior probabilities (one for each class)

        Perform a SVM classification using the model obtained via :py:meth:`~_ClassifyVect.train`::

          jim_sml=pj.classify.classify('svm', jim,model='/path/to/model.txt')

        Perform an ANN classification using the model obtained via :py:meth:`~_ClassifyVect.train`::

          jim_sml=pj.classify.classify('ann', jim,model='/path/to/model.txt')
        """
        kwargs.update({'method': method, 'model': model})
        self._jim_object._set(self._jim_object._jipjim.classify(kwargs))

    def sml(self, **kwargs):
        """Supervised classification of a Jim object, using symbolic machine learning algorithm SML.
        For training, one or more reference raster datasets with categorical values is expected as
        a JimList. The reference raster dataset is typically at a lower spatial
        resolution than the input raster dataset to be classified. Unlike
        the Jim.classify(), the training is performed not prior to the
        classification, but in the same process as the classification.

        Modifies the instance on which the method was called.

        :param kwargs: See below
        return: multi-band Jim object, where each band represents the probability for each class

        :keyword arguments:
            :band: Band index (starting from 0). The band order must correspond
                to the band names defined in the model. Leave empty to use all
                bands
            :classes: List of classes to extract from the reference. Leave empty
                to extract two classes only (1 against rest)

        Perform a SML classification using the model obtained via :py:meth:`~_Classify.trainSML`::

          jim_sml=pj.classify.sml(jim,model='/path/to/model.txt')

        The result is a multi-band :py:class`Jim` object where the number of bands equals the number of classes and each band represents the probability for the respective class. To create a discrete classification result, based on the maximum probability for each class::

          jim_sml._jipjim.band2plane()
          smlclassnp=np.argmax(jim_sml.np(),axis=0).astype(np.uint8)
          sml_class=pj.np2jim(smlclassnp)
          sml_class.properties.setProjection(jim.properties.getProjection())
          sml_class.properties.setGeoTransform(jim.properties.getGeoTransform())

        The result contains the indices in the range(0,number of classes). Use :py:meth:`~_Classify.reclass` to convert the indices to the actual class numbers.

        """
        kwargs.update({'method': 'sml'})
        if 'classes' in kwargs:
            kwargs.update({'class': kwargs.pop('classes')})
        self._jim_object._set(self._jim_object._jipjim.classify(kwargs))

    def reclass(self, classes, reclasses, otype=None):
        """Reclassify a Jim object, replacing all pixels in the set classes to the corresponding values in reclasses

        :param classes: list of source values that need to be replaced
        :param reclasses: list of target values to which the pixels should be replaced

        Reclass a Jim object, replacing all values in [0,1,2] to [250,251,252]::

          jim.classify.reclass(classes=[0,1,2],reclasses[250,251,252])
        """
        kwargs = {'class': classes, 'reclass': reclasses}
        if otype:
            kwargs.update({'otype': otype})
        self._jim_object._jipjim.d_reclass(kwargs)


class _ClassifyList():

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class _ClassifyVect():

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    def train(self, method, output, **kwargs):
        """Train a supervised classifier based on extracted data including
        label information (typically obtained via :py:func:`geometry:extractOgr`).

        :param method: Classification method (:ref:`svm <svm>` or :ref:`ann <ann>`)
        :param output: output filepath where model will be written
        :param kwargs: See below

        keyword arguments:

        ======== =====================================================================================================================
        label    Attribute name for class label in training vector file (default: 'label')
        bandname List of band names to use that correspond to the fields in the vector dataset. Leave empty to use all bands
        class    List of alpha numeric class names as defined in the label attribute (use only if labels contain not numerical values)
        reclass  List of numeric class values corresponding to the list defined by the class key
        ======== =====================================================================================================================

        .. note::
          To train a symbolic machine learning classifier, use :py:meth:`~_Classify.trainSML`

        **Balancing the training sample**

        Keys used to balance the training sample:

        ======== ================================================================================================
        balance  Balance the input data to this number of samples for each class
        random   Randomize training data for balancing
        min      Set to a value to not take classes into account with a sample size that is lower than this value
        ======== ================================================================================================


        .. _svm:

        **Support vector machine**

        The support vector machine (SVM) supervised classifier is described `here <http://dx.doi.org/10.1007/BF00994018>`_. The implementation in JIPlib is based on the open source `libsvm <https://www.csie.ntu.edu.tw/~cjlin/libsvm/>`_.

        Keys specific to the SVM:

        ========== ======================================================================
        svmtype    Type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)","C_SVC")
        kerneltype Type of kernel function (linear,polynomial,radial,sigmoid) ","radial")
        kd         Degree in kernel function",3)
        gamma      Gamma in kernel function",1.0)
        coef0      Coef0 in kernel function",0)
        ccost      The parameter C of C_SVC, epsilon_SVR, and nu_SVR",1000)
        nu         The parameter nu of nu_SVC, one_class SVM, and nu_SVR",0.5)
        eloss      The epsilon in loss function of epsilon_SVR",0.1)
        cache      Cache memory size in MB",100)
        etol       The tolerance of termination criterion",0.001)
        shrink     Whether to use the shrinking heuristics",false)
        probest    Whether to train a SVC or SVR model for probability estimates",true,2)
        ========== ======================================================================

        Extract training data from a sample containing labeled features (fieldname is 'label') and write the result to a vector in memory (using the method :py:meth:`~geometry._Geometry.extractOgr` on a :py:class:`Jim` object in module :py:mod:`geometry`). The field 'label' is copied from the sample that will be used by the training. Then use the extracted vector to train a SVM and write the model to a file.::

           training=jim.extractOgr(sample,'training','oformat':'Memory','copy':'label'})
           training.train({'method':'svm','label':'label','model':args.model})

        Use :py:meth:`~_Classify.classify` to perform the classification

        .. _ann:

        **Artificial neural network**

        The artificial neural network (ANN) supervised classifier is based on the back propagation model as introduced by D. E. Rumelhart, G. E. Hinton, and R. J. Williams (Nature, vol. 323, pp. 533-536, 1986). The implementation is based on the open source C++ library `fann <http://leenissen.dk/fann/wp/>`_.


        Keys specific to the ANN:

        ========== ==========================================================================
        nneuron    List defining the number of neurons in each hidden layer in the neural network
        connection Connection rate (default: 1.0 for a fully connected network
        learning   Learning rate (default: 0.7)
        weights    Weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)
        maxit      Maximum epochs used for training the neural network (default: 500)
        ========== ==========================================================================

        .. note::
          To define two hidden layers with 3 and 5 neurons respectively, define a list of two values for the key 'nneuron': [3, 5].

        Extract training data from a sample containing labeled features (fieldname is 'label') and write the result to a vector in memory (using the method :py:meth:`~geometry._Geometry.extractOgr` on a :py:class:`Jim` object in module :py:mod:`geometry`). The field 'label' is copied from the sample that will be used by the training. Then use the extracted vector to train an ANN and write the model to a file.::

           training=jim.extractOgr(sample,'training','oformat':'Memory','copy':'label'})
           training.train({'method':'ann','label':'label','model':args.model})

        Use :py:meth:`~_Classify.classify` to perform the classification
        """
        kwargs.update({'method': method, 'model': output})
        self._jim_vect._jipjimvect.train(kwargs)


    def classify(self, method, model, **kwargs):
        """Supervised classification of a raster dataset.

        The classifier must have been trained via the train() method.
        The classifier can be selected with the key 'method'.

        Modifies the instance on which the method was called.

        :param method: Classification method ('svm', 'ann')
        :param model: Model filename for trained classifier
        :param kwargs: See below

        :keyword arguments:
            :band: Band index (starting from 0). The band order must correspond
                to the band names defined in the model. Leave empty to use all
                bands
            :priors: list of prior probabilities (one for each class)
            :output: output filename of classified vector dataset
            :f: output filename of classified vector dataset
            :co: creation option for output file. Multiple options can be
                specified as list
            :copy: copy these fields from input to output vector dataset
        """
        kwargs.update({'method': method, 'model': model})
        self._jim_vect._set(self._jim_vect._jipjimvect.classify(kwargs))