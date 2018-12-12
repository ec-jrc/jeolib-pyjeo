"""Module for operations connected to classification."""

import pyjeo as _pj


def classify(jim_object, method, model, **kwargs):
    """Supervised classification of a raster dataset.

    The classifier must have been trained via the train() method.
    The classifier can be selected with the key 'method'.

    :param jim_object: a Jim object
    :param method: Classification method (svm or ann)
    :param model: Model filename for trained classifier
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
        :priors: list of prior probabilities (one for each class)
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
    """
    kwargs.update({'method': 'sml'})
    return _pj.Jim(jim_object.classify(kwargs))


def reclass(jim_object, classes, reclasses, otype=None):
    kwargs = {'class': classes, 'reclass': reclasses}
    if otype:
        kwargs.update({'otype': otype})
    retJim=_pj.Jim(jim_object)
    retJim.d_reclass(kwargs)
    return retJim
    # return _pj.Jim(jim_object.reclass(kwargs))



class _Classify():
    """Define all classification methods."""

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
        :param model: Model filename for trained classifier
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
        """
        kwargs.update({'method': 'sml'})
        self._jim_object._set(self._jim_object.classify(kwargs))

    def reclass(self, classes, reclasses, otype=None):
        kwargs = {'class': classes, 'reclass': reclasses}
        if otype:
            kwargs.update({'otype': otype})
        # self._jim_object._set(self._jim_object.reclass(kwargs))
        self._jim_object.d_reclass(kwargs)


class _ClassifyList():
    """Define all classification methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class _ClassifyVect():
    """Define all classification methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    def train(self, method, output, **kwargs):
        """Train a supervised classifier based on extracted data including
        label information (typically obtained via :py:func:`geometry:extractOgr`).

        :param method: classification method: 'svm' (support vector machine),
           'ann' (artificial neural network)
        :param output: output filepath where model will be written
        :param kwargs: See below

        keyword arguments:

        ======== =====================================================================================================================
        label    Attribute name for class label in training vector file (default: 'label')
        bandname List of band names to use that correspond to the fields in the vector dataset. Leave empty to use all bands
        class    List of alpha numeric class names as defined in the label attribute (use only if labels contain not numerical values)
        reclass  List of numeric class values corresponding to the list defined by the class key
        ======== =====================================================================================================================

        **Balancing the training sample**

        Keys used to balance the training sample:

        ======== ================================================================================================
        balance  Balance the input data to this number of samples for each class
        random   Randomize training data for balancing
        min      Set to a value to not take classes into account with a sample size that is lower than this value
        ======== ================================================================================================

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

        **Artificial neural network**

        The artificial neural network (ANN) supervised classifier is based on the back propagation model as introduced by D. E. Rumelhart, G. E. Hinton, and R. J. Williams (Nature, vol. 323, pp. 533-536, 1986). The implementation is based on the open source C++ library fann (http://leenissen.dk/fann/wp/).


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

        """
        kwargs.update({'method': method, 'model': output})
        self._jim_vect.train(kwargs)

    def classify(self, method, output, **kwargs):
        """Supervised classification of a raster dataset.

        The classifier must have been trained via the train() method.
        The classifier can be selected with the key 'method'.

        Modifies the instance on which the method was called.

        :param method: Classification method (svm or ann)
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
        kwargs.update({'method': method, 'model': output})
        self._jim_vect._set(self._jim_vect.classify(kwargs))
