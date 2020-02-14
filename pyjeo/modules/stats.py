"""Module for statistical functions and interpolations."""

import numpy as _np

import pyjeo as _pj


def getHisto1d(jim_object):
    """Compute the frequency distribution of the grey levels of im.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo1d())


def getHisto2d(jim_object, jim_object2):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo2d(jim_object2._jipjim))


def getHisto3d(jim_object, jim_object2, jim_object3):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :param jim_object3: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo3d(jim_object2._jipjim,
                                              jim_object3._jipjim))


def getStatProfile(jim_object, function, **kwargs):
    kwargs.update({'function': function})
    return _pj.Jim(jim_object._jipjim.statProfile(kwargs))


# cannot overload?
# def getStatProfile(jimlist, function, **kwargs):
#     kwargs.update({'function': function})
#     return _pj.Jim(jimlist._jim_list._jipjimlist.statProfile(kwargs))


def getStats(jim_object, function=['min', 'max', 'mean'], **kwargs):
    """Compute basic statistics on a JimList object.

    Similar to the :py:meth:`~._Stats.getStats` method from Jim when
    jim_object is an instance of Jim, similar to
    the :py:meth:`~._StatsList.getStats` method from JimList when
    jim_object is an instance of JimList. For functions requiring two
    datasets (e.g., regression), use the objects in the list instead of bands.

    :param jim_object: either a Jim object or a JimList object
    :param function: (list of) statistical function(s) to calculate
        (default is ['min', 'max', 'mean'])
    :return: a dictionary with requested statistics
    """
    return jim_object.stats.getStats(function, **kwargs)


# cannot overload?
# def getStats(jimlist, function=['min','max','mean'], **kwargs):
#     kwargs.update({'function': function})
#     return jimlist._jim_list._jipjimlist.getStats(kwargs)


def stretch(jim_object, **kwargs):
    return _pj.Jim(jim_object._jipjim.stretch(kwargs))


class _Stats(_pj.modules.JimModuleBase):
    """Define all statistical methods."""

    def getHisto1d(self):
        """Compute the frequency distribution of the grey levels of im.

        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo1d())

    def getHisto2d(self, jim_object2):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo2d(jim_object2._jipjim))

    def getHisto3d(self, jim_object2, jim_object3):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :param jim_object3: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo3d(jim_object2._jipjim,
                                                        jim_object3._jipjim))

    def getHistoCumulative(self):
        """Compute the cumulative frequency distribution of the grey levels.

        :return: a Jim object
        """
        if self._jim_object.properties.getDataType() not in ['UInt32',
                                                             'Int32']:
            raise TypeError('Object must be of type UInt32 or Int32')
        return _pj.Jim(self._jim_object._jipjim.histo1dCumulative())

    def getStats(self, function=['min', 'max', 'mean'], **kwargs):
        """Compute basic statistics on a Jim object.

        For functions requiring two datasets (e.g., regression), use
        a multi-band Jim object (or use the :py:meth:`~_StatsList.getStats`
        method from JimList.

        :param function: (list of) statistical function(s) to calculate
            (default is ['min', 'max', 'mean'])
        :return: a dictionary with requested statistics

        .. _statFunctions:

        :statistical functions:

        +---------------------+-----------------------------------------------+
        | function            | description                                   |
        +=====================+===============================================+
        | invalid             | Count number of invalid pixels (nodata)       |
        +---------------------+-----------------------------------------------+
        | valid               | Count number of valid pixels (nodata)         |
        +---------------------+-----------------------------------------------+
        | basic               | Calculate min, max, mean, and stdev           |
        +---------------------+-----------------------------------------------+
        | mean                | Calculate the mean                            |
        +---------------------+-----------------------------------------------+
        | median              | Calculate the median                          |
        +---------------------+-----------------------------------------------+
        | var                 | Calculate the variance                        |
        +---------------------+-----------------------------------------------+
        | skewness            | Calculate the skewness as measure of assym.   |
        +---------------------+-----------------------------------------------+
        | kurtosis            | Calculate the kurtosis as measure of outliers |
        +---------------------+-----------------------------------------------+
        | stdev               | Calculate the standard deviation              |
        +---------------------+-----------------------------------------------+
        | minmax              | Calculate min and max                         |
        +---------------------+-----------------------------------------------+
        | min                 | Calculate min                                 |
        +---------------------+-----------------------------------------------+
        | min                 | Calculate max                                 |
        +---------------------+-----------------------------------------------+
        | histogram           | Calculate the historgram                      |
        +---------------------+-----------------------------------------------+
        | histogram2d         | Calculate the historgram based on 2 bands     |
        +---------------------+-----------------------------------------------+
        | rmse                | Calculate the RMSE based on 2 bands           |
        +---------------------+-----------------------------------------------+
        | regression          | Calculate the regression based on 2 bands     |
        +---------------------+-----------------------------------------------+

        Examples:

        Calculate min and max::

            jim = pj.Jim('/path/to/raster.tif')
            statDict = jim.stats.getStats(['min', 'max'])
            print('max value is: {}".format(statDict['max']))
            print('min value is: {}".format(statDict['min']))

        Calculate the histogram
        (returning a dictionary with keys 'bin' and 'histogram')::

            jim = pj.Jim('/path/to/raster.tif')
            histDict = jim.stats.getStats('histogram')

        Calculate the histogram, using 10 bins. Start reading from value 0 and
        stop reading up to value 100. Do not consider values equal to 0::

            jim = pj.Jim('/path/to/raster.tif')
            histDict = jim.stats.getStats('histogram', nbnin=10, src_min=0,
                                          src_max=100, nodata=0)
            histDict
            {'bin': [5.0, 15.0, 25.0, 35.0, 45.0,
                     55.0, 65.0, 75.0, 85.0, 95.0],
            'histogram': [38543.0,
                          29991.0,
                          49251.0,
                          40006.0,
                          12945.0,
                          2520.0,
                          112.0,
                          0.0,
                          0.0,
                          0.0]}

        """
        if not isinstance(function, list):
            function = function.split(',')

        statDict = dict()

        forceJiplib = False
        constraints = ('nodata', 'src_min', 'src_max')
        for key in constraints:
            if key in kwargs.keys():
                forceJiplib = True
                break
        if not forceJiplib:
            bands = []
            if 'band' in kwargs:
                try:
                    bands.extend(kwargs['band'])
                except TypeError:
                    bands.append(kwargs['band'])
            else:
                bands = range(0, self._jim_object.properties.nrOfBand())

            for f in function:
                statDict.update({f: []})
            for band in bands:
                if 'min' in function or 'max' in function:
                    min_max = self._jim_object._jipjim.getMiaMinMax(band)
                    if 'min' in function:
                        statDict['min'].append(min_max[1])
                    if 'max' in function:
                        statDict['max'].append(min_max[2])
                if 'mean' in function:
                    statDict['mean'].append(
                        _np.mean(self._jim_object.np(band)).item())
                if 'median' in function:
                    statDict['median'].append(
                        _np.median(self._jim_object.np(band)).item())

            for f in statDict:
                if len(statDict[f]) == 1:
                    statDict[f] = statDict[f][0]

            function_list = list()
            for f in function:
                if f not in ['min', 'max', 'mean', 'median']:
                    function_list.append(f)

            if len(function_list) > 0:
                kwargs.update({'function': function_list})
            if 'function' not in kwargs:
                return statDict
        else:
            kwargs.update({'function': function})

        if kwargs:
            statDict.update(self._jim_object._jipjim.getStats(kwargs))

        return statDict

    def statProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        self._jim_object._set(self._jim_object._jipjim.statProfile(kwargs))

    def stretch(self, **kwargs):
        self._jim_object._set(self._jim_object._jipjim.stretch(kwargs))


class _StatsList(_pj.modules.JimListModuleBase):
    """Define all statistical methods for JimLists."""

    def getStatProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        return _pj.Jim(self._jim_list._jipjimlist.statProfile(kwargs))

    def getStats(self, function=['min', 'max', 'mean'], **kwargs):
        """Compute basic statistics on a JimList object.

        Similar to the :py:meth:`~._Stats.getStats` method from Jim.
        For functions requiring two datasets (e.g., regression), use
        the objects in the list instead of bands

        :param function: (list of) statistical function(s) to calculate
            (default is ['min', 'max', 'mean'])
        :return: a dictionary with requested statistics

        .. _statFunctions_list:

        :statistical functions:

        +---------------------+-----------------------------------------------+
        | function            | description                                   |
        +=====================+===============================================+
        | invalid             | Count number of invalid pixels (nodata)       |
        +---------------------+-----------------------------------------------+
        | valid               | Count number of valid pixels (nodata)         |
        +---------------------+-----------------------------------------------+
        | basic               | Calculate min, max, mean, and stdev           |
        +---------------------+-----------------------------------------------+
        | mean                | Calculate the mean                            |
        +---------------------+-----------------------------------------------+
        | median              | Calculate the median                          |
        +---------------------+-----------------------------------------------+
        | var                 | Calculate the variance                        |
        +---------------------+-----------------------------------------------+
        | skewness            | Calculate the skewness as measure of assym.   |
        +---------------------+-----------------------------------------------+
        | kurtosis            | Calculate the kurtosis as measure of outliers |
        +---------------------+-----------------------------------------------+
        | stdev               | Calculate the standard deviation              |
        +---------------------+-----------------------------------------------+
        | minmax              | Calculate min and max                         |
        +---------------------+-----------------------------------------------+
        | min                 | Calculate min                                 |
        +---------------------+-----------------------------------------------+
        | min                 | Calculate max                                 |
        +---------------------+-----------------------------------------------+
        | histogram           | Calculate the historgram                      |
        +---------------------+-----------------------------------------------+
        | histogram2d         | Calculate the historgram based on 2 bands     |
        +---------------------+-----------------------------------------------+
        | rmse                | Calculate the RMSE based on 2 bands           |
        +---------------------+-----------------------------------------------+
        | regression          | Calculate the regression based on 2 bands     |
        +---------------------+-----------------------------------------------+

        Examples:

        Calculate regression between two Jim objects::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            pj.JimList([jim0,jim1]).stats.getStats('regression)
            {'c0': 10.0102, 'c1': 0.633352, 'r2': 0.491198}

        Calculate root mean square error between two Jim objects::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            pj.JimList([jim0, jim1]).stats.getStats('rmse')
            {'rmse': 10.4638}

        """
        keys = list(kwargs.keys())

        if not isinstance(function, list):
            function = function.split(',')

        if len(function) > 0:
            keys.extend(function)

        statDict = {key: None for key in keys}

        forceJiplib = False
        constraints = ('nodata', 'src_min', 'src_max')
        for key in constraints:
            if key in keys:
                forceJiplib = True
                break

        if function not in [['min', 'max'], ['max', 'min'], 'min', 'max']:
            forceJiplib = True

        if not forceJiplib:
            for item in self._jim_list:
                new_statDict = item.stats.getStats(function, **kwargs)
                if 'min' in keys and (statDict['min'] is None or
                                      new_statDict['min'] < statDict['min']):
                    statDict['min'] = new_statDict['min']
                if 'max' in keys and (statDict['max'] is None or
                                      new_statDict['max'] > statDict['max']):
                    statDict['max'] = new_statDict['max']
        else:
            kwargs.update({'function': function})
            statDict.update(self._jim_list._jipjimlist.getStats(kwargs))

        return statDict


class _StatsVect(_pj.modules.JimVectModuleBase):
    """Define all statistical methods for JimVects."""

    pass
