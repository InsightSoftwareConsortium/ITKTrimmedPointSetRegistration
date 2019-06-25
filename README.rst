ITKTrimmedPointSetRegistration
=================================

.. image:: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_apis/build/status/itktrimmedpointsetregistration?branchName=master
    :target: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_build/latest?definitionId=8&branchName=master
    :alt:    Build Status

.. image:: https://img.shields.io/pypi/v/itk-trimmedpointsetregistration.svg
    :target: https://pypi.python.org/pypi/itk-trimmedpointsetregistration
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/InsightSoftwareConsortium/ITKTrimmedPointSetRegistration/blob/master/LICENSE)
    :alt: License

Overview
--------

This is a module for trimmed point set registration.

This module implmenets a trimmed euclidean point set metric that implements a an adpative trimming of the point set used during the optimization in the registration process.

The itkTrimmedEuclideanDistancePointSetToPointSetMetricv4 is multithreaded and signficantly improves perfromance with the patch '<https://github.com/InsightSoftwareConsortium/ITK/pull/1037>'_ .

A Simple Example
----------------

Comparing Jensen, Euclidean and trimmed Euclinean metric on 2d toy example

Green - Input target point set

Black - Input source point set

Red -  Transformed source point set

.. image:: Documentation/Figures/jensen.png
.. image:: Documentation/Figures/euclidean.png
.. image:: Documentation/Figures/trimmed-euclidean.png


TODO
----

  - Implement multithreading.
  - Add a real world example.
  - Investigate Jensen unexpected behaviour in simple example.
  - Currently the point set metrics are not symmetric since the iterations is only over one fo the point sets. Consider adding symmetric versions. Or document the assymetric behaviour well.
   - A decorator to the PointSetToPointSetMetricv4 that overrides the accumulation of the value and derivative computation to use a trimmed number of points is in the module as well, but is as of yet not functional. The issue is access to protected members in the base class that I think would need to be fixed in the base class.


