
PMT Waveform Analysis Library
=============================

WORK IN PROGRESS (USE AT YOUR OWN RISK)

This is a collection of analysis and simulation scripts developed to
facilitate the Short-Range Forces with Mossbauer Spectroscopy project 
at Stanford University, under the direction of Professor Giorgio 
Gratta. 

Currently, the software only has limited performance and can only
load data acquired by a specific digitizer (Caen DT5724B), but is 
still in active development

Install
-------

From sources
````````````

To install system-wide, noting the path to the src since no wheels
exist on PyPI, use::

   pip install ./pmt_waveform

If you intend to edit the code and want the import calls to reflect
those changes, install in developer mode::

   pip install -e pmt_waveform

If you don't want a global installation (i.e. if multiple users will
engage with and/or edit this library) and you don't want to use venv
or some equivalent::

   pip install -e pmt_waveform --user

where pip is pip3 for Python3 (tested on Python 3.6.9). Be careful 
NOT to use ``sudo``, as the latter two installations make a file
``easy-install.pth`` in either the global or the local directory
``lib/python3.X/site-packages/easy-install.pth``, and sudo will
mess up the permissions of this file such that uninstalling is very
complicated.


Uninstall
---------

If installed without ``sudo`` as instructed, uninstalling should be 
as easy as::

   pip uninstall pmt_waveform

If installed using ``sudo`` and with the ``-e`` and ``--user`` flags, 
the above uninstall will encounter an error.

Navigate to the file ``lib/python3.X/site-packages/easy-install.pth``, 
located either at  ``/usr/local/`` or ``~/.local`` and ensure there
is no entry for ``pmt_waveform``.


License
-------

The package is distributed under an open license (see LICENSE file for
information).


Authors
-------

Charles Blakemore (chas.blakemore@gmail.com),
Gautam Venugopalan