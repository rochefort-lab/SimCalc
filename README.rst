SimCalc
=======
SimCalc (Simulated Calcium) simulates calcium imaging data, with the option for different neural shapes and locations.

This code was used for testing our FISSA toolbox, for more information see the FISSA `repository <https://github.com/rochefort-lab/fissa>`__ and `companion paper <https://www.doi.org/10.1038/s41598-018-21640-2>`__.

Currently this code can be quite memory heavy, so on standard laptops it's best not to define too big field-of-views or too long time periods.

Usage
-----
A general tutorial on the use of SimCalc can be found at:
https://rochefort-lab.github.io/SimCalc/examples/Simulated_data_tutorial.html

The corresponding jupyter notebook can be found in the examples folder of this repository.

Installation
------------
After downloading this repository, use a terminal to navigate into the base folder. Then: 

::

    pip install -e ./

Should do the trick.


Citing SimCalc
------------

If you use SimCalc for your research, please cite the following paper in
any resulting publications:

S. W. Keemink, S. C. Lowe, J. M. P. Pakan, E. Dylda, M. C. W. van
Rossum, and N. L. Rochefort. FISSA: A neuropil decontamination toolbox
for calcium imaging signals, *Scientific Reports*, **8**\ (1):3493,
2018.
`DOI:10.1038/s41598-018-21640-2 <https://www.doi.org/10.1038/s41598-018-21640-2>`__.

For your convenience, the FISSA package ships with a copy of this
citation in bibtex format, available at
`citation.bib <https://raw.githubusercontent.com/rochefort-lab/fissa/master/citation.bib>`__.

License
-------

Unless otherwise stated in individual files, all code is Copyright (c)
2015, Sander Keemink, Scott Lowe, and Nathalie Rochefort. All rights
reserved.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see http://www.gnu.org/licenses/.
