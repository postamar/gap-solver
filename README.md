(C) 2012, Marius Posta (mariusposta@gmail.com)

Generalized Assignment Problem solver
=====================================

Licensing
---------

All files in this repo, with the exception of minknap.c and minknap.h, are
distributed under the terms stated in LICENSE.txt.

The files minknap.c and minknap.h are the work of David Pisinger and fall under
a different license. Details within.

Description
-----------

See my paper, available at http://dx.doi.org/10.1007/s10589-011-9432-0

Usage
-----

The interface is rather crude. The executable requires the following 6 arguments:
1. path to problem instance
2. known upper bound (-1 if none)
3. maximum node lower bound evaluation iterations (100 is a good value)
4. time limit in seconds (-1 if none)
5. detailed progress log entry interval in seconds (0 for complete log)

The detailed progress log is printed to stderr and the summary log to stdout.
The summary log appears as follows:

instance <instance name>
<initial global lower bound> <time> <# of node lower bound iterations> <# of branch-and-bound nodes>
<same as above + 1> ...
<same as above + 2> ... 
...

