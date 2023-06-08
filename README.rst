Malbrid - Linear Hybrid System/Automaton Verification Library for Python
========================================================================

The Malbrid library is a Python library for simulating linear hybrid systems modelled as linear linear hybrid automata. 

While no full documentation is yet available, the example Jupyter notebooks show how to model automata and how to simulate them.


Table of contents
-----------------

.. contents:: 



Installation
------------
While the library is not available to be installed via PyPI, it can be installed by checking out the repository and running:

.. code-block:: raw

   pip install -e .


Getting started
---------------
Have a look at the examples in the "examples" directory which run with Jupyter Notebook. They describe a bouncing ball example as well as a ball-with-a-paddle example. The latter also shows how to construct the product of two linear hybrid automata on-the-fly.


Notes on Soundness
------------------
This examples of the library  contains a couple of checks that test if numerical imprecision during the simulation can cause the wrong discrete transitions to be taken. In such a case, an exception is thrown, and the library contains an exception class so that the respective exception can be caught in your own program.
