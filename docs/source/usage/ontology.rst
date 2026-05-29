========================
Ontology
========================

*This document defines the terms used in this project (code, data, etc.).
It is important because different words have different meanings to different people in different
situations. Here we try to use each word in only one place so that meanings are unambiguous.*

.. glossary::

   Run
       A collection of Simulations that share the same Scenario and Mechanism(s).

       *Because of this definition, it is important to say that we 'execute code' or 'execute a
       program', not 'run code' or 'run a program'.*

   Experiment
       A collection of Runs with a range of Scenarios and/or Mechanisms
       that are compared and contrasted to elucidate cause and effect relationships.

   Scenario
       A set of numeric parameters, to be used in a run. These are usually derived from physical,
       measurable characteristics.

       E.g.,

       * The number of tPA molecules in the container
       * The number of rows in the fibrin grid

   Mechanism
       A set of rules governing how the model will execute. These are usually represented by
       different code logic.

       E.g.,

       * The presence of red blood cells in the fibrin grid
       * The 'into-and-along' unbinding hypothesis

   Simulation
       A single execution of the fibrinolysis model.
       This can either just be the microscale model, or it can be both,
       but the macroscale model cannot be run without input from the microscale model first.
       A simulation must have an associated Scenario,
       and an associated Mechanism for each model being executed.

   Edge grid
       The rectilinear grid that dictates location in the Macroscale model.
       This grid is made up of nodes and edges between them.
       Fibers may exist along edges in the grid.
       Locations are based on the edges between adjacent nodes, however nodes and edges are numbered
       differently. There are different co-ordinate systems used for this grid in Fortran and Python.
       See the documentation accompanying the Python EdgeGrid class for more information.

   Edge grid row
       Nodes, edges, and locations in the edge grid are all arranged in rows, numbered from bottom to
       top.

   Edge grid column
       Nodes in the edge grid are arranged in columns, numbered from left to right.

   Edge grid rank
       Fibers and edges in a given edge grid row are numbered by rank from left to right.
       Note that a location's rank does NOT correspond with its adjoining node's column.

   Dataset
       A single table of data. One Run will produce many different datasets.
       (This usage of the term is for consistency with the H5Py library which being used to store data
       for this project.)

   Data Collection
       A group of datasets resulting from a single simulation. Currently, the three data collections
       are:

       :Microscale Out: The datasets output by the microscale model
       :Macroscale In: Those datasets, derived from the microscale output, that the macroscale model
                       needs as input.
       :Macroscale Out: The datasets output by the macroscale model.
