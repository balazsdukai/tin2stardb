================================================
TIN import/export into a Star-schema in database
================================================

A command line application for reading a TIN from a file and converting the
TIN into a Star-schema. What is a *Star-schema*? Read in the references.

Installation and usage
----------------------

At the moment, use Pipenv and the provided Pipfile to install all required
packages in a virtual environment. Once install, activate the virtual
environment and use the software from there.

Navigate to directory with the source code and do:

.. code-block::

    pipenv --python 3.6
    pipenv install --dev


Then activate the virtual environment with ``pipenv shell``. The command for
running the program is ``tin``.


Features
--------

- Read/write to Wavefront OBJ
- Use PostgreSQL+PostGIS as database backend
- Subset a OBJ(s) with a bounding box and write them back to files

Usage
-----

General usage:

.. code-block::

    Usage: tin [OPTIONS] CONFIGURATION COMMAND [ARGS]...

      Read/Write a TIN to a database from files.

      CONFIGURATION is a yaml file that contains the database credentials and
      table configuration options.


    Options:
      -v, --verbose  Increase verbosity. You can increment the level by chaining
                     the argument, eg. -vvv
      -q, --quiet    Decrease verbosity.
      --help         Show this message and exit.

    Commands:
      delete  Delete a TIN database.
      export  Export a TIN from a Star TIN database model to a file.
      import  Import TIN from a file into a database, using a Star-structure.
      subset  Subset a TIN from a file and write it out to a file.

Importing a single, or a directory of OBJs:

.. code-block::

    Usage: tin import [OPTIONS] FILENAME

      Import TIN from a file into a database, using a Star-structure.

      FILENAME can be a path to single file, or a directory containing the files
      to be imported. The importer tries to guess the file format from its
      suffix.

      Supported file formats: Wavefront OBJ (.obj)

    Options:
      --bbox FLOAT...    Import only a subset of the TIN in the 2D bbox: (minx
                         miny maxx maxy).
      --bboxes FILENAME  A GeoJSON file containing polygons of bounding boxes.
                         Each Feature must have a "tile" property, which is used
                         for finding the TIN file with equivalent name. Only used
                         when a directory of TIN files are provided as input.
      --help             Show this message and exit.


The software requires a configuration file for instructions on how to access
the database. This is a YAML file, and has the structure below.

.. code-block::

    database:
        dbname: tin_test
        host: localhost
        port: 5432
        user: someuser
        password:

    tin_schema:
        schema: tin
        table: tin
        field:
            id: vid
            geom: geom
            star: star

    epsg: 7415


References
----------

 Ledoux, H. & Meijers, M. (2013). A star-based data structure to store
 efficiently 3D topography in a database, Geo-spatial Information Science, 16:4, 256-266, DOI: 10.1080/10095020.2013.866618

 Ledoux, H. (2015). Storing and analysing massive TINs in a DBMS with a
 star-based data structure (No. 2015.01). Delft, the Netherlands: 3D geoinformation, Delft University of Technology.

 Kumar, K., Ledoux, H., & Stoter, J. (2016). COMPARATIVE ANALYSIS OF DATA
 STRUCTURES FOR STORING MASSIVE TINS IN A DBMS. ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLI-B2, 123–130. https://doi.org/10.5194/isprs-archives-xli-b2-123-2016
