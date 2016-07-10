INSTALLATION AND DOCUMENTATION:
-Should consider running from a virtual environment created from the
    requirements.txt file. 
-To create an environment in anaconda from requirements.txt that satisfies
    said requirements, can type from command line in project directory:
    "conda create --name [your_environment_name] --file requirements.txt"
    and switch to it with "activate [your_environment_name]"
-Can run spyder within said environment as well, as spyder is part of
    said package list - "switch to environment on command line, then run
    spyder" may be simplest method.
-Command line "python setup.py develop" - create in current python environment
    a fake package that links to the current work folder. Uninstall with
    "python setup.py develop --uninstall"
-Command line "python setup.py test" - run tests of everything in the tests
    folder using pytest.
-Command line "python setup.py install" - actually create a real install
    with a snapshot of the current work folder copied into ./lib/sitepackages
    not recommended, confuses installed version with work-in-progress version
-Note install_requires line of the setuptools script is borked, preventing
    an easy auto-install of required packages. Needs fixing, sorry.
-Command line "make html" to rebuild documentation. When it exists.
-Version info stored in top level __init__.py file as __version__


VERSION CONTROL:
-Standalone scripts for personal use don't necessarily need to be version
    controlled at all, but any changes to the "experimentdataanalysis"
    modules (bugfixes, added models, etc.) should be version controlled
    both so we can benefit from each other's work and so we don't fragment
    our code base.
-Git can get complicated, so make sure you remember to checkout a new
    branch _before_ making code changes! Working off the main branch is
    not a good idea and you don't want to have to migrate your changes
    to a separate branch later on, though anything is possible in Git
    with enough internet research.
-Git repository structure to follow the pattern at the following URL
    (it's pretty standard):
    http://nvie.com/posts/a-successful-git-branching-model/
-"Main" branch is reserved for stable versions you would trust telling
    someone to use who has no idea how the code works. If this were
    experiment running code, this is what you run
-Custom work should go on a new "feature" branch, from which some or all
    changes can be merged into the working "develop" branch which 


TODO (no particular order): 

0. ***update databrowser to:*** (most important)
    A. not do any of this "primary field" bs
    B. disable and hide all fitting stuff for now
    C. make sure saving images works
    D. make sure saving scandata/scandatasets works
    E. implement viewing vs fit parameter

0.1. fix sort: should be able to handle numeric OR non-numeric arguments for both priorities

0.2.swap out directory choosing from simply looking for CSVs to a system that handles any input
    will have to swap out most dcparsing calls for new ones to a source-agnostic function
    which will go for JSON/pickle > csvs/tsvs, but can have specific reading type explicitly
    given as a parameter (e.g. read this weird file format)

0.3. switch from pickle to JSON for saving files

0.4. better documentation/naming schemes for sorting/coordinates in scandatasets/scandatasetanalyzer/scandatamodels



1.
add unsorted & unfiltered checkboxes under plot (NOT related to fit; those can be filtered too except for fit param graph).
for 2D, comes back to how to display multiple types: I think standard is to only show those conforming to selected dataseries
xval-by-xval? in which case sorting order and filters of selected dataseries are used.

2.
add proper axes labels to databrowser

3.
fix data label parsing for 1D scan files (done?)

4.
get rid of dataclassfitting, dataclassgraphing, curvefitting
replace with
dataseriesprocessing, scandataprocessing, ____graphing

5. decide how to integrate ScanDataSet analysis with databrowser - if it is to
be combined at all. Using external code is in a sense its own paradigm, and
things like filters don't carry over well to a GUI interface. So how should
their relationship work? Ideas:
    5a) databrowser simply loads a series of scandata, and can perform simple
        fits on them using user-tweaked parameters. This is a separate
        paradigm from the scandatamodel system, to be used to either fix
        bad scans or to help work on the model in the first place. In other
        words, minimalist databrowser, zero integration with ScanDataSets.
*** vvv (probably this one) vvv
    5b) like 5a, but add the ability to import a scandataset by some
        method. then scandatamodel can be used to pre-load the fitting
        function, sort order, key index, etc. from the scandatamodel. fits can
        be adjusted and examined by hand, and then a new scandataset can be
        generated from the databrowser window in a similar way, e.g. replacing
        the imported set in place or added to workspace with a "_modified"
        tag at the end of the new scandataset's name.
***
    5c) use scandatasets directly, such the fit window is simply an interface
        to use ScanDataSet analysis functions. Even better would be to
        redisgn the UI slightly to be a ScanDataSet analyzer. Instead of a
        list view, we would want something like the view you get when creating
        a new folder in windows - a directory hierarchy of scandatasets and
        the scandata they contain.

5.1. Okay, if I go with 5b, then what to do about the fit window? Perhaps just
a table of fit parameters that can be edited with a live preview showing the
starting fit? Can start with just graying out the section and finishing the
plot display stuff - the live-fitting-in-the-fly is less important atm. anyway.

6. smarter file importation. perhaps use numpy's csv importation system?

7. perhaps get rid of the entire "scandata field ordering" system in favor of
    a name-based system. 








tests to make:
-operations on dataseries/scandata work correctly, in sequence. Such as:
--making sure operations on scandata dataseries work properly, such as adding
    excluded intervals and shifting times over. should be able to change
    field indices affected correctly, and do many in sequence and have
    correct result





