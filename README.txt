INSTALLATION:
-Recommended python package manager is Anaconda for its simplicity.
    Python 3.5 is recommended; multiprocessing requires it!
-Install Git, one simple option is Git For Windows:
    https://git-scm.com/downloads
-After setting up Git, use "git clone" to pull a copy of directory onto
    your computer, automatically set up to track shared online version.

    In directory where you want to put project, open Git Bash and type:
        git clone https://github.com/vsihlab/experimentdataanalysis.git

-Should immediately switch to a new custom branch with your name and upload
    it. This keeps your changes local to your "branch" both on your computer
    and online, and lets the rest of us see any proposed changes to be merged
    in.

    First, switch to the active development branch:
        git checkout develop

    Second, create your new branch from command line (in project folder) with:
        git checkout -b "[yourname]_[branchname]" develop

    Third, set up your branch to allow you to host it online via:
        git push --set-upstream origin "[yourname]_[branchname]"

    e.g. "git checkout develop"
         "git checkout -b michael_incremental_updates develop"
         "git push --set-upstream origin michael_incremental_updates"

-Can swap between branches via "git checkout". "git fetch" and "git status"
    can be used liberally to keep track of changes both local and remote
    and stay up to date. "git push origin" allows you to push your branch
    online to be shared, and "git pull origin" updates your local copy to
    stay current with the online version as long as there are no conflicts.

-Note: if you have other git repositories on different accounts, you may
    need to set repository-specific credential storage if you use the
    Git Credentials Manager.

    This can be done with the Git Bash line:
        git config --global credential.useHttpPath 

-Please read the section on version control below!
    


VIRTUAL ENVIRONMENTS:
-If changing code, consider running from a virtual environment - an
    alternate, more carefully controlled sub-installation of Python - based
    on the requirements.txt file. This ensures compatibility of packages.
    Some, like scipy, update often enough it can be important to have
    a new enough version, and it ensures new code doesn't accidentally add
    new requirements.
-To create an environment in Anaconda from requirements.txt that satisfies
    said requirements, can type from command line in project directory:
    "conda create --name [your_environment_name] --file requirements.txt"
    and switch to it with "activate [your_environment_name]"
-Can run spyder within said environment as well, as spyder is part of
    said package list - "switch to environment on command line, then enter
    'spyder'" may be simplest method. Can create shortcuts as well, look at
    default Anaconda start menu shortcuts for examples.

PACKAGE INSTALLATION:
-Not necessary to "install" in order to run, but doing a develop install (see
    below) allows you to run your scripts from anywhere on the computer,
    as otherwise they have to be run from the project root directory in order
    for the script to see the project modules
-Command line "python setup.py develop" - create in current python environment
    a fake package that links to the current work folder. Uninstall with
    "python setup.py develop --uninstall"
-Command line "python setup.py test" - run tests of everything in the tests
    folder using pytest. Not for installation, but great for testing!
-Command line "python setup.py install" - actually create a real install
    with a snapshot of the current work folder copied into ./lib/sitepackages
    not recommended, confuses installed version with work-in-progress version
-Note "install_requires" line of the setuptools script is borked, preventing
    an easy auto-install of required packages. Needs fixing, sorry.

VERSION CONTROL:
-Standalone scripts for personal use don't necessarily NEED to be version
    controlled at all, but any changes to the "experimentdataanalysis"
    modules (bugfixes, added models, etc.) need to be version controlled
    both so we can benefit from each other's work and so we don't fragment
    our code base.
-Really, though, it can only help to keep everything under version control,
    especially for syncing between computers and sharing with others.
-Lots of git tutorials out there, I recommend learning the basics:
    "git clone", "git push/pull", "git checkout", "git commit -m 'message'",
    "git status", "git fetch/merge", "git add/rm"...it's a lot to take in,
    but extremely useful both here and in any future coding work.
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
-"Release" branches should use "python setup.py develop" to update version
    after changing top-level "__init__.py". Don't forget CHANGES.txt either.
-Etiquette to be determined, but probably don't merge your branch into
    develop/main without talking to others. Can push your branches online
    without merging and affecting others' code, and changes to shared
    develop/main branches affect other branches downstream!

DOCUMENTATION:
-Erm, none yet. Except this readme and function docstrings.
-Working on implementing auto-generated documentation from the "docstrings"
    at the beginning of each python function/class/etc, using Sphinx.
-Command line "make html" in project dir to rebuild documentation.
    Once it exists.
-Note that version number is stored in top level __init__.py file
    (experimentdataanalysis\__init__.py) as "__version__ = 'X.X.X'"
    Both Sphinx and setup.py can read this variable; only git does not
    automatically tag version numbers (or even push version tags automatically)

IMPLEMENTATION NOTES ON ITERATORS AND ITERABLES:
-Most functions are set to accept iterables instead of say, lists. This means
    you can send it a list, a tuple, or even an iterator.
-Returning iterators is good practice if you can lazily process it. However,
    this can get confusing for people unfamiliar with the practice, as code
    that seems self-explanatory can fail:
        a = fcn_that_returns_iterator_instead_of_list()
        #  a = list(a) <- would fix problem, puts a's values in permanent list
        for x in a:
            print(x)  # runs fine, prints everything in a
        for x in a:
            print(x)  # _does nothing_, a is exhausted by previous loop
-In python 3.x, important built-in functions return iterators, e.g. zip()
-However, I've been moving towards returning lists instead of iterators.
    If code uses "list(fcn_returning_iterator_or_list)" everywhere it doesn't
    matter anyway, and processes like filtering tend to ruin the speedup.
    For our purposes the speed gains are unlikely to be worth the confusion.
-EXCEPTION: Parsing modules should evaluate lazily, and return iterators
    instead of lists. This is because one may want to process large amounts
    of data (e.g. a hard drive's worth), and by putting all the data in a list
    we force Python to read and store all the data at once!
    Note databrowser.py and scandatasetprocessing.py keep all data loaded at
    once, but future programs and modules by no means need to.
-RECOMMENDATION: whenever accepting a list/iterator/whatever from a function,
    use something like "returned_list = list(returned_list)" or similar
    to ensure iterators, iterables, tuples, etc. all transformed into a list.
    This means if that function is changed and the output type is modified,
    your code doesn't break. HOWEVER, despite the loss of speed, I recommend
    just returning a list rather than using fancy "yield" syntax to return
    iterators. This can be up for debate in the future.


















TODO SCRATCH NOTES(no particular order or organization): 
    -add JSON scandata saving!

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


1.
add unsorted & unfiltered checkboxes under plot (NOT related to fit; those can
    be filtered too except for fit param graph). for 2D, comes back to how to
    display multiple types: I think standard is to only show those conforming
    to selected dataseries xval-by-xval? in which case sorting order and
    filters of selected dataseries are used.

2.
add proper axes labels to databrowser 2D view. without crashing.

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







tests to make:
-operations on dataseries/scandata work correctly, in sequence. Such as:
--making sure operations on scandata dataseries work properly, such as adding
    excluded intervals and shifting times over. should be able to change
    field indices affected correctly, and do many in sequence and have
    correct result. Empty scandata/dataseries should not crash anything.





