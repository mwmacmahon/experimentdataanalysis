TODO (no particular order): 

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





