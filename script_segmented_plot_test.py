# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 11:42:56 2017

@author: Michael
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

scandata = scandata_list[4]
#dataseries = scandata.dataseries_list[2]
#xvals, yvals = zip(*dataseries.datatuples())
xvals, yvals = scandata.xy
segment_xvals = []
segments = []


last_xval = -999999
count_at_xval = 0
current_segment_yvals = []
for xval, yval in zip(xvals, yvals):
    if abs(xval - last_xval) > 0.0001:  # new segment found
        if last_xval != -999999:  # break off old segment
            segment_xvals.append(last_xval)
            segments.append(np.array(current_segment_yvals))
        last_xval = xval
        current_segment_yvals = []
    current_segment_yvals.append(yval)

plt.figure()

plt.subplot(2,1,1)
plt.title('Blue to Red at each field')
for xval, yvals in zip(segment_xvals, segments):
    for y_ind, yval in enumerate(yvals):
        scale = 1.0 * y_ind / len(yvals)
        color = (1.0 * scale, 0, 1 - 1.0 * scale)
        plt.plot(xval, yval, 'o',
                 markerfacecolor=color, markeredgecolor=color)

plt.subplot(2,1,2)
plt.title('Chronological, field shifts highlighted')
current_start_index = 0
for segment_index, (xval, yvals) in enumerate(zip(segment_xvals, segments)):
    if segment_index % 2 == 0:
        segment_color = 'b'
    else:
        segment_color = 'r'
    indices = current_start_index + np.arange(len(yvals))
    current_start_index += len(yvals)
    plt.plot(indices, yvals, '-', color=segment_color)
#    plt.plot(xval * np.ones_like(yvals), yvals, 'd', color=segment_color)