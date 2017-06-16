# General imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_dataframe_XYZ_pivot_tables(df, data_column,
                                   x_values_column=None,
                                   y_values_column=None,
                                   fill_value=np.nan):
    """
    Returns X, Y, Z dataframes corresponding to a 2D pivot of the
    given dataframe along X, Y. Expects a 2+ level multi-index and/or
    for x_values_column and y_values_column parameters to be given.
    'X' corresponds to pivot table column, and 'Y' to pivot table
    index, so that returned dataframes' .values attributes are
    plottable 2D matrices.

    If dataframe is 2+ level indexed _and_ one or both columns given,
    the indices will be used to determine row/col of data in table,
    but the values in X/Y matrix will utilize the [x/y]_values_column.
    Otherwise, X/Y values will be sorted numerically, although extra
    multi-index levels above two will be retained (in terms of .values,
    equivalent to stacking all 2D sub-dataframes along the y-axis).
    """
    # flatten as pivot_table() handles errors much better than unstack()
    if df.index.nlevels == 1:
        if (x_values_column is None) or (y_values_column is None):
            raise ValueError("plot_dataframe_2d requires either a 2+ level " +
                             "multi-indexed dataframe or both x_- and y_- " +
                             "values_column parameters to be given.")
        excess_index_columns = []
        x_index_column = x_values_column + '_copy'
        y_index_column = y_values_column + '_copy'
        flatdf = df.reset_index()  # just in case index important
        flatdf[x_index_column] = flatdf[x_values_column]
        flatdf[y_index_column] = flatdf[y_values_column]
    else:
        excess_index_columns = list(df.index.names[:-2])
        y_index_column = df.index.names[-2]
        x_index_column = df.index.names[-1]
        flatdf = df.reset_index()
        if x_values_column is None:
            x_values_column = x_index_column + '_copy'
            flatdf[x_values_column] = flatdf[x_index_column]
        if y_values_column is None:
            y_values_column = y_index_column + '_copy'
            flatdf[y_values_column] = flatdf[y_index_column]
    # add back any extra columns above 2d to keep y-axis ordering.
    pivot_index_columns = excess_index_columns + [y_index_column]
    # just quick-check to make sure these calls don't fail:
    flatdf[data_column]
    flatdf[x_values_column]
    flatdf[y_values_column]
    flatdf[pivot_index_columns]
    flatdf[x_index_column]
    pivot_df = flatdf.pivot_table(values=[data_column, x_values_column, y_values_column],
                                  index=pivot_index_columns,
                                  columns=[x_index_column],
                                  aggfunc=lambda x: x.head(1),
                                  fill_value=fill_value)
    xvals_df = pivot_df[x_values_column]  # pd.DataFrames in meshgrid()-style 
    yvals_df = pivot_df[y_values_column]
    zvals_df = pivot_df[data_column]
    return xvals_df, yvals_df, zvals_df

def get_dataframe_2d_matrix_and_axes_vecs(df, data_column,
                                          x_values_column=None,
                                          y_values_column=None,
                                          fill_value=np.nan):
    X, Y, Z = get_dataframe_XYZ_pivot_tables(df, data_column,
                                             x_values_column,
                                             y_values_column,
                                             fill_value)
    x_s = X.mean()  # pd.Series w/ALL values spread across rows
    y_s = Y.T.mean()
    return x_s.values, y_s.values, Z.values

# helper fcn for labelling axes with nonconsecutive values:
def get_inflection_points(values):
    trend_sign = np.sign(values[1] - values[0])
    vals_iterator = enumerate(values)
    last_ind, last_val = next(vals_iterator)  # pop off and add first (ind, val) pair
    inflection_point_indices = [last_ind]
    inflection_point_values = [last_val]
    for ind, val in yvals_iterator:
        if np.sign(val - last_val) == -1 * trend_sign:
            trend_sign = -1 * trend_sign
            inflection_point_indices.append(last_ind)
            inflection_point_values.append(last_val)
        last_ind, last_val = ind, val
    inflection_point_indices.append(ind)  # add last (ind, val) pair, too
    inflection_point_values.append(val)
    return inflection_point_indices, inflection_point_values

def plot_dataframe_waterfall(df, data_column,
                             num_waterfall_plots=None,
                             x_values_column=None, y_values_column=None,
                             fill_value=np.nan,
                             xlabel=None, ylabel=None,
                             ax=None):
    xvec, yvec, Zmat = \
        get_dataframe_2d_matrix_and_axes_vecs(df, data_column,
                                              x_values_column,
                                              y_values_column,
                                              fill_value)
    if x_values_column is None:  # disregard original indices
        xvec = np.arange(len(xvec))
    if y_values_column is None:
        yvec = np.arange(len(yvec))
    if num_waterfall_plots is None:  # default 5 plots
        num_waterfall_plots = 5
    num_waterfall_plots = min(num_waterfall_plots, len(yvec))
    waterfall_indices = np.linspace(0, len(yvec) - 1,
                                    num_waterfall_plots)
    waterfall_indices = np.int64(np.trunc(waterfall_indices))
    z_offset = 0
    for y_ind in waterfall_indices:
        zvals = Zmat[y_ind, :]
        valid_indices = ~np.isnan(zvals)
        zvals = zvals[valid_indices]
        zvals -= zvals.max()
        ax.plot(xvec[valid_indices], zvals + z_offset,
                'd-', label=str(yvec[y_ind]))
        z_offset += zvals.min()
    ax.legend()
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

def plot_matrix_colorplot(matrix, xvec=None, yvec=None,
                          xlabel=None, ylabel=None,
                          ax=None, **imshow_kwargs):
    nx, ny = matrix.shape
    if xvec is None:
        xvec = np.arange(nx)
    if yvec is None:
        yvec = np.arange(ny)
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    if 'origin' not in imshow_kwargs.keys():
        imshow_kwargs['origin'] = 'upper'
    if 'extent' not in imshow_kwargs.keys():
        if imshow_kwargs['origin'] == 'upper':
            imshow_kwargs['extent'] = [min(xvec), max(xvec),
                                       max(yvec), min(yvec)]
        else:
            imshow_kwargs['extent'] = [min(xvec), max(xvec),
                                       min(yvec), max(yvec)]
    width = abs(imshow_kwargs['extent'][1] - imshow_kwargs['extent'][0])
    height = abs(imshow_kwargs['extent'][3] - imshow_kwargs['extent'][2])
    natural_aspect_ratio = width / height
    if 'aspect' in imshow_kwargs.keys():
        imshow_kwargs['aspect'] *= natural_aspect_ratio
    else:
        imshow_kwargs['aspect'] = 1.6 * natural_aspect_ratio
    if 'cmap' not in imshow_kwargs.keys():
        imshow_kwargs['cmap'] = 'jet'
    if 'interpolation' not in imshow_kwargs.keys():
        imshow_kwargs['interpolation'] = 'nearest'
    ax.imshow(matrix, **imshow_kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
#    plt.show()  # screws up other plots subsequent to this one

def plot_dataframe_colorplot(df, data_column,
                             x_values_column=None, y_values_column=None,
                             fill_value=np.nan,
                             xlabel=None, ylabel=None,
                             ax=None, **imshow_kwargs):
    xvec, yvec, Zmat = \
        get_dataframe_2d_matrix_and_axes_vecs(df, data_column,
                                              x_values_column,
                                              y_values_column,
                                              fill_value)
    if x_values_column is None:  # disregard original indices
        xvec = np.arange(len(xvec))
    if y_values_column is None:
        yvec = np.arange(len(yvec))
    plot_matrix_colorplot(Zmat, xvec, yvec,
                          xlabel=xlabel, ylabel=ylabel,
                          ax=ax, **imshow_kwargs)

#     def get_axis_ticks(axis_values, use_inflection_points=False):
#         if use_inflection_points:
#             label_indices, label_values = get_inflection_points(axis_data_series)
#             if len(label_indices) > 2:  # keep only if not just linear
#                 return label_indices, label_values
#         else:
#             if len(axis_data_series) < 6:
#                 return axis_data_series.index, axis_data_series.values
#         iloc_indices_to_use = np.linspace(0, len(axis_data_series) - 1,
#                                           num_evenly_spaced_indices)
#         iloc_indices_to_use = np.int64(np.trunc(iloc_indices_to_use))
#         label_indices = axis_data_series.index.values[iloc_indices_to_use]
#         label_values = axis_data_series.loc[label_indices]
#         return label_indices, label_values

#     plot_label_column_names = ['wavelength', 'pump_power']
#     plotlabel = "\n".join(['{}: {}'.format(colname, dataframe.loc[run_id][colname].iloc[0])
#                            for colname in plot_label_column_names])
#     ax.text(1.1, 0.9, plotlabel, verticalalignment='top', horizontalalignment='left',
#             transform=ax.transAxes, color='black', fontsize=16)

#     x_tick_label_indices, x_tick_label_values = get_axis_ticks(x_s,
#                                                                use_inflection_points=False)
#     x_tick_labels = ["{}".format(val)
#                      for val in x_tick_label_values]
#     y_tick_label_indices, y_tick_label_values = get_axis_ticks(y_s,
#                                                                use_inflection_points=False)
#     y_tick_labels = ["{}".format(val)
#                      for val in y_tick_label_values]
#     if x_tick_label_indices is not None:
#         plt.xticks(x_tick_label_indices, x_tick_labels)
#     if y_tick_label_indices is not None:
#         plt.yticks(y_tick_label_indices, y_tick_labels)

# def plot_2d_with_run_id_slider(dataframe, data_column,
#                                x_values_column=None, y_values_column=None):
#     """
#     Expects a 3-level-multi-indexed dataframe and a column name corresponding
#     to the value to be plotted. Optionally, columns containing values corresponding
#     to each axis can be provided, otherwise the labels in the dataframe's index
#     will be used.
#     """
#     def plot_2d_by_3rd_index(dataframe, third_index):
#         plot_dataframe_2d(dataframe.xs(third_index, level=-3), 
#                           data_column, x_values_column, y_values_column)

#     run_id_slider = widgets.IntSlider(min=0, max=dataframe.index.get_level_values('run_id').max(),
#                                       value=0, description='Run ID:')
#     widgets.interact(plot_2d_by_3rd_index, dataframe=widgets.fixed(dataframe), run_id=run_id_slider);
