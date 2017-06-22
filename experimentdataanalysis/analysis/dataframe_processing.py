import numpy as np
import pandas as pd
from lmfit import minimize

def df_extract_dataset_indexed_matrices(df, column_names):
    """
    Given a pandas dataframe and a list of column names, returns
    a list of 2D numpy arrays corresponding to the data in each
    named column. Shape of matrix is given by multi-index:
    (# outer indices, # rows / inner indices per outer index)
    where it is expected that each outer index labels the
    same number of rows. If the dataframe does not have a
    multi-index, the resulting matrix has one row.
    
    The number of multi-index levels should not exceed 2.
    This function is intended for use in plugging multiple
    datasets into an equation at once, e.g. for cross-linked
    fitting.

    This function is _exceptionally_ fast if the inner index
    values are the same for each outer index (allowing a 2D array
    to be generated as simply "df.unstack()[column_name].values")

    e.g.  input parameters:
          df =                     [data columns]
          index_2d  index_1d       name    value     ...
                 0         0        'A'        1     ...
                 0         1        'B'        2     ...
                 1         0        'C'        3     ...
                 1         1        'D'        3     ...
                 2         0        'E'        2     ...
                 2         1        'F'        4     ...
               ...       ...        ...      ...     ...
          
          column_names = ['name', 'value', ...]
          
          output: [np.array([['A', 'B'], ['C', 'D'], ['E', 'F']]),
                   np.array([[1, 2], [3, 3], [2, 4]]),
                   ...]
    """
    if df.index.nlevels > 2:
        raise ValueError("dataframe multiindex should not exceed 2 levels")
    elif df.index.nlevels == 2:
        dataset_indices = df.index.levels[-2]
        # construct by unstack() [requires repeating inner index]
        mats = [df.unstack()[colname].values.copy()
                for colname in column_names]
        # check shape of final matrix is correct, else reindex and repeat
        dataset_indexing = df.loc[dataset_indices[0]].index
        if mats[0].shape[-1] != len(dataset_indexing):
            reindex_1d = lambda x: x.set_index(dataset_indexing)
            safe_df = df.groupby(level=-2).apply(reindex_1d)
            mats = [safe_df.unstack()[colname].values.copy()
                    for colname in column_names]
        return mats
    else:  # just upconvert vectors to one-row 2D arrays
        return [df[colname].values[np.newaxis, :].copy()
                for colname in column_names]


def df_extract_vector_lists_by_dataset(df, column_names):
    """
    Given a pandas dataframe and a list of column names, returns
    a list containing lists of the form:
    [column_1_array, column_2_array, ...]
    where the list corresponds to an outer index value in the given
    dataframe, and each array corresponds to all values of a
    named column with that outer index value.

    Resulting lists are returned in the order corresponding to
    the 2nd level of multiindex, given by df.index.levels[-2]

    e.g.  input parameters:
          df =                     [data columns]
          index_2d  index_1d       name    value     ...
                 0         0        'A'        1     ...
                 0         1        'B'        2     ...
                 1         0        'C'        3     ...
                 1         1        'D'        3     ...
                 2         0        'E'        2     ...
                 2         1        'F'        4     ...
               ...       ...        ...      ...     ...
          
          column_names = ['name', 'value', ...]
          
          output: [[np.array(['A', 'B']), np.array([1, 2])],
                   [np.array(['C', 'D']), np.array([3, 3])],
                   [np.array(['E', 'F']), np.array([2, 4])],
                   ...]
    """
    if df.index.nlevels > 2:
        raise ValueError("dataframe multiindex should not exceed 2 levels")
    elif df.index.nlevels == 2:
        dataset_indices = df.index.levels[-2]
        dataset_vecs_list = []
        for dataset_index in dataset_indices:
            dataset_df = df.loc[dataset_index]
            vecs = [dataset_df[colname].values.copy()
                    for colname in column_names]
            dataset_vecs_list.append(vecs)
        return dataset_vecs_list
    else:  # return one-dataset list
        return [[df[colname].values.copy()
                 for colname in column_names]]


def df_transform_dataset_df_to_fit_row(df, group_fit_params_dict,
                                       fit_params_to_add,
                                       column_aggregation_dict={},
                                       keep_const_columns=True):
    """
    Function to be used with DataFrame.groupby(level=-2) on a
    2-level-indexed dataframe to consolidate each
    outer-index-grouped dataset into a single row describing
    a fit performed on that dataset. By default, all columns
    will be kept if and only if they contain a constant value
    or an aggregation function is provided to handle that column.
    This will remove non-const data columns from fit automatically,
    which is likely what is wanted. Non-const columns can be kept
    via providing an aggregation function, and more columns can
    easily be dropped via "del df[colname]" if desired.

    Note that currently the resulting 1-row df is squeeze()'d
    in order to properly recombine into one row of the resulting
    dataframe after a groupby().apply() operation. So this function
    actually returns a series indexed by ordered column names.

    e.g. input parameters:
         df =                      x    y  bfield  elapsed_time
            index_2d  index_1d
                   4         0   1.0  0.1    40.0           0.0
                   4         1   1.5  0.3    40.0          10.0
                   4         2   2.0  0.5    40.0          20.0

         group_fit_params_list   ={..., 4: Parameters({'slope': 0.2,
                                                       'offset': 0.1}), ...}
         fit_params_to_add       =['slope']
         column_aggregation_dict ={'elapsed_time': lambda x: x.max()-x.min()}

         output (before squeeze()):
                        slope  bfield  elapsed_time
            index_2d  
                   4      0.2    40.0          20.0
    """
    new_df = df.head(1)
    for colname in list(new_df):
        if colname in column_aggregation_dict.keys():
            aggfcn = column_aggregation_dict[colname]
            new_df[colname] = aggfcn(df[colname])
            continue
        elif keep_const_columns:  # if const, no change needed
            first_value = df[colname].iloc[0]
            if all(df[colname] == first_value):
                continue
        del new_df[colname]
    dataset_index = df.index.get_level_values(level=-2)[0]
    fit_params = group_fit_params_dict[dataset_index]
    for param_name in fit_params_to_add:
        if param_name in fit_params.keys():
            param = fit_params[param_name]
            new_df[param_name] = param.value
            if param.stderr is not None:
                if param.stderr != 0:
                    param_error_str = param_name + '_error'
                    new_df[param_error_str] = param.stderr
    # groupby().apply() + squeeze() = 1-row-df -> series -> row-in-new-df
    # (returning a 1-row-df to apply() -> get doubled index, obnoxiously)
    # might need to revisit this sometime
    return new_df.squeeze(axis=0)

def df_minimize_fcn_on_datasets(df, residuals_fcn, fit_params,
                                independent_vars_columns,
                                measured_data_column,
                                *res_args,
                                column_aggregation_dict={},  # KEYWORDS ONLY!
                                keep_const_columns=True,
                                **res_kwargs):
    """
    Residuals function expected to take parameters
    (params, xvector1, xvector2, ..., yvector, *res_args, **res_kwargs)
    By default, drops all non-const columns in each dataset
    and adds all fit params to dataframe.
    """
    all_cols = independent_vars_columns + [measured_data_column]
    dataset_vecs_list = df_extract_vector_lists_by_dataset(df, all_cols)
    dataset_results_list = []  # will be in order of multiindex' outer indexing
    dataset_fit_params_list = []
    for vecs in dataset_vecs_list:
        xvecs = vecs[:-1]
        yvec = vecs[-1]
        # TODO: add option to scalar-ize as many as all-but-one vectors if const?
        if yvec.size == 0:
            continue
        result = minimize(residuals_fcn, fit_params,
                          args=(*xvecs, yvec, *res_args),
                          kws=res_kwargs)
        dataset_results_list.append(result)
        dataset_fit_params_list.append(result.params)
    dataset_indices = df.index.levels[-2].values
    group_fit_params_dict = dict(zip(dataset_indices,
                                     dataset_fit_params_list))
    fit_params_to_add = list(result.params)
    dfgroups = df.groupby(level=-2)
    new_df = dfgroups.apply(df_transform_dataset_df_to_fit_row,
                            group_fit_params_dict,
                            fit_params_to_add,
                            column_aggregation_dict,
                            keep_const_columns)
    return dataset_results_list, new_df


def df_minimize_fcn_across_linked_datasets(df, residuals_fcn, fit_params,
                                           dataset_params_unpacking_fcn,
                                           independent_vars_columns,
                                           measured_data_column,
                                           *res_args,
                                           column_aggregation_dict={},  # KEYWORDS ONLY!
                                           keep_const_columns=True,
                                           **res_kwargs):
    """
    Residuals function expected to take parameters
    (params, xvector1, xvector2, ..., yvector, *res_args, **res_kwargs)
    By default, drops all non-const columns in each dataset
    and adds all fit params to dataframe.
    """
    all_cols = independent_vars_columns + [measured_data_column]
    mats_list = df_extract_dataset_indexed_matrices(df, all_cols)
    xmats = mats_list[:-1]
    ymat = mats_list[-1]
    # TODO: add option to scalar-ize as many as all-but-one vectors if const?
    result = minimize(residuals_fcn, fit_params,
                      args=(xmats, ymat, *res_args),
                      kws=res_kwargs)
    dataset_fit_params_list = []
    for dataset_index in df.index.levels[-2].values:
        dataset_params = dataset_params_unpacking_fcn(result.params, dataset_index)
        dataset_fit_params_list.append(dataset_params)
    dataset_indices = df.index.levels[-2].values
    group_fit_params_dict = dict(zip(dataset_indices,
                                     dataset_fit_params_list))
    fit_params_to_add = list(dataset_param_dict.keys())
    dfgroups = df.groupby(level=-2)
    new_df = dfgroups.apply(df_transform_dataset_df_to_fit_row,
                            group_fit_params_dict,
                            fit_params_to_add,
                            column_aggregation_dict,
                            keep_const_columns)
    return result, new_df
