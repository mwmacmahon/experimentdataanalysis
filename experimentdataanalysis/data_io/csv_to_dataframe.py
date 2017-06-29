import os

import numpy as np
import pandas as pd

import experimentdataanalysis.data_io.metadataparsing as metadataparsing

def parse_filepath_and_header(filepath, header_string=None,
                              parsing_keyword_lists=None,
                              existing_metadata_dict=None):
    if existing_metadata_dict is None:
        existing_metadata_dict = {}
    metadata_dict = metadataparsing.analyze_scan_filepath(filepath,
                                                          existing_metadata_dict,
                                                          parsing_keyword_lists)
    if header_string is not None:
        metadataparsing.analyze_string_for_dict_pairs(header_string, metadata_dict)
    return metadata_dict

def process_directory_csvs_to_dataframes(parent_dir,
                                         filename_key='.dat',
                                         run_criteria='directory',
                                         pandas_read_csv_kwargs={},
                                         parsing_keyword_lists=None,
                                         metadata_processing_fcns=[],
                                         metadata_filter_fcns=[],
                                         metadata_tag_to_column_list=[],
                                         dataframe_processing_fcns=[]):
    """Searches given directory for data files and parses them into
    pandas dataframes using pandas' read_csv function. Performs the
    following processing functions for each csv file found:
    1. Creates a metadata dict using the data_io.metadataparsing
       package and optionally the user-provided parsing_keyword_lists.
    2. Determines a "Run ID" key in the metadata based on the
       user-provided 'run_criteria', by default based on file's directory.
    2. Optionally calls all functions in user-provided
       metadata_processing_fcns list with metadata dict as parameter.
    3. Optionally calls all functions in user-provided
       metadata_filter_fcns list with metadata dict as parameter,
       and aborts processing the file if any function
       returns False.
    4. Creates a pandas dataframe from the contents of the file,
       optionally passing the user-provided pandas_read_csv_kwargs.
    5. Optionally looks for all tuples of form
       ('metadata tag', 'column name') in user-provided
       metadata_tag_to_column_list, and for each creates a column
       with the given name, with each row populated by the the 
       value corresponding to the metadata dict's associated tag.
    6. Finally, optionally calls all functions in user-provided
       dataframe_processing_fcns list with dataframe as parameter.
       Any values returned by function are unused; editing of
       dataframe generally requires use of pandas functions
       with the 'inplace' keyword argument set to True.
    The resulting lists of filenames, dataframes, and
    metadata dicts are returned, omitting rows filtered out
    via metadata_filter_fcns from all three lists.
    """
    # LOAD AND NOMINALLY FILTER DATA:
    unfiltered_filepath_list = []
    for dirpath, dirnames, filenames in os.walk(parent_dir):
        for filename in filenames:
            if filename_key in filename:
                unfiltered_filepath_list.append(os.path.join(dirpath, filename))
    if run_criteria == 'same':
        pass
    elif run_criteria == 'directory':
        last_dir = ''
        run_counter = -1
    else:
        raise NotImplementedError("currently only supported run " +
                                  "criteria is grouping-by-directory")
    # extract and remove # header rows, but don't change original
    pandas_read_csv_kwargs = pandas_read_csv_kwargs.copy()
    num_headerlines = pandas_read_csv_kwargs['skiprows']
    pandas_read_csv_kwargs['skiprows'] = 0
    filtered_file_list_index = 0
    filtered_filepath_list = []
    file_metadata_list = []
    file_dataframes_list = []
    for filepath in unfiltered_filepath_list:
        with open(filepath) as file_lines_iterator:
            try:
                header_lines = [next(file_lines_iterator)
                                for line in range(num_headerlines)]
                file_dataframe = \
                    pd.read_csv(filepath_or_buffer=file_lines_iterator,
                                **pandas_read_csv_kwargs)
            except pd.errors.ParserError:
                print("Pandas read_csv parser error, skipping file...")
                print("Filepath: {}".format(filepath))
                continue
            except pd.errors.EmptyDataError:
                print("No data found, skipping file...")
                print("Filepath: {}".format(filepath))
                continue
            except StopIteration:
                print("Problem encountered in process_directory_csvs_to_dataframes():")
                print("Tried to skip {} header lines, ".format(num_headerlines) +
                      "but file was too short. Skipping file...")
                print("Filepath: {}".format(filepath))
                continue
            # TODO: catch other exceptions for read_csv errors?
        # parse filepath, file header, and determine Run ID
        file_metadata = parse_filepath_and_header(filepath,
                                                  ''.join(header_lines),
                                                  parsing_keyword_lists)
        if run_criteria == 'same':
            file_metadata['Run ID'] = 0
        elif run_criteria == 'directory':
            current_dir = filepath.split('\\')[-2]
            if current_dir != last_dir:
                last_dir = current_dir
                run_counter += 1
            file_metadata['Run ID'] = run_counter

        # metadata processing and filtering
        for metadata_processing_fcn in metadata_processing_fcns:
            metadata_processing_fcn(file_metadata)
        if not all([metadata_filter_fcn(file_metadata)
                    for metadata_filter_fcn in metadata_filter_fcns]):
            continue

        # dataframe processing
        file_dataframe['file_index'] = filtered_file_list_index
        for metadata_tag, column_name in metadata_tag_to_column_list:
            if metadata_tag in file_metadata.keys():
                file_dataframe[column_name] = file_metadata[metadata_tag]
        for dataframe_processing_fcn in dataframe_processing_fcns:
            dataframe_processing_fcn(file_dataframe)

        filtered_filepath_list.append(filepath)
        file_dataframes_list.append(file_dataframe)
        file_metadata_list.append(file_metadata)
        filtered_file_list_index += 1
    return filtered_filepath_list, file_dataframes_list, file_metadata_list


