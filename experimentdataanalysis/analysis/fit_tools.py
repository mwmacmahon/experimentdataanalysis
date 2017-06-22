from inspect import signature

import lmfit

def create_residuals_fcn(fit_function, independent_var_names):
    """Given a fit function that accepts parameters
        (x_1, x_2, ..., x_m, param_1, param_2, ..., param_n)
    and the list of independent variable names
        ['x_1', 'x_2', ..., 'x_m'],
    creates a function that accepts parameters
        (lmfit.Parameters(), x_1, x_2, ..., measured_data=None)
    where Parameters() obj contains ['param_1', 'param_2', ...]
    If measured_data is None, the function returns the fit
    function's result with all independent variables and
    parameters plugged in. If measured_data is not None, it
    returns the difference between this result and the
    provided measured_data.
    """
    n_indvars = len(independent_var_names)
    sig = signature(fit_function)
    fcn_params_set = set(sig.parameters).difference(independent_var_names)
    def residuals_fcn(params, *args, measured_data=None, **kwargs):
        param_values = params.valuesdict()
        matched_fcn_param_values = {key: param_values[key]
                                    for key in fcn_params_set}
        arglist = list(args)  # going to parse through this
        args_index = 0
        indep_vars_list = []
        posargs_okay = True
        for ind, name in enumerate(independent_var_names):
            if name in kwargs.keys():  # try filling w/ kwarg first
#                 print('found ' + name + ' in kwargs')
#                 print('   args remaining: {}'.format(arglist[args_index:]))
                indep_vars_list.append(kwargs[name])
                posargs_okay = False  # by python rules, only kwargs now
            elif posargs_okay and len(arglist) > args_index:
#                 print('found ' + name + ' in args')
#                 print('   args remaining: {}'.format(arglist[args_index + 1:]))
                indep_vars_list.append(arglist[args_index])
                args_index += 1
            else:  # no more pos. args, no kwarg match found
                raise ValueError("residuals_fcn: failed to find all " +
                                 "{} independent ".format(n_indvars) +
                                 "variables, please make sure names and " +
                                 "ordering are correct.")
        if measured_data is None:
            if len(arglist) > args_index:
#                 print('found measured_data in args')
                measured_data = arglist[args_index]
                args_index += 1
        if len(arglist) > args_index:  # done parsing, is list exhausted?
            raise ValueError("residuals_fcn: too many inputs, only "
                             "expecting params, {} ".format(n_indvars) +
                             "independent variables, and " +
                             "optionally a 'measured_data'.")
        fit_function_results = fit_function(*indep_vars_list,
                                            **matched_fcn_param_values)
        if measured_data is None:
            return fit_function_results
        else:
            return measured_data - fit_function_results
    return residuals_fcn