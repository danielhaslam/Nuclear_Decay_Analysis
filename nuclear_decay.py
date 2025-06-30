# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 22:40:52 2022

@author: Daniel Haslam
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin


def create_data_array():
    '''
    Substitutes the inputted filenames and their corresponding delimiters to
    yield the containing data and produce a data set containing all the files'
    information.

    Parameters
    ----------
    None.

    Returns
    -------
    data_array_made : array (2D)

    '''

    data_array_made =  np.zeros((0,3))
    filenames = ["Nuclear_data_1.csv", "Nuclear_data_2.csv"]

    try:
        for index, filename in enumerate(filenames):
            data_array_made = np.vstack((data_array_made, read_in_data(
                filename, ",")))

    except OSError:
        print(f"The data file {filename} wasn't found. Please ensure your\
data file is in the same folder as this code.")

    data_array_made = data_array_made[data_array_made[:,0].argsort()]

    return data_array_made


def read_in_data(file_name, delimiter):
    '''
    Creates a 2-dimensional array containing the data in a data file.

    Parameters
    ----------
    file_name : string
    delimiter : string

    Returns
    -------
    file_data_array : array

    '''

    skipped_first_line_check = False
    is_nan_checker = False
    line_counter = 0
    erroneous_lines = []
    file_data_array = np.zeros((0,3))

    with open(file_name, encoding="utf-8") as input_file:
        for line in input_file:
            line_counter += 1

            if not skipped_first_line_check:
                skipped_first_line_check = True
            else:

                data_row = line.split(delimiter)

                if validate_data(data_row):
                    temporary_array = np.array([float(data_row[0]), float(data_row[1]),
                                     float(data_row[2])])

                    for value in range(0,2):
                        if np.isnan(temporary_array[value]):

                            is_nan_checker = True

                    if is_nan_checker:
                        is_nan_checker = False
                    else:
                        file_data_array = np.vstack((file_data_array, temporary_array))

                else:
                    line_counter = int(line_counter)
                    erroneous_lines = np.append(erroneous_lines, line_counter)

    erroneous_lines_display = ""
    for index in range(0, len(erroneous_lines-1)):
        erroneous_lines_display = erroneous_lines_display\
            + f"{erroneous_lines[index]}, "
    print("Lines ", erroneous_lines_display,
          f"couldn't be read in from {file_name} due to value error.")

    return file_data_array

def validate_data(data_row_given):
    '''
    Parameters
    ----------
    data_row_given : array (1D)

    Returns
    -------
    bool
        Determines whether the data row is eligible to be stacked onto the data
        array.

    '''

    for _, datum in enumerate(data_row_given):
        try:
            numerical_value = float(datum)
            if not np.isnan(numerical_value):
                if float(data_row_given[2]) > 0:
                    return True
                return False
            return False
        except ValueError:
            return False
        

def predicted_function_generator(sr_lambda_given, rb_lambda_given, time_value):
    '''
    Uses the exponential model of radioactive decay to produce and return a
    function representing the decay of Sr-79, using the decay constants
    substituted into it.

    Parameters
    ----------
    Sr_Lambda_given : float
    Rb_Lambda_given : float
    time_value : float

    Returns
    -------
    predicted_function : float

    '''

    initial_number_of_nuclei = (10**-6) * 6.022*10**23

    if(sr_lambda_given != 0 or rb_lambda_given != 0):
        predicted_function = 10**-12 * initial_number_of_nuclei * ((np.dot(sr_lambda_given,
        rb_lambda_given))/(rb_lambda_given - sr_lambda_given)) * (np.exp(
            -sr_lambda_given*time_value*3600) - np.exp(-rb_lambda_given*time_value*3600))

    else:
        return 0

    return predicted_function


def chi_squared_calculator(lambda_values):
    '''
    Calculates and returns the chi-squared value of the fit of a generated
    function, against the data collected from the data files.

    Parameters
    ----------
    lambda_values : array (1D)

    Returns
    -------
    total_chi_squared : float

    '''

    total_chi_squared = 0
    sr_lambda_given = lambda_values[0]
    rb_lambda_given = lambda_values[1]
    for i in range(0, len(data_array[:,0])-1):

        chi_squared_numerator = data_array[i,1] -\
            predicted_function_generator(sr_lambda_given,
                                         rb_lambda_given, data_array[i,0])
        chi_squared_denominator = data_array[i,2]
        chi_squared = (chi_squared_numerator)**2/(chi_squared_denominator)**2

        total_chi_squared += chi_squared

    return total_chi_squared


def chi_squared_minimisation():
    '''
    Finds the decay constants that produces the function that fits the data
    from the data array best - as judged by the decay constants that minimise
    the chi-squared parameter.

    Returns
    -------
    lambda_values : array (1D)

    '''

    chi_squared_min = fmin(chi_squared_calculator,\
                           (0,0), full_output=True, disp=False)
    lambda_values = chi_squared_min[0]

    return lambda_values


def remove_outliers(lambda_values_given, data_array_given):
    '''
    Considers if each data point lies within 3 standard deviations of an
    initial functional fit of the data. If not, then these are separated from
    the array of data given, and instead attributed to a new array,
    outliers_array_output.

    Parameters
    ----------
    lambda_values_given : array (1D)
    data_array_given : array (2D)

    Returns
    -------
    filtered_data_array : array (2D)
    outliers_array : array (2D)

    '''

    filtered_data_array = np.zeros((0,3))
    outliers_array = np.zeros((0,3))
    y_value = data_array_given[:,1]
    standard_deviation = data_array_given[:,2]

    for index in range(0, len(data_array_given[:,0])-1):

        function_value = predicted_function_generator(lambda_values_given[0],
        lambda_values_given[1], data_array_given[index,0])

        if float(function_value - 3*standard_deviation[index]) <= float(y_value[index])\
            <= float(function_value + 3*standard_deviation[index]):
            filtered_data_array = np.vstack((filtered_data_array,
                                             data_array_given[index,:]))
        else:
            outliers_array = np.vstack((outliers_array, data_array_given[index,:]))

    return filtered_data_array, outliers_array


def mesh_array_generator(final_lambda_values_given):
    '''
    Produces 2-dimensional "mesh" arrays for each decay constant,
    outlining a range of values within 5% of the final decay constant found.

    Parameters
    ----------
    final_lambda_values_given : array (1D)

    Returns
    -------
    Sr_lambda_mesh : array (2D)
    Rb_lambda_mesh : array (2D)

    '''

    sr_lambda_values = np.linspace(final_lambda_values_given[0] -
                                   0.05*final_lambda_values_given[0],
                                   final_lambda_values_given[0] +
                                   0.05*final_lambda_values_given[0], 100)

    rb_lambda_values = np.linspace(final_lambda_values_given[1] -
                                   0.05*final_lambda_values_given[1],
                                   final_lambda_values_given[1] +
                                   0.05*final_lambda_values_given[1], 100)

    sr_lambda_mesh, rb_lambda_mesh = np.meshgrid(sr_lambda_values,
                                                 rb_lambda_values)

    return sr_lambda_mesh, rb_lambda_mesh


def find_lambda_uncertainties(final_lambda_values_given,
                              min_chi_squared_given):
    '''
    Uses the lambda meshes generated in the function to find all combinations
    of lambda values that produce a chi-squared value of 1 more than the found
    minimum. It then identifies and returns the maximum of each lambda among
    all valid combinations, minus the lambda value found. This is the uncertainty
    on the value.

    Parameters
    ----------
    final_lambda_values_given : array (1D)
    min_chi_squared_given : float

    Returns
    -------
    sr_lambda_uncertainty : float
    rb_lambda_uncertainty : float

    '''

    sr_lambda_array = []
    rb_lambda_array = []

    sr_mesh, rb_mesh = mesh_array_generator(final_lambda_values_given)

    for x_coord in range(0, len(sr_mesh[:,0]-1)):
        for y_coord in range(0, len(sr_mesh[0,:]-1)):

            temp_lambda_array = np.array([sr_mesh[x_coord,y_coord],
                                          rb_mesh[x_coord,y_coord]])

            if np.abs(chi_squared_calculator(temp_lambda_array) -
                      (min_chi_squared_given+1)) < 0.05:
                sr_lambda_array = np.append(sr_lambda_array,
                                            sr_mesh[x_coord,y_coord])
                rb_lambda_array = np.append(rb_lambda_array,
                                            rb_mesh[x_coord,y_coord])

    sr_lambda_uncertainty = sr_lambda_array.max() - final_lambda_values_given[0]
    rb_lambda_uncertainty = rb_lambda_array.max() - final_lambda_values_given[1]

    return sr_lambda_uncertainty, rb_lambda_uncertainty


def plot_data(data_array_given, outliers_array_given, final_lambda_values_given):
    '''
    Plots the data on a graph of activity against time.

    Parameters
    ----------
    data_array_given : array (2D)
    outliers_array_given : array (2D)
    final_lambda_values_given : array (1D)

    Returns
    -------
    None.

    '''

    data_graph = plt.figure()
    data_plot = data_graph.add_subplot(111)

    data_plot.set_xlabel("Time (hours)")
    data_plot.set_ylabel("Activity (TBq)")
    data_plot.set_title("Activity of Decay of " + "$^{79}Rb$" + " With Time")

    data_plot.errorbar(data_array_given[:,0], data_array_given[:,1],
                        data_array_given[:,2], fmt='go', label="Valid data")

    data_plot.errorbar(outliers_array_given[:,0], outliers_array_given[:,1],
                       outliers_array_given[:,2], fmt='ko', label="Outliers")

    data_plot.errorbar(data_array_given[:,0],predicted_function_generator
                       (final_lambda_values_given[0], final_lambda_values_given[1],
                         data_array_given[:,0]), fmt='r-',
                       label="Predicted function")

    data_plot.legend()
    data_graph.savefig("nuclear_decay_plot.png", bbox_inches="tight")


def contour_plot_function(final_lambda_values_given,
                          min_chi_squared_given):
    '''
    Uses combinations of entries in the mesh arrays generated to substitute
    into the chi-squared function and produce combinations of the values of
    decay constants that generated a chi-squared one above the minimum
    chi-squared found. This gives the uncertainties on each parameter.

    Parameters
    ----------
    final_lambda_values_given : array (1D)
    chi_squared_value_given : float

    Returns
    -------
    None.

    '''

    lambda_meshes = mesh_array_generator(final_lambda_values_given)
    sr_lambda_mesh, rb_lambda_mesh = lambda_meshes

    chi_squared_mesh = np.apply_along_axis(chi_squared_calculator,
                                           0, lambda_meshes)

    contour_graph = plt.figure()
    contour_plot = contour_graph.add_subplot(111)

    contour_plot.set_xlabel("$λ_{Sr}$ "+ "$(s^{-1})$")
    contour_plot.set_ylabel("$λ_{Rb}$ "+ "$(s^{-1})$")
    contour_plot.set_title("A series of " + "$λ_{Sr}$" + " and " + "$λ_{Rb}$"
                           + " values producing " + "$χ_{min}^{_{2}}$" + " + 1")

    contour_plot.scatter(final_lambda_values_given[0],
                         final_lambda_values_given[1])

    chi_squared_contour = contour_plot.contour(sr_lambda_mesh, rb_lambda_mesh,
                                               chi_squared_mesh,
                                               levels=[min_chi_squared_given+1])

    sr_uncertainty, rb_uncertainty = find_lambda_uncertainties(
        final_lambda_values_given, min_chi_squared_given)

    sr_lambda_max = final_lambda_values_given[0] + sr_uncertainty
    rb_lambda_max = final_lambda_values_given[1] + rb_uncertainty

    contour_plot.axvline(sr_lambda_max, color='g', label = "Maximum "
                         + "$λ_{Sr}$" + " value")
    contour_plot.axhline(rb_lambda_max, color='c', label = "Maximum "
                         + "$λ_{Rb}$" + " value")

    contour_plot.clabel(chi_squared_contour, fontsize = 12, colors = 'r')
    contour_plot.legend()
    contour_graph.savefig("chi_squared_contour_plot.png", bbox_inches="tight")

def find_activity_uncertainty(final_lambda_values_given, time_value_given,
                              lambda_uncertainties_given):
    '''
    Uses the derivatives with respect to the decay constants of the predicted
    activity function format in order to calculate the total uncertainty of the
    activity with the uncertainties of the decay constants.

    Parameters
    ----------
    final_lambda_values_given : array (1D)
    time_value_given : float
    lambda_uncertainties_given : array (1D)

    Returns
    -------
    activity_uncertainty : float

    '''

    sr_lambda = final_lambda_values_given[0]
    rb_lambda = final_lambda_values_given[1]
    sr_uncertainty = lambda_uncertainties_given[0]
    rb_uncertainty = lambda_uncertainties_given[1]
    t_seconds = time_value_given*3600

    sr_lambda_derivative_term_1 = (rb_lambda*(np.exp(-sr_lambda*t_seconds)-
                                            np.exp(-rb_lambda*t_seconds)))/((
                                                rb_lambda-sr_lambda)**2)

    sr_lambda_derivative_term_2 = (np.exp(-sr_lambda*t_seconds)*sr_lambda*t_seconds)/(
        rb_lambda - sr_lambda)

    sr_derivative = 10**-12 * (10**-6) * 6.022*10**23 * rb_lambda*(sr_lambda_derivative_term_1 -
                                        sr_lambda_derivative_term_2)

    rb_lambda_derivative_term_1 = (sr_lambda*(np.exp(-sr_lambda*t_seconds)-
                                            np.exp(-rb_lambda*t_seconds)))/((
                                                rb_lambda-sr_lambda)**2)

    rb_lambda_derivative_term_2 = (np.exp(-rb_lambda*t_seconds)*rb_lambda*t_seconds)/(
        rb_lambda - sr_lambda)

    rb_derivative = 10**-12 * (10**-6) * 6.022*10**23 * sr_lambda*(-rb_lambda_derivative_term_1 +
                                        rb_lambda_derivative_term_2)

    activity_uncertainty = ((sr_derivative*sr_uncertainty)**2 +
                            (rb_derivative*rb_uncertainty)**2)**(1/2)

    return activity_uncertainty


def activity_request(final_lambda_values_given, lambda_uncertainties_given):
    '''
    Requests whether the user would like to see an activity value at a specific
    time. The user then inputs the time and the activity, along with its
    uncertainty, is generated.

    Parameters
    ----------
    final_lambda_values_given : array (1D)
    lambda_uncertainties_given : array (1D)

    Returns
    -------
    None.

    '''

    answer = ""
    answer = input("Would you like to see the activity value at a certain time? (y/n): ")
    while answer not in ("y", "n"):
        print("Invalid answer. Try again. ")
        activity_request(final_lambda_values_given, lambda_uncertainties_given)

    if answer == "y":

        try:
            time_input = float(input("Give the time value, in minutes, for which you would\
 like to see the activity: "))

            sr_lambda_given, rb_lambda_given = final_lambda_values_given

            activity_value = predicted_function_generator(sr_lambda_given,
                                                          rb_lambda_given, time_input/60)
            activity_uncertainty_value = find_activity_uncertainty(
                final_lambda_values_given, time_input/60, lambda_uncertainties_given)
        except ValueError:
            print("That's not a value. Try again.")
            activity_request(final_lambda_values_given, lambda_uncertainties_given)

    if answer == "n":
        sys.exit()

    print("The activity at {0} minutes is {1:.3g} ± {2:.1g} TBq.".format(time_input,
        activity_value, activity_uncertainty_value))
    

def final_data_presentation(data_array_given, outliers_array_given):
    '''
    Outputs the final key data produced using this program.

    Parameters
    ----------
    final_lambda_values_given : array (1D)
    outliers_array_given : array (2D)

    Returns
    -------
    None.

    '''

    elements = ["Strontium-79","Rubidium-79"]
    final_lambda_values = chi_squared_minimisation()
    min_chi_squared = chi_squared_calculator(final_lambda_values)
    lambda_uncertainties = find_lambda_uncertainties(
        final_lambda_values, min_chi_squared)

    half_lives = np.log(2)/final_lambda_values
    half_life_uncertainties = half_lives - np.log(2)/(final_lambda_values
                                                      + lambda_uncertainties)

    activity_at_90_minutes = predicted_function_generator(
        final_lambda_values[0], final_lambda_values[1], 1.5)
    activity_uncertainty_at_90_minutes = find_activity_uncertainty(
        final_lambda_values, 1.5, lambda_uncertainties)

    plot_data(data_array_given, outliers_array_given, final_lambda_values)
    contour_plot_function(final_lambda_values, min_chi_squared)

    for element in range(0,2):
        print("{0} has decay constant {1:.3g} ± {2:.3g} per second."\
              .format(elements[element], final_lambda_values[element],
                      lambda_uncertainties[element]))
        print("Thus, the half-life for {0} is {1:.3g} ± {2:.1g} minutes."
              .format(elements[element], half_lives[element]/60,
                      half_life_uncertainties[element]/60))

    print("The reduced χ\u00b2 value is {:.3g}"\
          .format(min_chi_squared/(len(data_array_given[:,0])-2)))

    print("The predicted activity at 90 minutes is {0:.3g} ± {1:.1g} TBq.".format(
        activity_at_90_minutes, activity_uncertainty_at_90_minutes))

    activity_request(final_lambda_values, lambda_uncertainties)


# --- Run Program ---


data_array = create_data_array()
outliers_array = np.zeros((0, 3))

corrected_lambda_values = chi_squared_minimisation()

data_array, outliers_array = remove_outliers(corrected_lambda_values, data_array)

final_data_presentation(data_array, outliers_array)

