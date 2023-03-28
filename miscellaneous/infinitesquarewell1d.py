"""" This file will contain a function that will plot the first 4 eigenstates of the 1-D PIB"""
# import the necessary libraries

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy import constants


def infinitesqw_1d(lenght, princqn, xvals):
    """
    This function takes
    :param: "lenght" which is the size of the box
    :param: princqn which is the principal quantum number of interests
    :param: xvals which is the list of values from 0 to L used to generate the function
    returns two lists and a constant:
        the values of the function, the probability density and the associated
    energy
    """
    # defining the infinite square well wavefunction
    function = np.sqrt(2 / lenght) * np.sin(princqn * np.pi * xvals / lenght)

    # calculating the associated energy of the state specified by princqn
    energy = ((princqn**2) * (constants.pi**2) * (constants.hbar**2)) / (
        2 * constants.m_e * lenght**2
    )

    # defining the probability density of the wavefunction
    probability_density = function**2
    # we should take a dot product of the complex conj of the function with itself
    # since we are dealing with real values we take the square of the function
    return (function, probability_density, energy)


if __name__ == "__main__":
    # assigning variables of interests
    #calling the function
    #generating function to update graph
    SIZE = 4
    xvalues = np.linspace(0, SIZE, 1000)
    PRINCIPLEQNUMBER = 1
    # function call
    funcvals, probadensity, ENERGY = infinitesqw_1d(SIZE, PRINCIPLEQNUMBER, xvalues)

    zeros = np.zeros(len(xvalues))
    # we generate the zeros list to make the nodes visible on the graph

    # creating plot and sliders

    # plotting the eigenstates
    fig, (ax1, ax2, ax3) = plt.subplots(3,figsize=(6,6))
    plt.tight_layout()
    fig.suptitle("Infinite Square Well")
    ax1.plot(xvalues, funcvals, "b")
    ax2.plot(xvalues, probadensity, "r")
    ax3.plot(PRINCIPLEQNUMBER, ENERGY, ".")

    plt.subplots_adjust(bottom=0.2,top=0.95)  # generating space under graph to add slider

    ax1.set_ylabel("Wavefunction")
    ax1.set_xlabel("Lenght")
    ax2.set_ylabel("Prob density")
    ax3.set_ylabel("Energy [J]")

    # Adding slider functionality to plot
    # xposition, yposition, width and height
    ax1.slide = plt.axes([0.15, 0.1, 0.65, 0.05])

    # Properties of the slider that will change the principal quantum number
    principleQnumber_slider = Slider(
        ax1.slide, label="n", valmin=1, valmax=10, valinit=1, valstep=1
    )

    # Making the function to update the plot
    def update(val):
        """This function will update the graphs every time the slider is used
        It will make a function call to infinitesqw_1d
        It will clear the previous plots and make new plots on the same fig objects
        It won't clear the energy plot (optional to user) uncomment if you want
        """
        unused_value=val
        #just to avoid the W0613 from pylint until I find a better solution

        current_val = principleQnumber_slider.val
        functionvals, probdensity, energ = infinitesqw_1d(SIZE, current_val, xvalues)
        ax1.cla()
        ax1.plot(xvalues, functionvals, "b")
        ax1.plot(xvalues, zeros, "g--")
        ax2.cla()
        ax2.plot(xvalues, probdensity, "r")
        # ax3.cla()
        ax3.plot(current_val, energ, ".")
        ax2.set_xlim(0, SIZE)
        ax1.set_xlim(0, SIZE)

    principleQnumber_slider.on_changed(update)
    plt.show()
