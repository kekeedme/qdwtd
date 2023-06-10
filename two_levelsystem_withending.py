"""This file contains a function to calculate the occupation probability of a two-level system"""
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def wavefunction(
    interaction, energydifference, dissipation, coefficient1, coefficient2
):
    """This function returns the occupation probability of states in a two-level system
    It takes as input:
    :param interaction:the interaction energy between states
    :param energydifference: the energy difference between states
    :param: coefficients 1 and 2: the coefficients specifying the initial population of each states
    :param: the energy dissipation magnitude
     it:
        -builds the wavefunction at time 0
        -generates the time-independent hamiltonian matrix
        -defines the basis states phi_1 and phi_2
        -generates wavefunction at all times->list
        -calculates the occupation probability of each states
        -returns the values as float for plotting
        -the population is phenomenologically modeled by adding an imaginary term
        to the off-diagonal term in the Hamiltonian -i*dissipation which is the dissipation"""

    # making the time data points
    time = np.linspace(0, 2 * np.pi, 1000)
    # defining the basis states as column vectors
    phi_1 = np.array([[1, 0]]).T
    phi_2 = np.array([[0, 1]]).T
    # Definining the wavefunction at t = 0
    psi_0 = (1 / np.sqrt(coefficient1**2 + coefficient2**2)) * np.array(
        [[coefficient1, coefficient2]]
    ).T
    # Defining the Hamiltonian matrix (making it non hermitian)
    hamiltonian = np.array(
        [[0, interaction], [interaction, energydifference - i * dissipation]]
    )
    # generating the wavefunction at all times (evolution operator is no longer unitary, so the norm decreases with time)
    psi_time = [expm(-t * hamiltonian * i).dot(psi_0) for t in time]
    # Defining probability for each state at each time point, and only taking the real part
    probability1 = np.real(
        [(phi_1.T.dot(Psi_t)).dot(Psi_t.conj().T.dot(phi_1)) for Psi_t in psi_time]
    )
    probability2 = np.real(
        [(phi_2.T.dot(Psi_t)).dot(Psi_t.conj().T.dot(phi_2)) for Psi_t in psi_time]
    )
    # The probability will be real but python returns complex where the Im part is zero
    # converting the list of array to a list of floats (needed for plotting)
    prob11 = [float(p1) for p1 in probability1]
    prob22 = [float(p2) for p2 in probability2]
    return prob11, prob22


if __name__ == "__main__":
    i = 1j  # complex i
    INTERACTION = 1  # interaction energy
    ENERGYDIFFERENCE = 1  # energy difference
    timepoints = np.linspace(0, 2 * np.pi, 1000)
    COEFF_1 = 0
    COEFF_2 = 1
    LOSS = 1
    prob1, prob2 = wavefunction(INTERACTION, ENERGYDIFFERENCE, LOSS, COEFF_1, COEFF_2)

    # creating figure object and subplot to add to the figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    # creating axes to add sliders
    COEFFICIENT1_slider_axis = fig.add_axes([0.15, 0.92, 0.65, 0.05])
    COEFFICIENT2_slider_axis = fig.add_axes([0.15, 0.87, 0.65, 0.05])
    INTERACTIONENERGY_slider_axis = fig.add_axes([0.15, 0.82, 0.65, 0.05])
    DELTA_slider_axis = fig.add_axes([0.15, 0.77, 0.65, 0.05])
    ENERGLOSS_slider_axis = fig.add_axes([0.15, 0.73, 0.65, 0.05])
    # creating the sliders
    COEFFICIENT1_slider = Slider(
        COEFFICIENT1_slider_axis,
        label="c1",
        valmin=-1,
        valinit=COEFF_1,
        valmax=1,
        valstep=0.1,
    )
    COEFFICIENT2_slider = Slider(
        COEFFICIENT2_slider_axis,
        label="c2",
        valmin=-1,
        valinit=COEFF_2,
        valmax=1,
        valstep=0.1,
    )
    DELTA_slider = Slider(
        DELTA_slider_axis, label="\u0394", valmin=-10, valmax=10, valinit=1, valstep=1
    )
    INTERACTIONENERGY_slider = Slider(
        INTERACTIONENERGY_slider_axis,
        label="V",
        valmin=-4,
        valmax=4,
        valinit=1,
        valstep=0.5,
    )
    ENERGLOSS_slider = Slider(
        ENERGLOSS_slider_axis,
        label="\u0393",
        valmin=0,
        valmax=6,
        valinit=1,
        valstep=0.5,
    )
    # making the plot
    ax.plot(timepoints, prob1,label=f'$P_1(t)$')
    ax.plot(timepoints, prob2,label=f'$P_2(t)$')
    ax.legend(bbox_to_anchor=(1.14, 0.6), loc='center right', frameon=False)
    ax.set_xlabel("time [a.u.]", fontsize=15)
    ax.set_ylabel("Probability", fontsize=15)
    plt.subplots_adjust(bottom=0.25, top=0.7)

    def update(val):
        """This function will update the graphs every time the slider is used
        It will make a function call to wavefunction
        It will clear the previous plots and make new plots on the same fig objects
        """
        unused_value = val  # avoiding pylint error of unused value
        currentcoefficient1 = COEFFICIENT1_slider.val
        currentcoefficient2 = COEFFICIENT2_slider.val
        currentinteraction = INTERACTIONENERGY_slider.val
        currentdifference = DELTA_slider.val
        currentloss = ENERGLOSS_slider.val
        proba1, proba2 = wavefunction(
            currentinteraction,
            currentdifference,
            currentloss,
            currentcoefficient1,
            currentcoefficient2,
        )
        ax.cla()  # clear the plot every time the function is called, before plotting again
        ax.plot(timepoints, proba1,label=f'$P_1(t)$')
        ax.plot(timepoints, proba2,label=f'$P_2(t)$')
        ax.legend(bbox_to_anchor=(1.14, 0.6), loc='center right', frameon=False)
        ax.set_xlabel("time [a.u.]", fontsize=15)
        ax.set_ylabel("Probability", fontsize=15)

    # When the slider values are changed, the update function is call
    COEFFICIENT1_slider.on_changed(update)
    COEFFICIENT2_slider.on_changed(update)
    INTERACTIONENERGY_slider.on_changed(update)
    DELTA_slider.on_changed(update)
    ENERGLOSS_slider.on_changed(update)

    plt.show()
