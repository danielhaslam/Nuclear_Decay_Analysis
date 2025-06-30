
# Nuclear Decay Analysis - ⁷⁹Rb and ⁷⁹Sr Isotopes

This project analyzes the decay of the radioactive isotope ⁷⁹Rb, which is part of the decay chains of ⁷⁹Zr and ⁸¹Nb. The goal is to determine the decay constants and half-lives of both ⁷⁹Sr and ⁷⁹Rb using experimental data.

## Background

A 10⁻⁶ mole sample of ⁷⁹Sr was collected and observed as it decayed via β⁺ emission to ⁷⁹Rb, which subsequently decayed to ⁷⁹Kr, also via β⁺ emission, accompanied by a detectable γ ray. The γ emissions from the second step were the only measurable quantity, allowing indirect analysis of both decays.

CSV files were provided, which contain the measured activity values and their uncertainties, taken once per minute for one hour. Alongside the following theory, this data was then used to compute the nuclear decay constant for ⁷⁹Rb. 

## Theory

The radioactive decay of Sr follows an exponential law:

$N_{Sr}(t) = N_{Sr}(0) · exp(−λ_{Sr}·t)$,

where $λ_{Sr}$ denotes the radioactive decay constant for ⁷⁹Sr, and $Nₛᵣ(t)$ the number of ⁷⁹Sr nuclei at time $t$. 

Then, the decay of ⁷⁹Rb can be solved according to the differential equation:

$\frac{dN_{Rb}(t)}{dt} = -\lambda_{Rb} N_{Rb} + \lambda_{Sr} N_{Sr}$,

since $^{79}Sr$ decays to $^{79}Rb$. The arising solution to this is then given as: 

$N_{Rb}(t) = N_{Sr}(0) \cdot \frac{\lambda_{Sr}}{\lambda_{Rb} - \lambda_{Sr}} \cdot \left[ \exp(-\lambda_{Sr} t) - \exp(-\lambda_{Rb} t) \right]$.

Then, the activity of decay is given as: 

$A_{Rb}(t) = λ_{Rb} · N_{Rb}(t) = N_{Sr}(0) \cdot \frac{\lambda_{Rb}\lambda_{Sr}}{\lambda_{Rb} - \lambda_{Sr}} \cdot \left[ \exp(-\lambda_{Sr} t) - \exp(-\lambda_{Rb} t) \right]$.

Using a parameter grid of $\lambda_{Rb}$ and $\lambda_{Sr}$ values, the a non-linear $\chi^2$ minimisation fit can be performed to approximate the decay constants of both $^{79}Rb$ and $^{79}Sr$.

Finally, the half-lives of each decay process may be calculated directly as: 

$t_{1/2} = \frac{\ln2}{\lambda}$.





