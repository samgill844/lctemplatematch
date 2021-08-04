---
title: 'lctemplatematch: A Python package for to identify transit events from archival data'
tags:
  - Python
  - astronomy
authors:
  - name: Samuel Gill^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-4259-0155
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Centre for Exoplanets and Habitability, University of Warwick, Gibbet Hill Road, Coventry CV4 7AL, UK
   index: 1
 - name: Department of Physics, University of Warwick, Gibbet Hill Road, Coventry CV4 7AL, UK
   index: 2
date: 3 August 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

hey

# Statement of need

The Transiting Exoplanet Survey Satellite, TESS [@2015JATIS1a4003R] is well into the extended mission with over 131 published planets and over 4000 TESS objects of interest \footnote{tess.mit.edu; accessed 2021-07-20}. The majority of these planets have mean oribtal period of 8 days due in part to the geometric probability of a transit occurring and the 27-day per sector observing strategy. A number objects are seen to exhibit only a single planet-like transit in a TESS sector permitting only lower bounds to the orbital period (typically $\gtrsim$ 15 days). These systems are particularly interesting as they extend the possibility of recovering a planet within the habitable zone beyond late M-dwarfs which are not in the continuous viewing zone given the 27-day observing strategy. Work by [@2018A&A619A175C; @2019A&A631A83C] estimate 339 single transit candidates from southern ecliptic plane of the TESS primary mission with approximately 80\% expected to transit again in the extended mission. These numbers are broadly consistent with those by [@2019AJ15784V] and represents a bounteous sample of potentially long-period transiting planets and so efforts must be made to determine their physical properties via recovering the orbital period. 

Archival photometry provides means to "pre-cover" either the orbital period or a discrete set of possible aliases. Recent work by [@2019AJ15737Y] used a box-fitting least squared algorithm with fixed depth, width, and ingress/egress duration with injected transits in Kilodegree Extremely Little Telescope KELT; [@2007PASP119923P] data to efficiently recover the orbital periods of up to 1 year for planets larger than Saturn. other work by [@2020MNRAS4911548G] employed a bespoke template matching algorithm in using archival photometry from the Wide Angle Search for Planets WASP [@2006PASP1181407P] to constrain the orbital period to a discrete set of 3 possible aliases for the TESS single transit candidate TIC-238855958. The basis of this algorithm continues to be used as part of an semi-automated algorithm to detect transit events for NGTS observations of TESS single transit candidates that are being  monitored. There are subtleties to this algorithm that were not explicitly discussed by [@2020MNRAS4911548G]  or [@2020MNRAS4952713G] such as optimal step sizes and injection tests to determine a suitable threshold. The aim of this work is to formalise this algorithm and apply it to various test cases and active single transit candidates



# ``lctemplatematch``

The template matching algorithm described by estimates the difference in log-likelihood of a transit model over a null model. Given a series of $n$ observations described time, $t_n$, flux, $F_{n}$, and flux errors, $\sigma_{n}$ we calculate the chi-squared for the null model,

\begin{equation}
    \chi^2_{\rm null} = \sum_{n} \frac{(F_{n} - w)^2}{\sigma_n^2},
\end{equation}

where $w$ is the weighted mean of the dataset, 
\begin{equation}
w = \frac{\sum\limits_{n} F_n \sigma_n}{\sum\limits_{n} \sigma_n}.
\end{equation}

The log-likliehood of the null model is $\log \mathcal{L}_{\rm null} = -\chi^2_{\rm null} / 2$ and yields a quantity which can be compared to a that of a more complex transit model template, $m$. In this work we use the methods described in to solve Kepler's equations for the projected planet-star separation and the analytical approximation presented by to describe an object eclipsing a star with limb-darkening described by the power-2 law. Obtaining $m$ requires fitting an observed transit event with a fixed orbital period, $P_{s}$, to find the scaled orbital separation, $R_1 / a$, ratio of radii, $k = R_2/R_1$, and impact parameter,$b$. The limb-darkening parameters associated with the power-2 law $c$ and $\alpha$ can either be modelled directly or interpolated using the spectroscopic atmospheric parameters of the host star. 

![$\Delta \log \mathcal{L}$ for corresponding values of $T_i$ for the NGTS data of NGTS-2b are shown in the upper panel with peaks corresponding to $\Delta \log \mathcal{L} > 500$ annotated. Lower panels show each peak with the transit model (red), null model (red-dashed), raw NGTS data (black), and 5-minute binned data (blue).\label{fig:Figure_1}](Figure_1.png){ scale=80% }


Once the transit model template is established we create a series of epochs, $T_{i}$, for which the transit model will be centred and compared to the data. The extent of $T_{i}$ must be chosen such that any value of $t_n$ must be enclosed by a values of $T_{i}$ that extend to at least half the transit width, $W$, where,
%
\begin{equation}
    W = \frac{P_{s}}{\pi} \times \arcsin \left(  \frac{R_1}{a} \times \sqrt{\frac{(1+k^2) - b^2}{1 - b^2\left( \frac{R_1}{a}^2 \right)} }  \right).
\end{equation}
%
The step size of $T_{i}$ must be chosen to sufficiently and precisely recover transit events and is a function of $W$ and geometry. Typically, we find that $T_{i+1} - T_{i} = W/30$ is sufficient to recover epochs to a reasonable precision, although this may not be optimal for all transit shapes. For each epoch in $T_{i}$ we calculate,


\begin{equation}
    \chi^2_{\rm transit, i} = \sum_{n} \frac{(F_{n} - m_n(T_{i}))^2}{\sigma_n^2},
\end{equation}
%
and corresponding values of $\log \mathcal{L}_{\rm transit,i} = -\chi^2_{\rm transit, i} / 2$. In these tests, we set $w=1.$ when calculating $\chi^2_{\rm null}$. Epochs of $T_{i}$ where $\Delta \log \mathcal{L} = \log \mathcal{L}_{\rm transit,i} - \log \mathcal{L}_{\rm null}$ exceed 0 are where the transit model is well-matched to the data. In real datasets, some peaks are artefacts of white and red-noise and require sufficient thresholds for $\Delta \log \mathcal{L}$ to minimise false-positive detections of transit events. This can be achieved by estimating the typical noise profile of the dataset and creating a synthetic time series, $T_j$, between $-W/2$ and $W/2$ at the same cadence of the data. We inject the transit model at trial epochs, $T_{j}$, spanning 0 to $W$ and we can obtain the maximum value of $\Delta \log \mathcal{L}$ as a function of the span of data in transit. This is exemplified in Figure \autoref{fig:Figure1} with the NGTS data of the inflated hot-Jupiter NGTS-2b. Figure \ref{fig:Figure1} shows that $\Delta \log \mathcal{L}$ is not sensitive to transit detection when less than 30-minutes of data are present for a given epoch. This is expected as the equivalent drop in the transit model after the first 30 minutes of ingress is is approximately equal to the scatter of the raw data and so a null model cannot be excluded. 

We decided to search for transit events of NGTS-2b in NGTS data using a threshold of $\Delta \log \mathcal{L} > 100$ which is equivalent to at least 30 minutes of in-transit data (Figure \autoref{fig:Figure_2}). In all, we recover 11 transit events with the smallest value of $\Delta \log \mathcal{L} = 1117$. Using the best-fitting orbital solution there is an additional epoch which contains the last 15 minutes of egress. This did correspond to a peak of $\Delta \log \mathcal{L} = 20$ but is consistent with other peaks arising from red-noise and systematic offsets.


![$\Delta \log \mathcal{L}$ for corresponding values of $T_i$ for the NGTS data of NGTS-2b are shown in the upper panel with peaks corresponding to $\Delta \log \mathcal{L} > 500$ annotated. Lower panels show each peak with the transit model (red), null model (red-dashed), raw NGTS data (black), and 5-minute binned data (blue).\label{fig:Figure_2}](Figure_2.png)

# Acknowledgements

SG has been supported by STFC through consolidated grants ST/L000733/1 and ST/P000495/1

# References