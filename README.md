Graphene-IR-Spectrum V5
=======================

In this version, site-resolved impurity averaging is taken into account, and 
the impurity scattering kernel, gamma, will be given an imaginary part independent 
of the width of the Fadeeva function, to smooth out the singularities on the real axis.
This will mimick the effects of the full self-consistent Born approximation, as in the
article of Castro-Neto et al.

Also, to reduce computational weight, q-points will only be sampled in a single wedge,
and H(q,w) H(-q,w) will be computed in this wedge. I don't have the proof, but H(q,w) H(-q,w)
is manifestly symmetric under all D6h operations.

Computation of anomalous IR response of grafted graphene, within a tight-binding model 

This code implements the theory presented in the paper 

"Phonon-disorder valley scattering modulation of optical transparency in graphene" ( arXiv:1407.8141 )

The goal of the project is to compute the contribution to the optical conductivity of 
disordered graphene of the Feynman diagram presented in the paper. 

