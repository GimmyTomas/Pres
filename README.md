# Pres
# Description
This is a repository containing the code used in <a href="https://arxiv.org/abs/2507.15110">arXiv:2507.15110</a> [1] to study gravitational atoms in binary systems. The code is public in order for other people to use it, check it and improve it.

(i) The file `horizon-fluxes.cpp` calculates the cloud mass loss rate and the resonance power ("P_res"). Compiling instructions are inside the file. The content of the main function should be self-explanatory.

(ii) The folder `100_Mc_break` contains a script, `100_Mc_break.py`, used to determine the bound on the cloud's mass due to the resonance with |100>. This script makes use of the other files in the folder.

(iii) The folder `Pres-results` contains some of the data shown in Figures 1, 2, and 4 of [1].

# Changelog
**November 2026.** Two corrections to the function `omega_Leaver` in `horizon-fluxes.cpp`:

1. The angular separation constant of the Kerr Klein-Gordon equation was evaluated with the wrong convention for the spheroidal eigenvalue. With c&sup2; = &atilde;&sup2;(&omega;&sup2; &minus; &alpha;&sup2;), the continued fraction needs the eigenvalue of the angular equation containing +c&sup2;cos&sup2;&theta;, that is `SpheroidalEigenvalue(l, m, -c²) - c²`, while the code used `SpheroidalEigenvalue(l, m, c²)`. The resulting error in the bound-state energies is of order &atilde;&sup2;&alpha;<sup>7</sup>, beyond the hyperfine splitting. For the parameters used in [1] (&alpha; = 0.2), it shifts the hyperfine splitting between |211&#10217; and |210&#10217; by 1.3%, the positions of all resonances on equatorial orbits by less than 0.03% (0.3% for |322&#10217;&rarr;|320&#10217;), and the decay widths by about 0.2% or less.
2. When the high-precision type `cpp_complex_100` is selected, b = (1 &minus; &atilde;&sup2;)<sup>1/2</sup> is now computed in high precision as well. It was previously a `double`, which limited the accuracy of the imaginary part of the frequency to a few times 10<sup>&minus;23</sup>, more than the decay width of states with &ell; &ge; 4. Results obtained in double precision are not affected.

The data in `Pres-results` and `100_Mc_break/Pres-data` have been regenerated with the corrected code.

# Questions
For any questions, feedback or suggestions, please send an email to <a href="mailto:gimmytomas@gmail.com">gimmytomas@gmail.com</a>.

# Citing
If you make use of this code, please consider citing the corresponding <a href="https://inspirehep.net/literature/2952263">paper</a>.
