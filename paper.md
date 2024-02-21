---
title: 'TAMA: Topological Amorphous Material Analysis'
tags:
  - Python
  - materials science
  - topological data analysis
  - persistent homology
  - chemisty
authors:
  - name: Yossi Bokor Bleile
    orcid: 0000-0002-4861-9174
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Mathematical Sciences, Aalborg Univeristy, Denmark
   index: 1
 - name: School of Mathematics and Statistics, The Univeristy of Sydney, Australia
   index: 2
date: ???
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
---

# Summary

Understanding the relationship between structural features and physical properties of materials is an important aspect of chemistry and materials science, in particular developing new materials. There are various *scales* at which we can examine the structure of materials: what are the very *local* structures between neighbouring atoms, how do these local structures fit together, and finally, what is the global or long range structure? These structures are examined with various degrees of difficulty, depending on the material. 

# Statement of need

`TAMA`, or `TopologicalAmorphousMaterialAnalysis` is a Python application for researchers in 
chemsitry and materials science, who wish to use persistent homology to analyse the structure 
of materials.  `TAMA` is designed to enable users with varying degrees of programming knowledge 
to leverage persistent homology for their research. It has both a GUI and CLI, enabling use in a 
wide range of computing environments. Through the GUI, users are able to perform exploratory 
analysis and visualise structural features. This removes the need for individual users to master 
the mathematical concepts required to perform such analysis, allowing users to focus on the 
analysis itself. While a certain level of knowledge is required to interpret the output and 
visualisation, these are easily understood in the materials science context.

The application is easily used by both students and experienced researchers, in chemistry
and materials science, as well as data scienctists wishing to explore the use cases of 
topological data analysis, in particular persistent homology. 

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References