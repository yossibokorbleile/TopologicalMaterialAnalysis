---
title: 'ToMA: Analyse the structures of materials using tools from TDA'
tags:
  - Python
  - materials science
  - topological data analysis
  - persistent homology
  - 
authors:
  - name: Yossi Bokor Bleile
    orcid: 0000-0002-4861-9174
	corresponding: true
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Department of Mathematical Sciences, Aalborg Univeristy, Denmark
   index: 1
 - name: School of Mathematics and Statistics, University of Sydney, Australia
   index: 2
   index: 3
date: 02 July 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
---

# Summary

The relationship between physical properties of materials and the physical/chemical 
structure is an important aspect of materials science. Geometric and topological 
data analysis is a framework which extracts geometric and topological features from 
complex datasets. These features can reveal a variety of information about the material. 
In `@zifs-aau`, persistent homology and the accumulated persistence function were used
to investigate the relationship between medium range order and the thermal stability of 
pores in zeolitic imidazolate frameworks.

# Statement of need

`Topological Material Analysis (ToMA) ` is a package with a Graphical User Interface (GUI) 
and Command Line Interface (CLI) for chemists, physicists, and materials scientists interested
in using persistent homology in their research. In particular, it accepts as input a variety of 
formats containing information about atom locations in materials from simulations.

`ToMA` has been designed to be used in both exploratory analysis, and in well-establised pipelines
which make use of persistent homology. It includes functionality to generate persistent diagrams and 
accumulated persistence functions, as well as display representative cycles, currently limited to 1-cycles.

# Mathematics

We consider the labelled atom locations $X \subset \mathbb{R}^3$ as a set of *weighted points*, 
where the weights are determined by the atom type, see [CITATION] for how these weights
can be chosen. `ToMA` uses weighted $\alpha$-complexes to construct a filtered simplicial complex 
on atom locations. The persistent homology of this filtration is then computed. 

For readers in the definitions of weighted $\alpha$-complexes and persistent homology see [CITATION]
and [CITATION]. For a general introduction to Topological Data Analysis see [CITATION].

# Example Uses

Using the data from `@zifs-aau`, we demonstrate some of the outputs `ToMA` can produce. In particular, 
persistence diagrams and accumulated persistence functions. 


# Acknowledgements


# References