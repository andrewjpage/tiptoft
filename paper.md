---
title: 'TipToft: detecting plasmids contained in uncorrected long read sequencing data'
tags:
  - bioinformatics
  - plasmid typing
  - long read sequencing
  - bacteria
authors:
 - name: Andrew J. Page
   orcid: 0000-0001-6919-6062
   affiliation: 1
 - name: Torsten Seemann
   orcid: 0000-0001-6046-610X
   affiliation: 2
affiliations:
 - name: Quadram Institute Bioscience, Norwich Research Park, Norwich, UK.
   index: 1
 - name: Melbourne Bioinformatics, The University of Melbourne, Parkville, Australia.
   index: 2
date: 1 October 2018
bibliography: paper.bib
---

# Summary
With rapidly falling costs, long-read DNA sequencing technology from Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT), are beginning to be used for outbreak investigations [@Faria2017; @Quick2015] and rapid infectious disease clinical diagnostics [@Votintseva2017]. ONT instruments can produce data within minutes, and PacBio within hours compared to short-read sequencing technologies which takes hours/days. By reducing the time from swab to an actionable answer, genomics can begin to directly influence clinical decisions, with the potential for a positive impact for patients [@Gardy2018]. Clinically important genes, like those conferring animicrobial resistance or encoding virulence factors, can be horizontally acquired from plasmids. With the increased speed afforded by long-read sequencing technologies comes increased base errors rates. The high error rates inherent in long-read sequencing reads require specialised tools to correct the reads [@Koren2017], however, these methods require substantial computational requirements, and often take longer to run than the original time to generate the sequencing data, and can result in the loss of small, clinically important plasmids. 

We present ``TipToft`` which uses raw uncorrected reads to predict which plasmids are present in the underlying raw data. This provides an independent method for validating the plasmid content of a *de novo* assembly. It is the only tool which can do this from uncorrected long reads. ``TipToft`` is fast and can accept streaming input data to provide results in a realtime manner. We tested the software on 1975 samples (https://www.sanger.ac.uk/resources/downloads/bacteria/nctc/) sequenced with using long read sequencing technologies from PacBio, predicting plasmids from de novo assemblies using abricate (https://github.com/tseemann/abricate). It identified 84 samples containing plasmids with a 100% match to a plasmid sequence, but where no corresponding plasmid was present in the de novo assembly. Taking all the plasmids identified in the assemblies with 100% match, Tiptoft identified 97% (n=326) of these, representing 95% (236) of the samples. The software is written in Python 3 and is available under the open source GNU GPLv3 licence from https://github.com/andrewjpage/tiptoft.

# Acknowledgements
This work was supported by the Quadram Institute Bioscience BBSRC funded Core Capability Grant (project number BB/CCG1860/1) and by the Wellcome Trust (grant WT 098051).

# References
