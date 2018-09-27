---
title: 'TipToft: detecting plasmids contained in uncorrected long read sequencing data'
tags:
  - bioinformatics
  - plasmid typing
  - long reads
  - bacteria
authors:
 - name: Andrew J. Page
   orcid: 0000-0001-6919-6062
   affiliation: Quadram Institute Bioscience, Norwich Research Park, Norwich, UK.
 - name: Torsten Seemann
   orcid: 0000-0001-6046-610X
   affiliation: Doherty Institute, University of Melbourne, Australia.
  
date: 1 Oct 2018
bibliography: paper.bib
---

# Summary
With rapidly falling costs, long-read sequencing technologies, such as from Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT), are beginning to be used for outbreak investigations [@Faria2017; @Quick2015] and for rapid clinical diagnostics [@Votintseva2017]. Long-read sequencers from Oxford Nanopore can produce sequence reads in a matter of minutes and sequencers from PacBio can produce sequences in a number of hours compared to short-read sequencing technologies which takes hours/days. By reducing the time from swab to an actionable answer, genomics can begin to directly influence clinical decisions, with the potential to make a real positive impact for patients [@Gardy2018]. Clinically important genes can be horizontally acquired from plasmids such as those conferring anti-microbial resistance or virulence. With the increased speed afforded by long-read sequencing technologies comes increased base errors rates. The high error rates inherent in long-read sequencing reads require specialised tools to correct the reads [@Koren2017], however, these methods have substantial computational resource requirements often taking longer to run than the original time to generate the sequencing data, and can result in the loss of small, clinically important plasmids. 

We present TipToft which uses raw uncorrected reads to predict which plasmids are present in the underlying raw data. This provides an independent method for validating the plasmid content of a de novo assembly. It is the only tool which can do this from uncorrected long reads. TipToft is fast and can stream read data as it is produced, providing results more rapidly. We tested the software on 1975 samples (https://www.sanger.ac.uk/resources/downloads/bacteria/nctc/) sequenced with using long read sequencing technologies from PacBio, predicting plasmids from de novo assemblies using abricate (https://github.com/tseemann/abricate). It identified 84 samples containing plasmids with a 100% match to a plasmid sequence, but where no corresponding plasmid was present in the de novo assembly. Taking all the plasmids identified in the assemblies with 100% match, Tiptoft identified 97% (n=326) of these, representing 95% (236) of the samples. The software is written in Python 3 and is available under the open source license GNU GPL version 3 from https://github.com/andrewjpage/tiptoft.

# Acknowledgements
This work was supported by the Quadram Institute Bioscience BBSRC funded Core Capability Grant (project number BB/CCG1860/1) and by the Wellcome Trust (grant WT 098051).

# References