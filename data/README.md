# Dataset of microsm experiment with three arthropod predators on a shared prey

authors: Florian Dirk Schneider<sup>1</sup>*, Stefan Scheu<sup>2</sup> , and Ulrich Brose<sup>2</sup>

## Description

Initial and final population and biomass densities of full factorial combinations of three arthropod predator populations on one basal springtail population. 

The experiment ran in 30 x 30 x 15 cm microcosms over a period of 48 days.

Details can be found in 

- [Schneider, Scheu and Brose 2012 Body mass constraints on feeding rates determine the consequences of predator loss, *Ecology Letters* 15:436-443](http://onlinelibrary.wiley.com/doi/10.1111/j.1461-0248.2012.01750.x/abstract)
- [Schneider and Brose 2013 Beyond diversity: how nested predator effects control ecosystem functions, *Journal of Animal Ecology*  82:64-71](http://onlinelibrary.wiley.com/doi/10.1111/1365-2656.12010/abstract)

## Invalid replicates

**The replicate #43 was affected by extraordinarily high water content** and was excluded from the analyses. 

## key to dataset fields

- **ID**: replicate ID  
- **treat**: treatment binary code of scheme `*.*.*.*` with 0 for absence and 1 for presence of centipedes (Lithobius forficatus), spiders (Pardosa lugubris), predatory mites (Hypoaspis sp.) and springtails (Heteromurus nitidus), respectively. 
- **treat_name**: treatment name one of "null" (no populations),"control" (only springtails), "full" (full community), "lith", "pard", "hypo" (monocultures of centipedes, spiders, mites), "ko_lith", "ko_pard", "ko_hypo" (knockout cultures of centipedes, spiders, mites).
- **num_pred**: number of predator species	

Initial (t0) and final (t1) population densities given in individuals per microcosm (= 0.09 m^2)

- **N0_het**: average initial springtail density at t0 was 912 (± 528SD, n = 5) as estimated from heat extractions of 5 replicates at t0.
- **N0_hypo**: Due to delayed availability of mites at t0 and during the first week of the experiment, only 250 mites were introduced initially. Another 100 individuals were added after one week.
- **N0_pard**: counted manually
- **N0_lith**: counted manually
- **N1_het**: counts from heat extraction applied to a quarter of the
microcosm content. 
- **N1_hypo**: counts from heat extraction applied to a quarter of the
microcosm content. 
- **N1_pard**: counted manually
- **N1_lith**: counted manually

Initial and final biomass densities given in g per microcosm (= 0.09 m^2) 

- **B0_het**: estimated from population densities by multiplying with mean individual body mass of springtails = 0.10 mg (± 0.02SD)
- **B0_hypo**: estimated from population densities by multiplying with mean individual body mass of mites = 0.16 mg (± 0.02SD)
- **B0_pard**: weighed individually 
- **B0_lith**: weighed individually 
- **B1_het**: estimated from population densities by multiplying with mean individual body mass of springtails = 0.10 mg (± 0.02SD)
- **B1_hypo**: estimated from population densities by multiplying with mean individual body mass of mites = 0.16 mg (± 0.02SD)
- **B1_pard**: weighed individually 
- **B1_lith**: weighed individually 

- **B1_miclitt**: final microbial biomass on the litter layer was estimated from a fresh sample (2.8 g) taken at the end of the experiment by measuring substrate induced
O2 consumption in an electrolytic microrespirometer (see Schneider & Brose 2013 Journal of Animal Ecology 82:64-71).


## License

The experimental data are made available under the [Open Database License](http://opendatacommons.org/licenses/odbl/1.0/). Any rights in individual contents of the database are licensed under the [Database Contents License](http://opendatacommons.org/licenses/dbcl/1.0/)

----

<sup>1</sup> Institut des Sciences de l’Evolution, CNRS, Université Montpellier 2 - CC065, Montpellier Cedex 05, 34095, France

<sup>2</sup> Georg August University Göttingen, J.F. Blumenbach Institute of Zoology and Anthropology, Berliner Str. 28, 37073 Göttingen, Germany

\* Correspondence E-mail: florian.schneider@univ-montp2.fr

