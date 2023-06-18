# AutoGenescan
**If you use this script, please cite our paper [McAllister et al., 2022, Nat Neurosci 25 446-457.](https://pubmed.ncbi.nlm.nih.gov/35379994/)**

Automated CAG calling pipeline for HTT genescan data. Uses the excellent Fragman package to detect peaks (Giovanny Covarrubias-Pazaran et al. 2016, BMC Genetics 17: 62). These peaks are then used to caculate graphs, modal CAG lengths and various instability/expansion indicies for downstream analyses.

See the associated .docx file for instructions on how to use the script.

Changelog

v1.5
Added an additional output graph without any marked peaks for presenations/publication.

v1.5.1
Removed data.table as compiling was an issue on Macs.

v1.6
Many thanks to Zach McLean for adding weighted CAG measurements (useful for mice and longer CAG models particularly). Also for adding easier control over peak height cutoffs and PCR size adjustment.
