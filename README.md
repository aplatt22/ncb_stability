# Computer-Aided Design of Stability-Enhanced Nicotinamide Cofactor Biomimetics
## Alexandra P. Platt, Heidi Klem, Sam J. B. Mallinson, Yannick J. Bomble, Robert S. Paton

This is a repository containing structures and analysis scripts used for the publication "Computer-Aided Design of Stability-Enhanced Nicotinamide Cofactor Biomimetics"

#### Scripts include:
* Analysis/calculate_spearman_correlation.ipynb
    * Calculates the Spearman rank correlation for the DFT and semi-empirical predictive models
    * Required: Stability values (predicted and calculated) for the NCB library
    * Generated: PANDAS dataframe for R<sup>2</sup> and p-values
* Analysis/extract_dft_descriptors.ipynb
    * Reads DFT files (optimizations, vibrational frequencies, single point energies, NBO calculations) and creates a CSV file of relevant/calculated descriptors
    * Required: Goodvibes output file for NCB library (oxidized, reduced, and decomposed), NBO calculation QM output files, Fukui index CSV for each NCB
    * Generated: NCB_reduced_dft_descriptors.csv (descriptors) & key_atom_indices.csv (list of key atom indices)
* Analysis/select_ncbs_for_ts1_calcs.ipynb
    * Uses binning technique based on 8 key DFT descriptors to isolate NCBs to be used for determining the kinetic/thermodynamic correlation for the NCB library
    * Required: DFT level descriptors (only the 8 most important from ROBERT)
    * Generated: List of representative NCBs based on those 8 descriptors
* NCB_Generation/NCB_generation_clean.ipynb
    * Creates SMILES representations of NCBs after combinitorial addition to redox active core ring
    * Required: N/A
    * Generated: generated_NCBs.csv (SMILES of Library) & censo_assist_NCBs.csv (required for high-throughput calculation workflow)
* Calculation_Workflow/cdft.py
    * Calculates/extracts condensed Fukui indices from NBO output files, called in batch_submit_cdft.sh
    * Required: NBO output files for the original, reduced, and oxidized structures
    * Generated: CSV files with calculations including condensed Fukui indices at each atom of each structure
* Calculation_Workflow/fukui_nbo_input.py
    * Generates oxidized and reduced NBO QM input files
    * Required: NBO output files for original NCB structure
    * Generated: QM input files for single point energy calculations of NCBs +/- one electron
* Calculation_Workflow/high_throughput_NCB_calcs.sh
    * Performs conformational sampling of NCBs starting from the CSV generated wtih NCB_Generation/NCB_generation_clean.ipynb
    * Required: NCBs conda environment (NCBs_env.yml), generated_NCBs.csv from NCB_Generation/NCB_generation_clean.ipynb
    * Generated: QM input files for geometry optimizations and vibrational frequency calculations for NCB conformer ensembles
* Calculation_Workflow/make_batch_submit_cdft.py
    * Creates a file to submit cdft.py calculations for large amounts of files
    * Required: DFT NBO outputs for original, reduced, and oxidized NCBs with consistent suffix naming scheme
    * Generated: batch_submit_cdft.sh (bash script to submit for cdft.py calculations)
* Calculation_Workflow/sp_calc_setup.sh
    * Creates single point energy correction QM input files in Gaussian from optimized structures
    * Required: Geometry optimization Gaussian output files
    * Generated: NBO single point energy correction Gaussian input files & condensed Fukui index Gaussian input files
* Figure_Creation/cluster_conformers_energy_plotting.ipynb
    * Creates stripplot figure showing conformer energies of C1a with different clustering techniques
    * Required: Goodvibes output files for each type of clustering (found in Data/Conf_Sampling_Protocol/C1a_GoodVibes_Energies/)
    * Generated: Stripplot figure for C1a conformer energies with and without the use of CENSO
* Figure_Creation/descriptor_correlations.ipynb
    * Creates scatterplot figure showing relationships between various descriptors and NCB stability
    * Required: DFT descriptors (found in Data/Predictive_Modeling_Inputs/)
    * Generated: Scatterplots for descriptor relationships
* Figure_Creation/kinetic_thermodynamic_correlation.ipynb
    * Creates scatterplot figure showing the relationship between kinetic barriers and thermodynamic free energy differences in NCB decomposition
    * Required: Barriers and thermochemistry information for selected NCBs (found in Data/Kinetic_Thermodynamic_Correlation/)
    * Generated: Scatterplot showing correlation and calculating R<sup>2</sup>
* Figure_Creation/pca_plotting.ipynb
    * Creates principal component analysis representation for NCBs based on DFT descriptors
    * Required: 8 most important DFT descriptors (found in Data/Predictive_Modeling_Inputs/)
    * Generated: Scatterplot showing PCA representation of the forst two principal componenets
* Figure_Creation/robert_graphing.ipynb
    * Updates ROBERT-generated parity plots and heatmap
    * Required: ROBERT files (found in Data/ROBERT_Results/)
    * Generated: Parity plot and DFT heatmap SVG image files

#### Required Packages:
* RDKit (v2022.09.5)
* Seaborn (v0.12.2)
* Matplotlib (v3.7.1)
* Pandas (v2.0.0)
* GoodVibes (v3.2)
* AQME (v1.3.0)
* SciKit-Learn (v1.2.2)
* SciPy (v1.10.1)
* NumPy (v1.22.3)
* Open Babel (v3.1.1)
* DISCO (v1)
* cclib (v1.7.2)
* DBStep (v1.1.0)



