## FRAME & NFC processes - stages and scripts to process emissions data and concentration data into and from FRAME-UK for the NFC

######################################################################################
#### **These scripts create emission input files for FRAME-UK (point and diffuse emissions), plus input files for FRAME-Europe.**
######################################################################################

#### for NFC 2018: FRAME is v. frame-9.15_17_1kb_2018 (@ 1km x 1km). Run on Rogue.
######################################################################################

**Process steps:**
----------------

1. Prepare emission files for FRAME-Europe and run FRAME-EUROPE

   * Using the EMEP 0.1&deg; x 0.1&deg; netcdfs with data for 1990 - 2018 (new in 2019)
   * https://www.ceip.at/new_emep-grid/01_grid_data
   * QAQC input totals, create summary sheet
   
2. Re-grid 50km output to 1km FRAME-UK boundary conditions

3. Create FRAME-UK emissions input files;

   * Update and check auxiliary data where necessary;
   
      * Rainfall (Eire and UK 1981 - 2010 mean @1km<sup>2</sup> resolution)
      * Land Cover (CEH Land Cover Map 2015 @1km<sup>2</sup> resolution)
      * Wind data (u & v wind rose data)


   * British National Grid. Input domain: xmin = -230000, xmax = 700000, ymin = 0, ymax = 1250000. 
   * Resolution = 1km<sup>2</sup>
   * UK point source emissions = https://naei.beis.gov.uk/data/map-large-source
   * UK diffuse emissions (with distribution) = https://naei.beis.gov.uk/data/map-uk-das
   * Eire point source emissions = https://prtr.eea.europa.eu/#/home
   * Eire diffuse totals = https://www.ceip.at/ms/ceip_home1/ceip_home/webdab_emepdatabase/
   * Eire diffuse distributions (latest 2016) = https://projects.au.dk/mapeire/spatial-results/
   * QAQC input totals, create summary sheet

4. run FRAME-UK via Nemesis

   * Sense check results
   * Evaluate via UKEAP measurement data = https://uk-air.defra.gov.uk/networks/network-info?view=ukeap
   * Obtain median bias for NH3 concentration correction
   * Output N-dep surfaces
   * Output corrected (and non-corrected) NH3 conc surfaces and inspect

5. Extract previous 2 years of corrected NH3 data and create a 3-year rolling mean for Critical Levels modelling

   * Use 'NFC_conc_cells.csv' for locations for CLe modelling

-----------------------------------------------------------------------------------------------------------------


_All scripts can be found in:_

\\nercbuctdb.ad.nerc.ac.uk\projects1\NEC03642_Mapping_Ag_Emissions_AC0112\NAEI_data_and_SNAPS\SNAPS_NFC_data

_All emissions data downloaded to:_

\\nercbuctdb.ad.nerc.ac.uk\projects1\NEC03642_Mapping_Ag_Emissions_AC0112\NAEI_data_and_SNAPS


