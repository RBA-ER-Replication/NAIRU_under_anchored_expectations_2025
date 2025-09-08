# NAIRU_under_anchored_expectations_2025
**Authors**: Alexander Ballantyne and Tom Cusbert

## Description

This repo contains data and replication files for Ballantyne and Cusbert (2025) The NAIRU under Anchored Inflation Expectations in the Australian Economic Review, September Policy Forum.


## License
This repository contains both code and data, each under different licenses:

Code: Licensed under the BSD 3-Clause License.
All datasets in the Data subfolder: Licensed under Creative Commons Attribution 4.0 International (CC-BY-4.0)
Please refer to the respective license files for full details.

## Coding Languages

**Eviews:** We used Eviews 14 to run the code. The software is available at https://www.eviews.com/home.html. 

Questions, comments and bug reports can be sent to `cusbertt`  at domain `rba.gov.au`.

## Folder Contents

- `README` is this read me file.
- `LICENSE` contains the license agreement that covers all code in this repo (Data subfolder is covered by separate agreement).
- `Data` subfolder contains the inmput data, data description and license document.
- `ReplicationCode` subfolder contains the eviews programs (.prg) to run generate the results in the paper and a worfile (.wf1) with the results done already.

### Subfolder: ReplicationCode 


To build the models and generate the results run `nairu_model_main_code.prg`.

It calls:
- `model_orig_revised.prg`
- `model_orig_precovid.prg`
- `Uncertainty/nairu_uncertainty_main.prg` (which calls `Uncertainty/uncertainty_orig_revised.prg`)
- `anchored_decay.prg` (NB. this file changes the unemployment data so do not run other code afterward)
- `check_sspace.prg` is a subroutine called by the other scripts
 
The code produces the tables (param_breaks, unrestricted_tests, and appendixtable) in the paper and Graphs 2,4,5,7,8,9, 10 and 11. These will be in the eviews workfile with naming g1_infpc, g2_ulcpc, etc.

Graph 3 and Graph 6 are made from different models and sources.

Tghe `Output` subfolder has the spreadsheets constructed with graph data.

The code will produce a workfile but the folder also includes the completed workfile `NAIRU_model_Output.wf1`.


## Code Structure

```
Repo
├── README.md
├── LICENSE
├── Data
│   ├── LICENSE-data.md (License agreement for data in this subfolder)
│   ├── NAIRU_update_data.csv
│   └── Input data details and sources.txt 
└── ReplicationCode
    ├── model_orig_revised.prg
    |── model_orig_precovid.prg
    |── anchored_decay.prg
    |── check_sspace.prg
    |── NAIRU_model_Output.wf1
    |── Uncertainty
    |    |──nairu_uncertainty_main.prg
    |    └──uncertainty_orig_revised.prg
    └── Output
```
