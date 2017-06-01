# paper_2018_curr_opin_biotechnol
This repository provides the code used to generate the figures and tables in the Karr Lab's 2018 whole-cell modeling review in *Current Opinion in Biotechnology*. 

The code performs the following functions:
* Downloads the [BioModels](http://www.ebi.ac.uk/biomodels-main) database (Release 30)
* Parses all of downloaded SBML files for both the curated and non-curated models
* Uses heuristics to try to determine the mathematical type (flux balance analysis, logical, stochastic) of each model
* Extracts all of the model-level annotations (e.g. taxon, pathway, reference) of each model
* Inserts all of the model metadata (id, name, mathematical type, annotations) into a sqlite database
* Exports all of the model mdetadata to an Excel file
* Analyze the metadata and generate the figures and tables for the review paper

## How to run this code

### Requirements
* Pip
* Python

### Run the code
```
pip install -r requirements.txt
python paper_2018_curr_opin_biotechnol/biomodels_analysis.py
```

### View the results
All of the results will be stored in the directory `paper_2018_curr_opin_biotechnol/data`.

## How to cite this code
Check back in 2018

## License
The code is released under the [MIT license](LICENSE).

## Development team
This code was developed by the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.
