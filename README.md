# pheval_semphen

semsimian + phenio.db for phenopacket disease prioritization<br>
Example of how to run using pheval<br>
 - Download http://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz<br>
 - Edit config file located in test_configs so that the path_to_phenio references the path to the file we just downloaded
 - pheval run -i test_configs/ -t path/to/phenopackets/ -r "semphenphevalrunner" -o test_output/
<br><br>

Example of how to run just using python executable<br>
 - Set up environment
   - git clone -b disease_and_gene_prioritization --single-branch git@github.com:monarch-initiative/pheval.semphen.git
   - cd pheval.semphen
   - poetry config virtualenvs.in-project true
   - poetry install
   - poetry shell
 - Download necessary monarch kg data
   - python src/pheval_semphen/semphen.py -p directory/to/store/monarch_kg_data
 - Run on directory of phenopackets (or single phenopacket) for disease/gene prioritization
   - python src/pheval_semphen/semphen.py -i path/to/phenopacket(s) -o path/to/results/directory -d path/to/monarch_kg_data -m disease
   - python src/pheval_semphen/semphen.py -i path/to/phenopacket(s) -o path/to/results/directory -d path/to/monarch_kg_data -m gene
 - Run with command line input of HP terms for disease/gene prioritization
   - python src/pheval_semphen/semphen.py -it HP:0003635,HP:0009002 -ot path/to/results_outfile.tsv -d path/to/monarch_kg_data -m disease
   - python src/pheval_semphen/semphen.py -it HP:0003635,HP:0009002 -ot path/to/results_outfile.tsv -d path/to/monarch_kg_data -m gene

# Acknowledgements

This [cookiecutter](https://cookiecutter.readthedocs.io/en/stable/README.html) project was developed from the [monarch-project-template](https://github.com/monarch-initiative/monarch-project-template) template and will be kept up-to-date using [cruft](https://cruft.github.io/cruft/).
