# General imports
import os
import argparse
import copy
import multiprocessing as mp
from pathlib import Path
from typing import List, Set
from collections import Counter
import pandas as pd
import numpy as np
from typing import List, Set, Optional
from pydantic import BaseModel

# Pheval phenopacket utils
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil

# Semsimian
from semsimian import Semsimian




#########################################
### Preprocessing and data structures ###
class Patient(BaseModel):
	sample_name: str
	phenopacket_path: str
	
	# Pulled from phenopacket
	phenotype_ids: list
	phenotype_count: int
	
	disease_name: str
	disease_id: str
	
	gene_symbol: str
	gene_id: str
	
	# For sssom mappings
	disease_mapped: Optional[str] = None
	gene_mapped: Optional[str] = None


def read_sssom_to_lookup(fpath, exact_match=True):
	"""
	Assumes monarch initiative header style of # and then first line is header
	"""
	
	sssom_map = {}
	with open(fpath, 'r') as infile:
		header = 0
		cols = {}
		for line in infile:

			
			line = line.strip('\r').strip('\n')
			
			# Weird thing happening with the monarch sssom files (new line character most likely as the final line)
			if len(line) == 0:
				break
				
			if line[0] == "#":
				continue
			
			header += 1
			cols = line.split('\t')
			if header == 1:
				col_inds = {v:i for i,v in enumerate(cols)}
				continue
			
			# Our actual data here
			map_key = cols[col_inds["object_id"]]
			map_val = cols[col_inds["subject_id"]]
			map_type = cols[col_inds["predicate_id"]]
			
			# Only deal with exact matches
			if map_type != "skos:exactMatch":
				continue
			
			sssom_map.update({map_key:map_val})
	
	val_count = len(set(list(sssom_map.values())))
	print("- {} unique terms mapping to {} unique terms...".format(format(len(sssom_map), ','), format(val_count, ',')))
	
	# Now add self referencing keys
	curr_keys = list(sssom_map.keys())
	for k in curr_keys:
		v = sssom_map[k]
		sssom_map.update({v:v})
	
	return sssom_map


def gather_input_output_info(input_path, output_path, results_suffix="_results.tsv"):
	"""
	Input is allowed to be a directory containing .json phenopacket files
	Or input is allowed to be a filepath specifying a .json phenopacket
	"""
	
	if os.path.isdir(input_path):
		process_data = [os.path.join(input_path, fname) for fname in os.listdir(input_path) if fname.endswith(".json")]
		
	elif os.path.isfile(input_path) and input_path.endswith(".json"):
		process_data = [input_path]
	
	else:
		return None, None
	
	# Pre-format output file names
	out_data = [os.path.join(output_path, pname.split('/')[-1].replace(".json", results_suffix)) for pname in process_data]
	
	# Convert to dictionaries
	process_data = {f.split("/")[-1].replace(".json", ""):f for f in process_data}
	out_data = {f.split("/")[-1].replace(results_suffix, ""):f for f in out_data}
	
	return process_data, out_data


def phenopacket_paths_to_data(sample_path_dict):
	
	# Return data structure and function variables
	pobjs = {}
	multi_gene, multi_dis = {}, {}
	processed, total_samples = 0, len(sample_path_dict)
	
	# Loop through all phenopackets, extract / map data and return Patient objects
	for sname, spath in sample_path_dict.items():
		
		phenopacket_util = PhenopacketUtil(phenopacket_reader(spath))
		observed_phenotypes = phenopacket_util.observed_phenotypic_features()
		phenotype_ids = [observed_phenotype.type.id for observed_phenotype in observed_phenotypes]
		phenotype_count = len(phenotype_ids)
		
		# Default is to take first term 
		# Will display how many multi terms there are at end... should be few if not none)
		dis_obj = phenopacket_util.diseases()[0]
		dis_name, dis_id = dis_obj.term.label, dis_obj.term.id
		
		gene_obj = phenopacket_util.diagnosed_genes()[0]
		gene_symbol, gene_id = gene_obj.gene_symbol, gene_obj.gene_identifier
		
		if len(phenopacket_util.diagnosed_genes()) > 1:
			multi_gene.update({sname:''})
		
		if len(phenopacket_util.diseases()) > 1:
			multi_dis.update({sname:''})
		
		#print(sname)
		#print(dis_name, dis_id)
		#print(gene_symbol, gene_id)

		pobjs.update({sname:Patient.model_validate({"sample_name":sname,
													"phenopacket_path":spath,
													"phenotype_ids":phenotype_ids,
													"phenotype_count":phenotype_count,
													"disease_name":dis_name,
													"disease_id":dis_id,
													"gene_symbol":gene_symbol,
													"gene_id":gene_id})})
		
		processed += 1
		if processed % 1_000 == 0:
			print("- {}/{} phenopackets read into memory...".format(format(processed, ','),
																	format(total_samples, ',')))
		
	
	print("- Multi gene diagnosis phenopackets found {}...".format(len(multi_gene)))
	print("- Multi disease diagnosis phenopackets found {}...".format(len(multi_dis)))
	print("- {} Phenopackets information read into memory...".format(format(len(pobjs), ',')))
	return pobjs


def filter_non_zero_data(phen_data, sub_sample=False):
	
	# Remove zero phenotype count samples
	removed = 0
	for k in list(phen_data.keys()):
		if phen_data[k].phenotype_count == 0:
			del phen_data[k]
			removed += 1       
	print("- {} samples removed with zero phenotypes...".format(removed))
	
	# Default is no subsampleing
	if sub_sample != False:
	
		# Subsample (to make full pipeline connection easier)
		sub_samp = 0
		removed = 0
		for k in list(phen_data.keys()):
			if sub_samp >= sub_sample:
				del phen_data[k]
				removed += 1

			sub_samp += 1
		print("- {} samples removed with zero phenotypes...".format(removed))
		print("- {} samples remaining...".format(len(phen_data)))
	
	return phen_data
	
   
def parse_input_arguments_to_process_data(input_dir=None, output_dir=None, input_terms=None, output_file=None, mode="disease"):
	
	# Deal with input/output for phenopacket(s)
	if (input_dir != None) and (output_dir != None):
		inpaths, outpaths = gather_input_output_info(args.input_dir, 
													 args.output_dir, 
													 results_suffix="_{}_results.tsv".format(args.mode))

		# Make our output directory if it doesn't exist
		if not os.path.isdir(args.output_dir) :
			os.makedirs(args.output_dir)
		
		# Now read all phenopacket(s) into memory
		phen_data = phenopacket_paths_to_data(inpaths)
		phen_data = filter_non_zero_data(phen_data, sub_sample=False)
	
	# For single sample consisting of only HP terms
	elif (input_terms != None) and (output_file != None):

		# Split out phenotype terms
		phenotype_terms = args.input_terms.split(",")
		phen_data = {"single_sample":Patient.model_validate({"sample_name":"single_sample",
															 "phenopacket_path":'',
															 "phenotype_ids":phenotype_terms,
															 "phenotype_count":len(phenotype_terms),
															 "disease_name":'',
															 "disease_id":'',
															 "gene_symbol":'',
															 "gene_id":''})}
		# Now mimic the output file structure
		outpaths = {"single_sample":args.output_file}
	
	else:
		print("- Must provide either input_dir and output_dir, or input_terms and output_file")
		exit(1)
	
	return phen_data, outpaths


# For divying up data into batches for parallel processing
def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


#########################
### Semsimain wrapper ###
def get_phenotype_associations(semsim, phenotype_ids, outfile, symbol_map, name_map, mode="disease"):
	"""
	This algorithm leverages Semsimian + Monarchs phenio ontology to find the disease(s) 
	that are most associated with a patients phenotypes. A single .json phenopacket file can be passed in or
	a directory containing multiple .json phenopacket files. The patients observed phenotype terms are pulled
	from the data and are used as input to semsimian. Disease information is returned with the top associated
	ids appearing first in the list. 

	Subject prefixes within the phenio db begginning with "MONDO:" are compared
	"""

	# Define our db prefix term
	if mode == "disease":
		subject_prefix = "MONDO:"
	elif mode == "gene":
		subject_prefix = "HGNC:"

	# Ensure phenotype_ids are of set type
	if type(phenotype_ids) != type(set()):
		phenotype_ids = set(phenotype_ids)

	# Perform search (results are sorted in order of best ranking to worst ranking)
	results =  semsim.associations_search(object_closure_predicate_terms={"biolink:has_phenotype"},
											object_terms=phenotype_ids, # Must be set
											include_similarity_object=False,
											subject_terms=None,
											subject_prefixes=[subject_prefix],
											
											search_type="full",
											#score_metric="ancestor_information_content",
											score_metric="phenodigm",
											limit=10000,
											direction="object_to_subject")

	###results = [[0,0,"A"], [1,1,"B"], [2,2,"C"]] Testing purposes

	# Results are originally in form of [[score, details, mondo_id], ...]
	results = np.asarray(results).T

	# Convert our non-human readable names to human readable names
	symbols, names = [], []
	for v in results[2]:
		symbols.append(symbol_map[v])
		names.append(name_map[v])

	# Convert to df
	results_df = pd.DataFrame({"{}_id".format(mode):results[2],
							   "{}_symbol".format(mode):symbols,
							   "{}_name".format(mode):names,
							   "score":np.round(results[0].astype(float), decimals=4)})

	# Must set back to float otherwise results are not sorted properly
	results_df['score'] = results_df['score'].astype(float)

	# Now "rank" our results based on the score (pre sorted by best first)
	results_df['rank'] = results_df['score'].rank(method='dense', ascending=False)

	# Write out our data
	results_df.to_csv(outfile, sep='\t', header=True, index=False)

	return results_df


# Allows us to bulk process multiple samples without having to instantiate a new semsimian object for each sample
# This is a sub function of the main function, designed to be called in parallel (or single core)	
def proccess_samples(phenio_path, mode, mode_symbols, mode_names, patients_to_process):

	# Load necessary data into memory for semsimian processing
	semsim = Semsimian(predicates=["rdfs:subClassOf"], 
					spo=None, 
					resource_path=phenio_path)
	print("- Semsimian object loaded...")

	tt = 0
	for p in patients_to_process:
		outpath = p[1]
		phenotype_ids = p[0].phenotype_ids

		# Perform search (results are sorted in order of best ranking to worst ranking)
		results = get_phenotype_associations(semsim, 
											 phenotype_ids, 
											 outpath,
											 mode_symbols,
											 mode_names, 
											 mode=mode)
		
		tt += 1
		# Progress statement
		if tt % 100 == 0:
			print("- Processed {}/{}".format(format(tt, ','), format(len(patients_to_process), ',')))

	return None



if __name__ == "__main__":
	#################
	### ARG PARSE ###
	#################
	def mm():
		parser = argparse.ArgumentParser(description='This algorithm leverages Semsimian + Monarchs phenio ontology \
														to find the disease(s) or gene(s) that are most associated with a \
														patients phenotypes. A single .json phenopacket file can be passed in \
														or a directory containing multiple .json phenopacket files. \
														The patients observed phenotype terms are pulled from the data \
														and are used as input to semsimian. Disease information is returned \
														with the top associated ids appearing first in the list. \
														- Subject prefixes within the phenio db begginning with "MONDO:" are compared for diseases \
														- Subject prefixes within the phenio db begginning with "HGNC:" are compared for genes')

		# For multipe (or single) samples consisting of phenopackets. Must provide input_dir, output_dir
		parser.add_argument("-i", "--input_dir", help="Path to directory containing phenopackets", required=False, type=str, default=None)
		parser.add_argument("-o", "--output_dir", help="Path to output directory. Will be created if doesn't exist already", required=False, type=str, default=None)

		# For single sample consisting of only HP terms. Must provide terms, and outputfile path
		parser.add_argument("-it", "--input_terms", help="Comma seperated list of HP terms. For example -it HP:0001,HP:0002", required=False, type=str, default=None)
		parser.add_argument("-ot", "--output_file", help="Path to output file ", required=False, type=str, default=None)

		# Must provide algorithm with path to directory containing phenio.db file, and sssom mapping files
		parser.add_argument("-d", "--data_dir", help="Directory containing phenio.db, gene and disease mapping sssom files", required=True, type=str)
		parser.add_argument("-m", "--mode", help="Prioritization mode... disease or gene are allowed", required=True, choices=["disease", "gene"], type=str)
		
		# For multiprocessing purposes
		parser.add_argument("-c", "--num_proc", help="Number of cores to use for parallel processing", required=False, type=int, default=1)
		
		return parser.parse_args()
	###########
	args = mm()

	###############
	### PROGRAM ###
	###############

	# Deal with input / output arguments for sample data
	phen_data, outpaths = parse_input_arguments_to_process_data(args.input_dir,
																args.output_dir,
																args.input_terms,
																args.output_file,
																args.mode)
	
	# Copy our patient information / output paths for parallel processing (or single core processing)
	phen_base_data = [[copy.copy(v),copy.copy(outpaths[k])] for k,v in phen_data.items()] ###[0:100] # For testing

	# Load necessary data into memory for semsimian processing
	semsim = Semsimian(predicates=["rdfs:subClassOf"], 
					   spo=None, 
					   resource_path=os.path.join(args.data_dir, "phenio.db"))
	print("- Semsimian object loaded...")

	# Gene name mappings
	hgnc_df = pd.read_csv(os.path.join(args.data_dir, "hgnc_complete_set.txt"), sep="\t", header=0, low_memory=False)
	hgnc_symbols = {k:v for k,v in zip(hgnc_df["hgnc_id"], hgnc_df["symbol"])}
	hgnc_names = {k:v for k,v in zip(hgnc_df["hgnc_id"], hgnc_df["name"])}
	print("- Gene mappings loaded...")

	# Disease name mappings
	mondo_df = pd.read_csv(os.path.join(args.data_dir, "mondo_nodes.tsv"), sep="\t", header=0, low_memory=False)
	mondo_symbols = {k:v for k,v in zip(mondo_df["id"], mondo_df["name"]) if k.startswith("MONDO:")}
	mondo_names = {k:v for k,v in zip(mondo_df["id"], mondo_df["description"]) if k.startswith("MONDO:")}
	print("- Disease mappings loaded...")

	# Unused (currently... mapping arbitrary name spaces of genes / disease is necessary for benchmarking purposes)
	##gene_map = read_sssom_to_lookup(os.path.join(args.data_dir, "mondo.sssom.tsv"))
	##disease_map = read_sssom_to_lookup(os.path.join(args.data_dir, "gene_mappings.sssom.tsv"))
	
	# Set our input namespace mappings for semsimian function / results
	if args.mode == "gene":
		entity_symbols = hgnc_symbols
		entity_names = hgnc_names

	elif args.mode == "disease":
		entity_symbols = mondo_symbols
		entity_names = mondo_names

	# For dev / testing purposes...
	# # Now process all samples
	# tt = 0
	# for k,v in phen_data.items():
		
	# 	# Perform search (results are sorted in order of best ranking to worst ranking)
	# 	results = get_phenotype_associations(semsim, 
	# 										 v.phenotype_ids, 
	# 										 outpaths[k],
	# 										 entity_symbols,
	# 										 entity_names, 
	# 										 mode=args.mode)

	# 	################################################################################
	# 	### For debugg / testing purposes. Allows us to display relevant information ###
	# 	# Patient input data
	# 	# patient_df = {"sample_name":[v.sample_name],
	# 	# 			 "sample_phenotype_count":[len(v.phenotype_ids)],
	# 	# 			 "disease_id":[disease_map[v.disease_id]],
	# 	# 			 "disease_name":[v.disease_name],
	# 	# 			 "gene_id":[v.gene_id],
	# 	# 			 "gene_name":[v.gene_symbol]}
		
	# 	# Display data (if we want)
	# 	#display(pd.DataFrame(patient_df))
	# 	#display(results[0:20])
	# 	#tt += 1
	# 	########
		
	# 	# Progress statement
	# 	if tt % 500 == 0:
	# 		print("- Processed {}/{}".format(format(tt, ','), format(len(phen_data), ',')))

	
	# Now deal with potential parallel processing (if cpu cores > 1)
	if args.num_proc > 1:
		
		# Divy up necessary input data for semsimian into parallel processing chunks
		div_entity_symbols = [copy.copy(entity_symbols) for i in range(0, args.num_proc)]
		div_entity_names = [copy.copy(entity_names) for i in range(0, args.num_proc)]
		div_phenio_paths = [os.path.join(args.data_dir, "phenio.db") for i in range(0, args.num_proc)]
		div_modes = [copy.copy(args.mode) for i in range(0, args.num_proc)]
		div_phen_data = divide_workload(phen_base_data, num_proc=args.num_proc)
		print("- Parallel processing with {} cores...".format(args.num_proc))
		
		# Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
		output = mp.Queue()
		pool = mp.Pool(processes=args.num_proc)
		results = [pool.apply_async(proccess_samples, args=(ph, m, msy, mn, pdata)) for ph, m, msy, mn, pdata in zip(div_phenio_paths, 
																													 div_modes, 
																													 div_entity_symbols, 
																													 div_entity_names, 
																													 div_phen_data)]
		output = [p.get() for p in results]
		pool.close()
		pool.join()
	
	# Single core processing (do not need to use multiprocessing overhead)
	else:
		print("- Single core processing...")
		proccess_samples(os.path.join(args.data_dir, "phenio.db"), 
						 args.mode, 
						 entity_symbols, 
						 entity_names, 
						 phen_base_data)
	
	print("- Done processing {} samples".format(format(len(phen_base_data), ',')))