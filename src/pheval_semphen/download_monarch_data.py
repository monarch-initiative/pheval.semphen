# General imports
import os
import argparse
import requests
import tarfile
import pickle
import gzip
import shutil


def download_file_url(url: str, outdir: str, extract_gz: bool = False, overwrite: bool = False, extract_zip: bool = False, extract_tar: bool = False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/"
    """
    
    # Download and write file
    filename = os.path.join(outdir, url.split("/")[-1])
    fname = filename.split("/")[-1]

    if overwrite == True:
        if os.path.isfile(filename):
            print("- Warning, file {} already exists... Set overwrite to True to download and replace")
            return

    with open(filename, "wb") as f:
        r = requests.get(url)
        f.write(r.content)
    
    # Extract tar
    if extract_tar != False:
        file = tarfile.open(filename)
        file.extractall(kg_dir_path)
        file.close()

    
    # Extract gzip
    elif extract_gz != False:
        outpath = os.path.join(outdir, fname.replace(".gz", ""))
        with gzip.open(filename, 'rb') as f_in:
            with open(outpath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    elif extract_zip != False:
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(outdir)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Downloads monarch data necessary for running pheval-semphen')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        parser.add_argument("-m", "--mondo_release_tag", help="Release tag for mondo to download", required=False, type=str, default="v2025-04-01")
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # KG download nodes, edges paths
    phenio_path = os.path.join(args.project_dir, "phenio.db")
    mondo_path = os.path.join(args.project_dir, "mondo.sssom.tsv")
    gene_path = os.path.join(args.project_dir, "gene_mappings.sssom.tsv")
    hgnc_path = os.path.join(args.project_dir, "hgnc_complete_set.txt")
    mondo_nodes_path = os.path.join(args.project_dir, "mondo_nodes.tsv")

    # Create base project directory
    if not os.path.isdir(args.project_dir):
        print("- Creating project directory at {}".format(args.project_dir))
        os.makedirs(args.project_dir, exist_ok=True)
    

    # Download and upack 
    if not os.path.isfile(phenio_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking phenio database to {}".format(phenio_path))
        URL = 'http://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz'
        download_file_url(URL, args.project_dir, extract_gz=True, extract_zip=False, overwrite=False)
        print("- Download and upacking of phenio.db succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(phenio_path))
    

    # Download mondo sssom file
    if not os.path.isfile(mondo_path):

        # Fetch monarch mondo mappings
        print("- Downloading and upacking mondo sssom to {}".format(mondo_path))
        URL = 'https://data.monarchinitiative.org/mappings/latest/mondo.sssom.tsv'
        download_file_url(URL, args.project_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of mondo succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(mondo_path))

    
    # Download gene sssom file
    if not os.path.isfile(gene_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking gene sssom to {}".format(gene_path))
        URL = 'https://data.monarchinitiative.org/mappings/latest/gene_mappings.sssom.tsv'
        download_file_url(URL, args.project_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of genes succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(gene_path))
    

    # Download gene sssom file
    if not os.path.isfile(hgnc_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking hgnc data to {}".format(hgnc_path))
        URL = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
        download_file_url(URL, args.project_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of hgnc succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(hgnc_path))
    

    # Download mondo nodes
    if not os.path.isfile(mondo_nodes_path):

        # Fetch mondo nodes file
        print("- Downloading and upacking mondo data to {}".format(mondo_nodes_path))
        URL = 'https://github.com/monarch-initiative/mondo/releases/download/{}/mondo_nodes.tsv'.format(args.mondo_release_tag)
        download_file_url(URL, args.project_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of mondo nodes succesfull...")
    else:
        print("- Skipping mondo nodes download... File already exists at {}".format(mondo_nodes_path))
