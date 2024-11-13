import GEOparse
import ftplib
import os

def download_ftp_file(ftp_url, local_filename):
    """
    Downloads a file from an FTP server and saves it locally.

    Parameters:
    ftp_url (str): URL of the file on the FTP server.
    local_filename (str): Name of the file to save locally.
    """
    # Parse the FTP URL
    #parse = ftplib.parse257
    ftp_host, ftp_path = ftp_url.replace("ftp://", "").split("/", 1)
    
    # Connect to the FTP server
    with ftplib.FTP(ftp_host) as ftp:
        ftp.login()  # Anonymous login
        ftp.cwd(os.path.dirname(ftp_path))
        filename = os.path.basename(ftp_path)

        # Download the file
        with open(local_filename, 'wb') as f:
            ftp.retrbinary('RETR ' + filename, f.write)

def download_geo_data(geo_accession):
    """
    Downloads and processes a GEO dataset including its supplementary files.

    Parameters:
    geo_accession (str): The accession number of the GEO dataset to download.
    """
    # Download the dataset using GEOparse
    gse = GEOparse.get_GEO(geo=geo_accession)

    # Print basic information about the dataset
    print(f"Title: {gse.metadata['title'][0]}")
    print(f"Summary: {gse.metadata['summary'][0]}")
    print(f"Overall design: {gse.metadata['overall_design'][0]}")

    # Iterate over each sample and print its name
    for gsm in gse.gsms.values():
        print(f"Processing sample: {gsm.name}")

    # Check if there are supplementary files and attempt to download them
    if 'supplementary_file' in gse.metadata:
        for file_url in gse.metadata['supplementary_file']:
            if file_url.startswith('ftp://'):
                filename = file_url.split('/')[-1]
                print(f"Downloading supplementary file: {filename}")
                download_ftp_file(file_url, filename)

# Example usage
geo_accession = "GSE244901"
download_geo_data(geo_accession)
