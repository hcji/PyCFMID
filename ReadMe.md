# PyCFMID
***

This package is a python warpper for [CFM-ID](https://sourceforge.net/projects/cfm-id/) which provides a method for accurately and efficiently identifying metabolites in spectra generated by electrospray tandem mass spectrometry (ESI-MS/MS). The program 
uses Competitive Fragmentation Modeling to produce a probabilistic generative model for the MS/MS fragmentation process and machine 
learning techniques to adapt the model parameters from data.  

## Note
    1. This is not an official wrapper.  
    2. Only Windows platform is currently supported because of a compile error in Linux (issue #1).   
    3. The CFM-ID version is 2.4, download from https://sourceforge.net/projects/cfm-id/.  
    4. CFM-ID 3.x version is only available with webserver http://cfmid3.wishartlab.com.   
    5. Official CFM-ID 4.x version is already available with docker at https://hub.docker.com/r/wishartlab/cfmid.  

### Install
    
    pip install git+git://github.com/hcji/PyCFMID@master

## Usage

### fraggraph_gen
This program produces a complete fragmentation graph or list of feasible fragments for an input molecule. It systematically breaks bonds within the molecule and checks for valid resulting fragments as described in section 2.1.1 of the above publication.
    
    from PyCFMID.PyCFMID import fraggraph_gen
    frags = fraggraph_gen(smiles, max_depth=2, ionization_mode='+', fullgraph=True, output_file=None)

**smiles**: The smiles strings for the input molecule to fragment.  
**max depth**: The depth to which the program should recurse when computing the tree. e.g. depth 1 would be just the original molecule and its immediate descendants, depth 2 would allow those descendants to break one more time, etc.    
**ionization_mode**: Whether to generate fragments using positive ESI or EI, or negative ESI ionization. + for positive mode ESI [M+H], - for negative mode ESI [M-H], * for positive mode EI [M+].    
**fullgraph**: This specifies the type of output. fullgraph (default) will also return a list of the connections between fragments and their corresponding neutral losses. otherwise, it will return a list of unique feasible fragments with their masses.    
**output file**: (optional) The name and path of a file to write the output to. If this argument is not provided, it will make a dir in the working path.    

### cfm_predict
This program predicts spectra for an input molecule given a pre-trained CFM model.

    from PyCFMID.PyCFMID import cfm_predict
    spectra = cfm_predict(smiles, prob_thresh=0.001, param_file='', config_file='', annotate_fragments=False, output_file=None, apply_postproc=True, suppress_exceptions=False)
    
**smiles**: The smiles strings for the input molecule to fragment.  
**prob_thresh**: The probability below which to prune unlikely fragmentations during fragmentation graph generation (default 0.001).    
**ion_source**: The ion source of mass. Usually, 'EI' for GC-MS, and 'ESI' for LC-MS/MS. Will not used if param_file is given.       
**param_file**: (optional) The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in the current directory). This file is the output of cfm-train.    
**config_file**: (optional) The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in the current directory). This needs to match the file passed to cfm-train during training.    
**annotate_fragments**:(optional) Whether to include fragment information in the output spectra.
**output file**: (optional) The name and path of a file to write the output to. If this argument is not provided, it will make a dir in the working path.    
**apply_postproc**: (optional) Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first). If turned off, will output a peak for every possible fragment of the input molecule, as long as the prob_thresh argument above is set to 0.0.   
**suppress_exceptions**: (optional) Suppress most exceptions so that the program returns normally even when it fails to produce a result.   

### cfm_id
Given an input spectrum and a list of candidate smiles strings, this program computes the predicted spectrum for each candidate and compares it to the input spectrum. It returns a ranking of the candidates according to how closely they match. The spectrum prediction is done using a pre-trained CFM model.

    from PyCFMID.PyCFMID import cfm_id
    result = cfm_id(spectrum_file, candidate_file, num_highest=-1, ppm_mass_tol=10, abs_mass_tol=0.01, prob_thresh=0.001, param_file='', config_file='', score_type='Jaccard', apply_postprocessing=True, output_file=None)

**spectrum_file**: The filename where the input the spectrum. see [example_spectra.txt](https://sourceforge.net/p/cfm-id/code/HEAD/tree/cfm/example_spec.txt) as an example.   
**candidate_file**: The filename where the input list of candidate structures can be found as line separated 'id smiles_or_inchi' pairs. see[example_candidates.txt](https://sourceforge.net/p/cfm-id/code/HEAD/tree/cfm/example_candidates.txt) as an example.   
**num_highest** (optional): The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order).      
**ppm_mass_tol**: (optional) The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10 ppm).   
**abs_mass_tol**: (optional) The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 0.01 Da).      
**prob_thresh**: The probability below which to prune unlikely fragmentations during fragmentation graph generation (default 0.001).    
**ion_source**: The ion source of mass. Usually, 'EI' for GC-MS, and 'ESI' for LC-MS/MS. Will not used if param_file is given.   
**param_file**: (optional) The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in the current directory). This file is the output of cfm-train.    
**config_file**: (optional) The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in the current directory). This needs to match the file passed to cfm-train during training.    
**score_type**: (optional) The type of scoring function to use when comparing spectra. Options: Jaccard (default), DotProduct.    
**output file**: (optional) The name and path of a file to write the output to. If this argument is not provided, it will make a dir in the working path.    
**apply_postproc**: (optional) Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first). If turned off, will output a peak for every possible fragment of the input molecule, as long as the prob_thresh argument above is set to 0.0.   

### cfm_id_database
Given an input spectrum, this program retrieves candidates automatically and computes the predicted spectrum for each candidate and compares it to the input spectrum. It returns a ranking of the candidates according to how closely they match. The spectrum prediction is done using a pre-trained CFM model.

    from PyCFMID.PyCFMID import cfm_id_database
    cfm_id_database(spectrum_dataframe, formula, energy_level='high', database='biodb', input_dir=None, num_highest=-1, ppm_mass_tol=10, abs_mass_tol=0.01, prob_thresh=0.001, param_file='', config_file='', score_type='Jaccard', apply_postprocessing=True, output_file=None)
    
**spectrum_dataframe**: A two-column dataframe with m/z and intensity of a spectrum.  
**formula**: The formula of the candidates.  
**energy_level**: The energy_level of the mass spectrometry.  Options: high (default), medium, low.  
**database**: 'biodb' for biological database, 'pubchem' for PubChem database, or a file path for a custom candidate list. 
other parameters are the same as **cfm_id**
