import os
import platform
import json
import requests
import pubchempy as pc
from bs4 import BeautifulSoup
import subprocess
import pandas as pd
import numpy as np
import ssl
import PyCFMID
ssl._create_default_https_context = ssl._create_unverified_context

package_path = PyCFMID.__path__[0]
cwd = os.path.join(package_path, 'Windows')

def check_output_file(output_file=None):
    if output_file is None:
        try:
            os.mkdir('Output')
        except:
            pass
        output_file = os.path.join(os.getcwd(), 'Output', 'output.txt')
    return output_file
    

def check_input_file(input_dir=None):
    if input_dir is None:
        try:
            os.mkdir('Input')
        except:
            pass
        input_dir = os.path.join(os.getcwd(), 'Input')
    return input_dir


def fraggraph_gen(smiles, max_depth=2, ionization_mode='+', fullgraph=True, output_file=None):
    output_file = check_output_file(output_file)
    program = os.path.join(package_path, platform.platform().split('-')[0], 'fraggraph-gen.exe')
    cmd = os.path.join(os.getcwd(), program)
    cmd += ' ' + str(smiles)
    cmd += ' ' + str(max_depth)
    cmd += ' ' + str(ionization_mode)
    if fullgraph:
        cmd += ' fullgraph'
    else:
        cmd += ' fragonly'
    cmd += ' ' + str(output_file)
    subprocess.call(cmd, cwd =cwd)
    return parser_fraggraph_gen(output_file)
    
    
def parser_fraggraph_gen(output_file):
    with open(output_file) as t:
        output = t.readlines()
    output = [s.replace('\n', '') for s in output]
    nfrags = int(output[0])
    frag_index = [output[i].split(' ')[0] for i in range(1, nfrags+1)] 
    frag_mass = [output[i].split(' ')[1] for i in range(1, nfrags+1)] 
    frag_smiles = [output[i].split(' ')[2] for i in range(1, nfrags+1)]
    loss_from = [output[i].split(' ')[0] for i in range(nfrags+2, len(output))] 
    loss_to = [output[i].split(' ')[0] for i in range(nfrags+2, len(output))]
    loss_smiles = [output[i].split(' ')[0] for i in range(nfrags+2, len(output))]
    fragments = pd.DataFrame({'index': frag_index, 'mass': frag_mass, 'smiles': frag_smiles})
    losses = pd.DataFrame({'from': loss_from, 'to': loss_to, 'smiles': loss_smiles})
    return {'fragments': fragments, 'losses': losses}


def cfm_predict(smiles, prob_thresh=0.001, ion_source='ESI', ionization_mode='+', param_file='', config_file='', annotate_fragments=False, output_file=None, apply_postproc=True, suppress_exceptions=False):
    output_file = check_output_file(output_file)
    if ion_source == 'ESI':
        if ionization_mode == '+':
            config = 'esi_config'
        else:
            config = 'esi_config_neg'
    else:
        config = 'ei_config'
    if param_file == '':
        param_file = os.path.join(package_path, config + '.log')
    if config_file == '':
        config_file = os.path.join(package_path, config + '.txt')
    program = os.path.join(package_path, platform.platform().split('-')[0], 'cfm-predict.exe')
    cmd = os.path.join(os.getcwd(), program)
    cmd += ' ' + smiles
    cmd += ' ' + str(prob_thresh)
    cmd += ' ' + param_file
    cmd += ' ' +  config_file
    if annotate_fragments:
        cmd += ' ' + str(1)
    else:
        cmd += ' ' + str(0)
    cmd += ' ' + output_file
    if apply_postproc:
        cmd += ' ' + str(1)
    else:
        cmd += ' ' + str(0)
    if suppress_exceptions:
        cmd += ' ' + str(1)
    else:
        cmd += ' ' + str(0)
    subprocess.call(cmd, cwd =cwd)
    return parser_cfm_predict(output_file)
    
    
def parser_cfm_predict(output_file):
    with open(output_file) as t:
        output = t.readlines()
    output = [s.replace('\n', '') for s in output]
    low_energy = pd.DataFrame(columns=['mz', 'intensity'])
    medium_energy = pd.DataFrame(columns=['mz', 'intensity'])
    high_energy = pd.DataFrame(columns=['mz', 'intensity'])
    energy_level = 0
    for i in output:
        if 'energy0' == i:
            energy_level = 0
        elif 'energy1' == i:
            energy_level = 1
        elif 'energy2' == i:
            energy_level = 2
        elif '' == i:
            continue
        else:
            i = i.split(' ')
            i = [float(j) for j in i]
            if energy_level == 0:
                low_energy.loc[len(low_energy)] = i
            elif energy_level == 1:
                medium_energy.loc[len(medium_energy)] = i
            else:
                high_energy.loc[len(high_energy)] = i
    return {'low_energy': low_energy, 'medium_energy': medium_energy, 'high_energy': high_energy}


def cfm_id(spectrum_file, candidate_file, num_highest=-1, ppm_mass_tol=10, abs_mass_tol=0.01, prob_thresh=0.001, ion_source='ESI', ionization_mode='+', param_file='', config_file='', score_type='Jaccard', apply_postprocessing=True, output_file=None):
    output_file = check_output_file(output_file)
    if ion_source == 'ESI':
        if ionization_mode == '+':
            config = 'esi_config'
        else:
            config = 'esi_config_neg'
    else:
        config = 'ei_config'
    if param_file == '':
        param_file = os.path.join(package_path, config + '.log')
    if config_file == '':
        config_file = os.path.join(package_path, config + '.txt')
    program = os.path.join(package_path, platform.platform().split('-')[0], 'cfm-id.exe')
    cmd = os.path.join(os.getcwd(), program)
    cmd += ' ' + spectrum_file
    cmd += ' ' + 'AN_ID'
    cmd += ' ' + candidate_file
    cmd += ' ' + str(num_highest)
    cmd += ' ' + str(ppm_mass_tol)
    cmd += ' ' + str(abs_mass_tol)
    cmd += ' ' + str(prob_thresh)
    cmd += ' ' + param_file
    cmd += ' ' + config_file
    cmd += ' ' + score_type
    if apply_postprocessing:
        cmd += ' ' + str(1)
    else:
        cmd += ' ' + str(0)   
    cmd += ' ' + output_file
    subprocess.call(cmd, cwd =cwd)
    return parser_cfm_id(output_file)


def cfm_id_database(spectrum_dataframe, formula, energy_level='high', database='biodb', input_dir=None, num_highest=-1, ppm_mass_tol=10, abs_mass_tol=0.01, prob_thresh=0.001, ion_source='ESI', ionization_mode='+', param_file='', config_file='', score_type='Jaccard', apply_postprocessing=True, output_file=None):
    input_dir = check_input_file(input_dir)
    output_file = check_output_file(output_file)
    if ion_source == 'ESI':
        if ionization_mode == '+':
            config = 'esi_config'
        else:
            config = 'esi_config_neg'
    else:
        config = 'ei_config'
    if param_file == '':
        param_file = os.path.join(package_path, config + '.log')
    if config_file == '':
        config_file = os.path.join(package_path, config + '.txt')
    spectrum_file = os.path.join(input_dir, 'spectrum.txt')
    candidate_file = os.path.join(input_dir, 'candidate.txt')
    spectrum_file = write_spectrum(spectrum_dataframe, spectrum_file, energy_level)
    if database == 'biodb':
        candidates = search_biodatabase(formula, candidate_file)
    elif database == 'pubchem':
        candidates = search_pubchem(formula, candidate_file)
    else:
        candidates = pd.read_csv(database)
        output = pd.DataFrame({'ID': candidates.index, 'Smiles': candidates['SMILES']})
        output.to_csv(candidate_file, header=False, index=False, sep=' ')
    result = cfm_id(spectrum_file, candidate_file, num_highest, ppm_mass_tol, abs_mass_tol, prob_thresh, ion_source, ionization_mode, param_file, config_file, score_type, apply_postprocessing, output_file)
    return {'candidates':candidates, 'result':result}
    
    
def write_spectrum(spectrum_dataframe, spectrum_file, energy_level='high'):
    with open(spectrum_file, 'w+') as t:
        t.write('energy0' + '\n')
        if energy_level == 'low':
            for s in range(len(spectrum_dataframe)):
                t.write(str(spectrum_dataframe.iloc[s,0]) + ' ' + str(spectrum_dataframe.iloc[s,1]) + '\n')
        t.write('energy1' + '\n')
        if energy_level == 'medium':
            for s in range(len(spectrum_dataframe)):
                t.write(str(spectrum_dataframe.iloc[s,0]) + ' ' + str(spectrum_dataframe.iloc[s,1]) + '\n')
        t.write('energy2' + '\n')
        if energy_level == 'high':
            for s in range(len(spectrum_dataframe)):
                t.write(str(spectrum_dataframe.iloc[s,0]) + ' ' + str(spectrum_dataframe.iloc[s,1]) + '\n')
    return spectrum_file


def parser_cfm_id(output_file):
    output = pd.read_table(output_file, delim_whitespace=True, header=None, index_col=0)
    output.columns = ['Score', 'ID', 'Smiles']
    return output
    
    
def search_biodatabase(formula, structureDB, output_file=None):
    output_file = check_output_file(output_file)
    result = structureDB[structureDB['Formula'] == formula]
    output = pd.DataFrame({'ID': result.index, 'Smiles': result['SMILES']})
    output.to_csv(output_file, header=False, index=False, sep=' ')
    return result
    

def search_pubchem(formula, output_file=None, timeout=999):
    output_file = check_output_file(output_file)
    # get pubchem cid based on formula
    cids = pc.get_cids(formula, 'formula', list_return='flat')
    idstring = ''
    smiles = []
    inchikey = []
    all_cids = []
    # search pubchem via formula with pug
    for i, cid in enumerate(cids):
        idstring += ',' + str(cid)
        if ((i%100==99) or (i==len(cids)-1)):
            url_i = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + idstring[1:(len(idstring))] + "/property/InChIKey,CanonicalSMILES/JSON"
            res_i = requests.get(url_i, timeout=timeout)
            soup_i = BeautifulSoup(res_i.content, "html.parser")
            str_i = str(soup_i)
            properties_i = json.loads(str_i)['PropertyTable']['Properties']
            idstring = ''
            for properties_ij in properties_i:
                smiles_ij = properties_ij['CanonicalSMILES']
                if smiles_ij not in smiles:
                    smiles.append(smiles_ij)
                    inchikey.append(properties_ij['InChIKey'])
                    all_cids.append(str(properties_ij['CID']))
                else:
                    wh = np.where(np.array(smiles)==smiles_ij)[0][0]
                    all_cids[wh] = all_cids[wh] + ', ' + str(properties_ij['CID'])
    result = pd.DataFrame({'InChIKey': inchikey, 'SMILES': smiles, 'PubChem': all_cids})
    output = pd.DataFrame({'ID': result.index, 'Smiles': result['SMILES']})
    output.to_csv(output_file, header=False, index=False, sep=' ')
    return result
    