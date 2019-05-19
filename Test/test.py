# -*- coding: utf-8 -*-
"""
Created on Sun May 19 17:06:01 2019

@author: hcji
"""

import pandas as pd
from PyCFMID.PyCFMID import fraggraph_gen, cfm_predict, cfm_id_database

fragments = fraggraph_gen('CCCCN')
pred_ms = cfm_predict('CCCCN')

spectra = pd.DataFrame({'mz': [223.106608, 251.101730], 'intensity':[100.000000, 40.722900]})
cfm_id_database(spectra, formula='C17H22FN3O4S', database='pubchem')