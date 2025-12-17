import pytest
import pandas as pd
import numpy as np
from rob_norm import rob_norm

def test_compare_with_r_reference():

    df_raw = pd.read_csv('./data/simulated_measurements.txt', index_col=0, sep='\t')
    
    results = rob_norm(df_raw, gamma_0=0.5)
    
    df_results = pd.read_csv('./data/simulated_measurements_normalized.txt', index_col=0, sep='\t')

    assert np.allclose(df_results, results['norm_data'])
