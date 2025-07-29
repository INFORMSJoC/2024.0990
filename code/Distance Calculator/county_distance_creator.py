#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 05:37:42 2025

@author: davidbergman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 20:18:57 2024

@author: davidbergman
"""

import pandas as pd
import matplotlib.pyplot as plt
#from statsmodels.tsa.statespace.sarimax import SARIMAX
import numpy as np
import sys
# import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.utils import resample
# from scipy.stats import multivariate_normal
from scipy.linalg import cholesky
from numpy import linalg as la

from uszipcode import SearchEngine
from geopy.distance import geodesic

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3

def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

from math import radians, sin, cos, sqrt, atan2

def haversine(lat1, lon1, lat2, lon2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    
    # Radius of the Earth in kilometers
    R = 6371
    
    # Calculate the distance
    distance = R * c
    return distance


class ElectionPrediction():
    
    def __init__(self):
        
        self.path_to_us_cities = "../data/uscities.csv"
                
    def get_county_distances(self):
        
        self.us_cities_locations = pd.read_csv(self.path_to_us_cities)
        # sys.exit(1)
        self.us_cities_locations['lat_time_pop'] = self.us_cities_locations['lat']*self.us_cities_locations['population']
        self.us_cities_locations['lng_time_pop'] = self.us_cities_locations['lng']*self.us_cities_locations['population']
        
        group_columns = ['state_id', 'state_name', 'county_fips','county_name']
        self.county_distances = self.us_cities_locations.groupby(group_columns)[['lat_time_pop','lng_time_pop','population']].sum().reset_index()
        self.county_distances['weighted_lat'] = self.county_distances['lat_time_pop']/self.county_distances['population']
        self.county_distances['weighted_lng'] = self.county_distances['lng_time_pop']/self.county_distances['population']
        self.county_distances.drop(columns=['lat_time_pop','lng_time_pop'],inplace=True)
        self.county_distances.rename(columns={'population':'population_from_city_source'},inplace=True)
        
        self.county_fips_from_city_source = list(self.county_distances['county_fips'].unique())
        self.n_counties_from_city_source = len(self.county_fips_from_city_source)
        
        self.pairwise_distances = {}
        n_done = 0
        for county_1_index in range(self.n_counties_from_city_source):
            n_done += 1
            if n_done % 10 == 0:
                print("Done with ", n_done," out of ", self.n_counties_from_city_source)
            # if n_done == 100:
            #     sys.exit(1)
            county_1 = self.county_fips_from_city_source[county_1_index]
            county_1_extract = self.county_distances[self.county_distances['county_fips'] == county_1]
            if len(county_1_extract) != 1:
                print("ERROR")
                sys.exit(1)
            county_1_extract = county_1_extract.iloc[0]
            lat_1 = county_1_extract['weighted_lat']
            lng_1 = county_1_extract['weighted_lng']
            for county_2_index in range(county_1_index+1,self.n_counties_from_city_source):
                county_2 = self.county_fips_from_city_source[county_2_index]
                county_2_extract = self.county_distances[self.county_distances['county_fips'] == county_2]
                if len(county_2_extract) != 1:
                    print("ERROR")
                    sys.exit(1)
                county_2_extract = county_2_extract.iloc[0]
                lat_2 = county_2_extract['weighted_lat']
                lng_2 = county_2_extract['weighted_lng']
                distance = haversine(lat_1, lng_1, lat_2, lng_2)
                key = str(county_1) + "_" + str(county_2)
                self.pairwise_distances[key] = distance
                key = str(county_2) + "_" + str(county_1)
                self.pairwise_distances[key] = distance
            # sys.exit(1)
                
        
        # self.county_distances = self.us_cities_locations.groupby(group_columns)[['lat','lng','population']].apply(lambda x: ((x[['lat', 'lng']] * x['population']).sum() / x['population'].sum()))
   
if __name__ == "__main__":
    
    epred = ElectionPrediction()
    epred.get_county_distances()


    
    
    
