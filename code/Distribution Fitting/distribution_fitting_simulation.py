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

        self.path_to_county_pop = "../data/cc-est2022-agesex-all.csv"
        self.path_to_election_results = "../data/countypres_2000-2020.csv"
        self.path_to_us_cities = "../data/uscities.csv"
    
    def read_data(self):
        
        self.county_pop_raw = pd.read_csv(self.path_to_county_pop)
        self.election_results_raw = pd.read_csv(self.path_to_election_results)


        # step 0: assign arbitrary FIPS  to the three "county" voter logs that aren't actually voters

        self.election_results_raw['county_fips'] = np.where(
            (self.election_results_raw['state'] == 'CONNECTICUT') & (self.election_results_raw['county_name'] == 'STATEWIDE WRITEIN'),
            111111111,
            self.election_results_raw['county_fips']
            )
        
        self.election_results_raw['county_fips'] = np.where(
            (self.election_results_raw['state'] == 'MAINE') & (self.election_results_raw['county_name'] == 'MAINE UOCAVA'),
            111111112,
            self.election_results_raw['county_fips']
            )
        
        self.election_results_raw['county_fips'] = np.where(
            (self.election_results_raw['state'] == 'RHODE ISLAND') & (self.election_results_raw['county_name'] == 'FEDERAL PRECINCT'),
            111111113,
            self.election_results_raw['county_fips']
            )

        ## 1. fix kansas city
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == 'KANSAS CITY')&
                (self.election_results_raw['state_po'] == 'MO')
                ),
            'JACKSON',
            self.election_results_raw['county_name']
            )
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == 'LA SALLE')&
                (self.election_results_raw['state_po'] == 'LA')
                ),
            'LASALLE',
            self.election_results_raw['county_name']
            )
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == 'SAINT LOUIS')&
                (self.election_results_raw['state_po'] == 'MN')
                ),
            'ST. LOUIS',
            self.election_results_raw['county_name']
            )
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == 'DONA ANA')&
                (self.election_results_raw['state_po'] == 'NM')
                ),
            'DOÑA ANA',
            self.election_results_raw['county_name']
            )
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == 'SHANNON')&
                (self.election_results_raw['state_po'] == 'SD')
                ),
            'OGLALA LAKOTA',
            self.election_results_raw['county_name']
            )
        self.election_results_raw['county_name'] = np.where(
            (
                (self.election_results_raw['county_name'] == "ST MARY'S")&
                (self.election_results_raw['state_po'] == 'MD')
                ),
            "ST. MARY'S",
            self.election_results_raw['county_name']
            )
        map_of_fips = {}
        map_of_fips[24005] = 24510
        map_of_fips[2938000] = 29095
        map_of_fips[29510] = 29189
        map_of_fips[46102] = 46113
        map_of_fips[51515] = 51019
        map_of_fips[51600] = 51059
        map_of_fips[51620] = 51067
        map_of_fips[51760] = 51159
        map_of_fips[51770] = 51161
        
        for key in map_of_fips.keys():
            self.election_results_raw['county_fips'] = np.where(
                self.election_results_raw['county_fips'] == key,
                map_of_fips[key],
                self.election_results_raw['county_fips']
                )

        self.election_results_raw['county_name'] = self.election_results_raw['county_name'].str.upper()
        self.election_results_raw['county_name'] = self.election_results_raw['county_name'].str.replace(" COUNTY","")
        self.election_results_raw['county_name'] = self.election_results_raw['county_name'].str.replace(" PARISH","")
        self.election_results_raw['county_name'] = self.election_results_raw['county_name'].str.replace(" CITY","")
       
        self.county_pop_raw['CTYNAME'] = self.county_pop_raw['CTYNAME'].str.upper()
        self.county_pop_raw['CTYNAME'] = self.county_pop_raw['CTYNAME'].str.replace(" COUNTY","")
        self.county_pop_raw['CTYNAME'] = self.county_pop_raw['CTYNAME'].str.replace(" PARISH","")
        self.county_pop_raw['CTYNAME'] = self.county_pop_raw['CTYNAME'].str.replace(" CITY","")

        
    def predict_county_pop(self):
        
        self.voters_per_county = self.county_pop_raw[['STATE', 'COUNTY', 'STNAME', 'CTYNAME', 'YEAR','AGE18PLUS_TOT']]
        self.voters_per_county = self.voters_per_county.groupby(by=['STATE', 'STNAME', 'CTYNAME', 'YEAR'])['AGE18PLUS_TOT'].sum().reset_index()
        self.voters_per_county = self.voters_per_county[self.voters_per_county['YEAR'].isin([2,3,4])]

        grouped = self.voters_per_county.groupby(['STNAME', 'CTYNAME'])
        
        state = []
        county = []
        pop2020 = []
        pop2021 = []
        pop2022 = []
        pop2023 = []
        pop2024 = []

        for name, group in grouped:
            
            year_2020 = group[group['YEAR'] == 2]['AGE18PLUS_TOT'].iloc[0]
            year_2021 = group[group['YEAR'] == 3]['AGE18PLUS_TOT'].iloc[0]
            year_2022 = group[group['YEAR'] == 4]['AGE18PLUS_TOT'].iloc[0]
            
            growth_rate_2020_to_2021 = np.log(year_2021 / year_2020)
            growth_rate_2021_to_2022 = np.log(year_2022 / year_2021)
        
            mean_growth_rate = np.exp((np.log(1 + growth_rate_2020_to_2021) + np.log(1 + growth_rate_2021_to_2022)) / 2) - 1
    
            year_2023 = year_2022 * (1 + mean_growth_rate)
            year_2024 = year_2023 * (1 + mean_growth_rate)
            
            state.append(name[0])
            county.append(name[1])
            pop2020.append(year_2020)
            pop2021.append(year_2021)
            pop2022.append(year_2022)
            pop2023.append(year_2023)
            pop2024.append(year_2024)

        self.eligible_voters_per_county_per_year = pd.DataFrame(
            data = {
                'state':state,
                'county':county,
                'pop2020':pop2020,
                'pop2021':pop2021,
                'pop2022':pop2022,
                'pop2023':pop2023,
                'pop2024':pop2024
                }
            )
     
    def create_voter_preds(self):

        self.election_results_raw = self.election_results_raw[self.election_results_raw['year'] > 2000]
        self.election_results_raw['party_grouped'] = self.election_results_raw['party']
        self.election_results_raw['party_grouped'] = np.where(
            ~self.election_results_raw['party_grouped'].isin(['DEMOCRAT', 'REPUBLICAN']),
            'OTHER',
            self.election_results_raw['party_grouped']
            )
        self.votes_per_party_per_county_per_year = self.election_results_raw.groupby(['state','state_po','county_name','county_fips','party_grouped','year'])['candidatevotes'].sum().reset_index(name='votes')
        self.votes_per_county_per_year = self.election_results_raw.groupby(['state','state_po','county_name','county_fips','year'])['candidatevotes'].sum().reset_index(name='votes')
        self.votes_per_party_per_state_per_year = self.election_results_raw.groupby(['state','party_grouped','year'])['candidatevotes'].sum().reset_index(name='votes')
        self.votes_per_party_per_year = self.election_results_raw.groupby(['party_grouped','year'])['candidatevotes'].sum().reset_index(name='votes')
        
        self.votes_per_party_per_county_per_year['county_party'] = self.votes_per_party_per_county_per_year['county_fips'].astype(int).astype(str) + "_" + self.votes_per_party_per_county_per_year['party_grouped'].astype(str)
        self.pivot_df = self.votes_per_party_per_county_per_year.pivot_table(index='county_party',columns=['year'],values='votes')
        self.pivot_df = self.pivot_df.fillna(0)
        self.pivot_df['delta_0'] = (self.pivot_df[2008]-self.pivot_df[2004])/self.pivot_df[2004]
        self.pivot_df['delta_1'] = (self.pivot_df[2012]-self.pivot_df[2008])/self.pivot_df[2008]
        self.pivot_df['delta_2'] = (self.pivot_df[2016]-self.pivot_df[2012])/self.pivot_df[2012]
        self.pivot_df['delta_3'] = (self.pivot_df[2020]-self.pivot_df[2016])/self.pivot_df[2016]
        self.pivot_df['mean_delta'] = np.nanmean(self.pivot_df[['delta_0','delta_1','delta_2','delta_3']].replace([np.nan, np.inf, -np.inf], np.nan), axis=1)
        
        self.pivot_df['2024_pred'] = (1+self.pivot_df['mean_delta'])*self.pivot_df[2020]
        self.covariance = np.cov(self.pivot_df[[2004, 2008, 2012, 2016,2020]])
        self.correlation = np.corrcoef(self.pivot_df[[2004, 2008, 2012, 2016,2020]])
        
        county_party_id = []
        preds = []
        sds = []
        votes_2020 = []
        votes_2004 = []
        votes_2008 = []
        votes_2012 = []
        votes_2016 = []
        for i in range(len(self.pivot_df)):
            county_party_id.append(self.pivot_df.index[i])
            votes_2004.append(self.pivot_df.iloc[i][2004])
            votes_2008.append(self.pivot_df.iloc[i][2008])
            votes_2012.append(self.pivot_df.iloc[i][2012])
            votes_2016.append(self.pivot_df.iloc[i][2016])
            votes_2020.append(self.pivot_df.iloc[i][2020])
            preds.append(self.pivot_df.iloc[i]['2024_pred'])
            sds.append(np.sqrt(self.covariance[i,i]))
            
        self.by_county_summary = pd.DataFrame(
            data = {
                'county_party_id':county_party_id,
                '2004':votes_2004,
                '2008':votes_2008,
                '2012':votes_2012,
                '2016':votes_2016,
                '2020':votes_2020,
                'preds':preds,
                'sds':sds
                }
            )
        self.by_county_summary['ratio'] = self.by_county_summary['sds']/self.by_county_summary['preds']
        self.by_county_summary[['fip','party']] = self.by_county_summary['county_party_id'].str.split("_",expand=True)
        self.by_county_summary.drop(columns=['county_party_id'],inplace=True)
        self.by_county_summary = self.by_county_summary[['fip','party'] + [col for col in self.by_county_summary.columns if col not in ['fip','party']]]
        
        self.n_sims = 1000
        self.simulations = np.zeros(shape=(len(self.by_county_summary),self.n_sims))
        
        self.n_rows_done = 0
        
        for county_fip_row_index in self.by_county_summary.index:
            self.county_fip_row = self.by_county_summary.loc[county_fip_row_index]
            self.pred = self.county_fip_row['preds']
            self.sd = self.county_fip_row['sds']

            if self.sd == 0 or pd.isna(self.sd) or self.pred==0 or pd.isna(self.pred):
                print("Decided not to generate ", self.county_fip_row["fip"])
            else:
                log_mean = np.log(self.pred**2 / np.sqrt(self.pred**2 + self.sd**2))
                log_std_dev = np.sqrt(np.log(1 + self.sd**2 / self.pred**2))
                
                self.samples_log_normal = np.random.lognormal(log_mean, log_std_dev, self.n_sims)
                self.simulations[self.n_rows_done,:] = self.samples_log_normal

            self.n_rows_done += 1
            if self.n_rows_done % 100 == 0:
                print("Dont with ", self.n_rows_done, "/", len(self.by_county_summary))

        column_names = [f'sim_{i}' for i in range(self.simulations.shape[1])]
        self.simulations_df = pd.DataFrame(data=self.simulations, columns=column_names)
        self.simulations_df = pd.concat([self.by_county_summary[['fip','party']],self.simulations_df],axis=1)

        print("Done with sims now printing...")

        self.simulations_df.to_csv("calibration_set_pre_correlated.csv")
        
        self.n_rows_to_sim = len(self.by_county_summary)
        
        print("Trying to find nearest PD ... ")
        self.cov_2 = nearestPD(self.covariance)
        
        self.lognormal_means = self.by_county_summary['preds']
        self.lognormal_sds = self.by_county_summary['sds']
        self.normal_means = np.log(self.lognormal_means**2/np.sqrt(self.lognormal_sds**2+self.lognormal_means**2))
        self.normal_sds = np.sqrt(np.log(1+(self.lognormal_sds**2/self.lognormal_means**2)))
        
        self.uncorrelated_normals = np.random.normal(self.n_rows_to_sim, self.n_sims)
        print("Simulating correlated normals ... ")
        self.correlated_normals = np.random.multivariate_normal(self.normal_means, self.cov_2, self.n_sims).T
        
        print("Fixing sims for sd .. ")
        for i in range(self.n_rows_to_sim):
            self.correlated_normals[i, :] = self.correlated_normals[i, :] * self.normal_sds[i] / np.std(self.correlated_normals[i, :])
            self.correlated_normals[i, :] = self.correlated_normals[i, :] + self.normal_means[i]
        self.lognormal_simulations = np.exp(self.correlated_normals)
        self.lognormal_simulations = np.nan_to_num(self.lognormal_simulations, nan=0.0000000001)
        
        self.corr_simulations_df = pd.DataFrame(data=self.lognormal_simulations, columns=column_names)
        self.corr_simulations_df = pd.concat([self.by_county_summary[['fip','party']],self.corr_simulations_df],axis=1)

        print("Done with correlated sims now printing...")

        self.corr_simulations_df.to_csv("calibration_set_post_correlated.csv")

    def get_counties_distances_us_city_data(self):
        self.us_cities_locations = pd.read_csv(self.path_to_us_cities)
        group_columns = ['state_id', 'state_name', 'county_fips','county_name']
        self.county_info_from_city_data = self.us_cities_locations.groupby(group_columns)[['population']].sum().reset_index()

if __name__ == "__main__":
    
    epred = ElectionPrediction()
    epred.get_counties_distances_us_city_data()
    epred.read_data()
    epred.predict_county_pop()
    epred.create_voter_preds()


    

    
    
    