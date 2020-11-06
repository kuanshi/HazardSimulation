# -*- coding: utf-8 -*-
#
# Copyright (c) 2018 Leland Stanford Junior University
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of the SimCenter Backend Applications
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# You should have received a copy of the BSD 3-Clause License along with
# this file. If not, see <http://www.opensource.org/licenses/>.
#
# Contributors:
# Kuanshi Zhong
#

import os
import subprocess
import json
import random
import numpy as np
import pandas as pd

def create_earthquake_scenarios(scenario_info, stations):

    # Number of scenarios
    source_num = scenario_info.get('Number', 1)
    # Directly defining earthquake ruptures
    if scenario_info['Generator'] == 'Prescription':
        scenario_data = dict()
        for i in range(source_num):
            try:
                source_type = scenario_info['EqRupture'][i]['Type']
                source_model = scenario_info['EqRupture'][i]['Model']
                source_ID = scenario_info['EqRupture'][i]['Source_ID']
                source_rupture = scenario_info['EqRupture'][i]['Rupture_ID']
                scenario_data.update({i: {
                    'Type': source_type,
                    'RuptureForecast': source_model,
                    'SourceIndex': source_ID,
                    'RuptureIndex': source_rupture
                }})
            except:
                print('Please check the defined scenario number.')
    # Searching earthquake ruptures that fulfill the request
    elif scenario_info['Generator'] == 'Selection':
        source_type = scenario_info['EqRupture']['Type']
        source_model = scenario_info['EqRupture']['Model']
        source_name = scenario_info['EqRupture'].get('Name',False)
        min_M = scenario_info['EqRupture'].get('min_Mag', 6.0)
        max_M = scenario_info['EqRupture'].get('max_Mag', 9.0)
        max_R = scenario_info['EqRupture'].get('max_Dist', 1000.0)
        # Collecting all possible earthquake scenarios
        lat = []
        lon = []
        for s in stations['Stations']:
            lat.append(s['Latitude'])
            lon.append(s['Longitude'])
        # Reference location
        lat = np.mean(lat)
        lon = np.mean(lon)
        tmp = dict()
        tmp['Site'] = {'Type': 'SingleLocation'}
        tmp['Site']['Location'] = {'Latitude': lat,
                                   'Longitude': lon}
        tmp['EqRupture'] = {'Type': source_type,
                            'RuptureForecast': source_model,
                            'ExportGeoJson': True,
                            'MaxDistance': max_R,
                            'MaxSources': 200}
        with open('tmp_input.json', 'w') as f:
            json.dump(tmp, f, indent = 2)
        # Calling EQHazard to search ruptures
        _ = subprocess.call(['java', '-jar', './lib/EQHazard.jar',
                             'tmp_input.json', 'tmp_output.json'])
        # Processing the selection results
        with open('tmp_output.json', 'r') as f:
            tmp = json.load(f)
        feat = tmp['features']
        tag = []
        for i, cur_f in enumerate(feat):
            if source_name and (source_name not in cur_f['properties']['Name']):
                continue
            if min_M > cur_f['properties']['Magnitude']:
                continue
            tag.append(i)
        # Abstracting desired ruptures
        s_tag = random.sample(tag, min(source_num, len(tag)))
        tmp['features'] = list(feat[i] for i in s_tag)
        scenario_data = dict()
        for i, rup in enumerate(tmp['features']):
            scenario_data.update({i: {
                'Type': source_type,
                'RuptureForecast': source_model,
                'SourceIndex': rup['properties']['Source'],
                'RuptureIndex': rup['properties']['Rupture']
            }})
        del tmp
        # Cleaning tmp outputs
        os.remove('tmp_input.json')
        os.remove('tmp_output.json')
    return scenario_data


def create_wind_scenarios(scenario_info, stations, data_dir):

    # Number of scenarios
    source_num = scenario_info.get('Number', 1)
    # Directly defining earthquake ruptures
    if scenario_info['Generator'] == 'Simulation':
        # Collecting site locations
        lat = []
        lon = []
        for s in stations['Stations']:
            lat.append(s['Latitude'])
            lon.append(s['Longitude'])
        # Save Stations.csv
        df = pd.DataFrame({
            'lat': lat,
            'lon': lon
        })
        df.to_csv(data_dir + 'Stations.csv', index = False, header = False)
        # Save Lat_w.csv
        lat_w = np.linspace(min(lat) - 0.5, max(lat) + 0.5, 100)
        df = pd.DataFrame({'lat_w': lat_w})
        df.to_csv(data_dir + 'Lat_w.csv', index = False, header = False)
        # Parsing Terrain info
        df = pd.read_csv(data_dir + scenario_info['Terrain']['Longitude'],
                         header = None, index_col = None)
        df.to_csv(data_dir + 'Long_wr.csv', header = False, index = False)
        df = pd.read_csv(data_dir + scenario_info['Terrain']['Latitude'],
                         header = None, index_col = None)
        df.to_csv(data_dir + 'Lat_wr.csv', header = False, index = False)
        df = pd.read_csv(data_dir + scenario_info['Terrain']['Size'],
                         header = None, index_col = None)
        df.to_csv(data_dir + 'wr_sizes.csv', header = False, index = False)
        df = pd.read_csv(data_dir + scenario_info['Terrain']['z0'],
                         header = None, index_col = None)
        df.to_csv(data_dir + 'z0r.csv', header = False, index = False)
        # Parsing storm properties
        param = []
        param.append(scenario_info['Storm']['Landfall']['Latitude'])
        param.append(scenario_info['Storm']['Landfall']['Longitude'])
        param.append(scenario_info['Storm']['LandingAngle'])
        param.append(scenario_info['Storm']['Pressure'])
        param.append(scenario_info['Storm']['Speed'])
        param.append(scenario_info['Storm']['Radius'])
        df = pd.DataFrame({'param': param})
        df.to_csv(data_dir + 'param.csv', index = False, header = False)
        df = pd.read_csv(data_dir + scenario_info['Storm']['Track'],
                         header = None, index_col = None)
        df.to_csv(data_dir + 'Track.csv', header = False, index = False)
        # Saving del_par.csv
        del_par = [0, 0, 0] # default
        df =pd.DataFrame({'del_par': del_par})
        df.to_csv(data_dir + 'del_par.csv', header = False, index = False)
        # Parsing resolution data
        delta_p = [1000., scenario_info['Resolution']['DivRad'], 1000000.]
        delta_p.extend([0., scenario_info['Resolution']['DivDeg'], 360.])
        delta_p.extend([scenario_info['MeasureHeight'], 10,
                       scenario_info['MeasureHeight']])
        df = pd.DataFrame({'delta_p': delta_p})
        df.to_csv(data_dir + 'delta_p.csv', header = False, index = False)
    else:
        print('Currently only supporting Simulation generator.')
