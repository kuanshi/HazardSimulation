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
from FetchOpenSHA import *


def create_earthquake_scenarios(scenario_info, stations):

    # Number of scenarios
    source_num = scenario_info.get('Number', 1)
    if source_num == 'All':
        # Large number to consider all sources in the ERF
        source_num = 10000000
    # Directly defining earthquake ruptures
    if scenario_info['Generator'] == 'Simulation':
        # TODO:
        print('Physics-based earthquake simulation is under development.')
        return 1
    # Searching earthquake ruptures that fulfill the request
    elif scenario_info['Generator'] == 'Selection':
        # Collecting all possible earthquake scenarios
        lat = []
        lon = []
        for s in stations['Stations']:
            lat.append(s['Latitude'])
            lon.append(s['Longitude'])
        # Reference location
        lat = np.mean(lat)
        lon = np.mean(lon)
        ref_station = [lat, lon]
        # Getting earthquake rupture forecast data
        source_type = scenario_info['EqRupture']['Type']
        if source_type == 'ERF':
            source_model = scenario_info['EqRupture']['Model']
            source_name = scenario_info['EqRupture'].get('Name', None)
            min_M = scenario_info['EqRupture'].get('min_Mag', 5.0)
            max_M = scenario_info['EqRupture'].get('max_Mag', 9.0)
            max_R = scenario_info['EqRupture'].get('max_Dist', 1000.0)
            eq_source = getERF(source_model, True)
            erf_data = export_to_json(eq_source, ref_station, outfile = None, \
                                      EqName = source_name, minMag = min_M, \
                                      maxMag = max_M, maxDistance = max_R, \
                                      maxSources = np.max([500, source_num]))
            # Parsing data
            feat = erf_data['features']
            tag = []
            for i, cur_f in enumerate(feat):
                if source_name and (source_name not in cur_f['properties']['Name']):
                    continue
                if min_M > cur_f['properties']['Magnitude']:
                    continue
                tag.append(i)
            # Abstracting desired ruptures
            s_tag = random.sample(tag, min(source_num, len(tag)))
            erf_data['features'] = list(feat[i] for i in s_tag)
            scenario_data = dict()
            for i, rup in enumerate(erf_data['features']):
                scenario_data.update({i: {
                    'Type': source_type,
                    'RuptureForecast': source_model,
                    'SourceIndex': rup['properties']['Source'],
                    'RuptureIndex': rup['properties']['Rupture']
                }})
            # Cleaning tmp outputs
            del erf_data
        elif source_type == 'PointSource':
            scenario_data = dict()
            try:
                magnitude = scenario_info['EqRupture']['Magnitude']
                location = scenario_info['EqRupture']['Location']
                average_rake = scenario_info['EqRupture']['AverageRake']
                average_dip = scenario_info['EqRupture']['AverageDip']
                scenario_data.update({0: {
                    'Type': source_type,
                    'Magnitude': magnitude,
                    'Location': location,
                    'AverageRake': average_rake,
                    'AverageDip': average_dip
                }})
            except:
                print('Please check point-source inputs.')
    # return
    return scenario_data


def create_wind_scenarios(scenario_info, event_info, stations, data_dir):

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
        # Station list
        station_list = {
            'Latitude': lat,
            'Longitude': lon
        }
        # Track data
        try:
            track_file = scenario_info['Storm'].get('Track')
            df = pd.read_csv(os.path.join(data_dir, track_file), header = None, index_col = None)
            track = {
                'Latitude': df.iloc[:, 0].values.tolist(),
                'Longitude': df.iloc[:, 1].values.tolist()
            }
        except:
            print('CreateScenario: no storm track provided or file format not accepted.')
        # Save Lat_w.csv
        track_simu_file = scenario_info['Storm'].get('TrackSimu', None)
        if track_simu_file:         
            df = pd.read_csv(os.path.join(data_dir, track_simu_file), header = None, index_col = None)
            track_simu = df.iloc[:, 0].values.tolist()
        else:
            track_simu = track['Latitude']
        # Reading Terrain info (if provided)
        terrain_file = scenario_info.get('Terrain', None)
        if terrain_file:
            with open(os.path.join(data_dir, terrain_file)) as f:
                terrain_data = json.load(f)
        else:
            terrain_data = []
        # Parsing storm properties
        param = []
        param.append(scenario_info['Storm']['Landfall']['Latitude'])
        param.append(scenario_info['Storm']['Landfall']['Longitude'])
        param.append(scenario_info['Storm']['LandingAngle'])
        param.append(scenario_info['Storm']['Pressure'])
        param.append(scenario_info['Storm']['Speed'])
        param.append(scenario_info['Storm']['Radius'])
        # Monte-Carlo
        #del_par = [0, 0, 0] # default
        # Parsing mesh configurations
        mesh_info = [1000., scenario_info['Mesh']['DivRad'], 1000000.]
        mesh_info.extend([0., scenario_info['Mesh']['DivDeg'], 360.])
        # Wind speed measuring height
        measure_height = event_info['IntensityMeasure']['MeasureHeight']
        # Saving results
        scenario_data = dict()
        for i in range(source_num):
            scenario_data.update({i: {
                'Type': 'Wind',
                'CycloneParam': param,
                'StormTrack': track,
                'StormMesh': mesh_info,
                'Terrain': terrain_data,
                'TrackSimu': track_simu,
                'StationList': station_list,
                'MeasureHeight': measure_height
            }})
        # return
        return scenario_data
    else:
        print('Currently only supporting Simulation generator.')
