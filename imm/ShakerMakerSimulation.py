import numpy as np
import pandas as pd
import json
import os
from shakermaker import shakermaker
from shakermaker import crustmodel
from shakermaker.cm_library import AbellThesis
from shakermaker.cm_library import LOH
from shakermaker.cm_library import SOCal_LF
from shakermaker import station
from shakermaker import stationlist
from shakermaker import faultsource
from shakermaker import pointsource
from shakermaker.stf_extensions import brune

class ShakerMakerModel:

    def __init__(self, config, site, input_dir):
        """
        __init__: initializing the ShakerMaker model
        - config: configurations of the source and path
        - site: site locations
        """

        # getting simulation information
        self.dt = config.get('dt', 0.01)
        self.nfft = config.get('nfft', 4096)
        self.tb = config.get('tb', 1000)
        self.dk = config.get('dk', 0.3)
        self.kc = config.get('kc', 15.0)
        
        # getting source information
        self.STF = config['Source'].get('STF', 'Dirac')
        sourcefile = config['Source'].get('SourceFile', None)
        if sourcefile:
            source_flag = self.__config_source(os.path.abspath(os.path.join(input_dir, sourcefile)), self.STF)
            if source_flag:
                print('ShakerMakerSimulation: source configuration failed.')
        else:
            print('ShakeMakerSimulation: please provide a SourceFile.')

        # getting path information
        self.PathType = config['Path'].get('PathType', "SCEC_LOH_1")
        if (self.PathType == 'Abell'):
            self.num_split = config['Path'].get('SplitNumber')
        elif (self.PathType == 'UserDefined'):
            pathfile = config['Path'].get('PathFile', None)
            if pathfile:
                path_flag = self.__config_path(os.path.abspath(os.path.join(input_dir, pathfile)))
                if path_flag:
                    print('ShakeMakerSimulation: path configuration failed.')
                if (self.path_tkn[-1] != 0):
                    print('ShakerMakerSimulation: warning - the top layer should have a thickness of 0')
            else:
                print('ShakerMakerSimulation: please provide a PathFile.')

        # getting station/site information
        self.site_lat = [x['Latitude'] for x in site]
        self.site_lon = [x['Longitude'] for x in site]
        self.num_site = len(self.site_lat)
        self.site_alt = [0.0] * self.num_site
        if site[0].get('Altitude', None):
            self.site_alt = [x['Altitude'] for x in site]

        # setting up the coordinate space for the simulation
        self.__config_coord()


    def __config_source(self, sourcefile, stf):
        """
        __config_source: configuring the earthquake fault
        - sourcefile: directory to the file about the earthquake fault
        - stf: source time function type (Dirac, Brune, or Discrete)
        """

        df = pd.read_csv(sourcefile, header=0, index_col=None)
        self.num_source = len(df.index)
        self.source_lat = df['Latitude'].tolist()
        self.source_lon = df['Longitude'].tolist()
        self.source_dep = df['Depth'].tolist()
        self.source_stk = df['Strike'].tolist()
        self.source_dip = df['Dip'].tolist()
        self.source_rak = df['Rake'].tolist()
        self.source_tt0 = df['TriggerTime'].tolist()
        # readin additional parameters for Brune time
        if (stf == 'Brune'):
            try:
                self.source_slp = df['Slip'].tolist()
                self.source_chz = df['CornerFreq'].tolist()
                self.source_sdp = df['StressDrop'].tolist()
                self.source_sm0 = df['M0'].tolist()
                self.source_vsl = df['VsLocal'].tolist()
                self.source_smt = df['Smoothed'].tolist()
            except:
                print('ShakerMakerSimulation: please provide the following parameters for Brune STF:')
                print('Slip, CornerFreq, StressDrop, M0, VsLocal, Smoothed')
        
        # return
        return 0


    def __config_path(self, pathfile):
        """
        __config_path: configuring the crustal path
        - pathfiel: directory to the file about the crustal layers
        """

        df = pd.read_csv(pathfile, header=0, index_col=None)
        self.num_pathlayer = len(df.index)
        self.path_tkn = df['Thickness'].tolist()
        self.path_vp = df['Vp'].tolist()
        self.path_vs = df['Vs'].tolist()
        self.path_rho = df['Density'].tolist()
        self.path_qp = df['Qp'].tolist()
        self.path_qs = df['Qs'].tolist()

        # return
        return 0


    def __config_coord(self):
        """
        __config_coord: configuring the coordinate space
        """

        # earth raius (km)
        R = 6371.0
        # reference length
        dL = (1.0 / 180.0) * np.pi * R
        # coordinate origin
        Lat0 = np.mean(self.site_lat)
        Lon0 = np.mean(self.site_lon)
        Alt0 = 0.0
        # computing site coordinates
        self.site_X = np.array(self.site_lat - Lat0) * dL
        self.site_Y = np.array(self.site_lon - Lon0) * dL
        self.site_Z = Alt0 - np.array(self.site_alt)
        # computing source coordinates
        self.source_X = np.array(self.source_lat - Lat0) * dL
        self.source_Y = np.array(self.source_lon - Lon0) * dL


    def model_configuration(self):
        """
        model_configuration: setting up the simulaiton
        """

        # creating source
        self.source = []
        for i in range(self.num_source):
            tmp_loc = [self.source_X[i], self.source_Y[i], self.source_dep[i]]
            tmp_ang = [self.source_stk[i], self.source_dip[i], self.source_rak[i]]
            if (self.STF == 'Dirac'):
                tmp_s = pointsource.PointSource(tmp_loc, tmp_ang, tt=self.source_tt0[i])
            elif (self.STF == 'Brune'):
                tmp_stf = brune.Brune(slip=self.source_slp[i], f0=self.source_chz[i],
                    dsigma=self.source_sdp[i], M0=self.source_sm0[i],
                    Vs=self.source_vsl[i], smoothed=bool(self.source_smt[i]))
                tmp_s = pointsource.PointSource(tmp_loc, tmp_ang, tt=self.source_tt0[i], stf=tmp_stf)
            else:
                print('ShakerMakerSimulation: Discrete source time function is under development.')
            self.source.append(tmp_s)
        # collecting point sources
        self.fault = faultsource.FaultSource(self.source, metadata={"name":"fault"})
        
        # creating path
        if (self.PathType == 'SCEC_LOH_1'):
            self.crust = LOH.SCEC_LOH_1()
        elif (self.PathType == 'SCEC_LOH_3'):
            self.crust = LOH.SCEC_LOH_3()
        elif (self.PathType == 'Abell'):
            self.crust = AbellThesis.AbellThesis(split=self.num_split)
        elif (self.PathType == 'UserDefined'):
            self.crust = crustmodel.CrustModel(self.num_pathlayer)
            for i in range(self.num_pathlayer):
                self.crust.add_layer(self.path_tkn[i], self.path_vp[i], 
                    self.path_vs[i], self.path_rho[i], self.path_qp[i], self.path_qs[i])
        else:
            print('ShakerMakerSimulation: please use the right PathType.')
        
        # creating sites
        self.sitelist = []
        for i in range(self.num_site):
            tmp_s = station.Station([self.site_X[i], self.site_Y[i], self.site_Z[i]])
            self.sitelist.append(tmp_s)
        # collecting stations
        self.stations = stationlist.StationList(self.sitelist, {"name":"test"})


    def run_simulation(self, output_dir, output_format):
        """
        run_simulation: running the model
        """

        model = shakermaker.ShakerMaker(self.crust, self.fault, self.stations)
        model.run(dt=self.dt, nfft=self.nfft, dk=self.dk,  tb=self.tb, kc=self.kc)
        self.__parse_result(output_dir, output_format)


    def __parse_result(self, output_dir, output_format):
        """
        parse_result: parsing the simulation responses
        """

        # output EventGrid
        station_name = ['site'+str(j)+'.csv' for j in range(self.num_site)]
        df = pd.DataFrame({
            'Station': station_name,
            'Latitude': self.site_lat,
            'Longitude': self.site_lon            
        })
        df.to_csv(os.path.join(output_dir, 'EventGrid.csv'), index = False)

        for i in range(self.num_site):
            cur_site = self.sitelist[i]
            # get ground velocity responses
            vel_ns, vel_ew, vel_ud, t = cur_site.get_response()
            # get acceleration responses
            kms2g = 1000.0 / 9.81
            acc_ns = np.diff(vel_ns) / np.diff(t) * kms2g
            acc_ew = np.diff(vel_ew) / np.diff(t) * kms2g
            acc_ud = np.diff(vel_ud) / np.diff(t) * kms2g
            # output records
            rsn = 'SMN' + str(i)
            if (output_format == 'SimCenterEvent'):
                tmp = {
                    "name": rsn,
                    "dT": t[1] - t[0],
                    "data_x": acc_ew.tolist(),
                    "data_y": acc_ns.tolist(),
                    "data_z": acc_ud.tolist(),
                    "PGA_x": max(abs(acc_ew)),
                    "PGA_y": max(abs(acc_ns)),
                    "PGA_z": max(abs(acc_ud))
                }
                with open(os.path.join(output_dir, rsn+'.json'), 'w') as f:
                    json.dump(tmp, f, indent = 2)
            else:
                print('Currently only supporting SimCenterEvent')
            # output site#.csv
            df = pd.DataFrame({
                'TH_file': ['SMN'+str(i)],
                'factor': [1.0]
            })
            df.to_csv(os.path.join(output_dir, 'site'+str(i)+'.csv'), index = False)
        

