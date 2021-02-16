/*********************************************************************************
**
** Copyright (c) 2021 University of California, Berkeley
** Copyright (c) 2021 Leland Stanford Junior University
**
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 1. Redistributions of source code must retain the above copyright notice, this
** list of conditions and the following disclaimer.
**
** 2. Redistributions in binary form must reproduce the above copyright notice, this
** list of conditions and the following disclaimer in the documentation and/or other
** materials provided with the distribution.
**
** 3. Neither the name of the copyright holder nor the names of its contributors may
** be used to endorse or promote products derived from this software without specific
** prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
** EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
** SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
** INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
** BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
** CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
** IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
** SUCH DAMAGE.
**
***********************************************************************************/

// Contributors:
// Ajay B Harish, Post-Doc @ SimCenter, UC Berkeley
// Kuanshi Zhong, SimCenter, Stanford
// Dr. Frank McKenna, CTO of SimCenter, UC Berkeley
// Prof. Sanjay Govindjee, Director of SimCenter, UC Berkeley

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <string.h>
#include <vector>
#include "include/Eigen/Dense"
#include "include/rapidjson/document.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/ostreamwrapper.h"
#include "include/rapidjson/prettywriter.h"
#include "WindFieldModel.h"


int WindFieldModel::ConfigSimu(std::string config_file, std::string stn_file, 
    std::string trk_file, std::string latw_file)
{
    rapidjson::Document config_doc;
    std::ifstream ifs(config_file.c_str());
    rapidjson::IStreamWrapper isw(ifs);
    config_doc.ParseStream(isw);

    rapidjson::Value::Object scen = config_doc["Scenario"].GetObject();
    rapidjson::Value::Object strm = scen["Storm"].GetObject();
    rapidjson::Value::Object mesh = scen["Mesh"].GetObject();
    rapidjson::Value::Object evnt = config_doc["Event"].GetObject();

    // initializing modeling parameters
    this->param = Eigen::ArrayXd::Zero(6);
    this->delta_p = Eigen::ArrayXd::Zero(9);
    this->del_par = Eigen::ArrayXd::Zero(3);
    this->dP = 0.0;
    this->dV = 0.0;
    this->dR = 0.0;

    // landfall data
    rapidjson::Value::Object tmpObj = strm["Landfall"].GetObject();
    this->param(0) = tmpObj["Latitude"].GetDouble();
    this->param(1) = tmpObj["Longitude"].GetDouble();
    // landing angle
    this->param(2) = strm["LandingAngle"].GetDouble();
    // central pressure difference
    this->param(3) = strm["Pressure"].GetDouble();
    // forward translation speed
    this->param(4) = strm["Speed"].GetDouble();
    // radius of the maximum wind
    this->param(5) = strm["Radius"].GetDouble();
    // cyclone radius meshing
    this->delta_p(0) = 1000.0;
    this->delta_p(1) = mesh["DivRad"].GetDouble();
    this->delta_p(2) = 1000000.0;
    // cyclone angle meshing
    this->delta_p(3) = 0.0;
    this->delta_p(4) = mesh["DivDeg"].GetDouble();
    this->delta_p(5) = 360.0;
    // wind speed measuring height
    tmpObj = evnt["IntensityMeasure"].GetObject();
    rapidjson::Value::Array tmpAry = tmpObj["MeasureHeight"].GetArray();
    this->delta_p(6) = tmpAry[0].GetDouble();
    this->delta_p(7) = tmpAry[1].GetDouble();
    this->delta_p(8) = tmpAry[2].GetDouble();

    // stations
    std::ifstream stn_data(stn_file);
    std::string line, h;
    getline(stn_data, line);
    std::stringstream s(line);
    std::vector<std::string> stn_header;
    int tag_Lat = 0;
    int tag_Long = 0;
    int tmp_tag = 0;
    while (getline(s, h, ','))
    {
        stn_header.push_back(h);
        if (h == "Latitude") {
            tag_Lat = tmp_tag;
        } else if (h == "Longitude")
        {
            tag_Long = tmp_tag;
        }
        tmp_tag = tmp_tag + 1;
    }
    std::vector<double> tmp_Lat;
    std::vector<double> tmp_Long;
    while (getline(stn_data, line))
    {
        std::stringstream r(line);
        int tmp_tag = 0;
        while (getline(r, h, ','))
        {
            if (tmp_tag == tag_Lat)
            {
                tmp_Lat.push_back(std::stod(h));
            } else if (tmp_tag == tag_Long)
            {
                tmp_Long.push_back(std::stod(h));
            }
            tmp_tag = tmp_tag + 1;
        }
    }
    this->num_station = tmp_Lat.size();
    this->Lat_wout = Eigen::ArrayXd::Zero(this->num_station);
    this->Long_wout = Lat_wout;
    for (int i = 0; i < this->num_station; i++)
    {
        this->Lat_wout(i) = tmp_Lat[i];
        this->Long_wout(i) = tmp_Long[i];
    }

    // track
    std::ifstream trk_data(trk_file);
    std::vector<double> tmp_trkLat;
    std::vector<double> tmp_trkLong;
    while (getline(trk_data, line))
    {
        std::stringstream r(line);
        int tmp_tag = 0;
        while (getline(r, h, ','))
        {
            if (tmp_tag == 0)
            {
                tmp_trkLat.push_back(std::stod(h));
            } else if (tmp_tag == 1)
            {
                tmp_trkLong.push_back(std::stod(h));
            }
            tmp_tag = tmp_tag + 1;
        }
    }
    int num_trkp = tmp_trkLat.size();
    this->Lat_track = Eigen::ArrayXd::Zero(num_trkp);
    this->Long_track = this->Lat_track;
    for (int i = 0; i < num_trkp; i++)
    {
        this->Lat_track(i) = tmp_trkLat[i];
        this->Long_track(i) = tmp_trkLong[i];
    }

    // refined track latitude
    std::ifstream latw_data(latw_file);
    std::vector<double> tmp_latw;
    while (getline(latw_data, line))
    {
        tmp_latw.push_back(std::stod(line));
    }
    int num_latw = tmp_latw.size();
    this->Lat_w = Eigen::ArrayXd::Zero(num_latw);
    this->Long_w = this->Lat_w;
    for (int i = 0; i < num_latw; i++)
    {
        this->Lat_w(i) = tmp_latw[i];
    }
    // Calculate longitude of refence values (Long_w)
    // This is done by interpolating
    int value;
    for (int ii = 0; ii < num_latw; ii++)
    {
        value = 0;
        for (int jj = 0; jj < this->Lat_track.size() - 1; jj++)
        {
            if ((this->Lat_w(ii) >= this->Lat_track(jj)) && (this->Lat_w(ii) < this->Lat_track(jj + 1)))
            {
                value = jj;
                break;
            }
        }
        if (value > 0)
        {
            this->Long_w(ii) = abs(this->Long_track(value) + ((this->Lat_w(ii) - this->Lat_track(value)) 
                * (this->Long_track(value + 1) - this->Long_track(value))) / (this->Lat_track(value + 1) 
                - this->Lat_track(value)));
        } else {
            if (this->Lat_w(ii) < this->Lat_track(0))
            {
                this->Long_w(ii) = abs(this->Long_track(0));
            } else {
                this->Long_w(ii) = abs(this->Long_track(num_latw - 1));
            }
        }
    }

    return 0;
}


int WindFieldModel::PertubPath(std::string dpath_file)
{
    rapidjson::Document dpath_doc;
    std::ifstream ifs(dpath_file.c_str());
    rapidjson::IStreamWrapper isw(ifs);
    dpath_doc.ParseStream(isw);

    // \delta Anlge
    this->del_par(2) = dpath_doc["dAngle"].GetDouble();
    // \delta Latitude
    this->del_par(0) = dpath_doc["dLatitude"].GetDouble();
    // \delta Longitude
    this->del_par(1) = dpath_doc["dLongitude"].GetDouble();
    // \delta Pressure
    this->dP = dpath_doc["dP"].GetDouble();
    // \delta Speed
    this->dV = dpath_doc["dV"].GetDouble();
    // \delta Radius
    this->dR = dpath_doc["dR"].GetDouble();

    return 0;
}


int WindFieldModel::DefineTern(std::string refz0_file)
{
    rapidjson::Document z0_doc;
    std::ifstream ifs(refz0_file.c_str());
    rapidjson::IStreamWrapper isw(ifs);
    z0_doc.ParseStream(isw);

    rapidjson::Value::Array feat = z0_doc["features"].GetArray();
    this->num_region = feat.Size();
    // initialzing z0
    this->z0r = Eigen::ArrayXd::Zero(this->num_region);
    this->Wr_sizes = Eigen::ArrayXd::Zero(this->num_region);
    // getting z0 values and sizes
    for (int i = 0; i < this->num_region; i++)
    {
        rapidjson::Value::Object tmpObj = feat[i]["properties"].GetObject();
        this->z0r(i) = tmpObj["z0"].GetDouble();
        rapidjson::Value::Object tmpGeo = feat[i]["geometry"].GetObject();
        rapidjson::Value::Array tmpAry = tmpGeo["coordinates"].GetArray();
        this->Wr_sizes(i) = tmpAry.Size();
    }
    // getting polygons
    this->Lat_wr = Eigen::MatrixXd::Zero(this->Wr_sizes.maxCoeff(), this->num_region);
    this->Long_wr = Lat_wr;
    for (int i = 0; i < this->num_region; i++)
    {
        int jmax = this->Wr_sizes(i);
        rapidjson::Value::Object tmpObj = feat[i]["geometry"].GetObject();
        rapidjson::Value::Array tmpAry = tmpObj["coordinates"].GetArray();
        for (int j = 0; j < jmax; j++)
        {
            rapidjson::Value::Array curAry = tmpAry[j].GetArray();
            this->Lat_wr(j, i) = curAry[0].GetDouble();
            this->Long_wr(j, i) = curAry[1].GetDouble();
        }
    }

    return 0;
}


int WindFieldModel::ComputeStationZ0(std::string dirOutput)
{
    // Calculating the surface roughness z0 for staitons
    // Mapping individual stations to the reference z0r values of polygons
    Eigen::ArrayXd z0_station = Eigen::ArrayXd::Zero(this->Lat_wout.size());
    Eigen::ArrayXd PolyX = this->Lat_wr.block(0, this->Wr_sizes.size() - 1, this->Wr_sizes(Wr_sizes.size() - 1), 1);
    Eigen::ArrayXd PolyY = this->Long_wr.block(0, this->Wr_sizes.size() - 1, this->Wr_sizes(Wr_sizes.size() - 1), 1);
    for (int ii = 0; ii < this->Lat_wout.size(); ii++)
    {
        double stationX = this->Lat_wout(ii);
        double stationY = this->Long_wout(ii);
        z0_station(ii) = 0.01;
        int aux = 0;
        int oi = 0;
        while ((aux == 0) && (oi < this->Wr_sizes.size() - 1))
        {
            Eigen::ArrayXd PolyX = this->Lat_wr.block(0, oi, this->Wr_sizes(oi), 1);
            Eigen::ArrayXd PolyY = this->Long_wr.block(0, oi, this->Wr_sizes(oi), 1);
            aux = WindFieldModel::inpolygon(PolyX, PolyY, this->Wr_sizes(oi), stationX, stationY);
            oi++;
        }
        if (aux == 1)
        {
            z0_station(ii) = z0r(oi);
        }        
    }
    // Writing out z0_station to a result file
    std::ofstream z0File(dirOutput + "/StationZ0.csv");
    z0File << std::setprecision(4) << z0_station << "\n";
    z0File.close();

    return 0;
}


int WindFieldModel::SimulateWind(std::string dirOutput)
{
    /* debug
    std::vector<double> dbg_Lat_w, dbg_Long_w;
    std::vector<double> dbg_Lat_wout, dbg_Long_wout;
    */

    // central pressure difference
    double del_p = (this->param(3) + this->dP) * 100.0;
    // Holland B parameter
    double B = 1.38 + 0.00184 * (del_p / 100.0) - 0.00309 * (this->param(5) + this->dR);
    // translation speed
    double c = ((this->param(4) + this->dV) * 1000.0) / 3600.0;
    // radius in meter
    double r_m = (this->param(5) + this->dR) * 1000.0;
    // air density
    double rho = this->AIR_DENSITY;
    // eddy viscocity
    double k_m = this->EDDY_VISCOCITY;
    // Earth raius
    double R = this->R;
    // eps
    double eps = this->EPS;
    // ra 
    double ra = 180.0 / this->PI;

    // Calculating heading
    Eigen::ArrayXd beta_c = Eigen::ArrayXd::Zero(this->Lat_w.size());
    Eigen::ArrayXd temp00 = Eigen::ArrayXd::Zero(this->Lat_w.size());
    for (int ii = 0; ii < this->Lat_w.size() - 1; ii++)
    {
        double Delta = abs(this->Long_w(ii + 1)) - abs(this->Long_w(ii)) + eps * eps;
        beta_c(ii) = -this->del_par(2) + 90 + ra * atan2(sin(Delta / ra) * cos((this->Lat_w(ii + 1)) / ra), 
            cos((this->Lat_w(ii)) / ra) * sin((this->Lat_w(ii + 1)) / ra) - sin((this->Lat_w(ii)) / ra) 
            * cos((this->Lat_w(ii + 1)) / ra) * cos(Delta / ra));

        if (beta_c(ii) < 0)
        {
            beta_c(ii) = beta_c(ii) + 360.0;
        }
    }
    beta_c(this->Lat_w.size() - 1) = beta_c(this->Lat_w.size() - 2);

    /* debug
    std::vector<double> dbg_beta_c;
    for (int tag = 0; tag < beta_c.size(); tag++)
    {
        dbg_beta_c.push_back(beta_c(tag));
    }
    */

    // Calculate r - theta - zp
    int size = 1 + (this->delta_p(2) - this->delta_p(0)) / this->delta_p(1);
    Eigen::ArrayXd r = Eigen::ArrayXd::Zero(size);
    for (int ii = 0; ii < size; ii++)
    {
        r(ii) = this->delta_p(0) + ii * this->delta_p(1);
    }
    size = 1 + (this->delta_p(5) - this->delta_p(3)) / this->delta_p(4);
    Eigen::ArrayXd theta = Eigen::ArrayXd::Zero(size);
    for (int ii = 0; ii < size; ii++)
    {
        theta(ii) = this->delta_p(3) + ii * this->delta_p(4);
    }
    size = 1 + (this->delta_p(8) - this->delta_p(6)) / this->delta_p(7);
    Eigen::ArrayXd zp = Eigen::ArrayXd::Zero(size);
    for (int ii = 0; ii < size; ii++)
    {
        zp(ii) = this->delta_p(6) + ii * this->delta_p(7);
    }

    // Writing out zp to a result file
    std::ofstream zpFile(dirOutput + "/MeasureHeight.csv");
    zpFile << std::setprecision(4) << zp << "\n";
    zpFile.close();

    // Initialize matrices for later usage
    Eigen::ArrayXd u_max_tmp = Eigen::ArrayXd::Zero(this->Lat_wout.size());
    Eigen::MatrixXd u_max = Eigen::MatrixXd::Zero(this->Lat_wout.size(), zp.size());
    Eigen::MatrixXd u_tmp = Eigen::MatrixXd::Zero(zp.size(), r.size());
    Eigen::MatrixXd v_tmp = u_tmp;
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(theta.size() * r.size(), zp.size());
    Eigen::MatrixXd v = u;
    Eigen::MatrixXd v_g1 = Eigen::MatrixXd::Zero(theta.size(), r.size());
    Eigen::ArrayXd z0 = Eigen::ArrayXd::Zero(r.size());
    Eigen::VectorXd Wind_speed = Eigen::VectorXd::Zero(this->Lat_wout.size());
    Eigen::MatrixXd jvar2 = Eigen::MatrixXd::Zero(this->num_station, this->Lat_w.size());
    Eigen::MatrixXd kvar2 = Eigen::MatrixXd::Zero(this->num_station, this->Lat_w.size());

    // Perform the loop for different reference points
    // i.e., different instance locations in the entire storm path
    double Lat, Long;
    for (int ii = 0; ii < this->Lat_w.size(); ii++)
    {
        std::cout << "Computing wind speed field - storm track location #" << ii + 1 << "/" << this->Lat_w.size() << std::endl;
        // location and deading, including perturbation for MC case
        Lat = this->Lat_w(ii); // + 0*del_par(0);
        Long = abs(this->Long_w(ii)) - 0.3 * this->del_par(1);
        double beta = beta_c(ii);
        // Coriolis
        double omega = 0.00007292;
        double f = 2.0 * omega * sin((Lat * this->PI) / 180.0);

        Eigen::ArrayXd dd = acos(cos(this->Lat_wout / ra) * cos(Lat / ra) * cos((abs(this->Long_wout) - Long) / ra) 
            + sin(this->Lat_wout / ra) * sin(Lat / ra));
        double temp100 = 6371.0 * 180.0 / this->PI / ra * 1000.0;
        dd = dd * temp100;
        Eigen::ArrayXd Delta2 = abs(this->Long_wout) - Long + pow(eps, 2);

        /* debug
        std::vector<double> dbg_dd, dbg_Delta2;
        for (int tag = 0; tag < dd.size(); tag++)
        {
            dbg_dd.push_back(dd(tag));
            dbg_Delta2.push_back(Delta2(tag));
        }
        */

        // Bearing angle in degrees
        Eigen::ArrayXd bearing = Eigen::ArrayXd::Zero(this->Lat_wout.size());
        Eigen::ArrayXd jvar = Eigen::ArrayXd::Zero(this->Lat_wout.size());
        Eigen::ArrayXd kvar = Eigen::ArrayXd::Zero(this->Lat_wout.size());
        for (int mm = 0; mm < this->Lat_wout.size(); mm++)
        {
            bearing(mm) = 90.0 + ra * atan2(sin(Delta2(mm) / ra) * cos(this->Lat_wout(mm) / ra),
                cos(Lat / ra) * sin(this->Lat_wout(mm) / ra) - sin(Lat / ra) * cos(this->Lat_wout(mm) / ra) * cos(Delta2(mm) / ra));
            if (bearing(mm) < 0)
            {
                bearing(mm) = bearing(mm) + 360.0;
            }
            jvar(mm) = trunc(bearing(mm) / this->delta_p(4)) + 1;
            kvar(mm) = min(trunc(dd(mm) / this->delta_p(1)), r.size() - 1) + 1;
        }

        jvar2.col(ii) = jvar;
        kvar2.col(ii) = kvar;

        /* debug
        std::vector<double> dbg_Long_w, dbg_Lat_w;
        for (int tag = 0; tag < this->Long_w.size(); tag ++)
        {
            dbg_Long_w.push_back(this->Long_w(tag));
            dbg_Lat_w.push_back(this->Lat_w(tag));
        }
        std::vector<double> dbg_bearing, dbg_jvar, dbg_kvar;
        for (int tag = 0; tag < bearing.size(); tag++)
        {
            dbg_bearing.push_back(bearing(tag));
            dbg_jvar.push_back(jvar(tag));
            dbg_kvar.push_back(kvar(tag));
        }
        */

        // Calculate wind field model, for loop for different polar coordinates
        // theta in for loop and r element-wise
        for (int jj = 0; jj < theta.size(); jj++) //theta.size()
        {

            double Ctheta = -c * sin((theta(jj) - beta) / ra);
            double THETA;
            if ((theta(jj) >= 0.0) && (theta(jj) <= 90.0))
                THETA = 90.0 - theta(jj);
            else
                THETA = 450.0 - theta(jj);

            Eigen::ArrayXd Lat_t = Eigen::ArrayXd::Zero(r.size());
            Eigen::ArrayXd Long_t = Eigen::ArrayXd::Zero(r.size());
            Lat_t = ra * asin(sin(Lat / ra) * cos(r / R) + cos(Lat / ra) * sin(r / R) * cos(THETA / ra));
            for (int kk = 0; kk < r.size(); kk++)
            {
                Long_t(kk) = Long + ra * atan2(sin(THETA / ra) * sin(r(kk) / R) * cos(Lat / ra), cos(r(kk) / R) 
                    - sin(Lat / ra) * sin(Lat_t(kk)));
                Long_t(kk) = wrapTo360(Long_t(kk) + 540.0) - 180.0;
            }

            // Check if point is inside the polygon
            for (int kk = 0; kk < r.size(); kk++)
            {
                // Create arrays
                Eigen::ArrayXd PolyX = this->Lat_wr.block(0, Wr_sizes.size() - 1, this->Wr_sizes(this->Wr_sizes.size() - 1), 1);
                Eigen::ArrayXd PolyY = this->Long_wr.block(0, Wr_sizes.size() - 1, this->Wr_sizes(this->Wr_sizes.size() - 1), 1);
                z0(kk) = 0.01;
                int aux = 0;
                int oi = 0;
                while ((aux == 0) && (oi < this->Wr_sizes.size() - 2))
                {
                    Eigen::ArrayXd PolyX = this->Lat_wr.block(0, oi, this->Wr_sizes(oi), 1);
                    Eigen::ArrayXd PolyY = this->Long_wr.block(0, oi, this->Wr_sizes(oi), 1);
                    aux = inpolygon(PolyX, PolyY, this->Wr_sizes(oi), Lat_t(kk), Long_t(kk));
                    oi++;
                }
                if (aux == 1)
                {
                    z0(kk) = z0r(oi);
                }                    
            }

            double z10 = 10.0;
            double A = 11.4;
            Eigen::ArrayXd h = A * (z0.array().pow(0.86));
            Eigen::ArrayXd d = 0.75 * h;
            double kappa = 0.40;
            Eigen::ArrayXd Cd = (z10 + h - d) * z0.inverse();
            Cd = ((Cd.log()).pow(2)).inverse();
            Cd = kappa * kappa * Cd;
            Eigen::ArrayXd temp = r_m * (r.inverse());
            Eigen::ArrayXd temp5 = -(temp.pow(B));
            Eigen::ArrayXd temp6 = temp5.exp();
            Eigen::ArrayXd der_p = B * (pow(r_m, B)) * del_p * (r.pow(-(B + 1))) * temp6;
            Eigen::ArrayXd der_p_2 = (-(B + 1.0) * r.inverse() + B * (pow(r_m, B)) * r.pow(-(B + 1.0))) * der_p;
            v_g1.row(jj) = (0.5) * (Ctheta - f * r) + (((0.5) * (Ctheta - f * r)).array().pow(2) 
                + (1.0 / rho) * (r.cwiseProduct(der_p))).array().pow(0.5);

            Eigen::ArrayXd der_v_g1_r = -0.5 * f + (0.5) * ((((0.5) * (Ctheta - f * r)).array().pow(2) 
                + (1 / rho) * (r * der_p)).array().pow(-0.5)) * (-(Ctheta - f * r) * f * 0.5 
                + (1.0 / rho) * der_p + (1.0 / rho) * r * der_p_2);
            Eigen::ArrayXd der_v_g1_theta = (0.25 * c) * (cos((theta(jj) - beta) / ra)) * (-Ctheta + f * r) 
                * ((0.5 * (-Ctheta + f * r)).array().pow(2) + (1 / rho) * r * der_p).array().pow(-0.5) 
                - (0.5 * c) * cos((theta(jj) - beta) / ra);
            Eigen::ArrayXd BB = ((2.0 * k_m * r).inverse()) * der_v_g1_theta;
            Eigen::ArrayXd Eta = (((0.5) * (Ctheta - f * r)).pow(2) + (1.0 / rho) * (r * der_p)).pow(0.5);
            temp = v_g1.row(jj);
            Eigen::ArrayXd ALPHA = (1.0 / (2.0 * k_m)) * (f + 2.0 * (temp * r.inverse()));
            Eigen::ArrayXd BETA = (1.0 / (2.0 * k_m)) * (f + temp * r.inverse() + der_v_g1_r);
            Eigen::ArrayXd GAMMA = (1.0 / (2.0 * k_m)) * (temp * r.inverse());
            Eigen::ArrayXcd ALPHA_C = Eigen::ArrayXcd::Zero(ALPHA.size());
            Eigen::ArrayXcd BETA_C = Eigen::ArrayXcd::Zero(BETA.size());
            Eigen::ArrayXcd GAMMA_C = Eigen::ArrayXcd::Zero(GAMMA.size());
            Eigen::ArrayXcd BB_C = Eigen::ArrayXcd::Zero(BB.size());
            ALPHA_C.real() = ALPHA;
            BETA_C.real() = BETA;
            GAMMA_C.real() = GAMMA;
            BB_C.real() = BB;
            Eigen::ArrayXd XXX = -(ALPHA * BETA).pow(0.25);
            Eigen::ArrayXd YYY = -(ALPHA * BETA).pow(0.25);

            /* debug
            std::vector<double> dbg_ALPHA, dbg_BETA, dbg_GAMMA;
            for (int tag = 0; tag < ALPHA.size(); tag++)
            {
                dbg_ALPHA.push_back(ALPHA(tag));
                dbg_BETA.push_back(BETA(tag));
                dbg_GAMMA.push_back(GAMMA(tag));
            }
            */

            Eigen::ArrayXcd PP_zero = Eigen::ArrayXcd::Zero(XXX.size());
            Eigen::ArrayXcd PP_one = Eigen::ArrayXcd::Zero(XXX.size());
            Eigen::ArrayXcd PP_minus_one = Eigen::ArrayXcd::Zero(XXX.size());
            PP_zero.real() = XXX;
            PP_zero.imag() = YYY;
            std::complex<double> sc(1.0, 1.0);
            PP_one = -sc * (GAMMA_C + (ALPHA_C * BETA_C).pow(0.5) - BB_C).pow(0.5);
            PP_minus_one = -sc * (-GAMMA_C + (ALPHA_C * BETA_C).pow(0.5) - BB_C).pow(0.5);

            Eigen::ArrayXcd X1 = PP_zero + (f / k_m) * r * Cd - (2.0 / k_m) * Eta * Cd -
                                 (c * c * Cd.pow(2)) * ((4.0 * k_m * k_m * (PP_one - PP_minus_one.conjugate())).inverse()) +
                                 (c * c * Cd.pow(2)) * ((4.0 * k_m * k_m * (PP_one.conjugate() - PP_minus_one)).inverse());

            Eigen::ArrayXcd X2 = -PP_zero.conjugate() - (f / k_m) * r * Cd + (2.0 / k_m) * Eta * Cd -
                                 (c * c * Cd.pow(2)) * ((4.0 * k_m * k_m * (PP_one - PP_minus_one.conjugate())).inverse()) +
                                 (c * c * Cd.pow(2)) * ((4.0 * k_m * k_m * (PP_one.conjugate() - PP_minus_one)).inverse());

            Eigen::ArrayXcd X3 = Eigen::ArrayXcd::Zero(r.size());
            X3.imag() = -(2.0 / k_m) * Cd * (Eta - (f / 2.0) * r).pow(2);

            Eigen::ArrayXcd X4 = -(-PP_zero - (f / (2.0 * k_m)) * r * Cd + (Eta / k_m) * Cd) *
                                 (-PP_zero.conjugate() - (f / (2.0 * k_m)) * r * Cd + (Eta / k_m) * Cd).inverse();

            /* debug
            std::vector<std::complex<double>> dbg_X1, dbg_X2, dbg_X3, dbg_X4;
            for (int tag = 0; tag < X1.size(); tag++)
            {
                dbg_X1.push_back(X1(tag));
                dbg_X2.push_back(X2(tag));
                dbg_X3.push_back(X3(tag));
                dbg_X4.push_back(X4(tag));
            }
            */

            Eigen::ArrayXcd A_zero = -X3 * ((X1 + (X2 * X4)).inverse());
            std::complex<double> temvar1(0, 1);
            std::complex<double> temvar2(0, -beta);
            Eigen::ArrayXcd A_one = Eigen::ArrayXcd::Zero(Cd.size());
            A_one = temvar1 * c * Cd * exp(temvar2) * (A_zero + A_zero.conjugate()) * (4.0 * k_m * (PP_one - PP_minus_one.conjugate())).inverse();
            Eigen::ArrayXcd A_minus_one = -A_one.conjugate();

            // Loop for different heights currently only done for z=0
            for (int kk = 0; kk < zp.size(); kk++)
            {
                //int kk = 0;
                Eigen::ArrayXd u_zero = sqrt(ALPHA * inverse(BETA)) * (A_zero * exp(PP_zero * zp(kk))).real();
                Eigen::ArrayXd v_zero = (A_zero * exp(PP_zero * zp(kk))).imag();

                Eigen::ArrayXd u_one = sqrt(ALPHA * inverse(BETA)) * (A_one * exp(PP_one * zp(kk) + temvar1 * theta(jj) * (PI / 180.0))).real();
                Eigen::ArrayXd v_one = (A_one * exp(PP_one * zp(kk) + temvar1 * theta(jj) * (PI / 180.0))).imag();
                Eigen::ArrayXd u_minus_one = sqrt(ALPHA * inverse(BETA)) * (A_minus_one * exp(PP_minus_one * zp(kk) 
                    - temvar1 * theta(jj) * (PI / 180.0))).real();
                Eigen::ArrayXd v_minus_one = (A_minus_one * exp(PP_minus_one * zp(kk) - temvar1 * theta(jj) * (PI / 180.0))).imag();

                /* debug
                std::vector<double> dbg_u_zero, dbg_u_one, dbg_u_minus_one, dbg_v_zero, dbg_v_one, dbg_v_minus_one;
                for (int tag = 0; tag < X1.size(); tag++)
                {
                    dbg_u_zero.push_back(u_zero(tag));
                    dbg_u_one.push_back(u_one(tag));
                    dbg_u_minus_one.push_back(u_minus_one(tag));
                    dbg_v_zero.push_back(v_zero(tag));
                    dbg_v_one.push_back(v_one(tag));
                    dbg_v_minus_one.push_back(v_minus_one(tag));
                }
                */

                u_tmp.row(kk) = u_zero + u_one + u_minus_one;
                v_tmp.row(kk) = v_zero + v_one + v_minus_one;
            }

            /* debug
            std::vector<double> dbg_u_tmp, dbg_v_tmp;
            for (int tag = 0; tag < u_tmp.size(); tag++)
            {
                dbg_u_tmp.push_back(u_tmp(tag));
                dbg_v_tmp.push_back(v_tmp(tag));
            }
            */

            // Saving u_tmp and v_tmp to u and v
            for (int pp = 0; pp < r.size(); pp++)
            {
                u.row(jj * r.size() + pp) = u_tmp.col(pp).array();
                v.row(jj * r.size() + pp) = v_tmp.col(pp).array();
            }
        }

        // Looping over each height in zp
        for (int pp = 0; pp < zp.size(); pp++)
        {
            //printarray("v.csv", v.col(pp).array());
            Eigen::Map<Eigen::MatrixXd> tmp_v(v.col(pp).data(), r.size(), theta.size());
            //printmatrix("tmp_v.csv", tmp_v);
            //std::cout << "tmp_v " << tmp_v.rows() << "," << tmp_v.cols() << std::endl;
            Eigen::Map<Eigen::MatrixXd> tmp_u(u.col(pp).data(), r.size(), theta.size());
            //std::cout << "tmp_u " << tmp_u.rows() << "," << tmp_u.cols() << std::endl;
            Eigen::MatrixXd v_1 = v_g1.transpose() + tmp_v;
            //std::cout << "v_1 " << v_1.rows() << "," << v_1.cols() << std::endl;
            Eigen::MatrixXd Uvel = (v_1.array().pow(2) + tmp_u.array().pow(2)).array().sqrt();
            //std::cout << "Uvel " << Uvel.rows() << "," << Uvel.cols() << std::endl;
            for (int mm = 0; mm < Lat_wout.size(); mm++)
            {
                Wind_speed(mm) = Uvel(kvar(mm) - 1, jvar(mm) - 1);
                u_max_tmp(mm) = std::max(u_max(mm), Wind_speed(mm));
            }

            /* debug
            std::vector<double> dbg_u_max_tmp;
            for (int tag = 0; tag < u_max_tmp.size(); tag++)
            {
                dbg_u_max_tmp.push_back(u_max_tmp(tag));
            }
            */

            u_max.col(pp) = u_max_tmp;
        }
    }

    // Writing out u_max to a result file
    std::ofstream resultFile(dirOutput + "/MaxWindSpeed.csv");
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    resultFile << u_max.format(CSVFormat);
    resultFile.close();

    return 0;
}


int WindFieldModel::inpolygon(Eigen::ArrayXd PolyX, Eigen::ArrayXd PolyY, int N, double px, double py)
{
    int counter = 0;
    int ii;
    double xinters;
    double P1x = PolyX(0);
    double P1y = PolyY(0);
    double P2x, P2y;

    for (int ii = 1; ii < N + 1; ii++)
    {
        if (ii < N)
        {
            P2x = PolyX(ii);
            P2y = PolyY(ii);
        }
        else
        {
            P2x = PolyX(0);
            P2y = PolyY(0);
        }
        if (py > std::min(P1y, P2y))
        {
            if (py <= std::max(P1y, P2y))
            {
                if (px <= std::max(P1x, P2x))
                {
                    if (P1y != P2y)
                    {
                        xinters = (py - P1y) * (P2x - P1x) / (P2y - P1y) + P1x;
                        if (P1x == P2x || px <= xinters)
                        {
                            counter++;
                        }
                    }
                }
            }
        }
        P1x = P2x;
        P1y = P2y;
    }
    if (counter % 2 == 0)
        return (0);
    else
        return (1);
}


double WindFieldModel::wrapTo360(double angle)
{
    angle = fmod(angle, 360);

    if (angle < 0)
        angle += 360;

    return angle;
}


int WindFieldModel::min(int a, int b)
{
    if (a > b)
        return b;
    else
        return a;
}


int main(int argc, char *argv[])
{
    // input files and directories
    std::string config_file, site_file, track_file, latw_file;
    std::string perturb_file, terrain_file;
    // output files
    std::string z0_dir, pws_dir;
    // loading arguments
    if (argc == 17)
    {
        config_file = argv[2];
        site_file = argv[4];
        track_file = argv[6];
        latw_file = argv[8];
        perturb_file = argv[10];
        terrain_file = argv[12];
        z0_dir = argv[14];
        pws_dir = argv[16];
    } else {
        std::cout << "WindFieldModel: please include --config {} --site {} --track {} --latw {} --pert {} --terrain {} --z0 {} output {}" << std::endl;
        return 1;
    }
    std::cout << "WindFieldModel: computaiton started" << std::endl;
    WindFieldModel a;
    a.ConfigSimu(config_file, site_file, track_file, latw_file);
    a.PertubPath(perturb_file);
    a.DefineTern(terrain_file);
    a.ComputeStationZ0(z0_dir);
    a.SimulateWind(pws_dir);
    return 0;
}