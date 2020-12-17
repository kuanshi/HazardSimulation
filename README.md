# HazardSimulation

Command:
python HazardSimulation.py --hazard_config "PATH/TO/HazardConfiguration_EQ2.json"

User input:
1. A configuration json file (e.g., [HazardConfiguration_EQ1.json](https://github.com/kuanshi/HazardSimulation/blob/main/test/HazardConfiguration_EQ1.json))
2. A csv file for site locations and Vs30 values (e.g., [SiteFile.csv](https://github.com/kuanshi/HazardSimulation/blob/main/test/eq_input/SiteFile.csv)
The configuration json file contains 4 sections:

(1) Site list - providing the csv filename and specifying the site ID range
```json
"Site": {
    "Type": "From_CSV",
    "input_file": "SiteFile.csv",
    "min_ID": 1,
    "max_ID": 2
}
```

(2) Scenario information - specifying scenario type, number of scenarios, generator type, and earthquake rupture. Currently, two types of EqRupture are supported: "ERF" (earthquake rupture forecast) and "PointSource"
"ERF" provides users to select historical earthquakes as the scenarios to run. Historical earthquake data are available from 4 supported "Model"s: "WGCEP (2007) UCERF2 - Single Branch", "Mean UCERF3", "Mean UCERF3 FM3.1", and "Mean UCERF3 FM3.2". Users can specify keywords in "Name" to search earthquakes. "min_Mag", "max_Mag", and "max_Dist" are optional to filter earthquakes.
```json
"Scenario": {
    "Type": "Earthquake",
    "Number": 1,
    "Generator": "Selection",
    "EqRupture": {
        "Type": "ERF",
        "Model": "WGCEP (2007) UCERF2 - Single Branch",
        "Name": "San Andreas",
        "min_Mag": 7.0,
        "max_Mag": 8.0,
        "max_Dist": 200.0
    }
}
```
"PointSource" lets users to directly prescribe the earthquake properties including magnitude, earthquake focus location, average rake, and average dip.
```json
"Scenario": {
    "Type": "Earthquake",
    "Number": 1,
    "Generator": "Selection",
    "EqRupture": {
        "Type": "PointSource",
        "Magnitude": 6.0,
        "Location": {
            "Latitude": 37.9,
            "Longitude": -122.3,
            "Depth": 0.0
        },
        "AverageRake": 0.0,
        "AverageDip": 90.0
    }
}
```

(3) Event configuration - configuring how many ground motion records per site in each scenario ("NumberPerSite"), which ground motion prediction equation to be used ("GMPE"), which correlation models to be used ("CorrelationModel"), which intensity measures to be computed ("IntensityMeasure"), how much records can be scaled to be fit the target intensity measures ("ScalingFactor"), whether to save simulated IM ("SaveIM"), and whether to download records from "Database".
```json
"Event": {
    "NumberPerSite": 10,
    "GMPE": {
        "Type": "Chiou & Youngs (2014)",
        "Parameters": {}
    },
    "CorrelationModel": {
        "SaInterEvent": "Baker & Jayaram (2008)",
        "SaIntraEvent": "Markhvida et al. (2017)"
    },
    "IntensityMeasure": {
        "Type": "SA",
        "Periods": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    },
    "ScalingFactor": {
        "Maximum": 20.0,
        "Minimum": 0.1
    },
    "SaveIM": true,
    "Database": "NGAWest2",
    "UserName": "simcenter_test@outlook.com",
    "UserPassword": "123456",
    "OutputFormat": "SimCenterEvent"
}
```
Currently supported GMPEs: "Abrahamson, Silva & Kamai (2014)", "Boore, Stewart, Seyhan & Atkinson (2014)", "Campbell & Bozorgnia (2014)", and "Chiou & Youngs (2014)". Currently supported correlation models: "Baker & Jayaram (2008)" (for inter-event), "Jayaram & Baker (2009)", "Loth & Baker (2013)", and "Markhvida et al. (2017)" (for intra-event).

(4) Directory - defining the working, input, and output directories.
```json
"Directory":{
    "Work": "PATH/TO/HazardSimulation/",
    "Input": "PATH/TO/HazardSimulation/test/eq_input/",
    "Output": "PATH/TO/HazardSimulation/test/records/"
}
```

Output files:
1. A "SiteIM.json" file (if "SaveIM" is true) for simulated response spectra
2. A "records" folder (if "UserName", "UserPassword", and "Output" directory are provided) of selected ground motion records for the user-defined earthquake scenario(s).

License:
1. OpenSHA-1.5.2 is subjected to Apache-2.0 License
2. chromedriver is subjected to Apache-2.0 License

Extra Python packages:
[JPype1](https://pypi.org/project/JPype1/)

[selenium](https://pypi.org/project/selenium/)

[tqdm](https://pypi.org/project/tqdm/2.2.3/)



