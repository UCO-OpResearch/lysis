#Need Import Statements Here
import numpy as np
import h5py
import src.python.lysis as lysis
import os

class h5builder:
    """
    A class to build HDF5 files for Lysis.
    """
    #Need to create an Init that reads from the original fortran output files
    def __init__(self , experiment_code):
        e = lysis.util.Experiment(os.path.join("data"), experiment_code)
        e.read_file()
        e.macro_params.total_molecules

        #Creates a variable f for our HDF5 file
        f = h5py.File('schema_view.h5', 'w')

        #Creates folder and dataset structure
        runs = f.create_group("1-PKd") 
        micro_data = runs.create_group("micro_data") 
        macro_data = runs.create_group("macro_data")


        #Micro data set initializations
        pli_first_time = micro_data.create_dataset("pli_first_time", (50000 , ), dtype= np.float64)
        tpa_final_num = micro_data.create_dataset("tpa_final_num", (50000 , ), dtype= np.uint8)
        fiber_degraded = micro_data.create_dataset("fiber_degraded", (50000 , ), dtype= np.bool_)
        sim_final_time = micro_data.create_dataset("sim_final_time", (50000 , ), dtype= np.float64)
        pli_generated_num = micro_data.create_dataset("pli_generated_num", (50000 , ), dtype= np.uint16)
        tpa_leaving_time = micro_data.create_dataset("tpa_leaving_time", (50000 , ), dtype= np.float64)
        tpa_unbound_by_pli = micro_data.create_dataset("tpa_unbound_by_pli", (50000 , ), dtype= np.bool_)
        tpa_unbound_kinetic = micro_data.create_dataset("tpa_unbound_kinetic", (50000 , ), dtype = np.bool_)

        #Macro data set initializations

        for i in range(0 , e.macro_params.simulations):
            simulation = macro_data.create_group(f"sim_{i:02}")
            fiber_degrade_time = simulation.create_dataset("fiber_degrade_time", (1 , 3) , dtype=np.float64, maxshape = (None, 3))
            tpa_bind_events = simulation.create_dataset("tpa_bind_events", (1 , 10) , dtype=np.float64 , maxshape = (None, 10)) #snapshot time = 
            snapshot_time = simulation.create_dataset("snapshot_time", (1 , 1) , dtype=np.float64 , maxshape = (None, 1))
            #tpa_location_snapshot = simulation.create_dataset("tpa_location_snapshot", (e.macro_params.total_molecules , 1) , dtype=np.int32 , maxshape = (e.macro_params.total_molecules, None)) #I uncapped the max number of rows so the data can fit
            tpa_location_snapshot = simulation.create_dataset("tpa_location_snapshot", (e.macro_params.total_molecules , 1) , dtype=np.int32 , maxshape = (None, None))
            tpa_transit_time = simulation.create_dataset("tpa_transit_time", (e.macro_params.total_molecules , 1) , dtype=np.float64) #I uncapped the max number of rows so the data can fit
        
        #Reads in Micro scale data into file system
        file_code = "PLG2_tPA01_TB-xiii"

        pli_first_time[:] = np.fromfile(
            os.path.join(e.os_path, f"firstPLi_{file_code}.dat"),
        )

        tpa_final_num[:] = np.fromfile(os.path.join(e.os_path, f"lasttPA_{file_code}.dat"), dtype=np.int32)

        fiber_degraded[:] = np.fromfile(
            os.path.join(e.os_path, f"lyscomplete_{file_code}.dat"), 
            dtype=np.int32
        ).astype(bool)

        sim_final_time[:] = np.fromfile(os.path.join(e.os_path, f"lysis_{file_code}.dat"))

        pli_generated_num[:] = np.fromfile(os.path.join(e.os_path, f"PLi_{file_code}.dat"), dtype=np.int32).astype(np.uint16)

        tpa_leaving_time[:] = np.fromfile(os.path.join(e.os_path, f"tPA_time_{file_code}.dat"))

        tpa_unbound_by_pli[:] = np.fromfile(
            os.path.join(e.os_path, f"tPAPLiunbd_{file_code}.dat"), 
            dtype=np.int32
        ).astype(bool)

        tpa_unbound_kinetic[:] = np.fromfile(
            os.path.join(e.os_path, f"tPAPLiunbd_{file_code}.dat"), 
            dtype=np.int32
        ).astype(bool)

        #Reads in Macro scale Fiber Degrade Time Data
        for i in range (0 , 10):
            file_reference = f[f"1-PKd/macro_data/sim_0{i}/fiber_degrade_time"]
            macro_file_code = f"TB-xiii__21_105_0{i}"
            data = np.loadtxt(os.path.join(e.os_path, f"0{i}\\f_deg_list_{macro_file_code}.dat") , delimiter = ",")
            data = np.reshape(data, (-1, 3))
            data[: , 1] = data[: , 1] - 1   #The data in column 1 is 1 indexed, so we need to convert it to 0 indexed
            file_reference.resize(data.shape)
            file_reference[:] = data
            #Adding Attributes to datasets
            file_reference.attrs["units"] = ["second" , "dimensionless" , "second"]

        #Reads in Macro scale TPA Bind Events Data
        for i in range (0 , 10):
            file_reference = f[f"1-PKd/macro_data/sim_0{i}/tpa_bind_events"]
            macro_file_code = f"TB-xiii__21_105_0{i}"
            data = np.loadtxt(os.path.join(e.os_path, f"0{i}\\m_bind_t_{macro_file_code}.dat") , delimiter = ",")
            data = np.reshape(data, (-1, 4))
            data[:,[1,3]] = data[:,[1,3]] - 1   #The data in column 1 and 3 is 1 indexed, so we need to convert it to 0 indexed
            file_reference.resize(data.shape)
            file_reference[:] = data

        #Reads in Macro scale snapshot time data
        for i in range (0 , 10):
            file_reference = f[f"1-PKd/macro_data/sim_0{i}/tpa_location_snapshot"]
            macro_file_code = f"TB-xiii__21_105_0{i}"
            data = np.fromfile(os.path.join(e.os_path, f"0{i}\\m_loc_{macro_file_code}.dat") , dtype = np.int32).reshape(-1, e.macro_params.total_molecules)
            data = data - 1 #The data is 1 indexed, so we need to convert it to 0 indexed
            file_reference.resize(data.shape)
            file_reference[:] = data

        #Reads in Macro scale TPA Transit Time Data
        for i in range (0 , 10):
            file_reference = f[f"1-PKd/macro_data/sim_0{i}/tpa_transit_time"]
            macro_file_code = f"TB-xiii__21_105_0{i}"
            data = np.fromfile(os.path.join(e.os_path, f"0{i}\\mfpt_{macro_file_code}.dat") , dtype = np.float64).reshape(-1,1)
            file_reference[:] = data

        #Reads in Macro scale Snapshot Time Data
        for i in range (0 , 10):
            file_reference = f[f"1-PKd/macro_data/sim_0{i}/snapshot_time"]
            macro_file_code = f"TB-xiii__21_105_0{i}"
            data = np.fromfile(os.path.join(e.os_path, f"0{i}\\tsave_{macro_file_code}.dat") , dtype = np.float64).reshape(-1,1)
            file_reference.resize(data.shape)
            file_reference[:] = data

        #Closes the file
        f.close()

    #Also need to create a function that can read from an existing HDF file