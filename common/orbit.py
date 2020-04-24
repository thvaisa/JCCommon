import numpy as np
import json
import orekit
import math
from pathlib import Path
from org.orekit.time import AbsoluteDate, TimeScalesFactory 
from org.orekit.orbits import KeplerianOrbit
from scipy.spatial.transform import Rotation as R
from org.orekit.frames import FramesFactory
from orekit.pyhelpers import setup_orekit_curdir
import org.orekit.utils as utils

from org.orekit.orbits import KeplerianOrbit, PositionAngle # pylint: disable=import-error
from org.orekit.attitudes import Attitude, FixedRate # pylint: disable=import-error
from org.orekit.propagation.analytical import KeplerianPropagator # pylint: disable=import-error
from org.orekit.time import AbsoluteDate
from org.hipparchus.geometry.euclidean.threed import Rotation, RotationConvention, Vector3D 






RAD_CONVERSION = np.pi/180.0

def solve_rotation_matrix(RA,DEC):
    h = RA
    d = DEC 

    r_x = R.from_euler('x', DEC, degrees=False)
    r_z = R.from_euler('z', RA, degrees=False)
    return r_z*r_x

class PlanetaryObject:


    def __init__(self,input_data):


        if(type(input_data) is str):
            with open(fname) as json_file:
                self.data = json.load(json_file)
        elif(type(input_data) is dict):
            self.data = input_data

        #################### orekit VM init ####################
        OREKIT_DATA_FILE = self.data["orekit_data_path"]
        OREKIT_VM = orekit.initVM() # pylint: disable=no-member
        setup_orekit_curdir(str(OREKIT_DATA_FILE))
        self.vm = orekit.initVM()    

        self.ref_frame = FramesFactory.getICRF()
        self.trj = {}
        data = self.data
        trj = self.trj
        trj["a"]   = data["orbit"]["a"]*1000
        trj["e"]   = data["orbit"]["e"]
        trj["i"]   = data["orbit"]["i"]
        trj["omega"] = data["orbit"]["omega"]
        trj["Omega"]   = data["orbit"]["ohm"]
        trj["M"] = data["orbit"]["M0"]
        mu = utils.Constants.IAU_2015_NOMINAL_SUN_GM
        RA      = math.radians(data["common"]["RA"])
        DEC     = math.radians(data["common"]["DEC"])
        date = data["date"]
        #date["year"] = data['comet']["common"]["year"]
        #date["month"] = data['comet']["commont"]["month"]
        #date["hour"] = data['comet']["common"]["hour"]
        #date["minutes"] = data['comet']["common"]["minutes"]
        #date["seconds"] = data['comet']["common"]["seconds"]
       # date = trj["date"]

        timescale = TimeScalesFactory.getTDB()
        trj_date = AbsoluteDate(int(date["year"]),
                                 int(date["month"]),
                                 int(date["day"]),
                                 int(date["hour"]),
                                 int(date["minutes"]),
                                 float(date["seconds"]),
                                 timescale)
      
        date = data["encounter_date"]
        self.encounter_date = AbsoluteDate(int(date["year"]),
                                 int(date["month"]),
                                 int(date["day"]),
                                 int(date["hour"]),
                                 int(date["minutes"]),
                                 float(date["seconds"]),
                                 timescale)



        trj["T"] = data["common"]["T_period"]
        self.rotation_rate = 2.0 * math.pi / (trj["T"] * 3600)

        trajectory = KeplerianOrbit(trj["a"],
                                         trj["e"],
                                         math.radians(trj["i"]),
                                         math.radians(trj["omega"]),
                                         math.radians(trj["Omega"]),
                                         math.radians(trj["M"]),
                                         PositionAngle.MEAN,
                                         self.ref_frame,
                                         trj_date, 
                                         mu)

       
        rot_conv = RotationConvention.valueOf("VECTOR_OPERATOR")


        self.R = solve_rotation_matrix(RA,DEC).as_matrix()
        rotation = utils.AngularCoordinates(Rotation.IDENTITY, 
                                            Vector3D(0., 0., self.rotation_rate))
        attitude = Attitude(trj_date, self.ref_frame, rotation)
        att_provider = FixedRate(attitude)

        # Create propagator
        self.propagator = KeplerianPropagator(trajectory, att_provider)
 




    def get_pos(self):
        prop = self.propagator.propagate(self.encounter_date)
        self.pos = prop.getPVCoordinates(self.ref_frame).getPosition()
        self.rot = prop.getAttitude().getRotation().getAngle()
        return self.pos





