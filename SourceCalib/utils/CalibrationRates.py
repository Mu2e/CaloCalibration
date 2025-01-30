{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "ec2f37cd-5f5b-4e6d-b17e-b3e63d2312bf",
   "metadata": {},
   "source": [
    "#### This notebook is to calculate the calibration rates.\n",
    "#### The equations and quanities are taken from DocDB-2904\n",
    "##### Eq (1) OR (2) are the same, but calculate different quanties (i.e # of calibration gammas per crystal given calibration time OR time needed to calibrate each crystal given number of calibration gammas per crystal, respectively)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "3aae1e6d-03f1-4a38-812f-1fd8ed8626c0",
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "import math"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "edfa08a9-102d-4836-bfb4-0caece3689f3",
   "metadata": {},
   "source": [
    "#### List of quantities:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 30,
   "id": "23e32c98-46cd-4fe7-bda2-a415023d7f2e",
   "metadata": {},
   "outputs": [],
   "source": [
    "sim_gam = 6.5e+7     #total number of photons simulated in GEANT4\n",
    "n_xtal = 1348        #total number of crystals \n",
    "V_bath = 0.021       #Volume of bath at dt generator \n",
    "#n_16N = 1.4e+8       #Number of 16N in bath  \n",
    "a = 0.40             #Attenuation from bath to furthest crystal\n",
    "tao = 10.3           #Lifetime of 16N\n",
    "n_gam =63000         #number of calibration gammas per calibration per crystal\n",
    "e_geom = 0.3         #(n_xtal)/ (7e-5 * sim_gam ) #Geometric efficiency to enter crystal\n",
    "f_16N = 0.3          #fractions of inelastic reactions producing 16N\n",
    "w_bath=0.16          #width of bath as seen by neutron\n",
    "mfp = 3.0            #mean free path of neutron in fluid\n",
    "e_bath = 0.88        # geoemtric efficiency of bath"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4ff67599-7d83-4c16-b4d5-dd4c6f5990ae",
   "metadata": {},
   "source": [
    "##### To account for changing neutron yield: change the neutron rate from the DT generator, R_DT, below and multiply by the reduced yield"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 50,
   "id": "759971b9-1617-464b-a221-18e7ad72fd11",
   "metadata": {},
   "outputs": [],
   "source": [
    "R_DT = 1e+09 * 0.7        #neutron rate from DT generator"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "675127a7-9d91-4699-8e89-3fca9f0987eb",
   "metadata": {},
   "source": [
    "#### Equations from paper:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 51,
   "id": "e0399149-4e7d-47b3-958d-e2c57bd84ab3",
   "metadata": {},
   "outputs": [],
   "source": [
    "T_prod_16N = f_16N*R_DT*(1-math.exp(-w_bath/mfp))*e_bath #rate to produce 16N in bath \n",
    "n_16N = T_prod_16N*tao                                   #Number of 16N in bath  \n",
    "density_16N = n_16N/ V_bath                              #Number density of 16N in bath\n",
    "V_xtal =  0.0019/n_xtal                                  #Volume of fluid at a crystal\n",
    "n_N_xtal = V_xtal*density_16N*a                          #Number of 16N at crystal"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "bac64167-ffc3-4765-8b03-633576e45f1b",
   "metadata": {},
   "source": [
    "#### (1) Equation to calculate number of calib gammas per crystal, given calibration time (T_cal)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 52,
   "id": "971d8360-ea64-40f2-9409-65d38c4d10d4",
   "metadata": {},
   "outputs": [],
   "source": [
    "T_cal = 10*60        #time in seconds\n",
    "n_gam = (e_geom*T_cal*n_N_xtal)/tao  #number of calibration gammas per calibration per crystal"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 53,
   "id": "62235575-1a93-4f3a-9d40-d53b1236bd3a",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Number of photons simulated in GEANT4:  6.50e+07\n",
      "Number of calibration gammas per crystal:  4.64e+04\n",
      "Time needed for calibration (in minutes):  10.0\n",
      "Time needed for calibration (in seconds):  600\n"
     ]
    }
   ],
   "source": [
    "print(\"Number of photons simulated in GEANT4: \",\"{:0.2e}\".format(sim_gam))\n",
    "print(\"Number of calibration gammas per crystal: \", \"{:0.2e}\".format(n_gam))\n",
    "print(\"Time needed for calibration (in minutes): \",\"{:0.3}\".format(T_cal/60))\n",
    "print(\"Time needed for calibration (in seconds): \",T_cal)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "dfc9fea5-2414-4a0b-9fbb-5a9e60dc5986",
   "metadata": {},
   "source": [
    "#### (2) Equation to calculate calibration time, given number of calibration gammas per crystal (calib_gammas)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 54,
   "id": "878d73eb-8175-4137-b2f8-1e4c7fb20036",
   "metadata": {},
   "outputs": [],
   "source": [
    "calib_gammas = 63000\n",
    "T_calib = (calib_gammas*tao)/(e_geom*n_N_xtal) #Time duration of calibration (in seconds)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 55,
   "id": "12819c6e-1707-41bb-b726-c6ca1aa6e854",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Time needed for calibration provided 63000 calibration gammas:  13.6 minutes\n",
      "Time needed for calibration provided 63000 calibration gammas:  814.9742933152714 seconds\n"
     ]
    }
   ],
   "source": [
    "print(\"Time needed for calibration provided\",calib_gammas, \"calibration gammas: \",\"{:0.3}\".format(T_calib/60), \"minutes\")\n",
    "print(\"Time needed for calibration provided\",calib_gammas, \"calibration gammas: \",T_calib,\"seconds\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d57fdbee-f2b1-46db-8f52-4ebef5fb5476",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1b34a1a8-7bb2-4829-8372-4ddf10643a68",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.9"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
