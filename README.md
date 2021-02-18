# ChemProps
Python script that creates and updates ChemProps DB in mongoDB. This is a one-stop script that includes downloading google spreadsheet as raw data, processing raw data, creating and updating ChemProps DB. Need to configure some environmental parameters before use.

By Bingyin Hu

### 1. System preparations

Required packages:

- requests
  - https://pypi.org/project/requests/
  - For downloading google spreadsheets.

- xlrd
  - https://github.com/python-excel/xlrd
  - Read the input Excel files.

- pymongo
  - https://pypi.org/project/pymongo/
  - MongoDB API.

- MechanicalSoup
  - https://mechanicalsoup.readthedocs.io/en/stable/index.html
  - For SMILES translation.

- beautifulsoup4
  - https://pypi.org/project/beautifulsoup4/
  - For SMILES translation.

- logging
  - Python default package

Open the command or terminal and run
```
pip install -r requirements.txt
```

**Python 3 must be 3.6 and later release to run SMILES translation flawlessly**

### 2. How to run

1. Note this script is tailored to the NanoMine workflow. It starts from getting data from a google spreadsheet all the way to updating the MongoDB.

2. **If you already have NanoMine installed on your system, then you can skip this step.**
Before running this script, make sure you have MongoDB configured and running on your system. Details see: https://www.mongodb.com/download-center/community
After mongod service starts on your system, you will need to create an admin user and set the following environment parameters: $NM_MONGO_PORT (port for MongoDB, usually 27017) $NM_MONGO_USER (username for the MongoDB user) $NM_MONGO_PWD (password for the MongoDB user).

3. Open a terminal and copy the commands in `mongoUpdateAdminUser` to grant the admin user readWrite role in ChemProps and unknown DBs. Otherwise, you will not be able to update info in those DBs but only create them.

4. Follow `gs.config.example` and `mongo.config.example` to create the configuration files for google spreadsheet as `gs.config` and MongoDB as `mongo.config`.

5. In python, run
```
from nmChemPropsPrepare import nmChemPropsPrepare
nm = nmChemPropsPrepare()
```
This command will kickoff the downloading and processing jobs. To create or update the MongoDB ChemProps and unknowns DB, run
```
nm.updateMongoDB()
```

6. Logs can be found in `nmChemPropsPrepare.log`.

7. MongoDB Compass can be useful during testing. It provides a GUI for your MongoDB. See: https://www.mongodb.com/products/compass