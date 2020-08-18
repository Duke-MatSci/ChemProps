#!/usr/bin/env python

# 06/04/2019 Bingyin Hu
# 05/14/2020 Bingyin Hu updateMongoDB fn updated
import requests
import xlrd
from pymongo import MongoClient
import logging
import string
import os

class nmChemPropsPrepare():
    def __init__(self):
        # load logger
        logging.basicConfig(filename='nmChemPropsPrepare.log',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level = logging.INFO
                           )
        self.loadGSconfig()
        self.loadMGconfig()
        # downloadGS
        self.downloadGS()
        # if we changed something in mongo, gsUpdate list will record these changes
        # admins need to manually address these changes in the google spreadsheet
        self.gsUpdate = []
        # prepare filler and polymer data
        self.filler = dict()
        self.polymer = dict()
        self.prepFiller()
        self.prepPolymer()
        # mongo init
        # self.client = MongoClient('mongodb://%s:%s@%s:%s/tracking?authSource=admin'
        if 'NM_MONGO_CHEMPROPS_URI' in self.env:
            self.client = MongoClient(self.env['NM_MONGO_CHEMPROPS_URI'])
        else:
            self.client = MongoClient('mongodb://%s:%s@%s:%s/%s'
                                      %(self.env['NM_MONGO_USER'],
                                        self.env['NM_MONGO_PWD'],
                                        self.env['NM_MONGO_HOST'],
                                        self.env['NM_MONGO_PORT'],
                                        self.env['NM_MONGO_DB']
                                        )
                                      )

    # load google spreadsheet configurations
    def loadGSconfig(self):
        # read gs.config for configurations
        with open("gs.config", "r") as f:
            confs = f.read().split('\n')
        self.url_format = confs[0]
        self.key = confs[1]
        self.format = "xlsx"
        # gids
        self.gids = dict()
        for i in range(2,len(confs)):
            kv = confs[i]
            k = kv.split(':')[0].strip()
            v = kv.split(':')[1].strip()
            self.gids[k] = v
    
    # load mongo configurations
    def loadMGconfig(self):
        self.env = dict()
        cpuri = os.environ.get('NM_MONGO_CHEMPROPS_URI', None)
        if cpuri:
            self.env['NM_MONGO_CHEMPROPS_URI'] = cpuri
        else:
            # read mongo.config for configurations
            with open("mongo.config", "r") as f:
                confs = f.read().split('\n')
            for i in range(len(confs)):
                kv = confs[i]
                k = kv.split(':')[0].strip()
                v = kv.split(':')[1].strip()
                self.env[k] = v

    # download google spreadsheets
    def downloadGS(self):
        for fname in self.gids:
            resp = requests.get(self.url_format %(self.key,
                                                  self.format,
                                                  self.gids[fname]
                                                 )
                               )
            with open(fname + ".xlsx", "wb") as f:
                f.write(resp.content)
                logging.info("%s sheet is downloaded as %s.xlsx" %(fname, fname))

    # prepare ChemProps.polymer data
    def prepPolymer(self):
        xlfile = xlrd.open_workbook("matrixRaw.xlsx") # change the filename if gs.config is changed
        sheet = xlfile.sheets()[0] # one sheet per xlsx file
        header = sheet.row_values(0) # SMILES;uSMILES;std_name;density(g/cm3);density_std_err(g/cm3);abbreviations;synonyms;tradenames
        # create a map for headers
        hmap = {}
        for i in range(len(header)):
            hmap[header[i]] = i
        # loop
        for row in range(1, sheet.nrows):
            rowdata = sheet.row_values(row)
            # skip the unfilled items
            if (len(rowdata[hmap['abbreviations']]) == 0 and
                len(rowdata[hmap['synonyms']]) == 0 and
                len(rowdata[hmap['tradenames']]) == 0):
                continue
            # otherwise save the data to self.polymer
            if rowdata[hmap['uSMILES']] not in self.polymer:
                self.polymer[rowdata[hmap['uSMILES']]] = {
                    "_id": rowdata[hmap['uSMILES']],
                    "_stdname": rowdata[hmap['std_name']],
                    "_abbreviations": [],
                    "_synonyms": [],
                    "_tradenames": [],
                    "_density": rowdata[hmap['density(g/cm3)']],
                    "_boc": [self.bagOfChar(rowdata[hmap['std_name']])]
                }
            # abbreviations
            if len(rowdata[hmap['abbreviations']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_abbreviations'] = self.striplist(rowdata[hmap['abbreviations']].split(';'))
                for abb in self.polymer[rowdata[hmap['uSMILES']]]['_abbreviations']:
                    self.polymer[rowdata[hmap['uSMILES']]]['_boc'].append(self.bagOfChar(abb))
            # synonyms
            if len(rowdata[hmap['synonyms']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_synonyms'] = self.striplist(rowdata[hmap['synonyms']].split(';'))
                for syn in self.polymer[rowdata[hmap['uSMILES']]]['_synonyms']:
                    self.polymer[rowdata[hmap['uSMILES']]]['_boc'].append(self.bagOfChar(syn))
            # tradenames
            if len(rowdata[hmap['tradenames']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_tradenames'] = self.striplist(rowdata[hmap['tradenames']].split(';'))
                for tra in self.polymer[rowdata[hmap['uSMILES']]]['_tradenames']:
                    self.polymer[rowdata[hmap['uSMILES']]]['_boc'].append(self.bagOfChar(tra))
        # log
        logging.info("Finish processing the polymer data.")


    # prepare ChemProps.filler data
    def prepFiller(self):
        xlfile = xlrd.open_workbook("fillerRaw.xlsx") # change the filename if gs.config is changed
        sheet = xlfile.sheets()[0] # one sheet per xlsx file
        header = sheet.row_values(0) # nm_entry/std_name/density_g_cm3
        # create a map for headers
        hmap = {}
        for i in range(len(header)):
            hmap[header[i]] = i
        # loop
        for row in range(1, sheet.nrows):
            rowdata = sheet.row_values(row)
            if rowdata[hmap['std_name']] not in self.filler:
                self.filler[rowdata[hmap['std_name']]] = {
                    "_id":rowdata[hmap['std_name']],
                    "_density": rowdata[hmap['density_g_cm3']],
                    "_alias":[],
                    "_boc": [self.bagOfChar(rowdata[hmap['std_name']])]
                }
            self.filler[rowdata[hmap['std_name']]]['_alias'].append(rowdata[hmap['nm_entry']])
            self.filler[rowdata[hmap['std_name']]]['_boc'].append(self.bagOfChar(rowdata[hmap['nm_entry']]))
        # log
        logging.info("Finish processing the filler data.")

    # update MongoDB
    def updateMongoDB(self):
        # initPolymer = False # a flag inidicating whether this is the first time creating the ChemProps.polymer collection
        # initFiller = False # a flag inidicating whether this is the first time creating the ChemProps.filler collection
        # if 'NM_MONGO_CHEMPROPS_URI' not in self.env :
        #     dbnames = self.client.list_database_names() # check if db exists (not strictly necessary.
        #                                               # checking for collections within db is sufficient and does not require
        #                                               # special privileges.
        #     if 'ChemProps' not in dbnames:
        #         initPolymer = True
        #         initFiller = True
        cp = self.client.ChemProps
        clctnames = cp.list_collection_names() # check if collection exists
        if 'polymer' in clctnames:
            cp.polymer.drop()
            logging.info("The polymer collection in the ChemProps DB is dropped.")
        if 'filler' in clctnames:
            cp.filler.drop()
            logging.info("The filler collection in the ChemProps DB is dropped.")
        # create polymer
        pol = cp.polymer
        posted_polymer = pol.insert_many(list(self.polymer.values()))
        logging.info("The polymer collection in the ChemProps DB is created for the first time.")
        # create filler
        fil = cp.filler
        posted_filler = fil.insert_many(list(self.filler.values()))
        logging.info("The filler collection in the ChemProps DB is created for the first time.")
        ## old way of updating
        # # if ChemProps exists
        # if not initPolymer and 'polymer' not in clctnames:
        #     initPolymer = True
        # if not initPolymer and 'filler' not in clctnames:
        #     initFiller = True
        # ## first creation cases (polymer)
        # if initPolymer:
        #     pol = cp.polymer
        #     posted_polymer = pol.insert_many(list(self.polymer.values()))
        #     logging.info("The polymer collection in the ChemProps DB is created for the first time.")
        # ## update cases (polymer)
        # else: 
        #     # first loop through the cp.polymer collection, remove documents that no longer exists in self.polymer dict   
        #     for mgUSMILES in cp.polymer.find().distinct("_id"):   
        #         if mgUSMILES not in self.polymer: 
        #             cp.polymer.delete_one({"_id": mgUSMILES}) 
        #             logging.info("Remove polymer with uSMILES (_id): %s from ChemProps." %(mgUSMILES))
        #     # loop through the items in the self.polymer dict, see if everything matches
        #     for uSMILES in self.polymer:
        #         gsData = self.polymer[uSMILES] # google spreadsheet data
        #         # if MongoDB doesn't have this uSMILES document
        #         if cp.polymer.find({"_id": uSMILES}).count() == 0:
        #             # insert it directly
        #             cp.polymer.insert(self.polymer[uSMILES])
        #             logging.info("Insert polymer with uSMILES (_id): %s to ChemProps." %(uSMILES))
        #             continue
        #         # otherwise, take the first and only result
        #         mgData = cp.polymer.find({"_id": uSMILES})[0] # mongo data, find by _id
        #         # continue if there is no difference between gsData and mgData
        #         if gsData == mgData:
        #             continue
        #         # otherwise, find the difference
        #         # if gsData is a superset of mgData, update mgData
        #         # if mgData is a superset of gsData, record the difference in self.gsUpdate
        #         # the structure of result see compareDict() instruction
        #         d1name = 'google sheet'
        #         d2name = 'mongo'
        #         result = self.compareDict(d1 = gsData,
        #                                   d1name = d1name,
        #                                   d2 = mgData,
        #                                   d2name = d2name,
        #                                   imtbKeys = {'_id',
        #                                               '_stdname',
        #                                               '_density'
        #                                              }
        #                                  )
        #         # required updates for gsData go to self.gsUpdate
        #         for change in result['google sheet']:
        #             self.gsUpdate.append(
        #                 "%s %s to %s of the polymer with uSMILES: %s."
        #                 %(change[0], change[2], change[1], uSMILES))
        #         # apply/update the changes
        #         for change in result[d2name]:
        #             cp.polymer.update(
        #                 {"_id": uSMILES},
        #                 {"%s" %(change[0]): { change[1]: change[2]}}
        #             )
        #             logging.info("Apply %s with value %s to %s of the polymer with uSMILES (_id): %s in ChemProps."
        #                          %(change[0], change[2], change[1], uSMILES)
        #                         )
        #             # rebuild _boc from ground up as well, "_stdname", "_abbreviations", "_synonyms"
        #             boc = []
        #             document = cp.polymer.find({"_id": uSMILES})[0]
        #             boc.append(self.bagOfChar(document["_stdname"]))
        #             for abb in document["_abbreviations"]:
        #                 boc.append(self.bagOfChar(abb))
        #             for syn in document["_synonyms"]:
        #                 boc.append(self.bagOfChar(syn))
        #             cp.polymer.update(
        #                 {"_id": uSMILES},
        #                 {"$set": {"_boc": boc}})
        #             logging.info("_boc updated for the polymer with uSMILES (_id): %s in ChemProps." %(uSMILES))
        #     # end of the loop
        # ## first creation cases (filler)
        # if initFiller:
        #     fil = cp.filler
        #     posted_filler = fil.insert_many(list(self.filler.values()))
        #     logging.info("The filler collection in the ChemProps DB is created for the first time.")
        # ## update cases (filler)
        # else:
        #     # loop through the items in the self.filler dict, see if everything matches
        #     for std_name in self.filler:
        #         # same structure as self.polymer
        #         gsData = self.filler[std_name]
        #         # if MongoDB doesn't have this std_name document
        #         if cp.filler.find({"_id": std_name}).count() == 0:
        #             # insert it directly
        #             cp.filler.insert(self.filler[std_name])
        #             logging.info("Insert filler with std_name (_id): %s to ChemProps." %(std_name))
        #             continue
        #         # otherwise, take the first and only result
        #         mgData = cp.filler.find({"_id": std_name})[0]
        #         if gsData == mgData:
        #             continue
        #         d1name = 'google sheet'
        #         d2name = 'mongo'
        #         result = self.compareDict(d1 = gsData,
        #                                   d1name = d1name,
        #                                   d2 = mgData,
        #                                   d2name = d2name,
        #                                   imtbKeys = {'_id',
        #                                               '_density'
        #                                              }
        #                                  )
        #         # required updates for gsData go to self.gsUpdate
        #         for change in result['google sheet']:
        #             self.gsUpdate.append(
        #                 "%s %s to %s of the filler with std_name: %s."
        #                 %(change[0], change[2], change[1], std_name))
        #         # apply/update the changes
        #         for change in result[d2name]:
        #             cp.filler.update(
        #                 {"_id": std_name},
        #                 {"%s" %(change[0]): { change[1]: change[2]}}
        #             )
        #             logging.info("Apply %s with value %s to %s of the filler with std_name (_id): %s in ChemProps."
        #                          %(change[0], change[2], change[1], std_name)
        #                         )
        #             # rebuild _boc from ground up as well, "_id", "_alias"
        #             boc = []
        #             document = cp.filler.find({"_id": std_name})[0]
        #             boc.append(self.bagOfChar(document["_id"]))
        #             for ali in document["_alias"]:
        #                 boc.append(self.bagOfChar(ali))
        #             cp.filler.update(
        #                 {"_id": std_name},
        #                 {"$set": {"_boc": boc}})
        #             logging.info("_boc updated for the filler with std_name (_id): %s in ChemProps." %(std_name))
        #     # end of the loop
        # ## append gsUpdate records as WARNING to the log
        # for rec in self.gsUpdate:
        #     logging.warn(rec)
        # self.gsUpdate = [] # reset gsUpdate

    # remove leading and trailing white spaces
    def striplist(self, mylist):
        for i in range(len(mylist)):
            mylist[i] = mylist[i].strip()
        return mylist

    # compare two dicts, need to specify the keys to the immutable objects,
    # the function returns a result dict that indicates objects that do not
    # exist in the current dict but exist in the other dict for each dict.
    # DO NOT SUPPORT NESTED DICTS
    # example:
    # d1 = {'k1': [1,2], 'k2': 'new'}
    # d2 = {'k1': [1,3], 'k2': 'old'}
    # result = {'d1': [('$addToSet', 'k1', 3)], 'd2': [('$addToSet', 'k1', 2), ('$set', 'k2', 'new')]}
    def compareDict(self, d1, d1name, d2, d2name, imtbKeys):
        result = {d1name: [], d2name: []} # init output dict
        # prepPolymer guarantees d1 and d2 will have the same keys set even if
        # some keys will have empty string or list (except for _boc)
        allKeys = set(d1.keys())
        for key in allKeys:
            # skip _boc
            if key == '_boc':
                continue
            # immutables always trust d1 has the latest version
            if key in imtbKeys:
                if d1[key] != d2[key]:
                    result[d2name].append(('$set', key, d1[key]))
            # non immutables
            else:
                # use set.difference() function to get the result
                # in d1[key] not in d2[key]
                d1subd2 = set(d1[key]).difference(set(d2[key]))
                # in d2[key] not in d1[key]
                d2subd1 = set(d2[key]).difference(set(d1[key]))
                # update result
                for addTod2 in d1subd2:
                    result[d2name].append(('$addToSet', key, addTod2))
                for addTod1 in d2subd1:
                    result[d1name].append(('$addToSet', key, addTod1))
        return result

    # a helper method to generate a bocStr for a string based on the occurrence
    # of the chars. Example: (a,b,c,...,y,z,0,1,2,3,4,5,6,7,8,9) to 101214...01
    def bagOfChar(self, myStr):
        bag = []
        for alphabet in string.ascii_lowercase:
            bag.append(str(myStr.lower().count(alphabet)))
        for number in string.digits:
            bag.append(str(myStr.lower().count(number)))
        return ''.join(bag)

if __name__ == '__main__':
    nm = nmChemPropsPrepare()
