#!/usr/bin/env python

# 06/04/2019 Bingyin Hu

from pymongo import MongoClient
import os
import logging
import re # for query reformat
import string # for query reformat
import fillerDensityModule as fDM # filler density module
from SMILEStrans import SMILEStrans # uSMILES translation module
from fillerDensityModule import removeNano # remove words like nanoparticles

class nmChemPropsAPI():
    def __init__(self, nmid):
        # input the NanoMine ID, example: L999_Someone_2020_S2
        self.nmid = nmid
        # load logger
        logging.basicConfig(filename='nmChemPropsInteract.log',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level = logging.INFO
                           )
        self.loadMGconfig()
        # mongo init
        # self.client = MongoClient('mongodb://%s:%s@localhost:27017/tracking?authSource=admin'
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
        # access to DBs
        self.cp = self.client.ChemProps
        self.__reset__()

    def __reset__(self):
        # logging the wf distribution
        self.polymer_wf = dict()
        self.filler_wf = dict()
        # logging the top three candidates and accumulated wfs
        self.top3polymers = [[],[],[]]
        self.top3fillers = [[],[],[]]

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

    # a helper method to generate a bocStr for a string based on the occurrence
    # of the chars. Example: (a,b,c,...,y,z,0,1,2,3,4,5,6,7,8,9) to 101214...01
    # The occurrence is capped by 9.
    def bagOfChar(self, myStr):
        bag = []
        for alphabet in string.ascii_lowercase:
            bag.append(str(min(9, myStr.lower().count(alphabet))))
        for number in string.digits:
            bag.append(str(min(9, myStr.lower().count(number))))
        return ''.join(bag)

    def analysisPolymers(self, keywords):
        self.searchPolymers(keywords)
        print("Keywords are: {}".format(keywords))
        print("=================")
        print("The top 3 matched polymers are:")
        for i in self.top3polymers:
            print(i)
        print("=================")
        print("Find details in .polymer_wf")

    # the main search function for polymer infos
    # will call six sub-methods for mapping, wf stands for weighting factor
    # 00 wf 1.6 apple to apple comparison for uSMILES (translate here by SMILEStrans), wf 5
    # 10 wf 1.0 apple to apple comparison for '_stdname' with rptname
    # 11 wf 1.0 apple to apple comparison for '_abbreviations' with rptname
    # 12 wf 1.0 apple to apple comparison for '_synonyms' with rptname
    # 20 wf 0.8 apple to apple comparison for '_abbreviations' with rptabbr
    # 21 wf 0.8 apple to apple comparison for '_boc' with rptabbrBOC
    # 30 wf 0.0 relaxed bag-of-word comparison for '_tradenames' with rpttrad
    # 31 wf 0.0 apple to apple comparison for '_boc' (alphabet only) with tradnameBOCalph
    # 40 wf 1.2 apple to apple comparison for '_boc' with rptnameBOC
    # 41 wf 1.2 apple to apple comparison for '_boc' (alphabet only) with rptnameBOCalph
    # 50 wf 1.0 relaxed bag-of-word comparison for '_stdname' with rptname
    # 51 wf 2.0 relaxed bag-of-word comparison for '_synonyms' with rptname
    # input format:
    #   {'ChemicalName': 'Poly(styrene)', 'Abbreviation': 'PS', 'TradeName': 'Dylite', 'SMILES': ''}
    #   NanoMine schema guarantees 'ChemicalName' has minimum occurrence of 1
    #   'Abbreviation', 'TradeName', and 'SMILES' are not required, since users might leave them blank
    # output format:
    # if there is a match:
    #   {'StandardName': _stdname, 'uSMILES': _id, 'density': _density}
    #   for multiple matches, return the one with highest cummulated wf 
    #   self.cp.ukpolymer.find({"_inputsmiles": rptuSMILES}) before exit,
    #   examine again whether the reported 'Abbreviation' and 'TradeName' are recorded in ChemProps, log them and manually check, if confirmed to be correct, add them to the google spreadsheet
    # if there is not a match:
    #   insert _inputname, _inputabbr, _inputsmiles, _nmid[] to unknowns.polymer
    def searchPolymers(self, keywords):
        # strip keywords
        for key in keywords:
            keywords[key] = keywords[key].strip()
        # init
        rptuSMILES = ''
        rptname = ''
        rptabbr = ''
        rpttrad = ''
        candidates = dict() # use '_id' as keys
        # 0) apple to apple comparison for uSMILES (translated here from SMILES by SMILEStrans), wf 5
        if 'SMILES' in keywords:
            rptuSMILES = keywords['SMILES']
            if len(rptuSMILES) > 0:
                try:
                    rptuSMILES = SMILEStrans(keywords['SMILES']).translate()
                except:
                    logging.warning("Error occurred during the SMILEStrans call for SMILES: %s" %(keywords['SMILES']))
                    pass
                for cand in self.cp.polymer.find({'_id': rptuSMILES}):
                    if cand['_id'] not in candidates:
                        candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                    candidates[cand['_id']]['wf'] += 1.6
                    candidates[cand['_id']]['features'].append('00')
        # 1) apple to apple comparison for polymer names (in _stdname, _abbreviations, _synonyms), wf 3
        rptname = re.escape(keywords['ChemicalName'])
        if len(rptname) > 1:
            # query for '_stdname' with rptname
            for cand in self.cp.polymer.find({'_stdname': {'$regex': '\\b%s\\b' %(rptname), '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1
                candidates[cand['_id']]['features'].append('10')
            # query for '_abbreviations' array
            for cand in self.cp.polymer.find({'_abbreviations': {'$regex': '\\b%s\\b' %(rptname), '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1
                candidates[cand['_id']]['features'].append('11')
            # query for '_synonyms' array
            for cand in self.cp.polymer.find({'_synonyms': {'$regex': '\\b%s\\b' %(rptname), '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1
                candidates[cand['_id']]['features'].append('12')
        # 2) apple to apple comparison for abbreviations (in _abbreviations), wf 2+1
        if 'Abbreviation' in keywords:
            rptabbr = re.escape(keywords['Abbreviation'])
            if len(rptabbr) > 0:
                # query for '_abbreviations' array, wf 2
                for cand in self.cp.polymer.find({'_abbreviations': {'$regex': '\\b%s\\b' %(rptabbr), '$options': 'i'}}):
                    if cand['_id'] not in candidates:
                        candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                    candidates[cand['_id']]['wf'] += 0.8
                    candidates[cand['_id']]['features'].append('20')
                # boc, wf 1
                rptabbrBOC = self.bagOfChar(rptabbr)
                for cand in self.cp.polymer.find({'_boc': {'$regex': rptabbrBOC}}):
                    if cand['_id'] not in candidates:
                        candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                    candidates[cand['_id']]['wf'] += 0.8
                    candidates[cand['_id']]['features'].append('21')
        # # 3) relaxed bag-of-word comparison for tradenames (in _tradenames), wf 1
        # if 'TradeName' in keywords:
        #     rpttrad = keywords['TradeName'].replace('-','\-').replace('(','\\(').replace(')','\\)').replace('[','\\[').replace(']','\\]')
        #     if len(rpttrad) > 0:
        #         # query for '_tradenames' array
        #         relBOW = self.containAllWords(rpttrad, '_tradenames', self.cp.polymer)
        #         for cand in relBOW:
        #             if cand['_id'] not in candidates:
        #                 candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
        #             candidates[cand['_id']]['wf'] += 1
        #             candidates[cand['_id']]['features'].append('30')
        #         # query for bag-of-character for only alphabets (use $regex '^0100020...')
        #         tradnameBOC = self.bagOfChar(rpttrad)
        #         tradnameBOCalph = tradnameBOC[:-10] # remove number's index
        #         for cand in self.cp.polymer.find({'_boc': {'$regex': '^%s' %(tradnameBOCalph)}}):
        #             if cand['_id'] not in candidates:
        #                 candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
        #             candidates[cand['_id']]['wf'] += 1
        #             candidates[cand['_id']]['features'].append('31')
        # 4) bag-of-character comparison for polymer names (in _boc), wf 2
        if len(rptname) > 0:
            rptnameBOC = self.bagOfChar(rptname)
            for cand in self.cp.polymer.find({'_boc': {'$regex': '%s' %(rptnameBOC)}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1.2
                candidates[cand['_id']]['features'].append('40')
            # only alphabets version, wf 1
            rptnameBOCalph = rptnameBOC[:-10] # remove number's index
            for cand in self.cp.polymer.find({'_boc': {'$regex': '^%s' %(rptnameBOCalph)}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1.2
                candidates[cand['_id']]['features'].append('41')
        # 5) relaxed bag-of-word comparison for polymer names (in _stdname, _synonyms), wf 2
        # query for '_stdname' array
            relBOWstd = self.containAllWords(rptname, '_stdname', self.cp.polymer)
            for cand in relBOWstd:
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 1
                candidates[cand['_id']]['features'].append('50')
            relBOWsyn = self.containAllWords(rptname, '_synonyms', self.cp.polymer)
            for cand in relBOWsyn:
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0, 'features': []}
                candidates[cand['_id']]['wf'] += 2
                candidates[cand['_id']]['features'].append('51')
        # end of the query part
        self.polymer_wf = candidates
        # if there is not a match:
        #   insert _inputname, _inputabbr, _inputsmiles, _nmid[] to unknowns.polymer
        if len(candidates) == 0:
            # prepare the insertion dict
            ukdict = {'_nmid': [],
                      '_inputname': [],
                      '_inputabbr': [],
                      '_inputtrade': [],
                      '_inputsmiles': ''}
            ukdict['_nmid'].append(self.nmid)
            ukdict['_inputname'].append(rptname)
            if len(rptabbr) > 0: ukdict['_inputabbr'].append(rptabbr)
            if len(rpttrad) > 0: ukdict['_inputtrade'].append(rpttrad)
            if len(rptuSMILES) > 0: # if smile reported, see if it's already in unknowns.polymer
                ukdict['_inputsmiles'] = rptuSMILES
                if self.cp.ukpolymer.find({"_inputsmiles": rptuSMILES}).count() == 0: # if not exist, create the document
                    # insert it directly
                    self.cp.ukpolymer.insert(ukdict)
                    logging.info("Insert unknown polymer with _inputsmiles: '%s' to unknowns." %(rptuSMILES))
                else:
                    # update the difference
                    ukdata = self.cp.ukpolymer.find({"_inputsmiles": rptuSMILES})[0]
                    if not self.lowerIn(rptname, ukdata['_inputname']):
                        self.cp.ukpolymer.update(
                            {"_inputsmiles": rptuSMILES},
                            {"$addToSet": { "_inputname": rptname}}
                        )
                        logging.warn("The document with _inputsmiles: '%s' has multiple _inputname! '%s' is newly added!" %(rptuSMILES, rptname))
                    if len(rptabbr) > 0 and not self.lowerIn(rptabbr, ukdata['_inputabbr']):
                        self.cp.ukpolymer.update(
                            {"_inputsmiles": rptuSMILES},
                            {"$addToSet": { "_inputabbr": rptabbr}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _inputabbr of the polymer with _inputsmiles: '%s' in unknowns."
                                 %(rptabbr, rptuSMILES)
                                )
                    if len(rpttrad) > 0 and not self.lowerIn(rpttrad, ukdata['_inputtrade']):
                        self.cp.ukpolymer.update(
                            {"_inputsmiles": rptuSMILES},
                            {"$addToSet": { "_inputtrade": rpttrad}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _inputtrade of the polymer with _inputsmiles: '%s' in unknowns."
                                 %(rpttrad, rptuSMILES)
                                )
                    if not self.lowerIn(self.nmid, ukdata['_nmid']):
                        self.cp.ukpolymer.update(
                            {"_inputsmiles": rptuSMILES},
                            {"$addToSet": { "_nmid": self.nmid}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _nmid of the polymer with _inputsmiles: '%s' in unknowns."
                                 %(self.nmid, rptuSMILES)
                                )
                    else:
                        logging.warn("The document with _inputsmiles: '%s' has duplicate in _nmid! Check '%s'!" %(rptuSMILES, self.nmid))
            else: # if smile not reported, see if the _inputname is already in unknowns.polymer
                if self.cp.ukpolymer.find({"_inputname": {'$regex': '%s' %(rptname)}}).count() == 0: # if not exist, create the document
                    # insert it directly
                    self.cp.ukpolymer.insert(ukdict)
                    logging.info("Insert unknown polymer with _inputname: '%s' to unknowns. No _inputsmiles reported." %(rptname))
                else:
                    # update the difference
                    ukdata = self.cp.ukpolymer.find({"_inputname": {'$regex': '%s' %(rptname)}})[0]
                    if len(rptabbr) > 0 and not self.lowerIn(rptabbr, ukdata['_inputabbr']):
                        self.cp.ukpolymer.update(
                            {"_inputname": {'$regex': '%s' %(rptname)}},
                            {"$addToSet": { "_inputabbr": rptabbr}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _inputabbr of the polymer with _inputsmiles: N/A, _inputname: '%s' in unknowns."
                                 %(rptabbr, rptname)
                                )
                    if len(rpttrad) > 0 and not self.lowerIn(rpttrad, ukdata['_inputtrade']):
                        self.cp.ukpolymer.update(
                            {"_inputname": {'$regex': '%s' %(rptname)}},
                            {"$addToSet": { "_inputtrade": rpttrad}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _inputtrade of the polymer with _inputsmiles: N/A, _inputname: '%s' in unknowns."
                                 %(rpttrad, rptname)
                                )
                    if not self.lowerIn(self.nmid, ukdata['_nmid']):
                        self.cp.ukpolymer.update(
                            {"_inputname": {'$regex': '%s' %(rptname)}},
                            {"$addToSet": { "_nmid": self.nmid}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _nmid of the polymer with _inputsmiles: N/A, _inputname: '%s' in unknowns."
                                 %(self.nmid, rptname)
                                )
                    else:
                        logging.warn("The document with _inputsmiles: N/A, _inputname: '%s' has duplicate nmid! Check '%s'!" %(rptname, self.nmid))
        # if there is a match:
        #   {'StandardName': _stdname, 'uSMILES': _id, 'density': _density}
        #   for multiple matches, return the one with highest cummulated wf
        #   before exit, examine again whether the reported 'Abbreviation' and 'TradeName' are recorded in ChemProps, log them and manually check, if confirmed to be correct, add them to the google spreadsheet
        else:
            # find the candidate with the highest wf
            wfs = [] # init a list to store wf
            wf_high = 0
            cand_high = []
            for cand in candidates:
                wfs.append(candidates[cand]['wf']) # append wf to wfs
                if candidates[cand]['wf'] > wf_high:
                    wf_high = candidates[cand]['wf']
                    cand_high = [candidates[cand]]
                elif candidates[cand]['wf'] == wf_high:
                    cand_high.append(candidates[cand])
            # set and sort wfs
            wfs = list(set(wfs))
            wfs.sort(reverse = True)
            # put top 3 polymers into self.top3polymers
            for i in range(min(len(wfs), 3)):
                for cand in candidates:
                    if candidates[cand]['wf'] == wfs[i]:
                        self.top3polymers[i].append(cand)
            # always return the first cand_high, but log if there's more than one cand
            if len(cand_high) > 1:
                tieWarning = "For the search package '%s', multiple winning matches found. Weighting factors tie at %d. They are:" %(str(keywords), wf_high)
                for candidate in cand_high:
                    tieWarning += "\n\t%s" %(candidate['data']['_stdname'])
                logging.warn(tieWarning)
            # warn admin if wf_high is no bigger than 2
            if wf_high <= 2:
                logging.warn("Careful! Low weighting factor %d for nmid '%s'! Mapped _inputname: '%s' with _stdname: '%s'." %(wf_high, self.nmid, rptname, cand_high[0]['data']['_stdname']))
            # now let's check whether the reported Abbreviation and Tradename are in the ChemProps
            # log them, if they are manually confirmed to be correct, add them
            # to the google sheet, and run nmChemPropsPrepare again. This way,
            # the boc will be updated during nmChemPropsPrepare as well.
            if len(rptabbr) > 0 and not self.lowerIn(rptabbr, cand_high[0]['data']['_abbreviations']):
                logging.warn("Admins please check whether '%s' is the abbreviation of polymer '%s'." %(rptabbr, cand_high[0]['data']['_stdname']))
            if len(rpttrad) > 0 and not self.lowerIn(rpttrad, cand_high[0]['data']['_tradenames']):
                logging.warn("Admins please check whether '%s' is the tradename of polymer '%s'." %(rpttrad, cand_high[0]['data']['_stdname']))
            return cand_high[0]['data']
        # otherwise, return None
        return None

    # the main search function for filler infos
    # will call six sub-methods for mapping, wf stands for weighting factor
    # 0) apple to apple comparison for filler names (in _id, _alias), wf 3
    # 1) bag-of-character comparison for filler names (in _boc), wf 2+1
    # 2) relaxed bag-of-word comparison for filler names (in _id, _alias), wf 1
    # input format:
    #   {'ChemicalName': 'silicon dioxide', 'Abbreviation': 'silica', 'TradeName': 'something'}
    #   NanoMine schema guarantees 'ChemicalName' has minimum occurrence of 1
    #   'Abbreviation' and 'TradeName'are not required, since users might leave them blank
    # output format:
    # if there is a match:
    #   {'StandardName': _id, 'density': _density}
    #   for multiple matches, return the one with highest cummulated wf
    # if there is not a match:
    #   call getFillerDensityGoogle(filler) i.e. filler density module see if there's a result
    #   if there is a result: returns (stdname, clean_result)
    #       log it as warning
    #       save it to ChemProps as
    #       {'_id': stdname, '_density': float(clean_result), '_alias': [reported name]}
    #   if there is not a result: return (stdname, -1)
    #       insert _inputname, _rawname, _nmid[] to unknowns.filler
    def searchFillers(self, keywords):
        # strip keywords
        for key in keywords:
            keywords[key] = keywords[key].strip()
        # init
        rptname = ''
        candidates = dict() # use '_id' as keys
        # 0) apple to apple comparison for filler names (in _id, _alias), wf 3
        rptname = re.escape(removeNano(keywords['ChemicalName']))
        if len(rptname) > 0:
            # query for '_id' with rptname
            for cand in self.cp.filler.find({'_id': {'$regex': '^'+rptname+'$', '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 4
            for cand in self.cp.filler.find({'_id': {'$regex': rptname, '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 3
            # query for '_alias' array
            for cand in self.cp.filler.find({'_alias': {'$regex': '^'+rptname+'$', '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 4
            for cand in self.cp.filler.find({'_alias': {'$regex': rptname, '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 3
            # 1) bag-of-character comparison for filler names (in _boc), wf 2
            rptnameBOC = self.bagOfChar(rptname)
            for cand in self.cp.filler.find({'_boc': {'$regex': '%s' %(rptnameBOC)}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 2
            # only alphabets version, wf 1
            rptnameBOCalph = rptnameBOC[:-10] # remove number's index
            for cand in self.cp.filler.find({'_boc': {'$regex': '^%s' %(rptnameBOCalph)}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 1
            # 2) relaxed bag-of-word comparison for filler names (in _id, _alias), wf 1
            # query for '_id'
            relBOWstd = self.containAllWords(rptname, '_id', self.cp.filler)
            for cand in relBOWstd:
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 1
            relBOWsyn = self.containAllWords(rptname, '_alias', self.cp.filler)
            for cand in relBOWsyn:
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 1
        # 3) apple to apple comparison for abbreviations (in _id, _alias), wf 2+1
        if 'Abbreviation' in keywords and len(keywords['Abbreviation']) > 0:
            rptabbr = re.escape(keywords['Abbreviation'])
            if len(rptabbr) > 0:
                # query for '_id' array, wf 2
                for cand in self.cp.filler.find({'_id': {'$regex': '^'+rptabbr+'$', '$options': 'i'}}):
                    if cand['_id'] not in candidates:
                        candidates[cand['_id']] = {'data': cand, 'wf': 0}
                    candidates[cand['_id']]['wf'] += 2
                # query for '_alias' array, wf 1
                for cand in self.cp.filler.find({'_alias': {'$regex': '^'+rptabbr+'$', '$options': 'i'}}):
                    if cand['_id'] not in candidates:
                        candidates[cand['_id']] = {'data': cand, 'wf': 0}
                    candidates[cand['_id']]['wf'] += 1
        # end of the query part
        # if there is not a match:
        #   call getFillerDensityGoogle(filler) i.e. filler density module see if there's a result
        if len(candidates) == 0:
            googleResult = fDM.getFillerDensityGoogle(rptname)
            stdname, clean_result = googleResult
            if len(stdname) == 0: # if no standard name found on google, call removeDescription and removeNano in fillerDensityModule
                stdname = fDM.removeDescription(fDM.removeNano(rptname)).lower().capitalize()
            #   if there is a result: returns (stdname, clean_result)
            #       log it as warning
            #       save it to ChemProps as
            #       {'_id': stdname, '_density': float(clean_result), '_alias': [reported name]}
            if type(clean_result) != int:
                logging.warn("Filler density module find stdname '%s' and density '%s' for reported ChemicalName '%s'. Please verify." %(stdname, clean_result, rptname))
                clean_result = float(clean_result)
                # double check the existance of stdname in ChemProps.filler
                if self.cp.filler.find({"_id": stdname}).count() == 0:
                    # if stdname does not exist in ChemProps, insert it
                    resultDict = {"_id": stdname,
                                  "_density": clean_result,
                                  "_alias": [rptname],
                                  "_boc": [self.bagOfChar(stdname), self.bagOfChar(rptname)]
                                 }
                    self.cp.filler.insert(resultDict)
                    logging.warn("Filler density module results are inserted into ChemProps with '_id': '%s'. Please double check!" %(stdname))
                    return resultDict
                else:
                    # if stdname exists in ChemProps, i.e. we fail to identify the equivalence of rptname and stdname, just add rptname to alias will be enough
                    self.cp.filler.update(
                        {"_id": stdname},
                        {"$addToSet": {"_alias": rptname, "_boc": self.bagOfChar(rptname)}}
                    )
                    logging.warn("Apply $addToSet with value '%s' to _alias and corresponding boc to _boc of the filler with _id: '%s' in ChemProps." %(rptname, stdname))
                    return self.cp.filler.find({"_id": stdname})[0]
            #   if there is not a result: return (stdname, -1)
            #       insert _inputname, _rawname, _nmid[] to unknowns.filler
            else:
                # prepare the insertion dict
                ukdict = {'_nmid': [],
                          '_inputname': stdname,
                          '_rawname': []}
                ukdict['_nmid'].append(self.nmid)
                ukdict['_rawname'].append(rptname)
                # see if the _inputname is already in unknowns.filler
                if self.cp.ukfiller.find({"_inputname": {'$regex': '%s' %(stdname), '$options': 'i'}}).count() == 0: # if not exist, create the document
                    # insert it directly
                    self.cp.ukfiller.insert(ukdict)
                    logging.info("Insert unknown filler with _inputname: '%s' to unknowns." %(stdname))
                else:
                    # update the difference
                    ukdata = self.cp.ukfiller.find({"_inputname": {'$regex': '%s' %(stdname), '$options': 'i'}})[0]
                    if not self.lowerIn(rptname, ukdata['_rawname']):
                        self.cp.ukfiller.update(
                            {"_inputname": {'$regex': '%s' %(stdname), '$options': 'i'}},
                            {"$addToSet": { "_rawname": rptname}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _rawname of the filler with _inputname: '%s' in unknowns."
                                 %(rptname, stdname)
                                )
                    if not self.lowerIn(self.nmid, ukdata['_nmid']):
                        self.cp.ukfiller.update(
                            {"_inputname": {'$regex': '%s' %(stdname), '$options': 'i'}},
                            {"$addToSet": { "_nmid": self.nmid}}
                        )
                        logging.info("Apply $addToSet with value '%s' to _nmid of the filler with _inputname: '%s' in unknowns."
                                 %(self.nmid, stdname)
                                )
                    else:
                        logging.warn("The document with _inputname: '%s' has duplicate in _nmid! Check '%s'!" %(stdname, self.nmid))
        # if there is a match:
        #   {'StandardName': _id, 'density': _density}
        #   for multiple matches, return the one with highest cummulated wf
        else:
            # find the candidate with the highest wf
            wf_high = 0
            cand_high = []
            for cand in candidates:
                if candidates[cand]['wf'] > wf_high:
                    wf_high = candidates[cand]['wf']
                    cand_high = [candidates[cand]]
                elif candidates[cand]['wf'] == wf_high:
                    cand_high.append(candidates[cand])
            # always return the first cand_high, but log if there's more than one cand
            if len(cand_high) > 1:
                tieWarning = "For the search package '%s', multiple winning matches found. Weighting factors tie at %d. They are:" %(str(keywords), wf_high)
                for candidate in cand_high:
                    tieWarning += "\n\t%s" %(candidate['data']['_id'])
                logging.warn(tieWarning)
            # warn admin if wf_high is no bigger than 2
            if wf_high <= 2:
                logging.warn("Careful! Low weighting factor %d for nmid '%s'! Mapped _inputname: '%s' with _stdname: '%s'." %(wf_high, self.nmid, rptname, cand_high[0]['data']['_id']))
            return cand_high[0]['data']       
        # otherwise, return None
        return None

    # this function querys collection.field array for the collections that
    # contain all of the alphabetic words in the query
    # input:
    #   query ------- a string for query, e.g. "Poly[(polytetrahydrofuran 1000)-alt-(4,4'-diphenylmethane diisocynate)]"
    #   field ------- the field to conduct search, e.g. "_tradenames"
    #   collection -- the collection that contains the field, e.g. "self.cp.polymer"
    # output:
    #   candidates -- a list of collection dicts that meet the criteria
    def containAllWords(self, query, field, collection):
        pattern = re.compile('[^a-zA-Z]', re.UNICODE)
        query = pattern.sub(' ', query)
        words = query.split()
        nWords = len(words)
        ids = dict()
        for word in words:
            if len(word) == 1: # skip the single characters
                nWords -= 1
                continue
            for result in collection.find({field: {'$regex': '\\b%s\\b' %(word), '$options':'i'}}):
                if result['_id'] not in ids:
                    ids[result['_id']] = {'data': result, 'freq': 0}
                ids[result['_id']]['freq'] += 1
        # check if the length of words equals any of the freq in ids
        output = []
        for cand in ids:
            if ids[cand]['freq'] == nWords:
                output.append(ids[cand]['data'])
        return output

    # this function checks whether the lower case of the given string myStr is
    # contained in a given list myList
    def lowerIn(self, myStr, myList):
        for s in myList:
            if myStr.lower() == s.lower():
                return True
        return False