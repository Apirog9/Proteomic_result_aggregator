# -*- coding: utf-8 -*-


import pyopenms as mzml
import pandas as pd
import plotnine as plt
import numpy as np
import pickle
from os import path
import re
import copy
import itertools as it
import re


"""
Created on Thu Jun 30 16:52:11 2022

@author: APirog


resultpath = r"D:\data\peptidomics_psmfile_summarization\06_2022_fractionations\fractions_60min\fragger_7-50\mfagger_FA"
rawfilepath = r"D:\data\peptidomics_psmfile_summarization\06_2022_fractionations\fractions_60min\mzml"
dataname = "FA/ACN HLB Fractionation"
samples = {"10":"FA_10_test.mzml","20":"FA_20_test.mzml","30":"FA_30_test.mzml","40":"FA_40_test.mzml","60":"FA_60_test.mzml","90":"FA_90_test.mzml"}
mzmltoomit = "badmzmls"
filtering = True
namelist = "namespace"
moddict = {"[147]":"(UniMod:35)","n[43]":"(UniMod:1)"}
"""


"""
contain data about identifications present in raw files.
sourcefile is a mzml file
resultcontainer is a dict of results (named by user)

"""
class IdentifiedData:
    def __init__(self):
        self.sourcefile = None
        self.resultcontainer = {}
        
    """
    Add result, check if the result is for the same file (critical!)
    
    """    
    def add(self,result,resultname):
        if self.sourcefile == None:
            self.sourcefile = result.datapath
            self.resultcontainer[resultname] = result
        elif self.sourcefile != result.datapath:         #check if the results you're adding are from the same file
            print ("Wrong results, wrong file path")
            print("Previous")
            print(self.sourcefile)
            print("Adding")
            print(result.datapath)
        else:
            self.sourcefile = result.datapath
            self.resultcontainer[resultname] = result
            
    """
    extracts identified spectra for all results in resultnames list
    writes mzml file with all MS1 spectra and identified MS2 spectra
    critically requires correct spectrum index     
    
    """
    def write_identified(self,resultnames):
        all_identified_scans = []
        for result in resultnames:
            result = self.resultcontainer[result]
            if type(result.rawdatatable) == None:
                result.load_mzmldata()
            else:
                pass
            identified = result.rawdatatable[result.rawdatatable["Identified"] == True]
            identified_scans = list(identified["Scan"])
            all_identified_scans = all_identified_scans + identified_scans
        all_identified_scans = set(all_identified_scans)
        print("identified:")
        print(len(all_identified_scans))
        data = mzml.MSExperiment()
        mzml.MzMLFile().load(self.sourcefile, data)
        identifiedspectra = []
        for spectrum in data.getSpectra():
            if spectrum.getMSLevel ==1:
                identifiedspectra.append(spectrum)
            else:
                if int(spectrum.getNativeID().split("=")[-1]) in all_identified_scans:
                    identifiedspectra.append(spectrum)
                else:
                    pass
        data.setSpectra(identifiedspectra)
        mzml.MzMLFile().store(".".join(self.sourcefile.split(".")[:-1])+"_identified.mzML", data)
                    
        return 0
    
    """
    extracts not identified spectra for all results in resultnames list
    writes mzml file with all MS1 spectra and identified MS2 spectra
    critically requires correct spectrum index     
    
    """
    def write_unidentified(self,resultnames):
        all_identified_scans = []
        for result in resultnames:
            result = self.resultcontainer[result]
            if type(result.rawdatatable) == None:
                result.load_mzmldata()
            else:
                pass
            identified = result.rawdatatable[result.rawdatatable["Identified"] == True]
            identified_scans = list(identified["Scan"])
            all_identified_scans = all_identified_scans + identified_scans
        all_identified_scans = set(all_identified_scans)
        print("not identified:")
        data = mzml.MSExperiment()
        mzml.MzMLFile().load(self.sourcefile, data)
        unidentifiedspectra = []
        for spectrum in data.getSpectra():
            if spectrum.getMSLevel ==1:
                unidentifiedspectra.append(spectrum)
            else:
                if int(spectrum.getNativeID().split("=")[-1]) not in all_identified_scans:

                    unidentifiedspectra.append(spectrum)
                else:
                    pass
        print(len(unidentifiedspectra))
        data.setSpectra(unidentifiedspectra)
        mzml.MzMLFile().store(".".join(self.sourcefile.split(".")[:-1])+"_not_identified.mzML", data)
                    
        return 0
    
    """
    slow and unclean
    add kind of conflict - mod location,sequence,mass
    """
    
    def check_type(self,dictionary,unmodified,modified,accuracy):
        mass = False
        sequence = False
        isoleucine = False
        modification_num = False
        modification_pos = False
        #remove No_ID
        modified = [dictionary[x] for x in modified]
        unmodified = [dictionary[x] for x in unmodified]
        modified = [x for x in modified if x != "No_ID"]
        unmodified = [x for x in unmodified if x != "No_ID"]
        #generate pairs
        pairs_mod = list(it.combinations(modified, 2))
        pairs_unmod = list(it.combinations(modified, 2))
        sequence = any([True for x in pairs_unmod if x[0] != x[1]])
        if sequence:
            pairs_unmod = [(x[0].replace("I","L"),x[1].replace("I","L")) for x in pairs_unmod]
            check = any([True for x in pairs_unmod if x[0] != x[1]])
            if not check:
                isoleucine = True
                sequence = False
        #here probably need to use re
        #here probably need to use re
        for item in pairs_mod:
            
            mass_a = mzml.AASequence.fromString(item[0]).getMonoWeight()
            mass_b = mzml.AASequence.fromString(item[1]).getMonoWeight()
            error = abs(((mass_a - mass_b)/mass_a)*1000000)
            if item[0] != item[1]:
                print(error)
                print(mass_a)
                print(mass_b)
                print(item[0])
                print(item[1])
            if error > accuracy:
                mass = True
        if not sequence:
            pairs_mod_extracted = [[re.findall("\(UniMod:[0-9]+\)",x[0]),re.findall("\(UniMod:[0-9]+\)",x[1])] for x in pairs_mod]
            pairs_mod_extracted = [[tuple(y) for y in x ] for x in pairs_mod_extracted]
            modification_num = any([True for x in pairs_mod_extracted if x[0] != x[1]])
        if not (sequence or modification_num):
            modification_pos = any([True for x in pairs_mod if x[0] != x[1]])
        if modification_pos:
            for item in pairs_mod:
                print(item)
        if [mass,sequence,isoleucine,modification_num,modification_pos] == [False, False, True, False,False] :
            for item in pairs_mod:
                print(item)
            print("dupa")
        strings = ["Mass","Sequence","Isoleucine/Leucine","Modification_Number","Modification_Position"]
        values = [mass,sequence,isoleucine,modification_num,modification_pos]
        i =  0
        returnval = []
        for item in values:
            if item:
                returnval.append(strings[i])
            i+=1
        if returnval == []:
            returnval = "None"
        else:
            returnval = "_".join(returnval)
        
        return returnval
        
        
        
    def write_spectra_result_table(self,resultnames,check_type=True):
        all_identified_scans = []
        for result in resultnames:
            result = self.resultcontainer[result]
            if result.rawdatatable == None:
                print("loading raw data")
                result.load_mzmldata()
            else:
                pass
            identified = result.rawdatatable[result.rawdatatable["Identified"] == True]
            identified_scans = list(identified["Scan"])
            all_identified_scans = all_identified_scans + identified_scans
        all_identified_scans = set(all_identified_scans)
        all_identified_scans = list(all_identified_scans)
        all_identified_scans.sort()
        dataframecontent = []

        for scan in all_identified_scans:
            line = {}
            line["SpectrumIndex"] = scan
            seqs = []
            pepseqs = []
            modpepseqs = []
            for result in resultnames:
                if scan in self.resultcontainer[result].datatable["SpectrumIndex"].unique():
                    dataline = self.resultcontainer[result].datatable[self.resultcontainer[result].datatable["SpectrumIndex"]==scan]
                    line["PeptideSequence"+"_"+result] = dataline["PeptideSequence"].iloc[0]
                    pepseqs.append("PeptideSequence"+"_"+result)
                    line["ModifiedPeptideSequence"+"_"+result] = dataline["ModifiedPeptideSequence"].iloc[0]
                    modpepseqs.append("ModifiedPeptideSequence"+"_"+result)
                    seqs.append(dataline["ModifiedPeptideSequence"].iloc[0])
                    if check_type:
                        line["ConflictTypes"] = self.check_type(line,pepseqs,modpepseqs,accuracy = 15)
                else:
                    line["PeptideSequence"+"_"+result] = "No_ID"
                    line["ModifiedPeptideSequence"+"_"+result] = "No_ID"
            if len(set(seqs)) == 1:
                line["Conflict"] = "No"
            else:
                line["Conflict"] = "Yes"
            dataframecontent.append(line)
        dataframe = pd.DataFrame(dataframecontent)
        
        self.mzml_summary = dataframe
        return dataframe
                
                    
            
                                    
        
        

"""      
Unified (as far as possible) identified result. 
-moddict dictionary to unify modification files
-mzmlfile raw source file path
-tsvoutput tsv output of search engine file path
-name result name chosen by user
-result_dataframe - actual contents of results of search engine
-datatable results in format as unified as possible
-rawdatatable unified table of raw data useful for plotting

"""  
        
class Result:
    def __init__(self,moddict,mzmlfile,tsvoutput,name):
        self.moddict = moddict                                         
        self.datapath = mzmlfile.replace("/","\\")                     
        self.resultpath = tsvoutput                                    
        self.name = name                                               
        self.result_dataframe = None                                   
        self.datatable = None    
        self.rawdatatable = None                                      

    """
    draws simple histogram 
    arguments
    what - column to use
    limits - dictionary of form {column:(minimum,maximum)} for all columns to be limited
    unique_by - drop duplicated valuse by column, e.g to obtain unique peptides or psms
    #TODO add possibility to use raw_datatable
    """
    def make_histogram_results(self,what = "PPMerror",limits = None,unique_by= None):
        if unique_by:
            newdatatable = self.datatable.drop_duplicates(subset = unique_by)
        else:
            newdatatable = copy.copy(self.datatable)
        if limits:
            for limitcolumn in limits.keys():
                newdatatable = newdatatable[(newdatatable[limitcolumn]>=limits[limitcolumn][0]) & (newdatatable[limitcolumn]<=limits[limitcolumn][1])]
        plot = plt.ggplot(newdatatable) + plt.aes(x=what) + plt.geom_histogram() + plt.ggtitle(self.name)
        plot.draw(show =True)
        
        return 0
    
    """
    arguments
    what - list of columns to use (2 element)
    limits - dictionary of form {column:(minimum,maximum)} for all columns to be limited
    unique_by - drop duplicated valuse by column, e.g to obtain unique peptides or psms
    #TODO add possibility to use raw_datatable
    """
    def make_scatter_results(self,what = ["PPMerror","RTinSeconds"],limits = None,unique_by= None):
        if unique_by:
            newdatatable = self.datatable.drop_duplicates(subset = unique_by)
        else:
            newdatatable = copy.copy(self.datatable)
        if limits:
            for limitcolumn in limits.keys():
                newdatatable = newdatatable[(newdatatable[limitcolumn]>=limits[limitcolumn][0]) & (newdatatable[limitcolumn]<=limits[limitcolumn][1])]       
        plot = plt.ggplot(newdatatable) + plt.aes(x=what[0],y=what[1]) + plt.geom_point() + plt.ggtitle(self.name)
        plot.draw(show =True)
        
        return 0
    
    """
    Reads some data from mzml file into tablular format.
    (RT,MZ,charge of all MS2 spectra, as well as identification status)
    Save the contents in python pickle (primitive)
    This function read raw file.
    """
    def load_mzmldata_fromfile(self):
        identifiedscans = self.datatable["SpectrumIndex"].unique()
        data = mzml.MSExperiment()
        mzml.MzMLFile().load(self.datapath, data)
        dataframe = []
        for spectrum in data.getSpectra():
            RT = spectrum.getRT()
            scan = int(spectrum.getNativeID().split("=")[-1])
            ion = {}
            ion["Retention"] = RT
            ion["Scan"] = scan
            ion["Identified"] = ion["Scan"] in identifiedscans
            for item in spectrum.getPrecursors():
                charge = item.getCharge()
                mz = item.getMZ()
                mass = mz*charge
                ion["M/Z"] = mz
                ion["Charge"] = charge
                ion["Mass"] = mass
                dataframe.append(ion)
        dataframe = pd.DataFrame(dataframe)
        self.rawdatatable = dataframe
        filename = self.datapath.split("\\")[-1]
        with open(self.name+"_"+filename+"_data_raw_table.pickle","wb") as outpickle:
            pickle.dump(self.rawdatatable,outpickle)
    
    """
    Checks if there is raw_datatable present in pickled form. If yes, loads it.
    If no, run load_mzml_fromfile and read actual raw file
    
    """
    def load_mzmldata(self):
        filename = self.datapath.split("\\")[-1]
        if path.exists(self.name+"_"+filename+"_data_raw_table.pickle"):
            self.rawdatatable = pickle.load(open(self.name+"_"+filename+"_data_raw_table.pickle","rb"))
        else:
            self.load_mzmldata_fromfile()
            
    def make_histogram_spectralfile(self,what = "Charge",limits = None):
        
        return 0
    
    
    """
    simply get a list of peptide sequences, do not differentiate differently modified sequences
    
    """
    def retrieve_sequences(self):
        sequences = self.datatable["PeptideSequence"].unique()
        sequences = list(sequences)
        
        return sequences
    
    """
    simply get a list of all sequences including modifications
    
    """
    def retrieve_sequences_modified(self):
        sequences_m = self.datatable["ModifiedPeptideSequence"].unique()
        sequences_m = list(sequences_m)
        
        return sequences_m
    
    

        
"""
---------------loader for MSFragger psm.tsv data file-------------------------------------

Unified Read Fragpipe psm.tsv as Fragpipe result. Attempt to unify data as far as possible
-moddict dictionary to unify modification files
-mzmlfile raw source file path
-tsvoutput tsv output of search engine file path
-name result name chosen by user
-result_dataframe - actual contents of results of search engine
-datatable results in format as unified as possible

"""  
class FraggerResult(Result):
    def __init__(self,moddict,mzmlfile,tsvoutput,name):
        self.moddict = moddict
        self.datapath = mzmlfile.replace("/","\\")                     #unify path slashes
        self.resultpath = tsvoutput.replace("/","\\") 
        self.name = name
        self.result_dataframe = None
        self.datatable = None
        self.rawdatatable = None

        if not path.exists(self.datapath):
            print("No raw data file!")
            print("in")
            print(self.name)
        if not path.exists(tsvoutput):
            print("No result file")
            print("in")
            print(self.name)
        if path.exists(self.datapath) and path.exists(self.resultpath):
            self.read_result_multiple()
            self.transform_tableresult()
            
    """
    loads psm.tsv file that contain data only for single mzml file
    ####TODO check if this is really single-file psm.tsv and contain data for right file
    
    """        
    def read_result_single(self):
        result_dataframe = pd.read_csv(self.resultpath, sep = "\t")
        self.result_dataframe = result_dataframe
        
    """
    Loads results for single file from multiple-file psm.tsv
    #### TODO more robust method to separate file names. maybe by spectrum?
    """     
    def read_result_multiple(self):
        result_dataframe = pd.read_csv(self.resultpath, sep = "\t")
        fragger_file_name = "interact-"+".".join(self.datapath.split("\\")[-1].split(".")[:-1])+ ".pep.xml"    #transform filename to spectrum name used in psm.tsv
        result_dataframe = result_dataframe[result_dataframe["Spectrum File"].str.endswith(fragger_file_name)]
        if result_dataframe.shape[0] <=1:
           print("No data for file " + self.datapath)
           print("in")
           print(self.name)
        self.result_dataframe = result_dataframe
    
    """
    Transform psm.tsv into more unified form, add columns, original columns remain in place or 
    are renamed without content change
    rename "Retention":"RTinSeconds"
    rename "Peptide":"PeptideSequence"
    rename "Peptide Length":"PeptideLength"
    rename "Calibrated Observed Mass":"MeasuredMassCalibrated"
    rename "Observed Mass":"MeasuredMass"
    add "ModifiedPeptideSequence" containing modified sequence in unified format
    add "SpectrumIndex" containing spectrum index as unique spectrum identifier
    add "PPMerror" containing PPM error of m/z before calibration
    add "PPMerrorAfterRecal" containing PPM error of m/z after calibration
    ####TODO add empty feature column, to use for feature based id checking and comparison
    
    """    
    def transform_tableresult(self):
        
        def makemodifiedsequence(serieslike,modifications):
            if not type(serieslike["Modified Peptide"]) == str:
                returnval = serieslike["PeptideSequence"]
            else:
                pepseq = serieslike["Modified Peptide"]
                for before,after in modifications.items():
                    pepseq = pepseq.replace(before,after)
                returnval = pepseq
            return returnval
        
        def makespectrum(titlestring):

            title = titlestring.split(".")
            scannum = title[1].lstrip("0")
            scannum = int(scannum)
            index = scannum -1
            return index
        new_datatable = copy.copy(self.result_dataframe)
        new_datatable.rename({"Retention":"RTinSeconds"},axis=1,inplace=True)
        new_datatable.rename({"Peptide":"PeptideSequence"},axis=1,inplace=True)
        new_datatable.rename({"Peptide Length":"PeptideLength"},axis=1,inplace=True)
        new_datatable.rename({"Calibrated Observed Mass":"MeasuredMassCalibrated"},axis=1,inplace=True)
        new_datatable.rename({"Observed Mass":"MeasuredMass"},axis=1,inplace=True)
        new_datatable["ModifiedPeptideSequence"] = new_datatable[["PeptideSequence","Modified Peptide"]].apply(makemodifiedsequence,axis=1,args = [self.moddict])
        new_datatable["SpectrumIndex"] = new_datatable["Spectrum"].apply(makespectrum)
        new_datatable["PPMerror"] = ((new_datatable["Observed M/Z"] - new_datatable["Calculated M/Z"])/new_datatable["Observed M/Z"])*1000000
        new_datatable["PPMerrorAfterRecal"] = ((new_datatable["Calibrated Observed M/Z"] - new_datatable["Calculated M/Z"])/new_datatable["Calibrated Observed M/Z"])*1000000
        self.datatable = new_datatable
    
    def read_pepxml_contents(self):
        
        return 0


"""
---------------loader for MSPathfinderT/Promex result psm.tsv data file-------------------------------------
requires for now
-tsv file with Target-Decoy analysis
-mzml file
-ms1ft file located in the same folder as mzml file
Will check for existence of all of them.

Contain:

-moddict dictionary to unify modification files
-mzmlfile raw source file path
-tsvoutput tsv output of search engine file path
-name result name chosen by user
-fdr to cutoff Target-Decoy analysed file. Usually should be 0.005 for reasonable results
-result_dataframe - actual contents of results of search engine
-datatable results in format as unified as possible

"""
class PathFinderResult(Result):
    def __init__(self,moddict,mzmlfile,tsvoutput,name,fdr):
        self.moddict = moddict
        self.datapath = mzmlfile.replace("/","\\")                     #unify path slashes
        self.resultpath = tsvoutput.replace("/","\\") 
        self.ms1ftpath = self.datapath.replace(".mzML",".ms1ft")
        self.name = name
        self.result_dataframe = None
        self.fdr = fdr
        self.datatable = None
        self.rawdatatable = None

        if not path.exists(self.datapath):
            print("No raw data file!")
        if not path.exists(self.datapath.replace(".mzml",".ms1ft")):
            print("No msft  data file!")
        if not path.exists(tsvoutput):
            print("No result file")
        if path.exists(self.datapath) and path.exists(self.resultpath):
            self.read_result()
            self.read_ms1ft()
            self.transform_tableresult()
    
    """
    Reads contents of tsv output file, filters by fdr, removes reverse hits
    
    """
    def read_result(self):
        result_dataframe = pd.read_csv(self.resultpath, sep = "\t")
        result_dataframe = result_dataframe[result_dataframe["PepQValue"] <= self.fdr]
        result_dataframe = result_dataframe[~result_dataframe["ProteinName"].str.startswith("XXX")]
        self.result_dataframe = result_dataframe
        
    """
        Reads contents of ms1ft file to dataframe
        
    """
    def read_ms1ft(self):
        ms1ft_dataframe = pd.read_csv(self.ms1ftpath, sep = "\t")
        self.ms1ft_dataframe = ms1ft_dataframe
        
        
    """
        Transform tsv result into more unified form, add columns, original columns remain in place or 
        are renamed without content change
        rename "Sequence":"PeptideSequence"
        add "RTinSeconds" read data from ms1ft, taking middle of elution of feature. Require filling missing data from raw data file for low-mass ions
        add "PeptideLength" simply calculate from "Start" and "End" columns
        add "ObservedMass" read data from ms1ft. Require filling missing data from raw data file for low-mass ions.
        add "CalculatedMass" calculate from composition
        add "PPMerror" calculate from mass, not m/z!!!! 
        add "ModifiedPeptideSequence" containing modified sequence in unified format
        
        
        
        ####TODO fill data from raw file (probably new function is better)

    """   
    def transform_tableresult(self):
        def get_rt(feature,self):
            if feature == 0:
                rt =0
            else:
                line = self.ms1ft_dataframe[self.ms1ft_dataframe["FeatureID"] == feature]

                rt = ((float(line['MaxElutionTime'].iloc[0]) + float(line['MinElutionTime'].iloc[0]))/2)*60
            
            return rt
        def get_mass(feature,self):
            if feature == 0:
                mass =0
            else:
                line = self.ms1ft_dataframe[self.ms1ft_dataframe["FeatureID"] == feature]

                mass = float(line['MonoMass'].iloc[0])
            return mass
        def calculate_mass(composition):
            formula = composition.replace("(","").replace(")","").replace(" ","")
            formula = mzml.EmpiricalFormula(formula)
            mass = formula.getMonoWeight()
            
            return mass
        def makemodifiedsequence(serieslike,moddict):
            
            if not type(serieslike["Modifications"]) == str :
                sequence = str(serieslike["PeptideSequence"])
            else:
                sequence = str(serieslike["PeptideSequence"])
                modifications = serieslike["Modifications"].split(",")
                locations = []
                
                for item in modifications:
                    mod =item.split(" ")[0]
                    loc =int(item.split(" ")[1])
                    mod = moddict[mod]
                    locations.append([mod,loc])
                
                i=0
                while i < len(locations):
                    
                    mod_current = locations[i][0]
                    f = lambda x : [x[0],x[1]+ len(mod_current)]
                    sequence = sequence[:locations[i][1]] + mod_current + sequence[locations[i][1]:]
                    locations = [f(x) for x in locations]
                    i+=1
                
                
                    
                
            return sequence
            
            
        
        new_datatable = copy.copy(self.result_dataframe)
        new_datatable["RTinSeconds"] = new_datatable["Ms1Features"].apply(get_rt, args = [self])   # the easiest rt, for not all of features
        new_datatable.rename({"Sequence":"PeptideSequence"},axis=1,inplace=True)
        new_datatable["PeptideLength"] = new_datatable[["Start","End"]].apply(lambda x: int(x["End"]) - int(x["Start"]), axis=1)
        new_datatable["MeasuredMass"] = new_datatable["Ms1Features"].apply(get_mass, args = [self])
        new_datatable["CalculatedMass"] = new_datatable["Composition"].apply(calculate_mass)  #check
        new_datatable["PPMerror"] = ((new_datatable["MeasuredMass"] - new_datatable["CalculatedMass"])/new_datatable["MeasuredMass"])*1000000
        new_datatable["ModifiedPeptideSequence"] = new_datatable[["PeptideSequence","Modifications"]].apply(makemodifiedsequence,axis=1,args = [self.moddict])
        new_datatable["SpectrumIndex"] = new_datatable["Scan"] - 1
        self.datatable = new_datatable
        
    def fill_data_from_raw(self):  # fill or replace missing rts and masses
        
        return 0
    def read_mzid_contents(self):
        return 0
    



        
        
    
    
        
    



