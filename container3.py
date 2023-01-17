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
TODO
- do not reload rawdatatable if not needed
- why modification_position always present
- develop filre result addition and subtraction
"""


class IdentifiedData:
    """
    contain data about identifications present in raw files.
    sourcefile  - path to source mzml file. Data itself are not stored, because of size
    resultcontainer -  is a dict of results (named by user) in a form {"Name": Result object}. used to
                       access individual Result objects."""
    def __init__(self):
        self.sourcefile = None
        self.resultcontainer = {}

    def add(self,result,resultname):
        """Add result, check if the result is for the same file (critical!)
           Arguments:
               result - container result of known type. Result object
               resultname - name to further identify the result. string
        """ 
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
            
    def write_identified(self,resultnames):
        """extracts identified spectra for all results in resultnames list
        writes mzml file with all MS1 spectra and identified MS2 spectra
        critically requires correct spectrum index. Writes a new mzML file with
        name suffix _identified.mzML. Returns file name.
        Arguments:
        resultnames - names of results for which to write identified spectra. list"""
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
        outfilename = ".".join(self.sourcefile.split(".")[:-1])+"_identified.mzML"
        mzml.MzMLFile().store(outfilename, data)
                    
        return outfilename
    
    def write_unidentified(self,resultnames):
        """extracts not identified spectra for all results in resultnames list
        writes mzml file with all MS1 spectra and identified MS2 spectra
        critically requires correct spectrum index. Writes a new mzML file with
        name suffix _not_identified.mzML. Returns file name.
        Arguments:
        resultnames - result names which will be used to remove spectra from all spectra. list"""
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
        outfilename = ".".join(self.sourcefile.split(".")[:-1])+"_not_identified.mzML"
        mzml.MzMLFile().store(outfilename, data)
                    
        return outfilename
    
    def check_type(self,dictionary,unmodified,modified,accuracy):
        """return kind of conflict between identifications of the same spectrum. Check conflict type for all pairs in a set,
            and mark kinds of conflicts that are present. Generally used for write_spectra_result_table, but can be used for anything
            provided proper argument creation.
            Kinds of conflicts:
            mass - calculated mass difference more than accuracy
            sequence - different unmodified sequence. I and L treted as identical
            isoleucine - sequence identical except leucine/isoleucine
            modification_num - number of modifications is different
            modification_pos - modifications are the same, but positions differ
            Arguments:
            dictionary - line in a form of dictionary {"column name":"column value"}. Dictionary
            unmodified - dictionary keys for unmodified sequences. list
            modified - dictionary keys modified sequences. list
            accuracy - maximum allowed difference in ppm. Roughly 2x mass accuracy setting used for searches. float
        """
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
        for item in pairs_mod:
            mass_a = mzml.AASequence.fromString(item[0]).getMonoWeight()
            mass_b = mzml.AASequence.fromString(item[1]).getMonoWeight()
            error = abs(((mass_a - mass_b)/mass_a)*1000000)
            if error > accuracy:
                mass = True
        if not sequence:
            pairs_mod_extracted = [[re.findall("\(UniMod:[0-9]+\)",x[0]),re.findall("\(UniMod:[0-9]+\)",x[1])] for x in pairs_mod]
            pairs_mod_extracted = [[tuple(y) for y in x ] for x in pairs_mod_extracted]
            modification_num = any([True for x in pairs_mod_extracted if x[0] != x[1]])
        if not (sequence or modification_num):
            modification_pos = any([True for x in pairs_mod if x[0] != x[1]])
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
        """Attempt to write a summary table of a mzml file and particular set of results
        For every spectrum index, peptide identifications are reported as well as kind of conflicts between them.
        Arguments:
        resultnames - set of results to use. list
        check_type - whether to check co0nflict type. bool, default True"""
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
                flag = True
                if self.resultcontainer[result].datatable is None:
                    flag = False
                if flag and scan in self.resultcontainer[result].datatable["SpectrumIndex"].unique():
                    dataline = self.resultcontainer[result].datatable[self.resultcontainer[result].datatable["SpectrumIndex"]==scan]
                    line["PeptideSequence"+"_"+result] = dataline["PeptideSequence"].iloc[0]
                    pepseqs.append("PeptideSequence"+"_"+result)
                    line["ModifiedPeptideSequence"+"_"+result] = dataline["ModifiedPeptideSequence"].iloc[0]
                    modpepseqs.append("ModifiedPeptideSequence"+"_"+result)
                    seqs.append(dataline["ModifiedPeptideSequence"].iloc[0])
                    if check_type:
                        line["ConflictTypes"] = self.check_type(line,pepseqs,modpepseqs,accuracy = 200000)
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
                
                    
            
class Result:
    """Unified (as far as possible) identified result. 
    This class contain methods that can be applied to any Result type"""
    def __init__(self,moddict,mzmlfile,tsvoutput,name):
        """ highly probably unnecesary here"""
        self.moddict = moddict                                         
        self.datapath = mzmlfile.replace("/","\\")                     
        self.resultpath = tsvoutput                                    
        self.name = name                                               
        self.result_dataframe = None                                   
        self.datatable = None    
        self.rawdatatable = None                                      

    def make_histogram_results(self,what = "PPMerror",limits = None,unique_by= None):
        """draws simple histogram of whatever parameter.
        Arguments
        what - column with data to use. string
        limits - dictionary of form {column:(minimum,maximum)} for all columns to be limited.
                As is, it can be used to select particular part of data, not only to limit histogram 
                range. dict
        unique_by - drop duplicated valuse by column, e.g to obtain unique peptides or psms. string
        #TODO add possibility to use raw_datatable"""
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
    
    def make_scatter_results(self,what = ["PPMerror","RTinSeconds"],limits = None,unique_by= None):  
        """Draw simple scatterplot of two columns with numeric data
        arguments
        what - list of columns to use (2 element). 2 element list or tuple of strings
        limits - dictionary of form {column:(minimum,maximum)} for all columns to be limited.  
                As is, it can be used to select particular part of data, not only to limit histogram 
                range. dict
        unique_by - drop duplicated valuse by column, e.g to obtain unique peptides or psms. string
        #TODO add possibility to use raw_datatable"""
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
    
    def load_mzmldata_fromfile(self):
        """Reads some data from mzml file into tablular format.
            -RT
            -M/Z
            -charge
            -spectrum index
            -mass (with hydrogen atoms!)
            -identfied status for spectrum
        Saves the contents in python pickle (primitive)
        This function read raw file. Normally, pickled dataframe will be used for faster execution time
        Returns dataframe and sets self.rawdatatable  as this dataframe"""
        if self.datatable is not None:
            identifiedscans = self.datatable["SpectrumIndex"].unique()
        else:
            identifiedscans = []                                         #in case 0 ids for file
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
        return self.rawdatatable
    
    def load_mzmldata(self):
        """Checks if there is rawdatatable is present in pickled form. If yes, loads it.
        If not, run load_mzml_fromfile and read actual raw file."""
        filename = self.datapath.split("\\")[-1]
        if path.exists(self.name+"_"+filename+"_data_raw_table.pickle"):
            self.rawdatatable = pickle.load(open(self.name+"_"+filename+"_data_raw_table.pickle","rb"))
            return self.rawdatatable
        else:
            self.load_mzmldata_fromfile()
            
    def make_histogram_spectralfile(self,what = "Charge",limits = None):
        """TODO"""
        return 0
    
    def retrieve_sequences(self):
        """simply get a list of peptide sequences, do not differentiate differently modified sequences"""
        sequences = self.datatable["PeptideSequence"].unique()
        sequences = list(sequences)
        
        return sequences
    
    def retrieve_sequences_modified(self):
        """simply get a list of all sequences including modifications"""
        sequences_m = self.datatable["ModifiedPeptideSequence"].unique()
        sequences_m = list(sequences_m)
        
        return sequences_m
    
    
class FraggerResult(Result):
    """---------------loader for MSFragger psm.tsv data file-------------------------------------

    Unified reader for Fragpipe psm.tsv as Fragpipe result. Attempt to unify data as far as possible.
    class attributes
    -moddict dictionary to unify modification names
    -datapath raw source file path
    -resultpath tsv output of search engine file path
    -name result name chosen by user
    -result_dataframe - actual contents of results of search engine
    -datatable results in format as unified as possible
    -rawdatatable - summary table for mzml file"""
    
    def __init__(self,moddict,mzmlfile,tsvoutput,name):
        """Arguments to create FraggerResult instance
            moddict - dictionary to unify modification names
            mzmlfile -  raw source file path
            tsvoutput - tsv output of search engine file path
            name - result name used to identify it. string"""
        self.moddict = moddict
        self.datapath = mzmlfile.replace("/","\\")                     #unify path slashes
        self.resultpath = tsvoutput.replace("/","\\") 
        self.name = name
        self.result_dataframe = None
        self.datatable = None
        self.rawdatatable = None
        self.result_for_file_exists = True  # flag false to not attempt to transform empty results

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
        if self.result_for_file_exists:
            self.transform_tableresult()
            
    def read_result_single(self):
        """Loads psm.tsv file that contain data only for single mzml file
        ####TODO check if this is really single-file psm.tsv and contain data for right file. 
        Probably useless"""     
        result_dataframe = pd.read_csv(self.resultpath, sep = "\t")
        self.result_dataframe = result_dataframe
        
        return self.result_dataframe
            
    def read_result_multiple(self):
        """Loads results for single file from multiple-file or single-file psm.tsv
        #### TODO more robust method to separate file names. maybe by spectrum?""" 
        def check_spectrum_file_match(specstring,filename):
            specstring = specstring.split('.')
            if filename == specstring[0]:
                print(filename)
                print(specstring)
                returnval= True
            else:
                returnval= False
            return returnval
        result_dataframe = pd.read_csv(self.resultpath, sep = "\t")
        fragger_file_name = "interact-"+".".join(self.datapath.split("\\")[-1].split(".")[:-1])+ ".pep.xml"    #transform filename to spectrum name used in psm.tsv
        #result_dataframe = result_dataframe[result_dataframe["Spectrum File"].str.endswith(fragger_file_name)]
        #print(result_dataframe.shape[0])
        #if result_dataframe.shape[0] <=1:
        filename = self.datapath.split("\\")[-1].rstrip(".mzML")
        result_dataframe['fileset'] = result_dataframe["Spectrum"].apply(check_spectrum_file_match, args = [filename])
        result_dataframe = result_dataframe[result_dataframe['fileset'] == True]
        
        if result_dataframe.shape[0] <=1:    
            print("No data for file " + self.datapath)
            print("in")
            print(self.name)
            self.result_for_file_exists = False
        print(filename)
        print(result_dataframe.shape[0])
        self.result_dataframe = result_dataframe
        
        return self.result_dataframe
    
    def transform_tableresult(self):
        """Transform psm.tsv into more unified form, add columns, original columns remain in place or 
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
        ####TODO add empty feature column, to use for feature based id checking and comparison""" 
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


class PathFinderResult(Result):
    """---------------loader for MSPathfinderT/Promex result psm.tsv data file-------------------------------------
    Unified reader for Informed Proteomics package output files. Attempt to unify data as far as possible.
    class attributes:
    -moddict dictionary to unify modification files
    -datapath raw source file path.
    
    

"""

    def __init__(self,moddict,mzmlfile,tsvoutput,name,fdr):
        """"Arguments:
            -moddict dictionary to unify modification files
            -mzmlfile raw source file path. .ms1ft file(Promex output) is also used, must have the same name except extension, 
             and be located in the same folder.
            -tsvoutput tsv output of search engine file path. Requires Target-Decoy analysis output!
            -name result name chosen by user
            -fdr to cutoff Target-Decoy analysed file. Usually should be 0.005 for reasonable results"""
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
    



        
        
    
    
        
    



