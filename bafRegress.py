#!/usr/bin/env python

import os
import gzip
import array
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import struct
import sys
from signal import signal, SIGPIPE, SIGINT, SIG_DFL 
from subprocess import call

__version__ = '0.9.3'

signal(SIGPIPE,SIG_DFL) #quiet broken pipe exists
signal(SIGINT,SIG_DFL) #quiet ctrl-c exits

class FilterBox:
    def __init__(self):
        self.sampleFilters = list()
        self.markerFilters = list()

    def addMarkerFilter(self, filter):
        self.markerFilters.append(filter)

    def addSampleFilter(self, filter):
        self.sampleFilters.append(filter)

    def shouldIncludeMarker(self, markername):
        shouldInclude = True
        for filter in self.markerFilters:
            shouldInclude = shouldInclude and filter.shouldIncludeMarker(markername)
        return shouldInclude

    def shouldIncludeSample(self, samplename):
        shouldInclude = True
        for filter in self.sampleFilters:
            shouldInclude = shouldInclude and filter.shouldIncludeSample(samplename)
        return shouldInclude

class OutputManager:
    def __init__(self, prefix, snpcolname="SNP Name", samplecolname="Sample Name"):
        self.out = list()
        self.filters = FilterBox()
        self.prefix = prefix
        self.snpCol = snpcolname
        self.sampleCol = samplecolname
        self.lastsample = ""
        self.markerindex = 0;
        self.samples = [];
        self.markers = [];
        self.currentfile = None

    def open(self):
        for of in self.out:
            of.openwrite(self.prefix)

    def nextinput(self, document):
        self.currentfile = document.filename
        for of in self.out:
            if "nextinput" in dir(of):
                of.nextinput(document)

    def process(self,data):
        marker = data[self.snpCol]
        if self.sampleCol == "NONE":
            sample = self.currentfile
        else:
            sample = data[self.sampleCol] 
        if self.shouldIncludeSample(sample):
            if self.shouldIncludeMarker(marker):
                if self.lastsample != sample:
                    self.markerindex = 0;
                    self.samples.append(sample)
                    self.lastsample = sample
                if self.markerindex < len(self.markers):
                    if self.markers[self.markerindex] != marker:
                        raise(NameError("Expecting marker " + self.markers[self.markerindex] + " found " + marker + ", " + str(self.markerindex)))
                else:
                    self.markers.append(marker)

                for of in self.out:
                    of.process(len(self.samples)-1,self.markerindex,data, self)

                self.markerindex += 1
            else:
                pass #skip marker
        else:
            pass #skip sample

    def close(self):
        for of in self.out:
            of.close()

    def addOutputFile(self, file):
        self.out.append(file)

    def addMarkerFilter(self, filter):
        self.filters.addMarkerFilter(filter)

    def addSampleFilter(self, filter):
        self.filters.addSampleFilter(filter)

    def shouldIncludeMarker(self, markername):
        return self.filters.shouldIncludeMarker(markername)

    def shouldIncludeSample(self, samplename):
        return self.filters.shouldIncludeSample(samplename)

class InputManager:
    def __init__(self):
        self.inputs = list()
        self.samples = list()
        self.markers = list()
        self.includeSampleName = False
        self.includeMarkerName = False

    def open(self, prefix):
        for file in self.inputs:
            file.openread(prefix)

    def header(self):
        values = list()
        if self.includeSampleName:
            values.add("SAMPLE")
        if self.includeMarkerName:
            values.add("MARKER")
        for file in self.inputs:
            values.extend(file.headers)
        return values
    

    def read(self, sampleindex, markerindex):
        values = list()
        if self.includeSampleName:
            values.add(self.samples[sampleindex])
        if self.includeMarkerName:
            values.add(self.markers[markerindex])
        for file in self.inputs:
            values.extend(file.read(sampleindex, markerindex, self))
        return values
    
    def close(self):
        for file in self.inputs:
            file.close()

    def addInputFile(self, file):
        self.inputs.append(file)


class BasicFilter:
    def __init__(self, filtertype, default, items=None):
        if filtertype == "sample":
            self.samplefilter = True
            self.markerfilter = False
        elif filtertype == "marker":
            self.markerfilter = True
            self.samplefilter = False
        else:
            raise NameError("Invalid filter type")
        self.default = default
        if items:
            self.items = dict.fromkeys(items, True)
        else:
            self.items = dict()

    def additem(self,val):
        self.items[val] = True

    def load(self, file, col=0, sep=None):
        f = open(file)
        for line in f:
            cols = line.rstrip().split(sep)
            self.items[cols[col]] = True
        f.close()

    def shouldIncludeSample(self, samplename):
        if not self.samplefilter:
            return True
        return samplename in self.items != self.default

    def shouldIncludeMarker(self, markername):
        if not self.markerfilter:
            return True
        return markername in self.items != self.default

    
class ReferenceTextOutputFile:
    def __init__(self, suffix, referenceValue):
        self.suffix = suffix
        self.referenceValue = referenceValue
        assert referenceValue=="sample" or referenceValue=="marker"
        self.highestIndex = -1

    def openwrite(self, prefix):
        self.file = open(prefix + self.suffix, 'w')

    def process(self, sampleindex, markerindex, data, manager):
        if self.referenceValue =="sample":
            index = sampleindex
            value = manager.samples[sampleindex]
        else:
            index = markerindex
            value = manager.markers[markerindex]

        if index > self.highestIndex:
            print(value, file=self.file)
            self.highestIndex = index

    def close(self):
        self.file.close()

class BinaryOutputFile:
    def __init__(self, magicNumber, suffix, cols, arrayType, width, parser):
        self.magicNumber = magicNumber
        self.suffix = suffix
        self.cols = cols
        self.arrayType = arrayType
        self.recordWidth = width
        self.parser = parser
        self.lastSampleIndex = -1
        self.highestMarkerIndex = -1
        self.dumped = True
        self.dataarray = array.array(self.arrayType)
    
    def openwrite(self, prefix):
        self.file = open(prefix + self.suffix, 'wb')
        self.file.write(self.magicNumber)
        self.file.write("\x00"*15)
        self.dumped = True

    def process(self, sampleindex, markerindex, data, manager):
        if sampleindex != self.lastSampleIndex:
            if not self.dumped:
                self.dump()
            self.lastSampleIndex = sampleindex
        newdata =  self.parser([data[x] for x in self.cols])
        assert len(newdata) == self.recordWidth
        if markerindex > self.highestMarkerIndex:
            if markerindex == self.highestMarkerIndex+1:
                self.highestMarkerIndex = markerindex
            else:
                raise(NameError("Missing marker index"))
            map(self.dataarray.append, newdata)
        else:
            for i,v in enumerate(newdata):
                self.dataarray[markerindex*self.recordWidth+i] = v
        self.dumped = False

    def dump(self):
        self.dataarray.tofile(self.file)
        self.dumped = True

    def close(self):
        if not self.dumped:
            self.dump()
        self.file.close()

class Binary2DInputFile:
    def __init__(self, magicNumber, suffix, arrayType, width, nmarkers, headers, decode=None):
        self.magicNumber = magicNumber
        self.suffix = suffix
        self.arrayType = arrayType
        self.width = width
        self.sampleIndex = -1
        self.recsize = nmarkers
        self.headers = headers
        self.decode = decode

    def openread(self, prefix):
        name = prefix+self.suffix
        self.file = open(name, 'rb')
        mn = self.file.read(1)
        if (mn != self.magicNumber):
            raise NameError("Error in file " + name + ", expected magic number" + str(self.magicNumber) +
                " but found " + str(mn))
        self.dataarray = array.array(self.arrayType)

    def read(self, sampleindex, markerindex, manager):
        if (sampleindex != self.sampleIndex):
            self.file.seek(16 + sampleindex * self.recsize * self.dataarray.itemsize)
            self.dataarray.fromfile(self.file, self.recsize)
            self.sampleIndex = sampleindex
        d = [self.dataarray[x] for x in range(markerindex*self.width,markerindex*self.width+self.width)]
        if self.decode:
            d = self.decode(d)
        return d

    def close(self):
        self.file.close()

class Binary1DInputFile:
    def __init__(self, magicNumber, suffix, arrayType, width, nmarkers, headers):
        self.magicNumber = magicNumber
        self.suffix = suffix
        self.arrayType = arrayType
        self.width = width
        self.sampleIndex = -1
        self.recsize = nmarkers
        self.headers = headers

    def openread(self, prefix):
        name = prefix+self.suffix
        self.file = open(name, 'rb')
        mn = self.file.read(1)
        if (mn != self.magicNumber):
            raise NameError("Error in file " + name + ", expected magic number" + str(self.magicNumber) +
                " but found " + str(mn))
        self.file.read(15)
        self.dataarray = array.array(self.arrayType)
        self.dataarray.fromfile(self.file, self.recsize)

    def read(self, sampleindex, markerindex, manager):
        return [self.dataarray[x] for x in range(markerindex*self.width,markerindex*self.width+self.width)]

    def close(self):
        self.file.close()

class PopulationFrequencyFile:
    def __init__(self, magicNumber, suffix, cols, sep=None):
        self.magicNumber = magicNumber
        self.suffix = suffix
        self.cols = cols
        self.highestIndex = -1
        self.mafs = dict()
        self.headers = ("MAF",)
        self.sep = sep

    def load(self, filename):
        data = open(filename)
        #skip header
        line = data.readline()
        for line in data:
            a = line.rstrip().split(self.sep)
            val = a[self.cols[1]]
            if val != "NA" and float(val)>0:
                val = float(val)
                if val>.5:
                    val = 1-val
                self.mafs[a[self.cols[0]]] = val 

    def openwrite(self, prefix):
        self.file = open(prefix + self.suffix, 'wb')
        self.file.write(self.magicNumber)
        self.file.write("\x00"*15)

    def openread(self, prefix):
        pass

    def process(self, sampleindex, markerindex, data, manager):
        if markerindex > self.highestIndex:
            maf = self.mafs[manager.markers[markerindex]]
            self.file.write(struct.pack('d',float(maf)))
            self.highestIndex = markerindex

    def read(self, sampleindex, markerindex, manager): 
        maf = self.mafs.get(manager.markers[markerindex], "NA")
        return (maf,)

    def close(self):
        self.file.close()

    def shouldIncludeMarker(self, markername):
        return markername in self.mafs

class FinalReportStream:
    def __init__(self, options):
        self.isHeaderWritten = False

        self.colsample = options.colsample
        self.colmarker = options.colmarker
        self.colnormx = options.colnormx
        self.colnormy = options.colnormy
        self.colrawx = options.colrawx
        self.colrawy = options.colrawy
        self.colab1 = options.colab1
        self.colab2 = options.colab2
        self.colbaf = options.colbaf

        self.hasbaf = False
        self.hasab = False
        self.hasfreq = False

        self.freq = None

    def openwrite(self, prefix):
        pass

    def nextinput(self, fp):
        self.hasbaf = self.colbaf in fp.colnames
        self.hasab = self.colab1 in fp.colnames and self.colab2 in fp.colnames
        self.hasfreq = not self.freq is None

        if not self.isHeaderWritten:
            header = list()
            if self.hasbaf:
                header.append("BAF")
            if self.hasab:
                header.append("ABGENO")
            if self.hasfreq:
                header.append("MAF")
            print("\t".join(header))
            self.isHeaderWritten = True

    def process(self, sampleindex, markerindex, data, manager):
        cols=list()
        if self.hasbaf:
            cols.append(data[self.colbaf])
        if self.hasab:
            cols.extend(genoParser((data[self.colab1], data[self.colab2])))
        if self.hasfreq:
            cols.extend(self.freq.read(sampleindex, markerindex, manager))
        if len(cols):
            print("\t".join([str(x) for x in cols]))

    def close(self):
        pass

class FinalReport:
    def __init__(self, filename, sep="\t"):
        self.filename=filename
        self.sep = sep 

    def open(self):
        self.file = open(self.filename, 'r')
        self.header = dict()
        #parse header
        line = self.file.readline()
        while(line!="" and not line.startswith("[Data]")):
            line = line.rstrip()
            if line != "[Header]":
                a = line.split(self.sep)
                if len(a)<2:
                    sep = self.guesssep(line)
                    a = line.split(sep)
                if len(a)>1:
                    self.header[a[0]] = a[1]
            line = self.file.readline()
        if not line.startswith("[Data]"):
            raise(NameError("[Data] not found in file"))
        line = self.file.readline().rstrip()
        a = line.split(self.sep)
        if len(a)<3:
            self.sep = self.guesssep(line, 3)
            a = line.split(self.sep)
        if len(a)<3:
            sys.exit("Error: Unable to find at least three columns in header: " + line + "\n")
        self.colnames = a

    def guesssep(self, line, min=1):
        if line.count("\t")>min:
            sep = "\t"
        elif line.count(",")>min:
            sep=","
        else:
            sep=None
        return sep

    def read(self):
        for line in self.file:
            a = line.rstrip().split(self.sep)
            yield(dict(zip(self.colnames, a)))

    def close(self):
        self.file.close()

def bafEncode(vals):
    return [int(float(x)*10000) if x!="NaN" else -1 for x in vals]
def bafDecode(vals):
    return [float(x)/10000 if float(x)>-1 else float('NaN') for x in vals]

def genoParser(alleles):
    aballeles = ''.join(alleles)
    if aballeles == "AA":
        ab=0
    elif aballeles == "AB":
        ab=1
    elif aballeles == "BB":
        ab=2
    elif aballeles == "--":
        ab=3
    else:
        raise NameError("bad genotype:" + aballeles)
    return (int(ab),)

def AddFinalReportColumnOptions(parser, asgroup=True):
    if asgroup:
        target = OptionGroup(parser, "FinalReport", "Names of Final Report Columns")
    else:
        target = parser
    target.add_option("--colsample", action="store", dest="colsample", 
        metavar="col", default="Sample Name", help="sample column [default: %default]")
    target.add_option("--colmarker", action="store", dest="colmarker", 
        metavar="col", default="SNP Name", help="marker column [default: %default]")
    target.add_option("--colnormx", action="store", dest="colnormx", 
        metavar="col", default="X", help="Normalized X intensity column [default: %default]")
    target.add_option("--colnormy", action="store", dest="colnormy", 
        metavar="col", default="Y", help="Normalized Y intensity column [default: %default]")
    target.add_option("--colrawx", action="store", dest="colrawx", 
        metavar="col", default="X Raw", help="Raw X intensity column [default: %default]")
    target.add_option("--colrawy", action="store", dest="colrawy", 
        metavar="col", default="Y Raw", help="Raw Y intensity column [default: %default]")
    target.add_option("--colbaf", action="store", dest="colbaf", 
        metavar="col", default="B Allele Freq", help="B allele frequency column [default: %default]")
    target.add_option("--colab1", action="store", dest="colab1", 
        metavar="col", default="Allele1 - AB", help="AB Allele 1 column [default: %default]")
    target.add_option("--colab2", action="store", dest="colab2", 
        metavar="col", default="Allele2 - AB", help="AB Allele 2 column [default: %default]")
    if asgroup:
        parser.add_option_group(target)

def AddFrequencyOptions(parser, asgroup=True):
    if asgroup:
        target = OptionGroup(parser, "Population Frequency", "File containing population frequency information")
    else:
        target = parser
    target.add_option("--freqfile", action="store", dest="freqfile",
        metavar="FILE", help="Text file with marker name and population frequency")
    target.add_option("--freqcols", action="store", dest="freqcols", default="1,2",
        metavar="COLS", help="Columns for marker/frequency [default %default]")
    target.add_option("--freqsep", action="store", dest="freqsep", default=None,
        metavar="SEP", help="Split FILE using this [default white space]")
    if asgroup:
        parser.add_option_group(target)

def AddSubsetOptions(parser, asgroup=True, markers=True, samples=True):
    if asgroup:
        target = OptionGroup(parser, "Subset","Subset data")
    else:
        target = parser
    if markers:
        target.add_option("--extract", action="store", dest="extract", 
            metavar="EXTRACT", help="Include only these markers")
        target.add_option("--extractcol", action="store", dest="extractcol",type="int", 
            metavar="col", help="Use this column of EXTRACT [default: %default]", default=1)
        target.add_option("--extractsep", action="store", dest="extractsep", 
            metavar="sep", help="Split EXTRACT using this [default: white space]", default=None)
        target.add_option("--exclude", action="store",dest="exclude", 
            metavar="EXCLUDE", help="Include all but these markers")
        target.add_option("--excludecol", action="store",dest="excludecol",type="int", 
            metavar="col", help="Use this column of EXCLUDE [default: %default]", default=1)
        target.add_option("--excludesep", action="store", dest="excludesep", 
            metavar="sep", help="Split EXCLUDE using this [default: white space]", default=None)
    if samples:
        target.add_option("--keep",action="store", dest="keep", 
            metavar="KEEP", help="Include only these individuals")
        target.add_option("--keepcol", action="store", dest="keepcol",type="int", 
            metavar="col", help="Use this column of KEEP [default: %default]", default=1)
        target.add_option("--keepsep", action="store", dest="keepsep", 
            metavar="sep", help="Split KEEP using this [default: white space]", default=None)
        target.add_option("--remove", action="store", dest="remove", 
            metavar="REMOVE", help="Include all but these individuals")
        target.add_option("--removecol", action="store", dest="removecol",type="int", 
            metavar="col", help="Use this column of REMOVE [default: %default]", default=1)
        target.add_option("--removesep", action="store", dest="removesep", 
            metavar="sep", help="Split REMOVE using this [default: white space]", default=None)
    if asgroup:
        parser.add_option_group(target)

def ApplySubsetOptions(fm, options):
    if options.extract:
        extract = BasicFilter("marker",False)
        extract.load(options.extract, col=options.extractcol-1, sep=options.extractsep)
        fm.addMarkerFilter(extract)
    if options.exclude:
        exclude = BasicFilter("marker",True)
        exclude.load(options.exclude, col=options.excludecol-1, sep=options.excludesep)
        fm.addMarkerFilter(exclude)
    if options.keep:
        keep = BasicFilter("sample",True)
        keep.load(options.keep, col=options.keepcol-1, sep=options.keepsep)
        fm.addSampleFilter(keep)
    if options.remove:
        remove = BasicFilter("sample",True)
        remove.load(options.remove, col=options.removecol-1, sep=options.removesep)
        fm.addSampleFilter(remove)

def AddNegatableOptions(target, flagnames, defaults, destinations):
    for (f, d, x) in zip(flagnames, defaults, destinations):
        target.add_option("--"+f, action="store_true", dest=x, default=d, help="[default: %default]")
        target.add_option("--no"+f, action="store_false", dest=x, help=SUPPRESS_HELP)

def ConvertToBinary(callargs):
    parser = OptionParser(usage="usage: %prog convert [options] FILE", version = "%prog " + __version__)
    parser.add_option("--outprefix", action="store", dest="prefix", default="data")
    outgroup = OptionGroup(parser, "Output","These flags determine which output files are created")
    optfiles = dict(markers=True,samples=True,bfreq=True,abgeno=True, norm=False, ad=False)
    AddNegatableOptions(outgroup, optfiles.keys(), optfiles.values(), ["create"+x for x in optfiles.keys()])
    parser.add_option_group(outgroup)
    AddSubsetOptions(parser)
    AddFinalReportColumnOptions(parser)
    AddFrequencyOptions(parser)
    (options, args) = parser.parse_args(callargs)

    if len(args)<1:
        exit("Error: no final report files specified")

    om = OutputManager(options.prefix, snpcolname=options.colmarker, samplecolname=options.colsample)
    if options.createad:
        om.addOutputFile(BinaryOutputFile("\x01",".ad.bin",(options.colrawx, options.colrawy), 
            "H", 2, lambda x: map(int, x)))
    if options.createnorm:
        om.addOutputFile(BinaryOutputFile("\x11",".norm.bin",(options.colnormx, options.colnormy), 
            "f", 2, lambda x: map(float, x)))
    if options.createbfreq:
        om.addOutputFile(BinaryOutputFile("\x21",".bfreq.bin",(options.colbaf,), 
            "i", 1, bafEncode))
    if options.createabgeno:
        om.addOutputFile(BinaryOutputFile("\x31",".abgeno.bin",(options.colab1,options.colab2), 
            "B", 1, genoParser))
    if options.createsamples:
        om.addOutputFile(ReferenceTextOutputFile(".samples.txt","sample"))
    if options.createmarkers:
        om.addOutputFile(ReferenceTextOutputFile(".markers.txt","marker"))
    if options.freqfile:
        freq = PopulationFrequencyFile("\x25",".maf.bin",[int(x)-1 for x in options.freqcols.split(",")], options.freqsep)
        freq.load(options.freqfile)
        om.addOutputFile(freq)
        om.addMarkerFilter(freq)
    ApplySubsetOptions(om, options)

    om.open()

    for i,file in enumerate(args):
        print("processing ", file, "(", str(i+1), "/", str(len(args)), ")")
        f=FinalReport(file)
        f.open()
        om.nextinput(f)
        for line in f.read():
            om.process(line)
        f.close()

    om.close()

def ListSamplesRaw(callargs):
    parser = OptionParser(usage="usage: %prog listsamplesraw [options] FILE", version = "%prog " + __version__)
    AddFinalReportColumnOptions(parser)
    (options, args) = parser.parse_args(callargs)

    if len(args)<1:
        exit("Error: no final report files specified")
    if len(args)>1:
        exit("Error: only one final report file may be specified")

    f = FinalReport(args[0])
    f.open()
    lastSample = None
    print("sample")
    for line in f.read():
        sample = line[options.colsample] 
        if sample != lastSample:
            print(sample)
            lastSample = sample
    f.close()

def ReadFromBinary(callargs):
    parser = OptionParser(usage="usage: %prog readbin [options] PREFIX", version = "%prog " + __version__)
    parser.add_option("--sample", action="store", dest="sample")
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser)
    (options, args) = parser.parse_args(callargs)

    filter = FilterBox()
    ApplySubsetOptions(filter, options)

    if len(args)==1:
        inf = InputManager()

        markerlist = open(args[0] + ".markers.txt")
        inf.markers = markerlist.read().splitlines()
        markerlist.close()

        samplelist = open(args[0] + ".samples.txt")
        inf.samples = samplelist.read().splitlines()
        samplelist.close()

        inf.addInputFile(Binary2DInputFile("\x21",".bfreq.bin","i", 1, len(inf.markers), ("BAF",), bafDecode))
        inf.addInputFile(Binary2DInputFile("\x31",".abgeno.bin","B", 1, len(inf.markers), ("ABGENO",)))
        if options.freqfile:
            ff = PopulationFrequencyFile(0, "", [int(x)-1 for x in options.freqcols.split(",")], options.freqsep)
            ff.load(options.freqfile)
            inf.addInputFile(ff)
        elif os.path.isfile(args[0]+".maf.bin"):
            inf.addInputFile(Binary1DInputFile("\x25",".maf.bin","d", 1, len(inf.markers), ("MAF",)))

        if options.sample:
            try:
                si = inf.samples.index(options.sample)
            except ValueError:
                sys.exit('Sample ' + options.sample + ' not found')
            readsamples = (si,);
        else:
            readsamples = list()
            for i,s in enumerate(inf.samples):
                if filters.shouldIncludeSample(s):
                    readsamples.append(i)
        if len(readsamples)>1:
            inf.includeSampleName = True

        inf.open(args[0])

        print("\t".join([str(x) for x in inf.header()]))
        for s in readsamples:
            for m in range(len(inf.markers)):
                vals = inf.read(s,m)
                print("\t".join([str(x) for x in vals]))

        inf.close()
    else:
        sys.exit("Error: Exactly one data prefix required")

def ReadFromRaw(callargs):
    parser = OptionParser(usage="usage: %prog readraw [options] FILE", version = "%prog " + __version__)
    parser.add_option("--sample", action="store", dest="sample") 
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser)
    AddFinalReportColumnOptions(parser)
    (options, args) = parser.parse_args(callargs)
    sentHeader = False

    om = OutputManager("", snpcolname=options.colmarker, samplecolname=options.colsample)
    if(options.sample):
       om.addSampleFilter(BasicFilter("sample",False,(options.sample,))) 
    fps = FinalReportStream(options)
    if options.freqfile:
        ff = PopulationFrequencyFile(0,"",[int(x)-1 for x in options.freqcols.split(",")], options.freqsep)
        ff.load(options.freqfile)
        om.addMarkerFilter(ff)
        fps.freq = ff
    om.addOutputFile(fps)
    ApplySubsetOptions(om, options)

    om.open()

    for file in args:
        f=FinalReport(file)
        f.open()
        om.nextinput(f)
        for line in f.read():
            om.process(line) 
        f.close()

    om.close()


def AnalyzeSamplesBinary(callargs):
    parser = OptionParser(usage="usage: %prog estimatebin [options] PREFIX", version = "%prog " + __version__)
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser)
    (options, args) = parser.parse_args(callargs)

    pyscript = os.path.realpath(__file__)
    rscript = os.path.splitext(pyscript)[0]+".R"

    cmd = ["Rscript", rscript, "testsamplesbinary"]
    if len(args)==1:
        cmd.append(args[0])
        cmd.extend(["--py",pyscript])
        for opt in options.__dict__.keys():
            if getattr(options, opt):
                cmd.append("--" + opt)
                cmd.append(str(getattr(options, opt)))
        call(cmd)
    else:
        sys.exit("Must specify exactly one data prefix")

    
def AnalyzeSamplesRaw(callargs):
    parser = OptionParser(usage="usage: %prog estimate [options] FILE", version = "%prog " + __version__)
    parser.add_option("--sample", action="store", dest="sample", 
        default=None, help="sample name to analyze [default: first sample in file]")
    parser.add_option("--stacked", action="store_true", dest="stacked", 
        default=False, help="use if more than one sample per final report")
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser, samples=False)
    AddFinalReportColumnOptions(parser)
    (options, args) = parser.parse_args(callargs)

    if not options.freqfile:
        sys.exit("Error: --freqfile required")
    if len(args)<1:
        sys.exit("Error: must specify a final report file")
    
    pyscript = os.path.realpath(__file__)
    rscript = os.path.splitext(pyscript)[0]+".R"
    isFirst = True
    for file in args:
        fpfile = args[0]
        if options.stacked and not options.sample:
            #stacked, do all
            cmd = ["Rscript", rscript, "testsamplesraw"]
            cmd.append(file)
        else:
            #not stacked
            cmd = ["Rscript", rscript, "testonesampleraw"]
            if not options.sample:
                if options.colsample == "NONE":
                    sample = file
                else:
                    fp = FinalReport(file)
                    fp.open()
                    sample = next(fp.read())[options.colsample]
                    fp.close()
            else:
                sample = options.sample
            cmd.extend((file, sample))
        cmd.extend(["--py",pyscript])
        for opt in options.__dict__.keys():
            if getattr(options, opt) and opt != "sample" and opt != "stacked":
                cmd.append("--" + opt)
                cmd.append(str(getattr(options, opt)))
        if not isFirst:
            cmd.extend(("--noheader","True"))
        call(cmd)
        isFirst = False

def PlotSampleBinary(callargs):
    parser = OptionParser(usage="usage: %prog plotbin [options] PREFIX", version = "%prog " + __version__)
    parser.add_option("--sample", action="store", dest="sample", 
        default=None, help="sample name to plot [default: first sample in file]")
    parser.add_option("--outprefix", action="store", dest="outprefix", 
        default=None, help="file name without extension [default: --sample]")
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser, samples=False)
    (options, args) = parser.parse_args(callargs)

    if len(args)!=1:
        sys.exit("Error: must specify exactly one data prefix")

    pyscript = os.path.realpath(__file__)
    rscript = os.path.splitext(pyscript)[0]+".R"

    if not options.sample:
        with open(args[0] + ".samples.txt") as f:
            sample = f.readline().rstrip()
    else:
        sample = options.sample
    if not options.outprefix:
        options.outprefix=sample

    cmd = ["Rscript", rscript, "plotsamplebinary"]
    cmd.extend((args[0], sample))
    cmd.extend(("--py",pyscript))
    for opt in options.__dict__.keys():
        if getattr(options, opt) and opt != "sample":
            cmd.append("--" + opt)
            cmd.append(str(getattr(options, opt)))
    call(cmd)

def PlotSampleRaw(callargs):
    parser = OptionParser(usage="usage: %prog plot [options] FILE", version = "%prog " + __version__)
    parser.add_option("--sample", action="store", dest="sample", 
        default=None, help="sample name to plot [default: first sample in file]")
    parser.add_option("--outprefix", action="store", dest="outprefix", 
        default=None, help="file name without extension [default: --sample]")
    AddFrequencyOptions(parser)
    AddSubsetOptions(parser, samples=False)
    AddFinalReportColumnOptions(parser)
    (options, args) = parser.parse_args(callargs)

    if not options.freqfile:
        sys.exit("Error: --freqfile required")
    if len(args)!=1:
        sys.exit("Error: must specify exactly one final report file")
    
    pyscript = os.path.realpath(__file__)
    rscript = os.path.splitext(pyscript)[0]+".R"

    if not options.sample:
        if options.colsample=="NONE":
            sample = file
        else:
            fp = FinalReport(file)
            fp.open()
            sample = fp.read().next()[options.colsample]
            fp.close()
    else:
        sample = options.sample
    if not options.outprefix:
        options.outprefix = sample

    cmd = ["Rscript", rscript, "plotsampleraw"]
    cmd.extend((args[0], sample))
    cmd.extend(("--py",pyscript))
    for opt in options.__dict__.keys():
        if getattr(options, opt) and opt != "sample":
            cmd.append("--" + opt)
            cmd.append(str(getattr(options, opt)))
    call(cmd)

if __name__ == "__main__":

    actions = dict(convert=ConvertToBinary, 
        readbin=ReadFromBinary,
        readraw=ReadFromRaw,
        estimatebin=AnalyzeSamplesBinary,
        estimate=AnalyzeSamplesRaw,
        plotbin=PlotSampleBinary,
        plot=PlotSampleRaw,
        listsamplesraw=ListSamplesRaw)

    if len(sys.argv) < 2:
        sys.exit('Usage: %s [%s]\nVersion %s' % (sys.argv[0],",".join(actions.keys()), __version__))

    action = sys.argv[1]
    if action in actions:
        actions[action](sys.argv[2:])
    else:
        sys.exit('Unrecognized action: %s. Expecting: %s\nVersion %s' % ( action, ",".join(actions.keys()), __version__) )
