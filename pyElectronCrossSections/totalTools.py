#!/usr/bin/env python
################################################################################
# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2007 Hendrix Demers"
__license__ = ""

# Subversion informations for the file.
__svnRevision__ = "$Revision: 2636 $"
__svnDate__ = "$Date: 2011-12-15 16:58:59 -0500 (Thu, 15 Dec 2011) $"
__svnId__ = "$Id: totalTools.py 2636 2011-12-15 21:58:59Z hdemers $"

import sys
import os
import logging

################################################################################
#    global variables.
g_modelList = ["RutherfordRelativisticScreened",
                            "NISTDatabase64",
                            "MottCzyzewski",
                            "MottBrowning",
                            "MottDrouin"]
g_modelNameList = {"RutherfordRelativisticScreened": "Rutherford",
                                    "NISTDatabase64": "NIST",
                                    "MottCzyzewski": "Czyzewski",
                                    "MottBrowning": "Browning",
                                    "MottDrouin": "Drouin"}
g_totalResultType = ["Total_3D", "Total_E", "Total_Z"]

################################################################################
def getResultPath():
    try:
        resultPath = sys.argv[1]
        if resultPath[-1] != os.path.sep:
            resultPath += os.path.sep

#        print resultPath
        return resultPath
    except IndexError:
        print("Usage is:" + sys.argv[0] + " resultPath")
        sys.exit()

################################################################################
def getModelPresent(fileList):
    modelList = []

    for filename in fileList:
        for model in g_modelList:
            if model in filename:
                if not model in modelList:
                    modelList.append(model)

    return modelList

################################################################################
def getResultType(fileList):
    resultTypeList = []
    for filename in fileList:
        for resultsType in g_totalResultType:
            if resultsType in filename:
                if not resultsType in resultTypeList:
                    resultTypeList.append(resultsType)

    return resultTypeList

################################################################################
def getTotalFileList(fileList):
    modelList = getModelPresent(fileList)
    resultTypeList = getResultType(fileList)

    totalFileList = {}
    for model in modelList:
        if not model in totalFileList:
            totalFileList[model] = {}
        for resultsType in resultTypeList:
            totalFileList[model][resultsType] = {}

#    print len(fileList)
#    print len(modelList)
#    print len(resultTypeList)

    for filename in fileList:
        for model in modelList:
            if model in filename:
                for resultsType in resultTypeList:
                    if resultsType in filename:
                        totalFileList[model][resultsType]["filename"] = filename
#                        print model, resultsType, filename

    return totalFileList

################################################################################
def read3DData(filename):
    dataDict = {}
    atomicNumberList = []
    energyList = []

    for line in open(filename, 'r').readlines():
        line = line.strip()
        if not '#' in line:
            (atomicNumber, energy, sigma) = line.split()

            atomicNumber = int(atomicNumber)
            energy = float(energy)
            sigma = float(sigma)

            if not atomicNumber in atomicNumberList:
                atomicNumberList.append(atomicNumber)
            if not energy in energyList:
                energyList.append(energy)

            if not atomicNumber in dataDict:
                dataDict[atomicNumber] = {}
            if not energy in dataDict[atomicNumber]:
                dataDict[atomicNumber][energy] = {}
            if not "sigma" in dataDict[atomicNumber][energy]:
                dataDict[atomicNumber][energy]["sigma"] = sigma

#        else:
#            print "Comment line."

    return dataDict, atomicNumberList, energyList

################################################################################
def readZData(filename):
    dataDict = {}
    for line in open(filename, 'r').readlines():
        line = line.strip()
        if not '#' in line:
            (atomicNumber, energy, sigma) = line.split()
            atomicNumber = int(atomicNumber)
            energy = float(energy)
            sigma = float(sigma)

            if not atomicNumber in dataDict:
                dataDict[atomicNumber] = {}
            if not energy in dataDict[atomicNumber]:
                dataDict[atomicNumber][energy] = {}
            if not "sigma" in dataDict[atomicNumber][energy]:
                dataDict[atomicNumber][energy]["sigma"] = sigma

        else:
            print("Comment line.")

    return dataDict

################################################################################
# TODO: Implement readEData(filename)
#def readEData(filename):
#    energyList = getEnergyList(filename)
#
#    dataDict = {}
#    for line in open(filename, 'r').readlines():
#        line = line.strip()
#        if not '#' in line:
#            sigma = []
#            (atomicNumber, sigma) = line.split()
#            atomicNumber = int(atomicNumber)
#            sigma = float(sigma)
#
#            if not atomicNumber in dataDict:
#                dataDict[atomicNumber] = {}
#            if not "sigma" in dataDict[atomicNumber][energy]:
#                dataDict[atomicNumber][energy]["sigma"] = sigma
#
#        else:
#            print "Comment line."
#
#    return dataDict

################################################################################
def getAtomicNumberList(filename):
    line = open(filename).readline()

    atomicNumberList = line.split('\t')

    atomicNumberList = atomicNumberList[1:-1]

    return atomicNumberList

################################################################################
def getEnergyList(filename):
    line = open(filename).readline()

    energyList = line.split('\t')

    energyList = energyList[1:-1]

    return energyList

################################################################################
def creataDataForZ(atomicNumber, totalFileList, modelRef):
    typeKey = "Total_3D"
    dataKey = "data"

    dataModel = {}

    for model in g_modelList:
        if model in totalFileList:
            dataModel[model] = []

            if typeKey in totalFileList[model]:
                if dataKey in totalFileList[model][typeKey]:
                    data = totalFileList[model][typeKey][dataKey]

                    keys = data[atomicNumber].keys()
                    keys.sort()
                    for energy in keys:
                        sigma = data[atomicNumber][energy]["sigma"]
                        dataModel[model].append([energy, sigma])

    listRef = dataModel[modelRef]

    for model in g_modelList:
        if model in dataModel:
            newList = []
            for (energyList, energyRefList) in zip(dataModel[model], listRef):
                energy = energyList[0]
                sigmaRef = energyRefList[1]
#                print energy, sigmaRef
                ratioRef = 1.0

                sigma = energyList[1]
                sigmaDiff = sigma - sigmaRef
                sigmaDiffRelative = sigmaDiff/sigmaRef

                ratio = sigma/sigmaRef
                ratioDiff = ratio - ratioRef
                ratioDiffRelative = ratioDiff/ratioRef

                newList.append([energy, sigma, sigmaDiff, sigmaDiffRelative, ratio, ratioDiff, ratioDiffRelative])

            dataModel[model] = newList

#    print len(dataModel)

    return dataModel

################################################################################
def creataDataForE(energy, totalFileList, modelRef):
    typeKey = "Total_3D"
    dataKey = "data"

    dataModel = {}

    for model in g_modelList:
        if model in totalFileList:
            dataModel[model] = []

            if typeKey in totalFileList[model]:
                if dataKey in totalFileList[model][typeKey]:
                    data = totalFileList[model][typeKey][dataKey]

                    keys = data.keys()
                    keys.sort()
                    for atomicNumber in keys:
                        sigma = data[atomicNumber][energy]["sigma"]
                        dataModel[model].append([atomicNumber, sigma])

    listRef = dataModel[modelRef]

    for model in g_modelList:
        if model in dataModel:
            newList = []
            for (atomicNumberList, atomicNumberRefList) in zip(dataModel[model], listRef):
                atomicNumber = atomicNumberList[0]
                sigmaRef = atomicNumberRefList[1]
#                print energy, sigmaRef
                ratioRef = 1.0

                sigma = atomicNumberList[1]
                sigmaDiff = sigma - sigmaRef
                sigmaDiffRelative = sigmaDiff/sigmaRef

                ratio = sigma/sigmaRef
                ratioDiff = ratio - ratioRef
                ratioDiffRelative = ratioDiff/ratioRef

                newList.append([atomicNumber,
                                                sigma,
                                                sigmaDiff,
                                                sigmaDiffRelative,
                                                ratio,
                                                ratioDiff,
                                                ratioDiffRelative])

            dataModel[model] = newList

#    print len(dataModel)

    return dataModel

################################################################################
def createPlot(dataModel, resultPath, name, yIndex, xLog=False, yLog=False):
    try:
        import pyx

        graphDataList = []

        for model in g_modelList:
    #        print model
            if model in dataModel:
                xGraphAxis = pyx.graph.axis.linear()
                yGraphAxis = pyx.graph.axis.linear()

                if xLog == True:
                    xGraphAxis = pyx.graph.axis.logarithmic()
                if yLog == True:
                    yGraphAxis = pyx.graph.axis.logarithmic()

                g = pyx.graph.graphxy(width=16,
                                      x=xGraphAxis,
                                      y=yGraphAxis,
                                      key=pyx.graph.key.key(pos="tr", dist=0.1))


                graphDataList.append(pyx.graph.data.list(dataModel[model], x=1, y=yIndex, title=g_modelNameList[model]))

    #    print len(graphDataList)
    #    print len(dataModel)
    #    colorList = [color.rgb.black, color.rgb.red, color.rgb.blue, color.rgb.green, color.cmyk.gray]
        colorList = pyx.color.palette.Rainbow

        graphStyleLine = pyx.graph.style.line([colorList])
        graphStyleLineList = [graphStyleLine]

        g.plot(graphDataList, graphStyleLineList)
        outputFilename = name

    #    print resultPath + outputFilename
        g.writeEPSfile(resultPath + outputFilename)
        g.writePDFfile(resultPath + outputFilename)
        print(outputFilename)
    except ImportError:
        logging.error("pyx module not found")

################################################################################
def createPlotZ(totalFileList, atomicNumberList, resultPath):
    print("Create plot Z.")

    for atomicNumber in atomicNumberList:
#    for atomicNumber in range(1,2):

        atomicNumberStr = ""

        if int(atomicNumber) < 100:
            atomicNumberStr += "0"
        if int(atomicNumber) < 10:
            atomicNumberStr += "0"

        atomicNumberStr += str(atomicNumber)

        modelRef = "RutherfordRelativisticScreened"

        dataModel = creataDataForZ(atomicNumber, totalFileList, modelRef)

        outputFilename = "TotalSigma_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 2, xLog=True, yLog=True)

        outputFilename = "TotalSigmaDiff_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 3, xLog=True)

        outputFilename = "TotalSigmaDiffRel_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 4, xLog=True)

        outputFilename = "TotalRatio_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 5, xLog=True)

        modelRef = "NISTDatabase64"

        dataModel = creataDataForZ(atomicNumber, totalFileList, modelRef)

        outputFilename = "NISTTotalSigmaDiff_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 3, xLog=True)

        outputFilename = "NISTTotalSigmaDiffRel_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 4, xLog=True)

        outputFilename = "NISTTotalRatio_Z" + atomicNumberStr
        createPlot(dataModel, resultPath, outputFilename, 5, xLog=True)
################################################################################
def createPlotE(totalFileList, energyList, resultPath):
    print("Create plot E.")

#    for energy in energyList[:1]:
    for energy in energyList:

        energyStr = ""

        if int(energy) < 10000:
            energyStr += "0"
        if int(energy) < 1000:
            energyStr += "0"
        if int(energy) < 100:
            energyStr += "0"
        if int(energy) < 10:
            energyStr += "0"

        energyStr += str(int(energy))
        energyStr += "eV"

        modelRef = "RutherfordRelativisticScreened"

        dataModel = creataDataForE(energy, totalFileList, modelRef)

        outputFilename = "TotalSigma_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 2, yLog=True)

        outputFilename = "TotalSigmaDiff_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 3)

        outputFilename = "TotalSigmaDiffRel_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 4)

        outputFilename = "TotalRatio_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 5)

        modelRef = "NISTDatabase64"

        dataModel = creataDataForE(energy, totalFileList, modelRef)

        outputFilename = "NISTTotalSigmaDiff_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 3)

        outputFilename = "NISTTotalSigmaDiffRel_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 4)

        outputFilename = "NISTTotalRatio_E" + energyStr
        createPlot(dataModel, resultPath, outputFilename, 5)

################################################################################
def intersect(array):
    result = []
    if len(array) > 0:
        seq1 = array[0]
        for x in seq1:
            for seq in array[1:]:
                if not x in seq:
                    break
            else:
                result.append(x)

    return result

################################################################################
def readData(totalFileList, resultPath):
    print("Read data from file.")
    typeKey = "Total_3D"
    filenameKey = "filename"
    dataKey = "data"

    atomicNumberArray = []
    energyArray = []

    for model in totalFileList:
        if typeKey in totalFileList[model]:
            if filenameKey in totalFileList[model][typeKey]:
                filename = totalFileList[model][typeKey][filenameKey]
#                print filename

                (data, ZList, EList) = read3DData(resultPath+filename)

                atomicNumberArray.append(ZList)
                energyArray.append(EList)

                totalFileList[model][typeKey][dataKey] = data

    atomicNumberList = intersect(atomicNumberArray)
    energyList = intersect(energyArray)

#    print len(atomicNumberList), atomicNumberList[0], atomicNumberList[-1]
#    print len(energyList), energyList[0], energyList[-1]

    return totalFileList, atomicNumberList, energyList

def run():
    import os

    resultPath = getResultPath()
    fileList = os.listdir(resultPath)

    totalFileList = getTotalFileList(fileList)

    (totalFileList, atomicNumberList, energyList) = readData(totalFileList, resultPath)

    createPlotZ(totalFileList, atomicNumberList, resultPath)

    createPlotE(totalFileList, energyList, resultPath)

if __name__ == "__main__":
    run()
