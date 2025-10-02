(* ::Package:: *)

loadData[baseDir_,ibexDataDir_,ibexDataName_,spinAxDir_,outputDir_]:=Module[{ibexFiles},
ibexFiles = Map[FileNameTake,FileNames["*.qABC", ToFileName[ibexDataDir,ToString[ibexDataName]]]];
Get[ToFileName[spinAxDir, "spinAxEclTabWarsaw.m"]];
Import[baseDir<>"/geometricPackage.wl"];
Import[baseDir<>"/mapCreationPackage.wl"];
ibexHi = Map[ fun`ibexDataRead[ToFileName[ibexDataDir,ibexDataName], #]&, ibexFiles]
]


run[ibexHi_, healpixDir_] := Module[{hpPath, result, mapCountsMain, mapSignalMain, mapExposuresMain},
hpPath = ToFileName[healpixDir, "testXYZ"];
init[hpPath, tesselation];
AbsoluteTiming[
result = ParallelMap[
Module[{orbitCount = #, exposuretime, counts, signal, backgroundRate, ibexLatitude, ibexLongitude, 
rotationAxisAng, visibilityRangePixels, oneOrbitMap, mapCounts, mapExposures, mapSignal,mapColValues, normalizationFactor},
exposuretime = ToExpression["expE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
exposuretime=exposuretime*10^-3*1.;
counts = ToExpression["ctsE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
backgroundRate=ToExpression["bkgE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLatitude  = ToExpression["eclatE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLongitude = ToExpression["eclonE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
rotationAxisAng = {spinAxEcLon /. ibexHi[[orbitCount]], spinAxEcLat /. ibexHi[[orbitCount]]};
visibilityRangePixels = choseRing[rotationAxisAng];
Print["Loading orbit ", orbNo /. ibexHi[[orbitCount]]];
oneOrbitMap = calcOneOrbit[ibexLatitude, ibexLongitude, exposuretime, counts, backgroundRate, visibilityRangePixels, healpixringxyz];
mapCounts = Total[Module[{element = #}, First[Last[#]] & /@ element] & /@ oneOrbitMap];
mapExposures = Total[Module[{element = #}, Last[#][[2]] & /@ element] & /@ oneOrbitMap];
mapSignal = Total[Module[{element = #}, Last[#][[3]] & /@ element] & /@ oneOrbitMap];

{mapCounts, mapExposures, mapSignal}] &, Range[Length[ibexHi]]];
mapCountsMain    = Total[result[[All, 1]]];
mapExposuresMain = Total[result[[All, 2]]];
mapSignalMain = Total[result[[All, 3]]];
{mapCountsMain, mapExposuresMain, mapSignalMain}]
]


exportData[mapCountsMain_, mapExposuresMain_,mapSignalMain_,tesselation_,energyStep_,outputDir_]:=Module[{mapRatesMainMathematica,mapRatesMainHealpy,mapENAFluxMathematica,mapENAFluxHealpy,datMathematica,datHealpy},
mapRatesMainMathematica=calculateRatesMathematica[mapCountsMain,mapExposuresMain];
mapENAFluxMathematica=calculateENAFluxMathematica[mapSignalMain, mapExposuresMain];
datMathematica=ExportString[Transpose[Append[{Range[Length[mapCountsMain]],mapCountsMain,mapExposuresMain,mapRatesMainMathematica,mapSignalMain},mapENAFluxMathematica]],"Table"];
Export[ToFileName[outputDir]<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",datMathematica];
mapRatesMainHealpy=calculateRatesHealpy[mapCountsMain,mapExposuresMain];
mapENAFluxHealpy=calculateENAFluxHealpy[mapSignalMain, mapExposuresMain];
datHealpy=ExportString[Transpose[Append[{Range[Length[mapCountsMain]],mapCountsMain,mapExposuresMain,mapRatesMainHealpy,mapSignalMain},mapENAFluxHealpy]],"Table"];
Export[ToFileName[outputDir]<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"Null.txt",datHealpy];
]
