(* ::Package:: *)

loadData[baseDir_,ibexDataDir_,ibexDataName_,spinAxDir_,outputDir_]:=Module[{ibexFiles},
ibexFiles = Map[FileNameTake,FileNames["*.qABC", ToFileName[ibexDataDir,ToString[ibexDataName]]]];
Get[ToFileName[spinAxDir, "spinAxEclTabWarsaw.m"]];
Import[baseDir<>"/geometricPackage.wl"];
Import[baseDir<>"/mapCreationPackage.wl"];
ibexHi = Map[ fun`ibexDataRead[ToFileName[ibexDataDir,ibexDataName], #]&, ibexFiles]
]


run[ibexHi_, healpixDir_] := Module[{hpPath, result, mapCountsMain, mapExposuresMain},
hpPath = ToFileName[healpixDir, "testXYZ"];
init[hpPath, tesselation];
AbsoluteTiming[
result = ParallelMap[
Module[{orbitCount = #, exposuretime, counts, ibexLatitude, ibexLongitude, 
rotationAxisAng, visibilityRangePixels, oneOrbitMap, mapCounts, mapExposures},
exposuretime = ToExpression["expE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
counts = ToExpression["ctsE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLatitude  = ToExpression["eclatE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLongitude = ToExpression["eclonE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
rotationAxisAng = {spinAxEcLon /. ibexHi[[orbitCount]], spinAxEcLat /. ibexHi[[orbitCount]]};
visibilityRangePixels = choseRing[rotationAxisAng];
Print["Loading orbit ", orbNo /. ibexHi[[orbitCount]]];
oneOrbitMap = calcOneOrbit[ibexLatitude, ibexLongitude, exposuretime, counts, visibilityRangePixels, healpixringxyz];
mapCounts = Total[Module[{element = #}, First[Last[#]] & /@ element] & /@ oneOrbitMap];
mapExposures = Total[Module[{element = #}, Last[#][[2]] & /@ element] & /@ oneOrbitMap];
{mapCounts, mapExposures}] &, Range[Length[ibexHi]]];
mapCountsMain    = Total[result[[All, 1]]];
mapExposuresMain = Total[result[[All, 2]]];
{mapCountsMain, mapExposuresMain}]
]


exportData[mapCountsMain_, mapExposuresMain_,tesselation_,energyStep_,outputDir_]:=Module[{mapRatesMainMathematica,mapRatesMainHealpy,datMathematica,datHealpy},
mapRatesMainMathematica=calculateRatesMathematica[mapCountsMain,mapExposuresMain];
datMathematica=ExportString[Transpose[Append[{Range[Length[mapCountsMain]],mapCountsMain,mapExposuresMain},mapRatesMainMathematica]],"Table"];
Export[ToFileName[outputDir]<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",datMathematica];
mapRatesMainHealpy=calculateRatesHealpy[mapCountsMain,mapExposuresMain];
datHealpy=ExportString[Transpose[Append[{Range[Length[mapCountsMain]],mapCountsMain,mapExposuresMain},mapRatesMainHealpy]],"Table"];
Export[ToFileName[outputDir]<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"Null.txt",datHealpy];
]
