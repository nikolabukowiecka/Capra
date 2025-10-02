(* ::Package:: *)

createMapMean[coordinates_,data_,resolution_]:=Module[{healpixPixelsClassified,healpixPixelsClassifiedWithout0,zeroList,pixels1,pixels2,ranges,dataNew,pixels1Val,pixels2Val,listZero,list},
(*data[[1;;15]]="bad";*)
zeroList=Flatten[Position[data,"bad"]];
healpixPixelsClassified=classifyHEALPixPixels[coordinates,resolution];
healpixPixelsClassifiedWithout0=Complement[DeleteCases[Flatten[healpixPixelsClassified,1],0],Partition[zeroList,1]];
pixels1=Select[healpixPixelsClassifiedWithout0,Length[#]==1&];
pixels2=Complement[healpixPixelsClassifiedWithout0,pixels1];
pixels2=Delete[pixels2,Flatten[(Position[pixels2,#]&/@zeroList),1]];
pixels1=Flatten[pixels1];
ranges=Length[pixels2[[#]]]&/@Range[Length[pixels2]];
dataNew=ConstantArray[0,Length[data]];
pixels1Val=Transpose[Append[{pixels1},(data[[#]]&/@pixels1)]];
pixels2Val=Transpose[Append[{Flatten[pixels2]},(data[[#]]&/@Flatten[pixels2])]];
listZero=Transpose[Append[{zeroList},ConstantArray[0,Length[zeroList]]]];
list=SortBy[Join[pixels1Val,pixels2Val,listZero],1];
Module[{idx1=#},dataNew[[pixels1[[idx1]]]]=list[[pixels1[[idx1]],2]]]&/@Range[Length[pixels1]];
Module[{idx1=#},Module[{idx2=#},dataNew[[pixels2[[idx1]][[idx2]]]]=Mean[list[[pixels2[[idx1]]]][[;;,2]]]]&/@Range[ranges[[idx1]]]]&/@Range[Length[pixels2]];
If[zeroList!={},
dataNew=Last[With[{x=zeroList[[#]]},dataNew=ReplacePart[dataNew,x->Null]]&/@Range[Length[zeroList]]]];
dataNew
]

createMapTotal[coordinates_,data_,resolution_]:=Module[{healpixPixelsClassified,healpixPixelsClassifiedWithout0,zeroList,pixels1,pixels2,ranges,dataNew,pixels1Val,pixels2Val,listZero,list},
(*data[[1;;15]]="bad";*)
zeroList=Flatten[Position[data,"bad"]];
healpixPixelsClassified=classifyHEALPixPixels[coordinates,resolution];
healpixPixelsClassifiedWithout0=Complement[DeleteCases[Flatten[healpixPixelsClassified,1],0],Partition[zeroList,1]];
pixels1=Select[healpixPixelsClassifiedWithout0,Length[#]==1&];
pixels2=Complement[healpixPixelsClassifiedWithout0,pixels1];
pixels2=Delete[pixels2,Flatten[(Position[pixels2,#]&/@zeroList),1]];
pixels1=Flatten[pixels1];
ranges=Length[pixels2[[#]]]&/@Range[Length[pixels2]];
dataNew=ConstantArray[0,Length[data]];
pixels1Val=Transpose[Append[{pixels1},(data[[#]]&/@pixels1)]];
pixels2Val=Transpose[Append[{Flatten[pixels2]},(data[[#]]&/@Flatten[pixels2])]];
listZero=Transpose[Append[{zeroList},ConstantArray[0,Length[zeroList]]]];
list=SortBy[Join[pixels1Val,pixels2Val,listZero],1];
Module[{idx1=#},dataNew[[pixels1[[idx1]]]]=list[[pixels1[[idx1]],2]]]&/@Range[Length[pixels1]];
Module[{idx1=#},Module[{idx2=#},dataNew[[pixels2[[idx1]][[idx2]]]]=Total[list[[pixels2[[idx1]]]][[;;,2]]]]&/@Range[ranges[[idx1]]]]&/@Range[Length[pixels2]];
If[zeroList!={},
dataNew=Last[With[{x=zeroList[[#]]},dataNew=ReplacePart[dataNew,x->Null]]&/@Range[Length[zeroList]]]];
dataNew
]

createMapDataRelease[coordinates_,data_,resolution_,zeroList_]:=
Module[{healpixPixelsClassified,healpixPixelsClassifiedWithout0,pixels1,pixels2,ranges,dataNew,map},
healpixPixelsClassified=classifyHEALPixPixels[coordinates,resolution];
healpixPixelsClassifiedWithout0=DeleteCases[Flatten[healpixPixelsClassified,1],0];
map=Transpose[Append[{healpixPixelsClassifiedWithout0},Flatten[data]]];
pixels1=Select[healpixPixelsClassifiedWithout0,Length[#]==1&];
pixels2=Select[healpixPixelsClassifiedWithout0,Length[#]>1&];
(*pixels2=Delete[pixels2,Flatten[(Position[pixels2,#]&/@zeroList),1]];*)
pixels1=Flatten[pixels1];
ranges=Length[pixels2[[#]]]&/@Range[Length[pixels2]];
dataNew=ConstantArray[0,Length[healpixringxyz]];
Module[{idx1=#,x},x=Position[map,idx1];dataNew[[map[[x[[1,1]]]][[1,1]]]]=map[[x[[1,1]]]][[2]]]&/@pixels1;
Module[{idx1=#},
Module[{idx2=#},
x=Position[map,pixels2[[idx1]]][[1,1]];
dataNew[[map[[x,1,idx2]]]]=map[[x,2]]]&/@Range[ranges[[idx1]]]]&/@Range[Length[pixels2]];
If[zeroList!={},
dataNew=Last[Module[{x=#},dataNew=ReplacePart[dataNew,x->Null]]&/@zeroList]];
dataNew]
oneElementArea[idx_,dtheta_,dfi_]:=Module[{start,end,thetas},
start=onePixelResolution[[1]]/2;
end=180-start;
thetas=Range[start,end,dtheta];
Sin[thetas[[idx]] Degree]*dtheta Degree*dfi Degree
]//N

calcWage[onePixelResolution_]:=Module[{oneRow,wages,normalizationFactor},
oneRow=oneElementArea[#,onePixelResolution[[1]],onePixelResolution[[2]]]&/@Range[180/onePixelResolution[[1]]];
wages=ConstantArray[oneRow[[#]],(360/onePixelResolution[[2]])]&/@Range[180/onePixelResolution[[1]]];
(*Total[Total[wages]]/(4*Pi)*)
wages
]

getHPcoord[healpixringxyz_]:=Module[{healpixringang,healpixringangCOLATI,healpixringangLONG,coordinates},
healpixringang=Delete[#,1]&/@Map[getSpherCoord[#]&,healpixringxyz];
healpixringangCOLATI = lati2colatiDeg[healpixringang[[All,2]]];
healpixringangLONG =healpixringang[[All,1]];
coordinates = Transpose[{healpixringangCOLATI,healpixringangLONG}]
]
fun`classify[x_, ile_, top_]:=If[x!=top, IntegerPart[x/top (ile)]+1,1];
classifyHEALPixPixels[coordinates_,resolution_]:=Module[{map,colati,long,rows,columns},
map = ConstantArray[0,{resolution[[1]],resolution[[2]]}];
colati=coordinates[[;;,1]];
long=coordinates[[;;,2]];
rows=Map[fun`classify[#, resolution[[1]], 180]&,colati];
columns =Map[fun`classify[#, resolution[[2]], 360]&,long];
Module[{i=#},
Which[
Length[map[[rows[[i]],columns[[i]]]]]==0, map[[rows[[i]],columns[[i]]]]={i},
Length[map[[rows[[i]],columns[[i]]]]]==1, map[[rows[[i]],columns[[i]]]]=Append[map[[rows[[i]],columns[[i]]]],i],
Length[map[[rows[[i]],columns[[i]]]]]>=2,map[[rows[[i]],columns[[i]]]]=Append[map[[rows[[i]],columns[[i]]]],i]];
]&/@Range[Length[rows]];
map
]
uniformTessellationChange[data_,coordinates_,coordinatesNew_,tesselationNew_]:=Module[{coordinatesNewWithPixels,pixelsOriginalInNewGrid,dataNew,zeroList},
zeroList=Sort[Join[Flatten[Position[data,Null]],Flatten[Position[data,"bad"]]]];
coordinatesNewWithPixels=Transpose[Append[{Range[Length[coordinatesNew]]},coordinatesNew]];
(*Original pixels in the new grid*)
pixelsOriginalInNewGrid=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];
dataNew=ConstantArray[0,Length[coordinates]];
Module[{idx1=#,list},
list=Flatten[Position[pixelsOriginalInNewGrid,idx1]];
Module[{idx2=#},
If[Complement[list,zeroList]!={},
dataNew[[list[[idx2]]]]=Mean[Module[{idx3=#},data[[idx3]]]&/@Complement[list,zeroList]],dataNew[[list[[idx2]]]]=Null];
dataNew[[list[[idx2]]]]]&/@Range[Length[list]];
]&/@pixelsOriginalInNewGrid;
dataNew
]
partialTessellationChange[coordinatesChanged_,data_,coordinates_,coordinatesNew_,tesselationNew_]:=Module[{zeroList,coordinatesNewWithPixels,pixelsChanged,pixelsOriginalInNewGrid,dataNew},
zeroList=Sort[Join[Flatten[Position[data,Null]],Flatten[Position[data,"bad"]]]];
coordinatesNewWithPixels=Transpose[Append[{Range[Length[coordinatesNew]]},coordinatesNew]];
pixelsChanged=healpix`fun`ang2pixRing[tesselationNew,coordinatesChanged[[#,1]] Degree,coordinatesChanged[[#,2]]Degree]+1&/@Range[Length[coordinatesChanged]];
pixelsOriginalInNewGrid=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];
dataNew=data;
Module[{idx1=#,list},
list=Flatten[Position[pixelsOriginalInNewGrid,idx1]];
Module[{idx2=#},
If[Complement[list,zeroList]!={},
dataNew[[list[[idx2]]]]=Mean[Module[{idx3=#},data[[idx3]]]&/@Complement[list,zeroList]],dataNew[[list[[idx2]]]]=Null];
dataNew[[list[[idx2]]]]]&/@Range[Length[list]];
]&/@pixelsChanged;
dataNew
]
uniformSizeChange[data_]:=Module[{hpPathNew,healpixringxyzNew,coordinatesNew,coordinatesNewWithPixels,pixelsChanged,pixelsOriginalInNewGrid,healpixPixelsClassified,healpixPixelsClassifiedWithout0,list,dataNew},
(*hpPathNew = ToFileName[healpixDir,"testXYZ"];
healpixringxyzNew=init[hpPathNew,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyzNew];
coordinatesNewWithPixels=Transpose[Append[{Range[Length[coordinatesNew]]},coordinatesNew]];
pixelsChanged=healpix`fun`ang2pixRing[tesselationNew,coordinatesNew[[#,1]] Degree,coordinatesNew[[#,2]]Degree]+1&/@Range[Length[coordinatesNew]];
pixelsChanged=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];*)
(*pixelsOriginalInNewGrid=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];*)
healpixPixelsClassified=classifyHEALPixPixels[coordinates,resolution];
healpixPixelsClassifiedWithout0=DeleteCases[Flatten[healpixPixelsClassified,1],0];
If[zeroList!={},Module[{idx=#,p},p=Position[healpixPixelsClassifiedWithout0,zeroList[[idx]]];
healpixPixelsClassifiedWithout0=Delete[healpixPixelsClassifiedWithout0,p];]&/@Range[Length[zeroList]]];
list=Select[healpixPixelsClassifiedWithout0,Length[#]!=1&];
dataNew=data;
Module[{idx=#},dataNew[[list[[#]]]]=Mean[data[[list[[#]]]]]]&/@Range[Length[list]];
If[zeroList!={},
dataNew=Last[With[{x=zeroList[[#]]},dataNew=ReplacePart[dataNew,x->Null]]&/@Range[Length[zeroList]]]];
dataNew
]
partialSizeChange[data_,tesselationNew_,coordinatesChanged_]:=Module[{hpPathNew,healpixringxyzNew,coordinatesNew,coordinatesNewWithPixels,pixelsChanged,pixelsOriginalInNewGrid,healpixPixelsClassified,healpixPixelsClassifiedWithout0,list,listIndex,dataNew},
hpPathNew = ToFileName[healpixDir,"testXYZ"];
healpixringxyzNew=init[hpPathNew,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyzNew];
coordinatesNewWithPixels=Transpose[Append[{Range[Length[coordinatesNew]]},coordinatesNew]];
pixelsChanged=healpix`fun`ang2pixRing[tesselationNew,coordinatesChanged[[#,1]] Degree,coordinatesChanged[[#,2]]Degree]+1&/@Range[Length[coordinatesChanged]];
pixelsOriginalInNewGrid=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];
healpixPixelsClassified=classifyHEALPixPixels[coordinatesChanged,resolution];
healpixPixelsClassifiedWithout0=DeleteCases[Flatten[healpixPixelsClassified,1],0];
If[zeroList!={},Module[{idx=#,p},p=Position[healpixPixelsClassifiedWithout0,zeroList[[idx]]];
healpixPixelsClassifiedWithout0=Delete[healpixPixelsClassifiedWithout0,p];]&/@Range[Length[zeroList]]];
list=Select[healpixPixelsClassifiedWithout0,Length[#]!=1&];
listIndex=DeleteDuplicates[(Flatten[Position[list,#]&/@pixelsChanged,1])[[;;,1]]];
list=list[[listIndex]];
dataNew=data;
Module[{idx=#},dataNew[[list[[#]]]]=Mean[data[[list[[#]]]]]]&/@Range[Length[list]];
If[zeroList!={},
dataNew=Last[With[{x=zeroList[[#]]},dataNew=ReplacePart[dataNew,x->Null]]&/@Range[Length[zeroList]]]];
dataNew
]
fluxHEALPix[data_,tesselation_]:=Module[{omega,flux},
omega = Pi/(3*tesselation^2)*1.;
flux=data/omega;
If[zeroList!={},
flux=Last[With[{x=zeroList[[#]]},flux=ReplacePart[flux,x->Null]]&/@Range[Length[zeroList]]]];
flux]

fluxENAHEALPix[signal_, exposuretime_, tesselation_, energyStep_]:=Module[{centralEnergies, centralEnergy, geometricFactorTriples, geometricFactor, omega,flux},
geometricFactorTriples = {0.00013,0.00037,0.00073,.0014,0.0025,0.0042};
centralEnergies = {0.45, 0.71, 1.10, 1.74, 2.73, 4.29};
geometricFactor = geometricFactorTriples[[energyStep - 3]];
centralEnergy   = centralEnergies[[energyStep - 3]];
flux=MapThread[If[#2!=0,#1/(#2*geometricFactor*centralEnergy),0]&,{signal,exposuretime}]
]

relativeFluxHEALPix[data_]:=Module[{meanData,relativeData},
meanData=Mean[DeleteCases[data,Null]];
relativeData=data/meanData;
If[zeroList!={},
relativeData=Last[With[{x=zeroList[[#]]},relativeData=ReplacePart[relativeData,x->Null]]&/@Range[Length[zeroList]]]];
relativeData
]

(*This is for all en steps (lr,br - centre of the Ribbon)*)
transformToRibbonCoord[vec_,lr_,br_]:=Module[{UPecl,up,lt,rot1,rot2},
UPecl=makeVec[255Degree,5.14Degree];
up=Ry[90Degree-br Degree] . (Rz[lr Degree] . UPecl);
lt=getSpherCoord[up][[2]];
rot1=Ry[90Degree-br Degree] . (Rz[lr Degree] . vec);
rot2=Rz[lt Degree] . rot1;
rot2
]

classifyHEALPixValues[data_,coordinates_,resolution_,values_]:=Module[{map,colati,long,rows,columns},
map = ConstantArray[0,{resolution[[1]],resolution[[2]]}];
colati=coordinates[[;;,1]];
long=coordinates[[;;,2]];
rows=Map[fun`classify[#, resolution[[1]], 180]&,colati];
columns =Map[fun`classify[#, resolution[[2]], 360]&,long];
Module[{i=#},
Which[
Length[map[[rows[[i]],columns[[i]]]]]==0, map[[rows[[i]],columns[[i]]]]={data[[i]]},
Length[Dimensions[map[[rows[[i]],columns[[i]]]]]]==1, map[[rows[[i]],columns[[i]]]]=Mean[Append[{map[[rows[[i]],columns[[i]]]]},{data[[i]]}]],
Length[Dimensions[map[[rows[[i]],columns[[i]]]]]]>=2,Mean[map[[rows[[i]],columns[[i]]]]=Append[map[[rows[[i]],columns[[i]]]],{data[[i]]}]]];
]&/@Range[Length[rows]];
map
]

createFluxMapHP[healpixDir_,outputDir_]:=Module[{dataPath,dataCts,dataExp,dataRat,zeroList,flux,dat},
dataPath = ToFileName[outputDir];
dataCts=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,2]];
dataExp=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,3]];
dataRat=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRat,"bad"]];
If[zeroList!={},
dataRat=Last[With[{x=zeroList[[#]]},dataRat=ReplacePart[dataRat,x->Null]]&/@Range[Length[zeroList]]]];
flux=fluxHEALPix[dataRat,tesselation];
dat=

String[Transpose[Append[{Range[Length[flux]]},flux]],"Table"];
Export[ToFileName[outputDir]<>"data_flux_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]
createRelativeFluxMapHP[healpixDir_,outputDir_]:=Module[{dataPath,dataCts,dataExp,dataRat,zeroList,flux,fluxNew,relativeFluxNew,meanRelFluxNew,dat},
dataPath = ToFileName[outputDir];
dataCts=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,2]];
dataExp=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,3]];
dataRat=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRat,"bad"]];
If[zeroList!={},
dataRat=Last[With[{x=zeroList[[#]]},dataRat=ReplacePart[dataRat,x->Null]]&/@Range[Length[zeroList]]]];
flux=fluxHEALPix[dataRat,tesselation];
zeroList=Flatten[Position[dataNew,Null]];
fluxNew=fluxHEALPix[dataNew,tesselation];
If[zeroList!={},
fluxNew=Last[With[{x=zeroList[[#]]},fluxNew=ReplacePart[fluxNew,x->Null]]&/@Range[Length[zeroList]]]];
relativeFluxNew=relativeFluxHEALPix[fluxNew];
meanRelFluxNew=Mean[DeleteCases[relativeFluxNew,Null]];
Print["Mean realative flux = "];
Print[meanRelFluxNew];(*Should be 1*)
Print["(Should be 1)"];
dat=ExportString[Transpose[Append[{Range[Length[relativeFluxNew]]},relativeFluxNew]],"Table"];
Export[ToFileName[outputDir]<>"data_relativeFlux_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]

uniformTessellationChangeHP[healpixDir_,outputDir_,tesselationNew_]:=Module[{hpPath,coordinates,dataPath,dataRat,hpPathNew,coordinatesNew,dataNew,dat},
hpPath = ToFileName[healpixDir,"testXYZ"];
init[hpPath,tesselation];
coordinates=getHPcoord[healpixringxyz];
dataPath = ToFileName[outputDir];
dataRat=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
hpPathNew = ToFileName[healpixDir,"testXYZ"];
init[hpPathNew,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyz];
dataNew=uniformTessellationChange[dataRat,coordinates,coordinatesNew,tesselationNew];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_UniformTessellationChange"<>ToString[tesselationNew]<>"_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]


partialTessellationChangeHP[healpixDir_,outputDir_,tesselationNew_,point1_,point2_,marginLat1_,marginLong1_,marginLat2_,marginLong2_]:=Module[{hpPath,coordinates,coordinatesNew,coordinatesBoundaries,coordinatesChanged1,coordinatesChanged2,coordinatesChanged,dataPath,dataRat,zeroList,dataNew,dat},
hpPath = ToFileName[healpixDir,"testXYZ"];
init[hpPath,tesselation];
coordinates=getHPcoord[healpixringxyz];
init[hpPath,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyz];
coordinatesBoundaries={{point1[[1]]-marginLat1,point1[[1]]+marginLat1},{point1[[2]]-marginLong1,point1[[2]]+marginLong1}};
coordinatesChanged1=Select[coordinates,#[[1]]>=coordinatesBoundaries[[1,1]]&&#[[1]]<=coordinatesBoundaries[[1,2]]&&#[[2]]>=coordinatesBoundaries[[2,1]]&&#[[2]]<=coordinatesBoundaries[[2,2]]&];
coordinatesBoundaries={{point2[[1]]-marginLat2,point2[[1]]+marginLat2},{point2[[2]]-marginLong2,point2[[2]]+marginLong2}};coordinatesChanged2=Select[coordinates,#[[1]]>=coordinatesBoundaries[[1,1]]&&#[[1]]<=coordinatesBoundaries[[1,2]]&&#[[2]]>=coordinatesBoundaries[[2,1]]&&#[[2]]<=coordinatesBoundaries[[2,2]]&];
coordinatesChanged=Sort[Join[coordinatesChanged1,coordinatesChanged2]];
dataPath = ToFileName[outputDir];
dataRat=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRat,"bad"]];
If[zeroList!={},
dataRat=Last[With[{x=zeroList[[#]]},dataRat=ReplacePart[dataRat,x->Null]]&/@Range[Length[zeroList]]]];
dataNew=partialTessellationChange[coordinatesChanged,dataRat,coordinates,coordinatesNew,tesselationNew];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_PartialTessellationChange"<>ToString[tesselationNew]<>"_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]


uniformPixelShapeChangeHP[healpixDir_,resolution_,outputDir_]:=Module[{dataPath,hpPath,dataRat,zeroList,dataNew,dat},
(*It doesn't actually use healpixDir and resolution here, but uniforSizeChange funtion does internally, so it's here to remind you what it does*)
dataPath = ToFileName[outputDir];
hpPath = ToFileName[healpixDir,"testXYZ"];
init[hpPath,tesselation];
coordinates=getHPcoord[healpixringxyz];
dataRat=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRat,"bad"]];
If[zeroList!={},
dataRat=Last[With[{x=zeroList[[#]]},dataRat=ReplacePart[dataRat,x->Null]]&/@Range[Length[zeroList]]]];
dataNew=uniformSizeChange[dataRat];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_PixelShapeChangeRes"<>ToString[resolution]<>"_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]


transformToRibbonCoordHP[healpixDir_,outputDir_,resolution_]:=Module[{dataPath,dataRatesOriginal,zeroList,coordinatesOriginal,healpixPixelsClassifiedOriginal,originalFields,healpixringxyzTransformed,coordinatesChanged,hpChangedValues,newFields,dataNew,dat},
(*WITHOUT INCLUDING NON OBSERVED PIXELS-TO BE ADDED (only use for 4,5 and 6 energy steps)*)
Which[
energyStep==4,lr=220.51;br=29.96,
energyStep==5,lr=218.08;br=38.44,
energyStep==6,lr=214.68;br=34.13];
dataPath=ToFileName[outputDir];
dataRatesOriginal=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRatesOriginal,"bad"]];
coordinatesOriginal=getHPcoord[healpixringxyz];
healpixPixelsClassifiedOriginal=classifyHEALPixPixels[coordinatesOriginal,resolution];originalFields=Flatten[healpixPixelsClassifiedOriginal,1];
healpixringxyzTransformed=transformToRibbonCoord[healpixringxyz[[#]],lr,br]&/@Range[Length[healpixringxyz]];
coordinatesChanged=getHPcoord[healpixringxyzTransformed];
hpChangedValues=classifyHEALPixValues[dataRatesOriginal,coordinatesChanged,resolution,dataRatesOriginal];
newFields=Flatten[hpChangedValues,1];
dataNew=ConstantArray[0,Length[healpixringxyz]];
Module[{idx1=#,list},
list=originalFields[[idx1]];
Module[{idx2=#},
dataNew[[list[[idx2]]]]=newFields[[idx1]]]&/@Range[Length[list]];
]&/@Range[Length[originalFields]];
dataNew=Flatten[dataNew];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_RibbonCentered_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]


classifyHealPixData[coordinatesDataGoodWithValues_,resolution_]:=Module[{colati,long,rows,columns,val,i,map},
long=coordinatesDataGoodWithValues[[;;,1]][[;;,2]];
colati=coordinatesDataGoodWithValues[[;;,1]][[;;,1]];
val=coordinatesDataGoodWithValues[[;;,2]];
rows=Map[fun`classify[#, resolution[[1]], 180]&,colati];
columns =Map[fun`classify[#, resolution[[2]], 360]&,long];
map = ConstantArray[0,{resolution[[1]],resolution[[2]]}];
For[i=1,i<Length[rows]+1,i++,
Which[
Length[map[[rows[[i]],columns[[i]]]]]==0, map[[rows[[i]],columns[[i]]]]={i,colati[[i]],long[[i]],val[[i]]},
Dimensions[map[[rows[[i]],columns[[i]]]]]=={4}, map[[rows[[i]],columns[[i]]]]=Append[{map[[rows[[i]],columns[[i]]]]},{i,colati[[i]],long[[i]],val[[i]]}],
Dimensions[map[[rows[[i]],columns[[i]]]]][[1]]>=2,map[[rows[[i]],columns[[i]]]]=Append[map[[rows[[i]],columns[[i]]]],{i,colati[[i]],long[[i]],val[[i]]}]];
];
map]


totalHealPixValues[map_,resolution_]:=Module[{totalValues,j,i},
totalValues=ConstantArray[0,{resolution[[1]],resolution[[2]]}];
j=1;
While[j<resolution[[1]]+1,
i=1;
While[i<resolution[[2]]+1,
Which[
Length[map[[j,i]]]==0,totalValues[[j,i]]={j,i,0},
Dimensions[map[[j,i]]]=={4},totalValues[[j,i]]={j,i,map[[j,i]][[4]]},
Dimensions[map[[j,i]]][[1]]>=2,totalValues[[j,i]]={j,i,Total[map[[j,i]][[All,4]]]}];
i++];
j++];
totalValues]


makeHpToIbMap[coordinatesDataGoodWithValues_,resolution_]:=Module[{coordinates,mapClassified},
mapClassified = classifyHealPixData[coordinatesDataGoodWithValues,resolution];
Print["Total data count: ",Total[coordinatesDataGoodWithValues[[;;,2]]]];
totalHealPixValues[mapClassified,resolution]
]
oneElementArea[idx_,dtheta_,dfi_]:=Module[{start,end,thetas},
start=onePixelResolution[[1]]/2;
end=180-start;
thetas=Range[start,end,dtheta];
Sin[thetas[[idx]] Degree]*dtheta Degree*dfi Degree
]//N
normalize[onePixelResolution_,map_]:=Module[{oneRow,wages,normalizationFactor,nFactor},
oneRow=oneElementArea[#,onePixelResolution[[1]],onePixelResolution[[2]]]&/@Range[180/onePixelResolution[[1]]];
wages=ConstantArray[oneRow[[#]],(360/onePixelResolution[[2]])]&/@Range[180/onePixelResolution[[1]]];
(*Total[Total[wages]]/(4*Pi)*)
normalizationFactor=Total[Total[map*wages]]/Total[Total[wages]];
Print[normalizationFactor];
nFactor=normalizationFactor;
map/normalizationFactor
]
(*Sum of the wages should be equal to 4Pi*)


changeGridHP[healpixDir_,resolution_,outputDir_]:=Module[{dataPath,dataCts,dataExp,coordinates,data,coordinatesWithValues,mapFinalCounts,mapFinalExp,rates,zeroList1,zeroList2,zeroList,nonZeroList,ratesHP,mapFinalHP2IB,dat},
dataPath = ToFileName[outputDir];
dataCts=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,2]];
dataExp=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,3]];
coordinates=getHPcoord[healpixringxyz];
(*Make counts map*)
data=dataCts;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalCounts=makeHpToIbMap[coordinatesWithValues,resolution];
(*Make exp map*)
data=dataExp;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalExp=makeHpToIbMap[coordinatesWithValues,resolution];
rates=Flatten[mapFinalCounts[[;;,;;,3]]];
zeroList1=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0]];
zeroList2=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0.]];
zeroList=Join[zeroList1,zeroList2];
nonZeroList=Complement[Range[Length[rates]],zeroList];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->0]]&/@zeroList]];
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->(Flatten[mapFinalCounts[[;;,;;,3]]][[x]]/Flatten[mapFinalExp[[;;,;;,3]]][[x]])]]&/@nonZeroList];
rates=Partition[rates,resolution[[2]]];
ratesHP=rates;
rates=normalize[onePixelResolution,rates];
rates=Flatten[rates];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->Null]]&/@zeroList]];
rates=Partition[rates,resolution[[2]]];
mapFinalHP2IB=rates;
dat=ExportString[mapFinalHP2IB,"Table"];
Export[ToFileName[outputDir]<>"normalizedHP2IB_"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"_"<>ToString[onePixelResolution[[1]]]<>"_"<>ToString[onePixelResolution[[2]]]<>".txt",dat];
]

changeGridHP2[dataFlux_,resolution_,outputDir_]:=Module[{coordinates,data,coordinatesWithValues,mapFinalFlux,mapFinalFluxNormalized,flux,mapFinalHP2IB,dat},
coordinates=getHPcoord[healpixringxyz];
data=dataFlux;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalFlux=makeHpToIbMap[coordinatesWithValues,resolution];
flux= Partition[Flatten[mapFinalFlux[[;;,;;,3]]],resolution[[2]]];
(*mapFinalFluxNormalized=normalize[onePixelResolution,flux];*)
mapFinalFluxNormalized=mapFinalFluxNormalized/. 0.->Null;
dat=ExportString[mapFinalHP2IB,"Table"];
Export[ToFileName[outputDir]<>"normalized2HP2IB_"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"_"<>ToString[onePixelResolution[[1]]]<>"_"<>ToString[onePixelResolution[[2]]]<>".txt",dat];
]

changeGridRibbonHP[healpixDir_,resolution_,outputDir_]:=Module[{dataPath,dataRatesOriginal,zeroList,coordinatesOriginal,healpixPixelsClassifiedOriginal,originalFields,healpixringxyzTransformed,coordinatesChanged,hpChangedValues,newFields,dataNew,data,coordinatesWithValues,mapFinalRates,rates,mapFinalHP2IB,dat},
(*WITHOUT INCLUDING NON OBSERVED PIXELS-TO BE ADDED (only use for 4,5 and 6 energy steps)*)
Which[
energyStep==4,lr=220.51;br=29.96,
energyStep==5,lr=218.08;br=38.44,
energyStep==6,lr=214.68;br=34.13];
dataPath=ToFileName[outputDir];
dataRatesOriginal=Import[dataPath<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",{"Data"}][[All,4]];
zeroList=Flatten[Position[dataRatesOriginal,"bad"]];
coordinatesOriginal=getHPcoord[healpixringxyz];
healpixPixelsClassifiedOriginal=classifyHEALPixPixels[coordinatesOriginal,resolution];originalFields=Flatten[healpixPixelsClassifiedOriginal,1];
healpixringxyzTransformed=transformToRibbonCoord[healpixringxyz[[#]],lr,br]&/@Range[Length[healpixringxyz]];
coordinatesChanged=getHPcoord[healpixringxyzTransformed];
hpChangedValues=classifyHEALPixValues[dataRatesOriginal,coordinatesChanged,resolution,dataRatesOriginal];
newFields=Flatten[hpChangedValues,1];
dataNew=ConstantArray[0,Length[healpixringxyz]];
Module[{idx1=#,list},
list=originalFields[[idx1]];
Module[{idx2=#},
dataNew[[list[[idx2]]]]=newFields[[idx1]]]&/@Range[Length[list]];
]&/@Range[Length[originalFields]];
dataNew=Flatten[dataNew];
data=dataNew;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinatesChanged},data[[;;,2]]]];
mapFinalRates=makeHpToIbMap[coordinatesWithValues,resolution];
rates=Flatten[mapFinalRates[[;;,;;,3]]];
rates=Partition[rates,resolution[[2]]];
rates=normalize[onePixelResolution,rates];
mapFinalHP2IB=rates;
dat=ExportString[mapFinalHP2IB,"Table"];
Export[ToFileName[outputDir]<>"normalizedHP_Ribbon_2IB_"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"_"<>ToString[onePixelResolution[[1]]]<>"_"<>ToString[onePixelResolution[[2]]]<>".txt",dat];
]

changeGridR1deg[ibexHi_]:=Module[{dataCts,dataExp,colati,long,coordinates,data,coordinatesWithValues,mapFinalCounts,mapFinalExp,rates,zeroList1,zeroList2,nonZeroList,mapFinalIB,dat},
dataCts=Flatten[ToExpression["ctsE"<>ToString[energyStep]]/.ibexHi];
dataExp=Flatten[ToExpression["expE"<>ToString[energyStep]]/.ibexHi];
colati =lati2colatiDeg[Flatten[ToExpression["eclatE"<>ToString[energyStep]]/.ibexHi]];
long =Flatten[ToExpression["eclonE"<>ToString[energyStep]]/.ibexHi];
coordinates=Transpose[Append[{colati},long]];
(*Make counts map*)
data=dataCts;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalCounts=makeHpToIbMap[coordinatesWithValues,resolution];
(*Make exp map*)
data=dataExp;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalExp=makeHpToIbMap[coordinatesWithValues,resolution];
rates=Flatten[mapFinalCounts[[;;,;;,3]]];
zeroList1=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0]];
zeroList2=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0.]];
zeroList=Join[zeroList1,zeroList2];
nonZeroList=Complement[Range[Length[rates]],zeroList];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->0]]&/@zeroList]];
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->1.*(Flatten[mapFinalCounts[[;;,;;,3]]][[x]]/Flatten[mapFinalExp[[;;,;;,3]]][[x]])]]&/@nonZeroList];
rates=Partition[rates,resolution[[2]]];
rates=normalize[onePixelResolution,rates];
rates=Flatten[rates];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->Null]]&/@zeroList]];
rates=Partition[rates,resolution[[2]]];
mapFinalIB=rates;
dat=ExportString[mapFinalIB,"Table"];
Export[ToFileName[outputDir]<>"normalizedIB_"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"_"<>ToString[onePixelResolution[[1]]]<>"_"<>ToString[onePixelResolution[[2]]]<>".txt",dat];
]

bin1to6[list_]:=Module[{i=#},list[[i;;i+7]]]&/@Range[1,14400,8]

changeGridR6deg[ibexHi_]:=Module[{dataCts,dataExp,cts,exp,colati,long,coordinates,data,coordinatesWithValues,mapFinalCounts,mapFinalExp,rates,zeroList1,zeroList2,nonZeroList,mapFinalIB2, dataIB1deg,dataIB1degSorted,long6deg,colat6deg,cts6deg,exp6deg,dat},
dataCts=Flatten[ToExpression["ctsE"<>ToString[energyStep]]/.ibexHi];
dataExp=Transpose[Append[{Range[Length[dataCts]]},Flatten[ToExpression["expE"<>ToString[energyStep]]/.ibexHi]]];
dataCts=Transpose[Append[{Range[Length[dataCts]]},Flatten[ToExpression["ctsE"<>ToString[energyStep]]/.ibexHi]]];
colati =lati2colatiDeg[Flatten[ToExpression["eclatE"<>ToString[energyStep]]/.ibexHi]];
long =Flatten[ToExpression["eclonE"<>ToString[energyStep]]/.ibexHi];
cts=dataCts[[;;,2]];
exp=dataExp[[;;,2]];
dataIB1deg=Transpose[Append[{long,colati,cts},exp]];
dataIB1degSorted=SortBy[dataIB1deg,First];
long6deg=Flatten[bin1to6[dataIB1degSorted[[;;,1]]]];
colat6deg=Flatten[bin1to6[dataIB1degSorted[[;;,2]]]];
cts6deg=Mean[bin1to6[dataIB1degSorted[[;;,3]]][[#]]]&/@Range[Length[bin1to6[dataIB1degSorted[[;;,3]]]]];
cts6deg=Flatten[ConstantArray[cts6deg[[#]],8]&/@Range[Length[cts6deg]]];
exp6deg=Mean[bin1to6[dataIB1degSorted[[;;,4]]][[#]]]&/@Range[Length[bin1to6[dataIB1degSorted[[;;,4]]]]];
exp6deg=Flatten[ConstantArray[exp6deg[[#]],8]&/@Range[Length[exp6deg]]];
dataCts=cts6deg;
dataExp=exp6deg;
coordinates=Transpose[Append[{colat6deg},long6deg]];
(*Make counts map*)
data=dataCts;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalCounts=makeHpToIbMap[coordinatesWithValues,resolution];
(*Make exp map*)
data=dataExp;
data=Transpose[Append[{Range[Length[data]]},data]];
coordinatesWithValues=Transpose[Append[{coordinates},data[[;;,2]]]];
mapFinalExp=makeHpToIbMap[coordinatesWithValues,resolution];
rates=Flatten[mapFinalCounts[[;;,;;,3]]];
zeroList1=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0]];
zeroList2=Flatten[Position[Flatten[mapFinalExp[[;;,;;,3]]],0.]];
zeroList=Join[zeroList1,zeroList2];
nonZeroList=Complement[Range[Length[rates]],zeroList];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->0]]&/@zeroList]];
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->1.*(Flatten[mapFinalCounts[[;;,;;,3]]][[x]]/Flatten[mapFinalExp[[;;,;;,3]]][[x]])]]&/@nonZeroList];
rates=Partition[rates,resolution[[2]]];
rates=normalize[onePixelResolution,rates];
rates=Flatten[rates];
If[zeroList!={},
rates=Last[Module[{x=#},rates=ReplacePart[rates,x->Null]]&/@zeroList]];
rates=Partition[rates,resolution[[2]]];
mapFinalIB2=rates;
dat=ExportString[mapFinalIB2,"Table"];
Export[ToFileName[outputDir]<>"normalizedIB2_"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"_"<>ToString[onePixelResolution[[1]]]<>"_"<>ToString[onePixelResolution[[2]]]<>".txt",dat];
]
