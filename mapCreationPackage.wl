(* ::Package:: *)

loadData[baseDir_,ibexDataDir_,ibexDataName_,spinAxDir_,outputDir_]:=Module[{ibexFiles},
ibexFiles = Map[FileNameTake,FileNames["*.qABC", ToFileName[ibexDataDir,ToString[ibexDataName]]]];
Get[ToFileName[spinAxDir, "spinAxEclTabWarsaw.m"]];
Import[baseDir<>"/geometricPackage.wl"];
Import[baseDir<>"/mapCreationPackage.wl"];
ibexHi = Map[ fun`ibexDataRead[ToFileName[ibexDataDir,ibexDataName], #]&, ibexFiles]
]


Options[run] = {"smoothing" -> Null};

run[ibexHi_, healpixDir_, OptionsPattern[]] :=
Module[{smoothing, domega, hpPath, result, mapCountsMain, mapSignalMain, mapExposuresMain, mapRatesMain, mapENAFluxMain, geometricFactorTriples, centralEnergies, geometricFactor, centralEnergy},
smoothing = If[OptionValue["smoothing"] === "Gaussian", "Gaussian", Null];

hpPath = ToFileName[healpixDir, "testXYZ"];
domega = Pi/(3*tesselation^2)*1.;
init[hpPath, tesselation];

AbsoluteTiming[
With[{sm = smoothing},
result =
ParallelMap[
Module[{
orbitCount = #, exposuretime, counts, signal, backgroundRate, ibexLatitude, ibexLongitude, rotationAxisAng, visibilityRangePixels, oneOrbitMap, mapCounts, mapExposures, mapSignal},
exposuretime = ToExpression["expE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
exposuretime = exposuretime*10^-3*1.; (* seconds *)
counts = ToExpression["ctsE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
backgroundRate = ToExpression["bkgE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLatitude = ToExpression["eclatE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
ibexLongitude = ToExpression["eclonE" <> ToString[energyStep]] /. ibexHi[[orbitCount]];
rotationAxisAng = { spinAxEcLon /. ibexHi[[orbitCount]], spinAxEcLat /. ibexHi[[orbitCount]]};
visibilityRangePixels = choseRing[rotationAxisAng];
Print["Loading orbit ", orbNo /. ibexHi[[orbitCount]]];

oneOrbitMap = calcOneOrbit[ibexLatitude, ibexLongitude, exposuretime, counts, backgroundRate,visibilityRangePixels, healpixringxyz, sm];
mapCounts = Total[(Module[{element = #}, First[Last[#]] & /@ element] & /@ oneOrbitMap)];
mapExposures =Total[(Module[{element = #}, Last[#][[2]] & /@ element] & /@ oneOrbitMap)];
mapSignal = Total[(Module[{element = #}, Last[#][[3]] & /@ element] & /@ oneOrbitMap)];
{mapCounts, mapExposures, mapSignal}] &,Range[Length[ibexHi]]
	];
		];	
	mapCountsMain = Total[result[[All, 1]]];
	mapExposuresMain = Total[result[[All, 2]]];
	mapSignalMain = Total[result[[All, 3]]];
	mapSignalMain = mapSignalMain /. x_?NumericQ /; x < 0 -> 0;
	mapCountsMain = (mapCountsMain/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
	mapExposuresMain = (mapExposuresMain/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
	mapSignalMain = (mapSignalMain/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
	
	mapRatesMain = If[mapExposuresMain[[#]]==0||mapExposuresMain[[#]]===Null,Null,mapSignalMain[[#]]/(mapExposuresMain[[#]])]&/@Range[Length[mapExposuresMain]];
	geometricFactorTriples = {0.00013, 0.00041, 0.00075, 0.0013, 0.0024, 0.0045}; (*source: https://ibex.princeton.edu/sites/g/files/toruqf1596/files/documents/IBEX_CMAD_signed_final.pdf*)
	centralEnergies = {0.45, 0.71, 1.08, 1.85, 2.70, 4.09}; (*source: above*)
	geometricFactor = geometricFactorTriples[[energyStep]];
	centralEnergy   = centralEnergies[[energyStep]];
	mapENAFluxMain = If[mapExposuresMain[[#]]==0||mapExposuresMain[[#]]===Null,Null,mapSignalMain[[#]]/(mapExposuresMain[[#]]*geometricFactor*centralEnergy)]&/@Range[Length[mapExposuresMain]];
	mapRatesMain = (mapRatesMain/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
	mapENAFluxMain = (mapENAFluxMain/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
   
   {mapCountsMain, mapExposuresMain, mapSignalMain, mapRatesMain, mapENAFluxMain}
   ]
  ];


exportData[mapCountsMain_, mapExposuresMain_,mapSignalMain_,mapRatesMain_,mapENAFluxMain_,tesselation_,energyStep_,outputDir_]:=Module[{datHealpy},
datHealpy=ExportString[Transpose[Append[{Range[Length[mapCountsMain]],mapCountsMain,mapExposuresMain,mapSignalMain,mapRatesMain},mapENAFluxMain]],"Table"];
Export[ToFileName[outputDir]<>"data_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>"Null.txt",datHealpy];
]


fun`selectGoodData[data_, indran_, pwords_]:= Module[{ddd,bbb},
ddd = Transpose[data];
(*select spin angle, bin positions in the sky, and data relevant for the given energy step*)
bbb=Transpose[Join[{Last[ddd]},ddd[[1;;2]], ddd[[indran[[1]];;indran[[2]]]]]];
bbb=Transpose[Select[bbb, #[[4]] >= 0 && #[[5]] >=0&]];
Table[pwords[[i]]-> bbb[[i]], {i, Length[bbb]}]
]


fun`ibexDataRead[dir_, fname_]:= Module[{data0, data, spinax,
indranE1 = {3, 6}, indranE2, indranE3, indranE4, indranE5, indranE6, res},
(*define data index ranges fo individual energy steps*)
indranE2 = Last[indranE1]+{1,4 };
indranE3 = Last[indranE2]+{1,4 };
indranE4 = Last[indranE3]+{1,4 };
indranE5 = Last[indranE4]+{1,4 };
indranE6 = Last[indranE5]+{1,4 };
(*read in the data file*)
data0 = Import[ToFileName[dir, fname], "Table"];
(*toss the header*)
data = Rest[data0];
(*add spin angle, make it the last column*)
data = Transpose[Append[Transpose[data], (Range[360]-1)0.5]];
(*import spin angle pointing and time of observations*)
spinax=Select[Drop[spinAxEclTab`waw,2], ( "0"<>#[[1]]) ==(FileBaseName[fname])&]//First;
(*call function to extract data for individual energy steps*)
res=Map[fun`selectGoodData[data,#[[1]], #[[2]]]&, {{indranE1, {spinAngE1,eclatE1, eclonE1, ctsE1, expE1,  bkgE1, bkgErrE1}},
{indranE2, {spinAngE2,eclatE2, eclonE2, ctsE2, expE2, bkgE2, bkgErrE2}},
{indranE3, {spinAngE3,eclatE3, eclonE3, ctsE3, expE3, bkgE3, bkgErrE3}},
{indranE4, {spinAngE4,eclatE4, eclonE4, ctsE4, expE4,  bkgE4, bkgErrE4}},
{indranE5, {spinAngE5,eclatE5, eclonE5, ctsE5, expE5,  bkgE5, bkgErrE5}},
{indranE6, {spinAngE6,eclatE6, eclonE6, ctsE6, expE6,  bkgE6, bkgErrE6}}}];
(*Collect together information universal for all energy steps and the data*)
Flatten[Join[{orbNo-> "o0"<>spinax[[1]],(*orbit label*)
filename-> fname, (*parent data file name*)
dirname-> dir,(*parent data directory (as a list)*)
obsTime-> Mean[spinax[[2;;3]]],(*decimal year of the observation time*)
spinAxEcLon-> spinax[[4]],(*ecliptic longitude of the spin axis for the IBEX orbit*)
spinAxEcLat-> spinax[[5]] (*ecliptic latitude of the spin axis for the IBEX orbit*)
}, res],1]
]


init[path_,tesselation_]:=Module[{z},
healpixringxyz=Import[path<>ToString[tesselation]<>".dat","Table"];
mapCountsMain=0;
mapExposuresMain=0;
mapRatesMain=0;
mapSignalMain=0;
mapENAFluxMain=0;
]


collHi[R_,alpha_,Ctrans_]:= (*function specific to IBEX-Hi collimator, due to Marzena Kubiak*)
Module[{Beta,Gamma,rest,aCorn,aBase,yCorn,yBase,R0,Betaint, Betamod,yR},
Beta=alpha -15.0;
Betaint=IntegerPart[Beta];
Betamod=Mod[IntegerPart[Beta],60];
rest=Beta-Betaint;
Gamma=Abs[Betamod+rest];
If[Gamma>30.0, Gamma=Gamma-60.0];
R0 =R;
aBase=R0*Cos[Gamma Degree];
Which[
aBase>7.35,yR=0.0,
aBase==0,yR=1,
aBase<7.35,
yBase=1.00014+(-0.180305+0.0060478*aBase)*aBase;
aCorn=R0*2.0*Cos[Gamma Degree]/Sqrt[3.0];
yCorn=0.993855+aCorn*(-0.136622+(-0.0110679+0.00158511*aCorn)*aCorn);
yR=R0*((yCorn-yBase)/(aCorn-aBase))+yCorn-aCorn*((yCorn-yBase)/(aCorn-aBase))];
yR=yR/Ctrans
]


(*collBkg[R_]:= (*function specific to IBEX background monitor, currently not used - have to implemented separate PSF*)
Module[{h2d}, 
h2d = 6.67178; (*height to deimeter ratio of the holes in the collimator*)
(2/Pi)*( ArcCos[ h2d * Tan[R Degree] ] - h2d * Tan[R Degree] * Sqrt[1 - (h2d*Tan[R Degree])^2])
] (*https://link.springer.com/article/10.1007/s11214-008-9439-8*)*)


coll[tesselation_,healpixringxyz_,angle1_,angle2_,colPixelsDistances_,smoothing_]:=Module[{collResult,collCentrePix,collCentreVector,collPixelVectors,alphas,colim1,domega,a1,a1bkg,n1,a2},
(*collCentreVector=makeVec[angle2*1. Degree,angle1*1. Degree];*)
(*collPixelVectors=healpixringxyz[[colPixelsDistances[[;;,1]]]];*)
(*alphas=angleVec[collCentreVector,#]&/@collPixelVectors;
colim1=collHi[colPixelsDistances[[#]][[2]]/Degree,alphas[[#]]/Degree,1]&/@Range[Length[alphas]];*)
Which[
ToString[smoothing] === "Gaussian",
colim1 = GaussianWeight[colPixelsDistances[[;;,2]],Cos[colRadius*1. Degree]];
domega = Pi/(3*tesselation^2)*1.;
colim1=colim1*domega;
a1=Total[colim1*domega];
n1=Total[colim1/a1];
a2=a1*n1;
collResult=Transpose[Append[{colPixelsDistances[[;;,1]]},(colim1/a2)]];
,
TrueQ[smoothing === Null],
collCentrePix=SortBy[colPixelsDistances,Min][[1]][[1]];
collResult=Transpose[{colPixelsDistances[[;;,1]],ConstantArray[0,Length[colPixelsDistances[[;;,1]]]]}];
collResult[[Position[collResult[[;;,1]],collCentrePix][[1,1]]]][[2]]=1;
];

(*colimBkg=colimBkg*domega; (*part of the background monitor implementation*)
colimBkg=collBkg[colPixelsDistances[[#]][[2]]]&/@Range[Length[colPixelsDistances]];
a1bkg=Total[colimBkg*domega];
n1bkg=Total[colimBkg/a1bkg];
a2bkg=a1bkg*n1bkg; 
{Transpose[Append[{colPixelsDistances[[;;,1]]},(colim1/a2)]], Transpose[Append[{colPixelsDistances[[;;,1]]},(colimBkg/a2bkg)]]}*)

collResult
]


Options[GaussianWeight]={"Units"->"Radians"};
GaussianWeight[alpha_,width_,OptionsPattern[]]:=Module[{},
If[OptionValue["Units"]==="Degrees",alpha=alpha Degree;width=width Degree;];
Exp[-(alpha^2)/(2 width^2)]]


lati2colatiDeg[\[Phi]_]:= 90-\[Phi]
lati2colatiRad[\[Phi]_]:= Pi/2-\[Phi]


choseRing[rotationAxisAng_]:=Module[{deltaPix, oneRingSpare, rotationAxis, maxVisibilityRange, minVisibilityRange, deltaVisiblityRange, angLengthsFromRotAxis, visibilityRangePixels},
rotationAxis=makeVec[rotationAxisAng[[1]]*1. Degree,rotationAxisAng[[2]]*1. Degree];
oneRingSpare  = Pi/(3. tesselation^2);        (* sr, Healpix pixel area *)
deltaPix  = Sqrt[oneRingSpare/Pi];           (* rad, approximate pixel radius *)
maxVisibilityRange =
  Cos[(scanRadius - colRadius) Degree - deltaPix];
minVisibilityRange =
  Cos[(scanRadius + colRadius) Degree + deltaPix];
visibilityRangePixels=Select[Transpose[Prepend[
{Dot[rotationAxis,#]&/@healpixringxyz},Range[Length[healpixringxyz]]]],#[[2]]<=maxVisibilityRange&&#[[2]]>=minVisibilityRange&][[All,1]]
]


calculateLengthsForColPixels[angle1_,angle2_,visibilityRangePixels_]:=Module[{vectorLatLong,angLengths,colPixels,lengthCosine},
vectorLatLong = makeVec[angle2*1. Degree,angle1*1. Degree];
angLengths=Transpose[Prepend[
{vectorLatLong . healpixringxyz[[#]]&/@visibilityRangePixels},visibilityRangePixels]];
colPixels=Select[angLengths,#[[2]]>=Cos[colRadius*1. Degree]&];
Transpose[Append[{colPixels[[;;,1]]},If[colPixels[[#,2]]>1.,colPixels[[#,2]]=1;ArcCos[colPixels[[#,2]]],ArcCos[colPixels[[#,2]]]]&/@Range[Length[colPixels]]]]
]

 calcOneOrbit[ibexLatitude_,ibexLongitude_,exposuretime_,counts_,backgroundRate_,visibilityRangePixels_,healpixringxyz_, smoothing_]:=Module[{measurementIndex=#,angle1,angle2,exposuretimeValue,countValue,backgroundRateValue,domega,colPixelsDistances,colValues,nonColPixelsWithZeros,collimatorLevel,colPixelsWithValues},
domega = Pi/(3*tesselation^2)*1.;
angle1 = ibexLatitude[[measurementIndex]];
angle2=ibexLongitude[[measurementIndex]];
countValue=counts[[measurementIndex]];
exposuretimeValue=exposuretime[[measurementIndex]];
backgroundRateValue=backgroundRate[[measurementIndex]];
colPixelsDistances=calculateLengthsForColPixels[angle1,angle2,visibilityRangePixels];
nonColPixelsWithZeros=#->{0,0,0}&/@Complement[Range[Length[healpixringxyz]],colPixelsDistances[[;;,1]]];

(*{colValues,colValuesBkg} = coll[tesselation,healpixringxyz,angle1,angle2,colPixelsDistances];*)
colValues = coll[tesselation,healpixringxyz,angle1,angle2,colPixelsDistances,smoothing];
colPixelsWithValues=MapThread[(Module[{idx=#1},
idx[[1]]->{If[exposuretimeValue==0,Null,idx[[2]]*countValue],
 If[exposuretimeValue==0,Null,idx[[2]]*exposuretimeValue], 
If[exposuretimeValue==0,Null, (idx[[2]]*countValue - (idx[[2]]*backgroundRateValue*exposuretimeValue))]
}])&,{colValues}];


collimatorLevel=Sort[Join[nonColPixelsWithZeros,colPixelsWithValues]]
]&/@Range[Length[ibexLatitude]];
