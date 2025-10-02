(* ::Package:: *)

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
]


collHi[R_,alpha_,Ctrans_]:=
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


collBkg[R_]:=
Module[{h2d},
h2d = 6.67178; (*height to deimeter ratio of the holes in the collimator*)
(2/Pi)*( ArcCos[ h2d * Tan[R Degree] ] - h2d * Tan[R Degree] * Sqrt[1 - (h2d*Tan[R Degree])^2])
]


coll[tesselation_,healpixringxyz_,angle1_,angle2_,colPixelsDistances_]:=Module[{collCentreVector,collPixelVectors,alphas,alphas2,colim1,colim2,colimBkg,colimBkg2,domega,colimatorBeforeFilter,colimatorAfterFilter1,collPixelVectors2,a1,a1bkg,n1,n1bkg,a2,a2bkg,n2,n2bkg,colimatorAfter2},
collCentreVector=makeVec[angle2*1. Degree,angle1*1. Degree];
collPixelVectors=healpixringxyz[[colPixelsDistances[[;;,1]]]];
alphas=VectorAngle[collCentreVector,#]&/@collPixelVectors;
colim1=collHi[colPixelsDistances[[#]][[2]]/Degree,alphas[[#]]/Degree,1]&/@Range[Length[alphas]];
colimBkg=collBkg[colPixelsDistances[[#]][[2]]]&/@Range[Length[colPixelsDistances]];
(*colim1=Transpose[Append[{colPixelsDistances[[;;,1]]},colim1]];
colimatorBeforeFilter=Transpose[Append[{colPixelsDistances[[;;,1]],colPixelsDistances[[;;,2]]},colim1[[;;,2]]]];
colimatorAfterFilter1=Select[colimatorBeforeFilter,#[[3]]>0&];
domega = Pi/(3*tesselation^2)*1./Degree;
a1=1.Total[Map[colimatorAfterFilter1[[#,3]]*domega&,Range[Length[colimatorAfterFilter1]]]];
n1=1.Total[Map[colimatorAfterFilter1[[#,3]]/a1&,Range[Length[colimatorAfterFilter1]]]];
a2=a1*n1;
collPixelVectors2=healpixringxyz[[colimatorAfterFilter1[[;;,1]]]];
alphas2=VectorAngle[collCentreVector,#]&/@collPixelVectors2;
colimatorAfter2=collHi[colimatorAfterFilter1[[#]][[2]]/Degree,alphas2[[#]]/Degree,a2]&/@Range[Length[alphas2]];
colimatorAfter2=Transpose[Append[{colimatorAfterFilter1[[;;,1]],colimatorAfterFilter1[[;;,2]]},colimatorAfter2]];
n2=1.Total[Map[colimatorAfter2[[#,3]]&,Range[Length[colimatorAfter2]]]];
colimatorAfter2[[;;,{1,3}]]*)
domega = Pi/(3*tesselation^2)*1./Degree;

a1=Total[colim1*domega];
n1=Total[colim1/a1];
a2=a1*n1;

a1bkg=Total[colimBkg*domega];
n1bkg=Total[colimBkg/a1bkg];
a2bkg=a1bkg*n1bkg;

{Transpose[Append[{colPixelsDistances[[;;,1]]},(colim1/a2)]], Transpose[Append[{colPixelsDistances[[;;,1]]},(colimBkg/a2bkg)]]}
]


lati2colatiDeg[\[Phi]_]:= 90-\[Phi]
lati2colatiRad[\[Phi]_]:= Pi/2-\[Phi]


choseRing[rotationAxisAng_]:=Module[{oneRingSpare, rotationAxis, maxVisibilityRange, minVisibilityRange, deltaVisiblityRange, angLengthsFromRotAxis, visibilityRangePixels},
rotationAxis=makeVec[rotationAxisAng[[1]]*1. Degree,rotationAxisAng[[2]]*1. Degree];
oneRingSpare=(Pi/(3*tesselation))*1.;
maxVisibilityRange = Cos[(scanRadius-colRadius)Degree]*1.+oneRingSpare;
minVisibilityRange=Cos[(scanRadius+colRadius)Degree]*1.-oneRingSpare;
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

 calcOneOrbit[ibexLatitude_,ibexLongitude_,exposuretime_,counts_,backgroundRate_,visibilityRangePixels_,healpixringxyz_]:=Module[{measurementIndex=#,angle1,angle2,exposuretimeValue,countValue,backgroundRateValue,colPixelsDistances,exposuretimeValueNormalized,colValues,colValuesBkg,nonColPixelsWithZeros,collimatorLevel,colPixelsWithValues},
angle1 = ibexLatitude[[measurementIndex]];
angle2=ibexLongitude[[measurementIndex]];
exposuretimeValue=exposuretime[[measurementIndex]];
countValue=counts[[measurementIndex]];
backgroundRateValue=backgroundRate[[measurementIndex]];
colPixelsDistances=calculateLengthsForColPixels[angle1,angle2,visibilityRangePixels];
{colValues,colValuesBkg} = coll[tesselation,healpixringxyz,angle1,angle2,colPixelsDistances];
If[exposuretimeValue!= 0,
exposuretimeValueNormalized=exposuretimeValue/((Total[ConstantArray[exposuretimeValue,Length[colValues]]])/exposuretimeValue),exposuretimeValueNormalized=exposuretimeValue];
nonColPixelsWithZeros=#->{0,0,0}&/@Complement[Range[Length[healpixringxyz]],colValues[[;;,1]]];
colPixelsWithValues=MapThread[(Module[{idx=#1,idx2=#2},idx[[1]]->{idx[[2]]*countValue, exposuretimeValueNormalized, (idx[[2]]*countValue - idx2[[2]]*(backgroundRateValue*exposuretimeValue))}])&,{colValues, colValuesBkg}];
(*Print["Signal, so count - (background*exp) from measurement: ", countValue-backgroundRateValue*exposuretimeValue];*)
collimatorLevel=Sort[Join[nonColPixelsWithZeros,colPixelsWithValues]]

(*Print["Count value: ", countValue, " , exposure time  ", exposuretimeValue, " and : count - (background*exp) from measurement: ", countValue-backgroundRateValue*exposuretimeValue];
Print["Count value and exposure time value on the collimator: ",collimatorLevel[[;;,2]][[;;,1;;2]]//Total, ", count - (background * exposure): ",(collimatorLevel[[;;,2]][[;;,3]]//Total) ];*)
]&/@Range[Length[ibexLatitude]];


calculateRatesMathematica[counts_,exposures_]:=Module[{idx=#,rates},If[exposures[[idx]]!=0,rates=counts[[idx]]/exposures[[idx]],rates="bad"]]&/@Range[Length[counts]];


calculateENAFluxMathematica[signal_, exposuretime_]:=Module[{centralEnergies, centralEnergy, geometricFactorTriples, geometricFactor, omega,flux},
geometricFactorTriples = {0.00013,0.00037,0.00073,.0014,0.0025,0.0042};
centralEnergies = {0.45, 0.71, 1.10, 1.74, 2.73, 4.29};
geometricFactor = geometricFactorTriples[[energyStep - 3]];
centralEnergy   = centralEnergies[[energyStep - 3]];
flux=MapThread[If[#2!=0,#1/(#2*geometricFactor*centralEnergy),"bad"]&,{signal,exposuretime}]
]


calculateRatesHealpy[counts_,exposures_]:=Module[{idx=#,rates},If[exposures[[idx]]!=0,rates=counts[[idx]]/exposures[[idx]],rates=Null]]&/@Range[Length[counts]];


calculateENAFluxHealpy[signal_, exposuretime_]:=Module[{centralEnergies, centralEnergy, geometricFactorTriples, geometricFactor, omega,flux},
geometricFactorTriples = {0.00013,0.00037,0.00073,.0014,0.0025,0.0042};
centralEnergies = {0.45, 0.71, 1.10, 1.74, 2.73, 4.29};
geometricFactor = geometricFactorTriples[[energyStep - 3]];
centralEnergy   = centralEnergies[[energyStep - 3]];
flux=MapThread[If[#2!=0,#1/(#2*geometricFactor*centralEnergy),0]&,{signal,exposuretime}]
]
