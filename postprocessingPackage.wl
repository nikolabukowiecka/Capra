(* ::Package:: *)

normalizeRelative[data_, tesselation_, energyStep_] :=
Module[{\[CapitalOmega], w, dataMasked, validMask, meanFluxW, rel, dat},
\[CapitalOmega] = Pi/(3 tesselation^2);
w = ConstantArray[\[CapitalOmega], Length[data]];
dataMasked = data /. Null -> Missing["NoData"];
validMask = NumericQ /@ dataMasked;

(*Weighted mean using valid pixels only *)
meanFluxW =
Total[ Pick[dataMasked * w, validMask] ]/
Total[ Pick[w, validMask] ];

(*Relative flux *)
rel = dataMasked/meanFluxW;

(*Missing back to Null for output*)
rel = rel /. Missing["NoData"] -> Null;

(*Diagnostic: mean should be about 1 *)
validMask = Map[NumericQ, rel];
Print["Weighted mean (should be 1): ", Total@Pick[rel*w, validMask] / Total@Pick[w, validMask]];

rel=rel/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a,Plus[a_?NumericQ,Null]:>a};

dat = ExportString[Transpose[{Range[Length[rel]], rel}],"Table"];
Export[ToFileName[outputDir] <>"data_normalizedRelative_t" <> ToString[tesselation] <>"_" <> ToString[energyStep] <> ".txt",dat];
rel
]


normalizeGDF[data_, dataExposure_, tesselation_, energyStep_] :=
Module[{\[CapitalOmega], fluxMasked, expMasked, flux, exp,
offMask1D, offMask, valid, selector,
gdfSample, expSample, baseline, usedWeights,
relFluxGDF, offRelMean, ribbonPeak, allSkyMean},

\[CapitalOmega] = Pi/(3 tesselation^2);
fluxMasked = Flatten[data] /. Null -> Missing["NoData"];
expMasked  = Flatten[dataExposure] /. Null -> Missing["NoData"];

flux = fluxMasked;
exp  = expMasked;

(*Off-ribbon mask *)
offMask1D = calcOffMask[data, dataExposure];

(* Trim in case of tiny length mismatches *)
offMask =
If[Length[offMask1D] >= Length[flux],
Take[offMask1D, Length[flux]],
PadRight[offMask1D, Length[flux], False]
];

(* Valid numeric + exposed pixels *)
valid = MapThread[NumericQ[#1] && NumericQ[#2] && #2 > 0 &, {flux, exp}];

(* Combine masks correctly: elementwise AND, not Times *)
selector = MapThread[And, {offMask, valid}];

(* Exposure-weighted median baseline *)
gdfSample = Pick[flux, selector, True];
expSample = Pick[exp,  selector, True];

If[gdfSample === {} || expSample === {} || Total[expSample] == 0,
baseline   = Median@Select[flux, NumericQ];
usedWeights = False,
baseline   = weightedMedian[gdfSample, expSample];
usedWeights = True
];

(* Normalize, keeping Nulls in output *)
relFluxGDF = fluxMasked / baseline;
relFluxGDF = relFluxGDF /. Missing["NoData"] -> Null;

(* Diagnostics *)
allSkyMean = Mean@Select[flux, NumericQ];
offRelMean = Mean@Pick[relFluxGDF, selector, True];
ribbonPeak = Max@Select[relFluxGDF, NumericQ];

Print["Baseline method: ", If[usedWeights, "Exposure-weighted median", "Simple median"]];
Print["Baseline value = ", baseline];
Print["GDF / all-sky mean = ", baseline/allSkyMean];
Print["Off-ribbon samples used = ", Length[gdfSample]];
Print["Off-ribbon mean (should be ~1) = ", offRelMean];
Print["Ribbon peak = ", ribbonPeak];

relFluxGDF
]


calcOffMaskGrid[data_, lonGrid_, latGrid_] := 
 Module[{flux, lon, lat, toRad, vec, R, q85, ribbonCandIdx, Rr, cov, vals, vecs, nHat, delPsi, lonNose, latNose, noseVec, tailVec, noseAng, tailAng,
noseMask, tailMask, psiCut, offRibbonArc, qhi, fluxGate, gt, le, pos, isnum, masks, offMask},

flux = Flatten@data;         
lon  = Flatten@lonGrid;    
lat  = Flatten@latGrid; 

toRad[x_] := x Degree;
vec[la_, be_] := Module[{th = toRad[90 - be], ph = toRad[la]}, {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}];

(* Unit look vectors*)
R = MapThread[vec, {lon, lat}];

(*Bright ribbon pixels for great-circle fit *)
q85 = Quantile[Select[flux, NumericQ], .85];
ribbonCandIdx = Flatten@Position[flux, _?(NumericQ[#] && # >= q85 &)];
Rr = R[[ribbonCandIdx]];

(* Great-circle pole -> smallest-variance direction *)
cov  = N@Covariance[Rr];
{vals, vecs} = Eigensystem[cov];
nHat = Normalize@vecs[[-1]];
(* Angular distance from ribbon great circle *)
delPsi = N@(ArcSin[Abs[R . nHat]] / Degree);

(* Nose / tail exclusions *)
lonNose = 255.7; latNose = 5.1;
noseVec = vec[lonNose, latNose];
tailVec = -noseVec;
noseAng = N@(ArcCos[Clip[R . noseVec, {-1, 1}]] / Degree);
tailAng = N@(ArcCos[Clip[R . tailVec, {-1, 1}]] / Degree);
noseMask = noseAng > 30;
tailMask = tailAng > 30;
(* Off-ribbon arc *)
psiCut = 45.;(* degrees from fitted ribbon arc *)
offRibbonArc = delPsi > psiCut;
(*Robust flux gating*)
qhi = Quantile[Select[flux, NumericQ], .70];
fluxGate = flux <= qhi;

gt[v_, t_] := Map[TrueQ, Thread[v > t]];
le[v_, t_] := Map[TrueQ, Thread[v <= t]];
pos[v_]    := Map[TrueQ, Thread[v > 0]];
isnum[v_]  := Map[TrueQ @* NumericQ, v];

(* mask *)
offRibbonArc = gt[delPsi, psiCut];
noseMask     = gt[noseAng, 30];
tailMask     = gt[tailAng, 30];
validQ       = isnum[flux];
fluxGate     = le[flux, qhi];
masks = {offRibbonArc, noseMask, tailMask, validQ, fluxGate};
If[ !And @@ (VectorQ[#, BooleanQ] & /@ masks) || 
Not@SameQ @@ (Length /@ masks),
Print["Mask shape/type mismatch"]; Abort[]
];
offMask = MapThread[And, masks];
Partition[offMask, Length[latGrid[[1]]]]
]


calcOffMask[data_,dataExposure_]:=Module[{flux,coordinates,lon,lat,validQ,expoQ,toRad,vec,R,q85,lonNose,latNose,ribbonCandIdx,Rr,cov,vals,vecs,nHat,delPsi,noseVec,tailVec,noseAng,tailAng,noseMask,tailMask,psiCut,offRibbonArc,qhi,fluxGate,gt,le,pos,isnum,masks,offMask},
flux=Flatten@data; 
coordinates=getHPcoord[healpixringxyz];
lon=Flatten@coordinates[[;;,2]]; 
lat=Flatten@coordinates[[;;,1]];
validQ=NumericQ/@flux; 
expoQ=(Flatten@dataExposure)>0;

toRad[x_]:=x Degree;
vec[la_,be_]:=Module[{th=toRad[90-be],ph=toRad[la]},{Sin[th] Cos[ph],Sin[th] Sin[ph],Cos[th]}];

(*Unit look vectors*)
R=MapThread[vec,{lon,lat}];

(*bright candidates for the ribbon circle pole*)
q85=Quantile[Select[flux,NumericQ],.85];
ribbonCandIdx=Flatten@Position[flux,_?(NumericQ[#]&&#>=q85&)];
Rr=R[[ribbonCandIdx]];

cov=N@Covariance[Rr];
{vals,vecs}=Eigensystem[cov];
nHat=Normalize@vecs[[-1]];   (*eigenvector with smallest eigenvalue*)

(*Angular distance from great circle*)
delPsi=N@(ArcSin[Abs[R . nHat]]/Degree);

(*Optional nose/tail exclusions*)
lonNose=255.7;latNose=5.1;
noseVec=vec[lonNose,latNose];
tailVec=-noseVec;
noseAng=N@(ArcCos[Clip[R . noseVec,{-1,1}]]/Degree);
tailAng=N@(ArcCos[Clip[R . tailVec,{-1,1}]]/Degree);
noseMask=noseAng>30; (*exclude within 30deg of nose*)
tailMask=tailAng>30; (*exclude within 30deg of tail*)

(*Off-ribbon mask by distance from ribbon*arc**)
psiCut=45.;
offRibbonArc=delPsi>psiCut;

(*gating to suppress residual structures in "background"*)
qhi=Quantile[Select[flux,NumericQ],.80];
fluxGate=flux<=qhi;


gt[v_,t_]:=Map[TrueQ,Thread[v>t]];
le[v_,t_]:=Map[TrueQ,Thread[v<=t]];
pos[v_]:=Map[TrueQ,Thread[v>0]]; 
isnum[v_]:=Map[TrueQ@*NumericQ,v];


psiCut=45.;
offRibbonArc=gt[delPsi,psiCut];
noseMask=gt[noseAng,30.];
tailMask=gt[tailAng,30.];
validQ=isnum[flux];                               
expoQ=pos[dataExposure/. Null->0/. Missing[_]->0];

qhi=Quantile[Select[flux,NumericQ],.70];
fluxGate=le[flux,qhi]; (*suppress bright tails in GDF*)

(*Sanity:same length+Boolean*)
masks={offRibbonArc,noseMask,tailMask,validQ,expoQ,fluxGate};
If[!And@@(VectorQ[#,BooleanQ]&/@masks)||Not@SameQ@@(Length/@masks),Print["Mask shape/type mismatch:"," heads=",Head/@masks," lengths=",Length/@masks];
Abort[]];

offMask=MapThread[And,masks]
]


getHPcoord[healpixringxyz_]:=Module[{healpixringang,healpixringangCOLATI,healpixringangLONG,coordinates},
healpixringang=Delete[#,1]&/@Map[getSpherCoord[#]&,healpixringxyz];
healpixringangCOLATI = lati2colatiDeg[healpixringang[[All,2]]];
healpixringangLONG =healpixringang[[All,1]];
coordinates = Transpose[{healpixringangCOLATI,healpixringangLONG}]
]


weightedMedian[val_List,w_List]:=Module[{pairs,sorted,cum,total},pairs=Select[Transpose[{Flatten[val],Flatten[w]}],NumericQ[#[[1]]]&&NumericQ[#[[2]]]&&#[[2]]>0&];
If[pairs==={}||Total[pairs[[All,2]]]==0,Missing["NotAvailable"],sorted=SortBy[pairs,First];
cum=Accumulate[sorted[[All,2]]];
total=Last[cum];
FirstCase[Thread[{cum/total,sorted[[All,1]]}],{x_,y_}/;x>=0.5:>y,Median[val]]]]
(*From the list of {x, y} pairs, find the first pair where the first element x is \[GreaterEqual] 0.5.
Then, return the second element y from that pair*)
(*x - is the cumulative exposure fraction C_i/Ctot
y - is the flux F_i
So it finds the first flux value where the cumulative exposure reaches 50% of the total - i.e. the wighted meadian sESA*)


Options[makeGridFromBinnedDataWeighted]={"Normalize"->"relative",(*"relative" or "none"*)
"Op"->"mean",(*"mean" for flux/rates;"sum" for exposure*)
"Tag"->"ibex"              (*or "theseus" for filenames*)};

makeGridFromBinnedDataWeighted[colatitudes_,longitudes_,data_,energyStep_,outputDir_,resolution_,OptionsPattern[]]:=Module[{vals,lats,lons,valid,nRows,nCols,dtheta,dphi,thetaCenters,weights,weightAssoc,rowIdx,colIdx,pairs,rules,assoc,mapGrid,areaSum,fourPiCheck,sanity,dat,mapOut,op,norm,tag,sanityPost,mapWeighted},
op=OptionValue["Op"];(*"mean"|"sum"*)
norm=OptionValue["Normalize"];(*"relative"|"none"*)
tag=OptionValue["Tag"];(*Flatten consistently*)
lats=Flatten[colatitudes]//N;
lons=Flatten[longitudes]//N;
vals=Flatten[data]//N;
valid=MapThread[NumericQ[#1]&&NumericQ[#2]&&NumericQ[#3]&,{lats,lons,vals}];
lats=Pick[lats,valid,True];
lons=Pick[lons,valid,True];
vals=Pick[vals,valid,True];
nRows=resolution[[1]];
nCols=resolution[[2]];
dtheta=180./nRows;dphi=360./nCols;
thetaCenters=Range[dtheta/2,180-dtheta/2,dtheta];
(*cell solid-angle weights (for diagnostics/relative normalization)*)
weights=Table[Sin[thetaCenters[[r]] Degree]*dtheta Degree*dphi Degree,{r,nRows},{c,nCols}];
weightAssoc=Association@Table[{r,c}->weights[[r,c]],{r,nRows},{c,nCols}];
(*Assign samples to cells*)rowIdx=Clip[Floor[(lats-1)/(180/nRows)+1],{1,nRows}];
colIdx=Clip[Floor[(lons-1)/(360/nCols)+1],{1,nCols}];
(*Build {row,col}->list of values*)pairs=Transpose[{rowIdx,colIdx,vals}];
rules=MapThread[Rule,{pairs[[All,{1,2}]],pairs[[All,3]]}];
assoc=GroupBy[rules,First->Last];
(*Cell aggregation:mean for flux/rates,sum for exposure*)
mapGrid=Table[With[{cellVals=Lookup[assoc,Key[{r,c}],{}]},Which[Length[cellVals]==0,Null,op==="sum",Total[cellVals],True,Mean[cellVals]  (*default=mean*)]],{r,1,nRows},{c,1,nCols}];
(*Diagnostics (only meaningful for flux/rate maps)*)
areaSum=Total[Flatten[weights]];
fourPiCheck=areaSum/(4 Pi);
sanity=(Total[Flatten[((mapGrid/. Null->0)*weights)]]/areaSum);
Print["Check area sum (sr): ",N[areaSum]];
Print["4\[Pi] normalization check: ",N[fourPiCheck]];
Print["Sanity check (global weighted mean): ",N[sanity]];
(*Optional relative normalization (for Part 1)*)
mapOut=If[norm==="relative"&&op=!="sum",mapWeighted=mapGrid/sanity;sanityPost=(Total[Flatten[(mapWeighted/. Null->0)*weights]]/areaSum);
Print["Sanity check (global weighted mean, post-normalization): ",N[sanityPost]];mapWeighted,(*make mean=1 relative map for flux/rates*) mapGrid 
 (*keep physical scale (flux) or totals (exposure)*)];
(*Export*)
(*dat=ExportString[mapOut/. {(_?NumericQ) Null->Null},"Table"];
Export[ToFileName[outputDir]<>If[tag==="theseus","theseus","normalizedIB"]<>"_"<>If[norm==="relative","relFlux_","phys_"]<>ToString[energyStep]<>"_"<>ToString[resolution[[1]]]<>"_"<>ToString[resolution[[2]]]<>".txt",dat];*)
mapOut];


(*The area weighting accounts for unequal pixel solid angles at different latitudes (since HEALPix pixels are equal-area on the sphere, but your 2D rectangular grid cells are not)*)
areaWeightedMean[assoc_, resolution_] := 
Module[{nRows, nCols, dtheta, dphi, thetaCenters, weights, weightAssoc, map, areaSum, fourPiCheck, sanity},
  
{nRows, nCols} = resolution;
dtheta = 180./nRows;
dphi   = 360./nCols;
thetaCenters = Range[dtheta/2, 180 - dtheta/2, dtheta];
  
(* Solid-angle weights sin \[Theta] \[CapitalDelta]\[Theta] \[CapitalDelta]\[CurlyPhi], same per row *)
weights = Table[Sin[thetaCenters[[r]] Degree] * dtheta Degree * dphi Degree, {r, nRows}, {c, nCols}];
weightAssoc = Association@Table[{r, c} -> weights[[r, c]], {r, nRows}, {c, nCols}];
(* Weighted mean per bin *)
map = Table[
Module[{vals = Lookup[assoc, Key[{r, c}], {}],
ws = ConstantArray[weightAssoc[{r, c}], Length@Lookup[assoc, Key[{r, c}], {}]]},
If[Length[vals] == 0, Null, Total[vals*ws]/Total[ws]]],
{r, 1, nRows}, {c, 1, nCols}];
(* Diagnostic prints *)
areaSum = Total[Flatten[weights]];
fourPiCheck = areaSum/(4 Pi);
sanity = ((Total[Flatten[(map*weights)]]/. Null->0)/areaSum);
Print["Check area sum (sr): ", N[areaSum]];
Print["4\[Pi] normalization check: ", N[fourPiCheck]];
Print["Sanity check (global weighted mean): ", N[sanity]];
map
];


Options[changeGridHPtoLayout] = {"weighted" -> False};
changeGridHPtoLayout[data_, resolution_, OptionsPattern[]] := 
Module[{mapMean, coordinates, longitudes, colatitudes, val, nRows, 
nCols, colatBins, longBins, rowIdx, colIdx, pairs, rules, assoc, 
mapTotalParitioned, dat, weightedQ},

weightedQ = TrueQ[OptionValue["weighted"]];

coordinates = getHPcoord[healpixringxyz];  (* {colat, lon} in deg *)
longitudes   = coordinates[[All, 2]];
colatitudes  = coordinates[[All, 1]];
val = Flatten[data];

nRows = resolution[[1]];
nCols = resolution[[2]];
(*bin edges *)
colatBins = Subdivide[0, 180, nRows];
longBins  = Subdivide[0, 360, nCols];

rowIdx = Clip[Floor[(colatitudes - 1)/(180/nRows) + 1], {1, nRows}];
colIdx = Clip[Floor[(longitudes - 1)/(360/nCols) + 1], {1, nCols}];

pairs = Transpose[{rowIdx, colIdx, val}];
rules = MapThread[Rule, {pairs[[All, {1, 2}]], pairs[[All, 3]]}];
assoc = GroupBy[rules, First -> Last];

(*mode: area-weighted OR simple mean*)
mapMean = If[weightedQ, 
areaWeightedMean[assoc, resolution],
Table[ With[{vals = Select[Lookup[assoc, Key[{r, c}], {}] /. {Null -> Sequence[], Missing[_] -> Sequence[]}, NumericQ]},
If[vals === {}, Null, Mean[vals]]],
{r, 1, nRows}, {c, 1, nCols}]
];
mapTotalParitioned = mapMean/.{(_?NumericQ) Null -> Null, 
Plus[Null, a_?NumericQ] :> a, Plus[a_?NumericQ, Null] :> a};
mapTotalParitioned=mapTotalParitioned/.{(_?NumericQ) Null -> Null, 
Plus[Null, a_?NumericQ] :> a, Plus[a_?NumericQ, Null] :> a};
(*dat = ExportString[mapTotalParitioned, "Table"];
  Export[
   ToFileName[outputDir] <> 
    "normalizedHP2IB_" <> ToString[tesselation] <> "_" <> 
     ToString[energyStep] <> "_" <> ToString[resolution[[1]]] <> "_" <> 
     ToString[resolution[[2]]] <> ".txt", dat];*)
 mapTotalParitioned
];


uniformTessellationChange[data_,coordinates_,coordinatesNew_,tesselationNew_]:=Module[{coordinatesNewWithPixels,pixelsOriginalInNewGrid,dataNew,zeroList},
zeroList=Flatten[Position[data,Null]];
coordinatesNewWithPixels=Transpose[Append[{Range[Length[coordinatesNew]]},coordinatesNew]];
(*Original pixels in the new grid*)
pixelsOriginalInNewGrid=healpix`fun`ang2pixRing[tesselationNew,coordinates[[#,1]] Degree,coordinates[[#,2]]Degree]+1&/@Range[Length[coordinates]];
dataNew=ConstantArray[0,Length[coordinates]];
Module[{idx1=#,list},
	list=Flatten[Position[pixelsOriginalInNewGrid,idx1]];
	Module[{idx2=#}, If[Complement[list,zeroList]!={},
		dataNew[[list[[idx2]]]]=Mean[Module[{idx3=#},data[[idx3]]]&/@Complement[list,zeroList]],dataNew[[list[[idx2]]]]=Null];
		dataNew[[list[[idx2]]]]]&/@Range[Length[list]];
		]&/@pixelsOriginalInNewGrid;
	dataNew
]

partialTessellationChange[coordinatesChanged_,data_,coordinates_,coordinatesNew_,tesselationNew_]:=Module[{zeroList,coordinatesNewWithPixels,pixelsChanged,pixelsOriginalInNewGrid,dataNew},
zeroList=Flatten[Position[data,Null]];
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

uniformTessellationChangeHP[data_,outputDir_,tesselationNew_]:=Module[{hpPath,coordinates,dataPath,dataRat,hpPathNew,coordinatesNew,dataNew,dat},
hpPath = ToFileName[healpixDir,"testXYZ"];
init[hpPath,tesselation];
coordinates=getHPcoord[healpixringxyz];
hpPathNew = ToFileName[healpixDir,"testXYZ"];
init[hpPathNew,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyz];
dataNew=uniformTessellationChange[data,coordinates,coordinatesNew,tesselationNew];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_UniformTessellationChange"<>ToString[tesselationNew]<>"_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]

partialTessellationChangeHP[data_,outputDir_,tesselationNew_,point1_,point2_,marginLat1_,marginLong1_,marginLat2_,marginLong2_]:=Module[{hpPath,coordinates,coordinatesNew,coordinatesBoundaries,coordinatesChanged1,coordinatesChanged2,coordinatesChanged,dataPath,dataRat,zeroList,dataNew,dat},
hpPath = ToFileName[healpixDir,"testXYZ"];
init[hpPath,tesselation];
coordinates=getHPcoord[healpixringxyz];
init[hpPath,tesselationNew];
coordinatesNew=getHPcoord[healpixringxyz];
coordinatesBoundaries={{point1[[1]]-marginLat1,point1[[1]]+marginLat1},{point1[[2]]-marginLong1,point1[[2]]+marginLong1}};
coordinatesChanged1=Select[coordinates,#[[1]]>=coordinatesBoundaries[[1,1]]&&#[[1]]<=coordinatesBoundaries[[1,2]]&&#[[2]]>=coordinatesBoundaries[[2,1]]&&#[[2]]<=coordinatesBoundaries[[2,2]]&];
coordinatesBoundaries={{point2[[1]]-marginLat2,point2[[1]]+marginLat2},{point2[[2]]-marginLong2,point2[[2]]+marginLong2}};coordinatesChanged2=Select[coordinates,#[[1]]>=coordinatesBoundaries[[1,1]]&&#[[1]]<=coordinatesBoundaries[[1,2]]&&#[[2]]>=coordinatesBoundaries[[2,1]]&&#[[2]]<=coordinatesBoundaries[[2,2]]&];
coordinatesChanged=Sort[Join[coordinatesChanged1,coordinatesChanged2]];
(*data=data/. 0.->Null;*)
dataNew=partialTessellationChange[coordinatesChanged,data,coordinates,coordinatesNew,tesselationNew];
dat=ExportString[Transpose[Append[{Range[Length[dataNew]]},dataNew]],"Table"];
Export[ToFileName[outputDir]<>"data_PartialTessellationChange"<>ToString[tesselationNew]<>"_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
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

transformToRibbonCoordinates[data_,resolution_,mapType_,outputDir_]:=Module[{lr,br,coordinatesOriginal,healpixringxyzTransformed,coordinatesChanged,pixelNew,longitudes,colatitudes,nRows,nCols,colatBins,longBins,rowIdx,colIdx,mapHP,pairs,rules,assoc,map,mapTotal,mapTotalParitioned,mapHp,dat},
lr=221.;br=39.; (*Average Ribbon centre*)
(*Source: Variability in the Position of the IBEX Ribbon over Nine Years: More Observational
Evidence for a Secondary ENA Source - M. A. Dayeh*)
coordinatesOriginal=getHPcoord[healpixringxyz];
healpixringxyzTransformed=transformToRibbonCoord[healpixringxyz[[#]],lr,br]&/@Range[Length[healpixringxyz]];
coordinatesChanged=getHPcoord[healpixringxyzTransformed];
pixelNew=MapThread[healpix`fun`ang2pixRing[tesselation,#1,#2]&,{coordinatesChanged[[;;,1]]Degree,coordinatesChanged[[;;,2]]Degree}];
mapHP = data[[pixelNew]];
Which[
mapType=="HP",
dat=ExportString[Transpose[Append[{Range[Length[mapHP]]},mapHP]],"Table"];
Export[ToFileName[outputDir]<>"data_RibbonCentered_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];,

mapType=="layout",
longitudes=coordinatesChanged[[;;,2]];
colatitudes=coordinatesChanged[[;;,1]];

nRows=resolution[[1]];
nCols=resolution[[2]];
(*Compute bin edges*)
colatBins=Subdivide[0,180,nRows];
longBins=Subdivide[0,360,nCols];

rowIdx=Clip[Floor[(colatitudes-1)/(180/nRows)+1],{1,nRows}];
colIdx=Clip[Floor[(longitudes-1)/(360/nCols)+1],{1,nCols}];

(*Combine all into triplets*)
pairs=Transpose[{rowIdx,colIdx,mapHP}];

rules=MapThread[Rule,{pairs[[All,{1,2}]],pairs[[All,3]]}];
(*Group all pixel numbers belonging to the same cell*)
assoc=GroupBy[rules,First->Last];

(*Fill the 2D map with lists of pixel numbers*)
map=Table[Lookup[assoc,Key[{r,c}],{}],{r,1,nRows},{c,1,nCols}];

mapTotal=Total[#]&/@Flatten[map,1];
mapTotalParitioned=Partition[mapTotal,resolution[[2]]];
mapTotalParitioned=(mapTotalParitioned/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a})/. (Plus[a_?NumericQ,Null]:>a);
dat=ExportString[mapTotalParitioned,"Table"];
Export[ToFileName[outputDir]<>"data_RibbonCentered_t"<>ToString[tesselation]<>"_"<>ToString[energyStep]<>".txt",dat];
]
]


bin1to6[list_]:=Module[{i=#},list[[i;;i+5]]]&/@Range[1,Length[list],6]
bin2to6[list_]:=Module[{i=#},list[[i;;i+3]]]&/@Range[1,Length[list],4]
