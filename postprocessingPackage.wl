(* ::Package:: *)

exposureWeightedGlobalFlux[fluxMap_, expMap_] :=
 Module[{flux, exp, valid},
  flux = Flatten[fluxMap];
  exp = Flatten[expMap];

  valid =
   MapThread[
    NumericQ[#1] && NumericQ[#2] && #2 > 0 &,
    {flux, exp}
   ];

  Total[Pick[flux*exp, valid, True]]/
   Total[Pick[exp, valid, True]]
]


circularMeanDeg[x_List]:=Module[{ang,s,c,a},ang=N[x] Degree;
s=Mean[Sin[ang]];
c=Mean[Cos[ang]];
a=Mod[ArcTan[c,s]/Degree,360];
a]


gdfNormalizeGrid[mapFluxGrid_,mapExpGrid_,offMaskGrid_,label_:""]:=Module[{fluxFlat,expFlat,maskFlat,valid,selector,offFlux,offExp,baseline,mapGDF,diagMedian,diagMean,diagWeightedMean},fluxFlat=Flatten[mapFluxGrid];
expFlat=Flatten[mapExpGrid];
maskFlat=Flatten[offMaskGrid];
valid=MapThread[NumericQ[#1]&&NumericQ[#2]&&#2>0&,{fluxFlat,expFlat}];
selector=MapThread[And,{maskFlat,valid}];
offFlux=Pick[fluxFlat,selector,True];
offExp=Pick[expFlat,selector,True];
If[offFlux==={}||offExp==={}||Total[offExp]==0,Print[label," GDF baseline failed: no valid off-ribbon exposed samples."];
Return[$Failed];];
baseline=weightedMedian[offFlux,offExp];
mapGDF=mapFluxGrid/baseline;
mapGDF=mapGDF/. {(_?NumericQ) Null->Null,Plus[Null,a_?NumericQ]:>a,Plus[a_?NumericQ,Null]:>a,ComplexInfinity->Null,Indeterminate->Null};
diagMedian=weightedMedian[offFlux/baseline,offExp];
diagMean=Mean[offFlux/baseline];
diagWeightedMean=Total[(offFlux/baseline)*offExp]/Total[offExp];
Print[label," GDF baseline = ",N[baseline]];
Print[label," off-ribbon exposure-weighted median after GDF = ",N[diagMedian]];
Print["Means don't need to be 1. GDF is median normalized"];
Print[label," off-ribbon unweighted mean after GDF = ",N[diagMean]];
Print[label," off-ribbon exposure-weighted mean after GDF (should be 1) = ",N[diagWeightedMean]];
Print[label," off-ribbon samples used = ",Length[offFlux]];
mapGDF]


Options[normalizeGDF] = Join[
  Options[calcOffMask],
  {
   "PrintDiagnostics" -> True
  }
];

normalizeGDF[data_, dataExposure_, tesselation_, energyStep_, OptionsPattern[]] :=
Module[{fluxMasked, expMasked, flux, exp,
offMask1D, offMask, valid, selector,
gdfSample, expSample, baseline, usedWeights,
relFluxGDF, offRelMean, ribbonPeak, allSkyMean,
offRel, offExp, offRelMedian, offRelWeightedMean},

fluxMasked = Flatten[data] /. Null -> Missing["NoData"];
expMasked  = Flatten[dataExposure] /. Null -> Missing["NoData"];

flux = fluxMasked;
exp  = expMasked;

(*Off-ribbon mask *)
offMask1D =
   calcOffMask[
    data,
    dataExposure,
    "PsiCut" -> OptionValue["PsiCut"],
    "NoseTailCut" -> OptionValue["NoseTailCut"],
    "FluxPercentile" -> OptionValue["FluxPercentile"],
    "RibbonFitPercentile" -> OptionValue["RibbonFitPercentile"],
    "NoseLon" -> OptionValue["NoseLon"],
    "NoseLat" -> OptionValue["NoseLat"]
   ];

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

allSkyMean = Mean@Select[flux, NumericQ];
ribbonPeak = Max@Select[relFluxGDF, NumericQ];

(* Diagnostics *)
offRel = Pick[relFluxGDF, selector, True];
  offExp = Pick[exp, selector, True];

  offRelMedian =
   If[usedWeights,
    weightedMedian[offRel, offExp],
    Median[offRel]
   ];

  offRelMean = Mean[offRel];

  offRelWeightedMean =
   If[usedWeights,
    Total[offRel*offExp]/Total[offExp],
    Missing["NoWeights"]
   ];

  If[TrueQ[OptionValue["PrintDiagnostics"]],
   Print["Baseline method: ",
    If[usedWeights, "Exposure-weighted median", "Simple median"]];
   Print["Baseline value = ", N[baseline]];
   Print[
    "Off-ribbon exposure-weighted median after GDF (should be 1) = ",
    N[offRelMedian]
   ];
   Print["Means don't need to be 1. GDF is median normalized"];
   Print["Off-ribbon unweighted mean after GDF = ", N[offRelMean]];
   Print[
    "Off-ribbon exposure-weighted mean after GDF = ",
    N[offRelWeightedMean]
   ];
   Print["GDF / all-sky mean = ", N[baseline/allSkyMean]];
   Print["Off-ribbon samples used = ", Length[gdfSample]];
   Print["Ribbon peak = ", N[ribbonPeak]];
  ];

  relFluxGDF
]


Options[calcOffMaskGrid] = {
  "PsiCut" -> 45.,
  "NoseTailCut" -> 30.,
  "FluxPercentile" -> .70,
  "RibbonFitPercentile" -> .85,
  "NoseLon" -> 255.7,
  "NoseLat" -> 5.1
};

calcOffMaskGrid[data_, lonGrid_, latGrid_, OptionsPattern[]] :=
 Module[
  {
   flux, lon, lat,
   psiCut, noseTailCut, fluxPercentile, ribbonFitPercentile,
   lonNose, latNose,
   toRad, vec, R,
   validFlux, validForFit,
   qRibbon, ribbonCandIdx, Rr, cov, vals, vecs, nHat,
   delPsi, noseVec, tailVec, noseAng, tailAng,
   offRibbonArc, noseMask, tailMask,
   qFlux, fluxGate, masks, offMaskFlat,
   nRows, nCols
  },

  flux = Flatten[data];
  lon = Flatten[lonGrid];
  lat = Flatten[latGrid];

  nRows = Length[latGrid];
  nCols = Length[latGrid[[1]]];

  psiCut = OptionValue["PsiCut"];
  noseTailCut = OptionValue["NoseTailCut"];
  fluxPercentile = OptionValue["FluxPercentile"];
  ribbonFitPercentile = OptionValue["RibbonFitPercentile"];
  lonNose = OptionValue["NoseLon"];
  latNose = OptionValue["NoseLat"];

  toRad[x_] := x Degree;

  vec[la_, be_] :=
   Module[{th = toRad[90 - be], ph = toRad[la]},
    {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}
   ];

  If[! SameQ[Length /@ {flux, lon, lat}],
   Print["Length mismatch in calcOffMaskGrid: ",
    Length /@ {flux, lon, lat}];
   Abort[];
  ];
(* Unit look vectors*)
  R = MapThread[vec, {lon, lat}];

  validFlux = Map[NumericQ, flux];
  validForFit = Pick[flux, validFlux, True];

  If[validForFit === {},
   Print["calcOffMaskGrid failed: no numeric flux values."];
   Abort[];
  ];

  (* Fit ribbon great circle using bright valid pixels. *)
  qRibbon = Quantile[validForFit, ribbonFitPercentile];

  ribbonCandIdx =
   Flatten@Position[
     flux,
     _?(NumericQ[#] && # >= qRibbon &)
   ];

  If[Length[ribbonCandIdx] < 3,
   Print["calcOffMaskGrid failed: too few ribbon-candidate pixels."];
   Abort[];
  ];

  Rr = R[[ribbonCandIdx]];
(* Great-circle pole -> smallest-variance direction *)
  cov = N@Covariance[Rr];
  {vals, vecs} = Eigensystem[cov];
  nHat = Normalize[vecs[[-1]]];
(* Angular distance from ribbon great circle *)
  delPsi = N@(ArcSin[Abs[R . nHat]]/Degree);

  noseVec = vec[lonNose, latNose];
  tailVec = -noseVec;

  noseAng = N@(ArcCos[Clip[R . noseVec, {-1, 1}]]/Degree);
  tailAng = N@(ArcCos[Clip[R . tailVec, {-1, 1}]]/Degree);

  offRibbonArc = Map[TrueQ[# > psiCut] &, delPsi];
  noseMask = Map[TrueQ[# > noseTailCut] &, noseAng];
  tailMask = Map[TrueQ[# > noseTailCut] &, tailAng];

  qFlux = Quantile[validForFit, fluxPercentile];

  fluxGate =
   Map[
    If[NumericQ[#], TrueQ[# <= qFlux], False] &,
    flux
   ];

  masks = {
    offRibbonArc,
    noseMask,
    tailMask,
    validFlux,
    fluxGate
  };

  If[
   ! And @@ (VectorQ[#, BooleanQ] & /@ masks) ||
    Not@SameQ @@ (Length /@ masks),
   Print[
    "Mask shape/type mismatch in calcOffMaskGrid: heads = ",
    Head /@ masks,
    " lengths = ",
    Length /@ masks
   ];
   Abort[];
  ];

  offMaskFlat = MapThread[And, masks];

  Partition[offMaskFlat, nCols]
]


Options[calcOffMask] = {
  "PsiCut" -> 45.,
  "NoseTailCut" -> 30.,
  "FluxPercentile" -> .70,
  "RibbonFitPercentile" -> .85,
  "NoseLon" -> 255.7,
  "NoseLat" -> 5.1
};

calcOffMask[data_, dataExposure_, OptionsPattern[]] :=
 Module[
  {
   flux, exp, coordinates, lon, lat,
   psiCut, noseTailCut, fluxPercentile, ribbonFitPercentile,
   lonNose, latNose,
   toRad, vec, R,
   validFlux, validExp, validForFit,
   qRibbon, ribbonCandIdx, Rr, cov, vals, vecs, nHat,
   delPsi, noseVec, tailVec, noseAng, tailAng,
   offRibbonArc, noseMask, tailMask,
   qFlux, fluxGate, masks, offMask
  },

  flux = Flatten[data];
  exp = Flatten[dataExposure] /. Null -> 0 /. Missing[_] -> 0;

  coordinates = getHPcoord[healpixringxyz];
  lon = Flatten[coordinates[[All, 2]]];
  lat = Flatten[coordinates[[All, 1]]];

  psiCut = OptionValue["PsiCut"];
  noseTailCut = OptionValue["NoseTailCut"];
  fluxPercentile = OptionValue["FluxPercentile"];
  ribbonFitPercentile = OptionValue["RibbonFitPercentile"];
  lonNose = OptionValue["NoseLon"];
  latNose = OptionValue["NoseLat"];

  toRad[x_] := x Degree;

  vec[la_, be_] :=
   Module[{th = toRad[90 - be], ph = toRad[la]},
    {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}
   ];


(*Unit look vectors*)
  R = MapThread[vec, {lon, lat}];

  validFlux = Map[NumericQ, flux];
  validExp = Map[TrueQ[# > 0] &, exp];

  If[! SameQ[Length /@ {flux, exp, lon, lat, R}],
   Print["Length mismatch in calcOffMask: ",
    Length /@ {flux, exp, lon, lat, R}];
   Abort[];
  ];

  (* Fit ribbon great circle using bright valid pixels. *)
  validForFit = Pick[flux, validFlux, True];

  If[validForFit === {},
   Print["calcOffMask failed: no numeric flux values."];
   Abort[];
  ];

(*bright candidates for the ribbon circle pole*)
  qRibbon = Quantile[validForFit, ribbonFitPercentile];

  ribbonCandIdx =
   Flatten@Position[
     flux,
     _?(NumericQ[#] && # >= qRibbon &)
   ];

  If[Length[ribbonCandIdx] < 3,
   Print["calcOffMask failed: too few ribbon-candidate pixels."];
   Abort[];
  ];

  Rr = R[[ribbonCandIdx]];

  cov = N@Covariance[Rr];
  {vals, vecs} = Eigensystem[cov];
  nHat = Normalize[vecs[[-1]]]; (*eigenvector with smallest eigenvalue*)

(*Angular distance from great circle*)
  delPsi = N@(ArcSin[Abs[R . nHat]]/Degree); 
(*Optional nose/tail exclusions*)
  noseVec = vec[lonNose, latNose];
  tailVec = -noseVec;

  noseAng = N@(ArcCos[Clip[R . noseVec, {-1, 1}]]/Degree);
  tailAng = N@(ArcCos[Clip[R . tailVec, {-1, 1}]]/Degree);
(*Off-ribbon mask by distance from ribbon*arc**)
  offRibbonArc = Map[TrueQ[# > psiCut] &, delPsi];
  noseMask = Map[TrueQ[# > noseTailCut] &, noseAng]; (*exclude within  x deg of nose*)
  tailMask = Map[TrueQ[# > noseTailCut] &, tailAng];
 (*suppress bright tails in GDF*)
  qFlux = Quantile[validForFit, fluxPercentile];
  fluxGate =
   Map[
    If[NumericQ[#], TrueQ[# <= qFlux], False] &,
    flux
   ];

  masks = {
    offRibbonArc,
    noseMask,
    tailMask,
    validFlux,
    validExp,
    fluxGate
  };

  If[
   ! And @@ (VectorQ[#, BooleanQ] & /@ masks) ||
    Not@SameQ @@ (Length /@ masks),
   Print[
    "Mask shape/type mismatch in calcOffMask: heads = ",
    Head /@ masks,
    " lengths = ",
    Length /@ masks
   ];
   Abort[];
  ];

  offMask = MapThread[And, masks];

  offMask
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


Options[makeGridFromBinnedDataWeighted] = {
   (*"Normalize" -> "relative",   (* "relative" or "none" *)*)
   "Op" -> "mean",              (* "mean" for flux/rates; "sum" for exposure, "exposureWeightedMean" for flux reconstructed from counts/exp *)
   "Tag" -> "ibex",
   "SourceDTheta" -> 6.,        (* source-bin size in colat, deg *)
   "SourceDPhi" -> 6.,          (* source-bin size in lon, deg *)
   "PrintDiagnostics" -> True
};

binSolidAngle[thetaCenter_, dtheta_, dphi_] :=
 Module[{th1, th2},
  th1 = Max[0., thetaCenter - dtheta/2] Degree;
  th2 = Min[180., thetaCenter + dtheta/2] Degree;
  dphi Degree * (Cos[th1] - Cos[th2])
 ];

makeGridFromBinnedDataWeighted[
   colatitudes_, longitudes_, data_, energyStep_, outputDir_, resolution_,
   OptionsPattern[]
] :=
 Module[
  {
   validMaskIn, observedArea, vals, lats, lons, valid,
   nRows, nCols, dthetaT, dphiT, thetaCentersT,
   targetWeights, rowIdx, colIdx,
   op, norm, tag, dthetaS, dphiS,
   srcWeights, records, assoc,
   mapGrid, assignedAreaGrid, validMask,
   nativeMean, remapMeanSource, remapMeanPainted,
   mapOut, sanityPost, areaSum, fourPiCheck, den1, den2
  },

  op = OptionValue["Op"];
  (*norm = OptionValue["Normalize"];*)
  tag = OptionValue["Tag"];
  dthetaS = OptionValue["SourceDTheta"];
  dphiS = OptionValue["SourceDPhi"];

  (* Flatten inputs *)
  lats = Flatten[colatitudes] // N;
  lons = Flatten[longitudes] // N;
  vals = Flatten[data] // N;

  valid = MapThread[
    NumericQ[#1] && NumericQ[#2] && NumericQ[#3] &,
    {lats, lons, vals}
  ];

  lats = Pick[lats, valid, True];
  lons = Pick[lons, valid, True];
  vals = Pick[vals, valid, True];

  nRows = resolution[[1]];
  nCols = resolution[[2]];

  (* target-grid geometry *)
  dthetaT = 180./nRows;
  dphiT = 360./nCols;
  thetaCentersT = Range[dthetaT/2, 180 - dthetaT/2, dthetaT];

  targetWeights = Table[
    Sin[thetaCentersT[[r]] Degree] * dthetaT Degree * dphiT Degree,
    {r, nRows}, {c, nCols}
  ];

  (* source-bin solid angles *)
  srcWeights = binSolidAngle[#, dthetaS, dphiS] & /@ lats;

  (* assign source bins to target cells *)
  rowIdx = Clip[Floor[(lats - 1)/(180./nRows) + 1], {1, nRows}];
  colIdx = Clip[Floor[(lons - 1)/(360./nCols) + 1], {1, nCols}];

  (* records = {row, col, value, sourceArea} *)
  records = Transpose[{rowIdx, colIdx, vals, srcWeights}];

  (* group records by target cell *)
  assoc = GroupBy[records, #[[1 ;; 2]] &];

  (* aggregate per cell *)
  mapGrid = Table[
    With[{cell = Lookup[assoc, Key[{r, c}], {}]},
      Which[
        cell === {}, Null,
        op === "sum",
          Total[cell[[All, 3]]],
        True,
          Total[cell[[All, 3]]*cell[[All, 4]]] / Total[cell[[All, 4]]]
      ]
    ],
    {r, 1, nRows}, {c, 1, nCols}
  ];

  (* total source area assigned into each target cell *)
  assignedAreaGrid = Table[
    With[{cell = Lookup[assoc, Key[{r, c}], {}]},
      If[cell === {}, 0., Total[cell[[All, 4]]]]
    ],
    {r, 1, nRows}, {c, 1, nCols}
  ];

  areaSum = Total[Flatten[targetWeights]];
  fourPiCheck = areaSum/(4 Pi);

  validMask = Map[NumericQ, mapGrid, {2}];
  observedArea = Total[
    Pick[Flatten[targetWeights], Flatten[validMask], True]
  ];

  (* source-conserved native mean *)
  nativeMean = Total[vals*srcWeights]/Total[srcWeights];

  (* source-conserved remap mean: should match nativeMean *)
  den1 = Total[Flatten[assignedAreaGrid]];
  remapMeanSource =
    If[den1 == 0, Missing["NoData"],
      Total[Flatten[Replace[mapGrid, Null -> 0, {2}] * assignedAreaGrid]] / den1
    ];

  (* painted-target mean: diagnostic only *)
  den2 = observedArea;
  remapMeanPainted =
    If[den2 == 0, Missing["NoPaintedArea"],
      Total[
        Pick[
          Flatten[Replace[mapGrid, Null -> 0, {2}] * targetWeights],
          Flatten[validMask],
          True
        ]
      ] / den2
    ];

  If[TrueQ[OptionValue["PrintDiagnostics"]],
    Print["Check area sum target (sr): ", N[areaSum]];
    Print["4\[Pi] normalization check: ", N[fourPiCheck]];
    Print["Observed area (target cells with data) (sr): ", N[observedArea]];
    Print["Native mean (source-area weighted): ", N[nativeMean]];
    Print["Remap mean (source-conserved):      ", N[remapMeanSource]];
    Print["Remap mean (painted target):        ", N[remapMeanPainted]];
  ];

  (*mapOut =
    If[norm === "relative" && op =!= "sum",
      Module[{mapWeighted},
        mapWeighted = mapGrid / remapMeanSource;
        sanityPost =
          Total[
            Flatten[Replace[mapWeighted, Null -> 0, {2}] * assignedAreaGrid]
          ] / den1;
        Print["Post-normalization mean (source-conserved): ", N[sanityPost]];
        mapWeighted
      ],
      mapGrid
    ];*)

  mapGrid
]


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


Options[changeGridHPtoLayout] = {
  "Op" -> "mean",              (* "mean" for flux/rate, "sum" for exposure/counts *)
  "PrintDiagnostics" -> True
};

changeGridHPtoLayout[data_, resolution_, OptionsPattern[]] := 
 Module[
  {
   coordinates, longitudes, colatitudes, val, nRows, nCols,
   rowIdx, colIdx, triplets, validTriplets, rules,
   assocVals, assocCnt, mapMean, countGrid,
   dtheta, dphi, thetaCenters, targetWeights, validMask,
   nativeMean, remapMeanSource, remapMeanPainted,
   den1, den2,op
  },

  coordinates = getHPcoord[healpixringxyz];
  longitudes = coordinates[[All, 2]];
  colatitudes = coordinates[[All, 1]];
  val = Flatten[data];

  nRows = resolution[[1]];
  nCols = resolution[[2]];

  rowIdx = Clip[Floor[(colatitudes)/(180./nRows) + 1], {1, nRows}];
  colIdx = Clip[Floor[(longitudes)/(360./nCols) + 1], {1, nCols}];

  triplets = Transpose[{rowIdx, colIdx, val}];
  validTriplets = Select[triplets, NumericQ[#[[3]]] &];

  rules = ({#[[1]], #[[2]]} -> #[[3]]) & /@ validTriplets;

  assocVals = GroupBy[rules, First -> Last];
  assocCnt = Counts[First /@ rules];

  (*mapMean = Table[
    With[{vals = Lookup[assocVals, Key[{r, c}], {}]},
      If[vals === {}, Null, Mean[vals]]
    ],
    {r, 1, nRows}, {c, 1, nCols}
  ];*)
  op = OptionValue["Op"];

	mapMean = Table[
	  With[{vals = Lookup[assocVals, Key[{r, c}], {}]},
	    Which[
	      vals === {}, Null,
	      op === "sum", Total[vals],
	      True, Mean[vals]
	    ]
	  ],
	  {r, 1, nRows}, {c, 1, nCols}
	];

  countGrid = Table[
    Lookup[assocCnt, Key[{r, c}], 0],
    {r, 1, nRows}, {c, 1, nCols}
  ];

  dtheta = 180./nRows;
  dphi = 360./nCols;
  thetaCenters = Range[dtheta/2, 180 - dtheta/2, dtheta];

  targetWeights = Table[
    Sin[thetaCenters[[r]] Degree] * dtheta Degree * dphi Degree,
    {r, nRows}, {c, nCols}
  ];

  validMask = Map[NumericQ, mapMean, {2}];

  nativeMean = Mean[Cases[val, _?NumericQ]];

  den1 = Total[Flatten[countGrid]];
  remapMeanSource =
    If[den1 == 0, Missing["NoData"],
      Total[Flatten[Replace[mapMean, Null -> 0, {2}] * countGrid]]/den1
    ];

  den2 = Total[
    Pick[Flatten[targetWeights], Flatten[validMask], True]
  ];
  remapMeanPainted =
    If[den2 == 0, Missing["NoPaintedArea"],
      Total[
        Pick[
          Flatten[Replace[mapMean, Null -> 0, {2}] * targetWeights],
          Flatten[validMask],
          True
        ]
      ]/den2
    ];

  If[TrueQ[OptionValue["PrintDiagnostics"]],
    Print["HP native mean (equal-area source): ", N[nativeMean]];
    Print["HP remap mean (source-conserved):   ", N[remapMeanSource]];
    Print["HP remap mean (painted target):     ", N[remapMeanPainted]];
  ];
Print["HP target cells total: ", nRows*nCols];
Print["HP target cells with >=1 HP pixel center: ", Count[Flatten[countGrid], _?(# > 0 &)]];
Print["HP target cells with numeric output: ", Count[Flatten[mapMean], _?NumericQ]];
Print["HP empty target cells: ", Count[Flatten[countGrid], 0]];
Print["HP Null target cells: ", Count[Flatten[mapMean], Null]];
Print["HP min/max HP pixels per painted cell: ",
  MinMax[Select[Flatten[countGrid], # > 0 &]]
];
  mapMean
]


areaWeightedObservedMeanGrid[map_, resolution_] :=
 Module[
  {nRows, nCols, dtheta, dphi, thetaCenters, weights,
   validMask, observedArea, mean},

  {nRows, nCols} = resolution;

  dtheta = 180./nRows;
  dphi = 360./nCols;

  thetaCenters = Range[dtheta/2, 180 - dtheta/2, dtheta];

  weights = Table[
    Sin[thetaCenters[[r]] Degree] * dtheta Degree * dphi Degree,
    {r, nRows}, {c, nCols}
  ];

  validMask = Map[NumericQ, map, {2}];

  observedArea =
    Total[Pick[Flatten[weights], Flatten[validMask], True]];

  mean =
    Total[
      Pick[
        Flatten[Replace[map, Null -> 0, {2}] * weights],
        Flatten[validMask],
        True
      ]
    ]/observedArea;

  mean
]


relativeNormalizeGrid[map_, resolution_, label_:""] :=
 Module[{mean},
  mean = areaWeightedObservedMeanGrid[map, resolution];

  If[label =!= "",
    Print[label, " relative denominator: ", N[mean]]
  ];

  map/mean
]


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
bin2to6[list_]:=Module[{i=#},list[[i;;i+2]]]&/@Range[1,Length[list],3]
