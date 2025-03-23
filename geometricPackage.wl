(* ::Package:: *)

(*healpix`version = "healPix_2017-01-11"*)

healpix`fun`ang2pixRing[nSide_, \[Theta]_, \[Phi]_] := 
    Module[{piover2 = Pi/2, twopi = 2*Pi, z0 = 2./3., z, za, tt, tp, tmp, ir, 
      ip, kshift, nl2, nl4, ncap, npix, jp, jm, jn, ipix1, phi}, 
     z = Cos[\[Theta]]; za = Abs[z]; phi = Which[\[Phi] >= twopi, 
        \[Phi] - twopi, \[Phi] < 0, \[Phi] + twopi, True, \[Phi]]; 
      tt = \[Phi]/piover2; nl2 = 2*nSide; nl4 = 4*nSide; 
      ncap = nl2*(nSide - 1); npix = 12*nSide*nSide; 
      ipix1 = If[za <= z0, jp = Floor[nSide*(0.5 + tt - 0.75*z)]; 
         jm = Floor[nSide*(0.5 + tt + 0.75*z)]; ir = nSide + 1 + jp - jm; 
         kshift = If[EvenQ[ir], 1, 0]; 
         ip = IntegerPart[Floor[(jp + jm - nSide + kshift + 1)/2]] + 1; 
         ip = If[ip > nl4, ip - nl4, ip]; ncap + nl4*(ir - 1) + ip, 
        tp = tt - Floor[tt]; tmp = Sqrt[3.*(1. - za)]; 
         jp = Floor[nSide*tp*tmp]; jm = Floor[nSide*(1 - tp)*tmp]; 
         ir = jp + jm + 1; ip = Floor[tt*ir] + 1; 
         ip = If[ip > 4*ir, ip - 4*ir, ip]; ipix1 = 2*ir*(ir - 1) + ip; 
         If[z <= 0., npix - 2*ir*(ir + 1) + ip, ipix1]]; ipix1 - 1]
 
\[Phi] = 0.09023352232810683
 
healpix`fun`nPixels[nSide_] := 12*nSide^2

(*timGeomCoordLibr`notebookName = "timeGeomCoordLibr.01.02.nb"
timGeomCoordLibr`revision = "2021-10-29"
timGeomCoordLibrDate = "2021-10-29T13:13:41"
timGeomCoordLibrRevision = "2021-10-29"
timGeomCoordLibr`author = "M. Bzowski, bzowski at cbk.waw.pl"
timGeomCoordLibr`LGTM = "M. Kubiak, M. Strumik, N. Bukowiecka"*)
 
Rx[a_] := With[{ca = Cos[a], sa = Sin[a]}, {{1, 0, 0}, {0, ca, sa}, 
      {0, -sa, ca}}]
 
Ry[b_] := With[{cb = Cos[b], sb = Sin[b]}, {{cb, 0, -sb}, {0, 1, 0}, 
      {sb, 0, cb}}]
 
Rz[g_] := With[{cg = Cos[g], sg = Sin[g]}, {{cg, sg, 0}, {-sg, cg, 0}, 
      {0, 0, 1}}]
 
lati2colatiRad[\[Phi]_] := Pi/2 - \[Phi]
 
lati2colatiDeg[\[Phi]_] := 90 - \[Phi]
 
colati2latiRad[\[Phi]_] := -\[Phi] + Pi/2
 
colati2latiDeg[\[Phi]_] := -\[Phi] + 90
 
makeVec[l_, b_] := {Cos[l]*Cos[b], Sin[l]*Cos[b], Sin[b]}
 
makeCoVec[l_, b_] := {Cos[l]*Sin[b], Sin[l]*Sin[b], Cos[b]}
 
getSpherCoordRad[vec_] := Module[{r, l, b}, 
     {r, l, b} = getSpherCoord\[Pi][vec]; If[l < 0, l = l + 2*Pi, l]; 
      {r, l, b}]
 
getSpherCoord\[Pi][vec_] := Module[{r, l, b}, r = Sqrt[vec . vec]; 
      If[r == 0. || r == 0, Return[{0, 0, 0}]]; b = ArcSin[vec[[3]]/r]; 
      l = If[(vec[[1]] == 0. || vec[[1]] == 0) && (vec[[2]] == 0 || 
          vec[[2]] == 0.), 0, ArcTan[vec[[1]], vec[[2]]]]; {r, l, b}]
 
getSpherCoord[vec_] := With[{qq = getSpherCoordRad[vec]}, 
     {qq[[1]], qq[[2]]/Degree, qq[[3]]/Degree}]
 
getSpherCoord180[vec_] := Module[{tmp}, tmp = getSpherCoord\[Pi][vec]; 
      tmp/{1, Degree, Degree}]
 
getSpherCoCoord[vec_] := Module[{r, l, b}, {r, l, b} = getSpherCoord180[vec]; 
      {r, l, lati2colatiDeg[b]}]
 
angleCoord[l1_, b1_, l2_, b2_] := Module[{vec1, vec2}, 
     vec1 = makeVec[l1, b1]; vec2 = makeVec[l2, b2]; angleVec[vec1, vec2]]
 
angleVec[vec1_, vec2_] := ArcCos[vec1 . vec2/Sqrt[vec1 . vec1*vec2 . vec2]]
 
angleCoordDegree[l1_, b1_, l2_, b2_] := Module[{vec1, vec2}, 
     angleCoord[l1*Degree, b1*Degree, l2*Degree, b2*Degree]]
 
vecLength[vec_] := Sqrt[vec . vec]
 
vtrFun[lB_, bB_, th_, fi_] := With[{cl = Cos[lB], sl = Sin[lB], cb = Cos[bB], 
      sb = Sin[bB], cth = Cos[th], sth = Sin[th], cfi = Cos[fi], 
      sfi = Sin[fi]}, Chop[{cb*cl*cth + cfi*cl*sb*sth - sfi*sl*sth, 
       cb*cth*sl + cl*sfi*sth + cfi*sb*sl*sth, cth*sb - cb*cfi*sth}]]
 
makeMatr2vec[r1_, r2_] := Module[{r3, xd, yd, zd, M}, 
     r3 = Cross[r1, r2]; xd = r1/vecLength[r1]; 
      yd = With[{q = Cross[r3, r1]}, q/vecLength[q]]; zd = r3/vecLength[r3]; 
      M = {xd, yd, zd}]
 
makeMatr2spherDeg[r1`ang_, r2`ang_] := 
    With[{r1 = makeVec[r1`ang[[1]]*Degree, r1`ang[[2]]*Degree], 
      r2 = makeVec[r2`ang[[1]]*Degree, r2`ang[[2]]*Degree]}, 
     makeMatr2vec[r1, r2]]
 
makeMatrPolVecspherDeg[r1`ang_, r2`ang_] := 
    With[{r1 = makeVec[r1`ang[[1]]*Degree, r1`ang[[2]]*Degree], 
      r2 = makeVec[r2`ang[[1]]*Degree, r2`ang[[2]]*Degree]}, 
     makeMatrPolVec[r1, r2]]
 
makeMatrPolVec[rz_, rzx_] := Module[{xd, yd, zd}, 
     zd = rz; yd = Cross[rz, rzx]; xd = Cross[yd, zd]; 
      {xd/vecLength[xd], yd/vecLength[yd], zd/vecLength[zd]}]
 
positionAngle[\[Alpha]1_, \[Delta]1_, \[Alpha]2_, \[Delta]2_] := 
    ArcTan[Cos[\[Delta]2]*Tan[\[Delta]1] - Sin[\[Delta]2]*
       Cos[\[Alpha]1 - \[Alpha]2], Sin[\[Alpha]1 - \[Alpha]2]]
