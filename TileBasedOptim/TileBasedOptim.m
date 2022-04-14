(* fibo-hilbert.m
   V.O. version 2002/12/14
*)
 
(****************** params *******************)
SetOptions[Graphics, ImageSize -> { 1024, Automatic}];

mf := MatrixForm
T :=  Transpose
PI = Pi//N;
known := ValueQ
i2s[n_,len_:6] := ToString[NumberForm[n, len-1, NumberPadding -> "0"]]
r2s[n_,len_:6,dec_:5] := ToString[NumberForm[n, {len,dec}]]
tab2s[tab_] := StringJoin @ Table[" "<>ToString[tab[[i]]]<>" ",{i,Length[tab]}]
tab2snosep[tab_] := StringJoin @ Table[ ToString[tab[[i]]] ,{i,Length[tab]}]
tab2ssep[tab_] := StringJoin @ Table["_"<>ToString[tab[[i]]],{i,Length[tab]}]
tab2ssepComma[tab_] := StringJoin[ToString[tab[[1]]], Table[","<>ToString[tab[[i]]],{i,Min[2,Length[tab]], Length[tab]}] ]
digits2str[digits_] := StringJoin[ToString /@ digits] 
str2n[str_] := FromDigits[#, 2] & @ Table[StringTake[str, {i}] // ToExpression, {i, StringLength[str]}]

(*------------------------- new code 2022-04 vo, based on fibo-hilbert.m (version 2002/12/14) -------------------------*)
(*------------------------- constants -------------------------*)
phi = N[GoldenRatio,32];
oneoverphi = 1/phi;
oneoverphi2 = oneoverphi^2;

mxRot0 =   {{1,0}, {0,1}};
mxRot90 =  {{0, -1}, {1, 0}};
mxRot180 = {{-1,0}, {0,-1}};
mxRot270 = {{0, 1}, {-1, 0}};

type1 = 1;
type2 = 2;
type3 = 3;
type4 = 4;
type5 = 5;
type6 = 6;

(* debugging flags *)
showTileType = 1;
showTilefcode = 2;
showSFC = 4;
showValue = 8;

(*------------------------- end of constants -------------------------*)
FIBOF[symbols_] := With[ {s = Reverse[symbols]}, Total@Table[Fibonacci[i+1 ] s[[i]], {i, Length[s]}] ]
FIBOFinv[symbols_] := Total@Table[Fibonacci[i+1 ] symbols[[i]], {i, Length[symbols]}]
FIBOFxy[symbols_] := FIBOF/@symbols

FIBOphitab = Table[phi^-i, {i, 32}] // N;
FIBOPhi[s_] := Sum[FIBOphitab[[i]] s[[i]], {i, Length[s]}] 


demoFiboSFC[niters_:15] :=
    Module[ {},
    	dbg = False;
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
		Graphics[ getFiboSFCTilesGL[tlst] ]//Print;
		Do[
			tlst = subdivFiboSFCTiles @ tlst;
			Graphics[ getFiboSFCTilesGL[tlst], PlotLabel-> iter ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,niters}];
		Graphics[ getFiboSFCTilesGL[tlst,showSFC] ]//Print;
		Graphics[ getFiboSFCTilesGL[tlst,showValue] ]//Print;
]


subdivFiboSFCTiles[tlst_] :=
    Block[ {res={}, tileType,refPt,v1,v2,samplingPt,fcode },
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,fcode} = tlst[[ind]];
            Switch[tileType
              ,1, If[fcode == {} || Last[fcode] == 0,
	                  AppendTo[res,{type4,refPt,{oneoverphi v1,v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type5,refPt+oneoverphi v1+v2,{ (mxRot270.v1),oneoverphi2 (mxRot270.v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
              ,2, If[Last[fcode] == 0,
	                  AppendTo[res,{type3,refPt,{ v1, oneoverphi v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type6,refPt+oneoverphi v2+v1,{oneoverphi2 (mxRot90.v1), (mxRot90.v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
              ,3, If[Last[fcode] == 0,
	                  AppendTo[res,{type1,refPt,{oneoverphi v1, v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type4,refPt+oneoverphi v1,{ (oneoverphi2 v1), (v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
              ,4, If[Last[fcode] == 0,
	                  AppendTo[res,{type2,refPt,{v1,oneoverphi v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type3,refPt+oneoverphi v2,{ (v1), (oneoverphi2 v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
              ,5, If[Last[fcode] == 0,
	                  AppendTo[res,{type3,refPt,{oneoverphi v1, v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type2,refPt+oneoverphi v1 + v2,		{ oneoverphi2 (mxRot270.v1),  (mxRot270.v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
              ,6, If[Last[fcode] == 0,
	                  AppendTo[res,{type4,refPt,{v1,oneoverphi v2},samplingPt,Append[fcode,0]}];
    	              AppendTo[res,{type1,refPt+oneoverphi v2 + v1,{ (mxRot90.v1),oneoverphi2 (mxRot90.v2)},samplingPt,Append[fcode,1]}];
              	  ,(*ELSE: Last[fcode] == 1 *)
	                  AppendTo[res,{tileType,refPt,{v1,v2},samplingPt,Append[fcode,0]}];
              	  ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivFiboSFCTiles *)

getsfcFiboSFC[tlst_] :=
    Block[ {sfc={}, tileType,refPt,v1,v2,samplingPt,fcode,delta=7},
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,fcode} = tlst[[ind]];
   			AppendTo[sfc,refPt + (v1/Norm[v1]+v2/Norm[v2])/1/phi^((Length[fcode]+delta)/2) ];
            Switch[tileType
              ,1, 
  				AppendTo[sfc,refPt + v1 + (v2/Norm[v2]-v1/Norm[v1])/1/phi^((Length[fcode]+delta)/2) ];
              ,2, 
  				AppendTo[sfc,refPt + v2 + (v1/Norm[v1]-v2/Norm[v2])/1/phi^((Length[fcode]+delta)/2) ];
              ,_, 
  				AppendTo[sfc,refPt + v1 + v2 - (v1/Norm[v1]+v2/Norm[v2])/1/phi^((Length[fcode]+delta)/2) ];
            ];
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcFiboSFC *)

getFiboSFCTilesGL[tlst_,params_:showValue+showSFC] :=
    Block[ {gl={},tileType,refPt,v1,v2,samplingPt,fcode,cont,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,fcode} = tlst[[ind]];
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
    		If[BitAnd[params,showValue] > 0, AppendTo[gl,{GrayLevel[FIBOPhi[Reverse@fcode]],Polygon@cont}] ];
			AppendTo[gl,Flatten[#,1]& @ {bortedStyle,Line[cont] } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{-1,-1}]} ] ];		
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
    	,{ind,Length[tlst]}];
    	sfc = getsfcFiboSFC[tlst];
    	If[BitAnd[params,showSFC] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}] ];
    	Return[gl]
    ] (* getFiboSFCTilesGL *)
    

