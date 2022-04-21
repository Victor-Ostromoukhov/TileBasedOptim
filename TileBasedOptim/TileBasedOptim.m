(* TileBasedOptim.m
   2022-04 vo, based on fibo-hilbert.m (version 2002/12/14)
*)
 
(****************** globals *******************)
SetDirectory[ToFileName[$HomeDirectory,"TileBasedOptim/"]];
SetOptions[Graphics, ImageSize -> {1024, Automatic}];

mf := MatrixForm
T :=  Transpose
PI = Pi//N;
known := ValueQ
i2s[n_,len_:6] := ToString[NumberForm[n, len-1, NumberPadding -> "0"]]
r2s[n_,len_:6,dec_:5] := ToString[NumberForm[n, {len,dec}, NumberPadding -> "0"]]
tab2s[tab_] := StringJoin @ Table[" "<>ToString[tab[[i]]]<>" ",{i,Length[tab]}]
tab2snosep[tab_] := StringJoin @ Table[ ToString[tab[[i]]] ,{i,Length[tab]}]
tab2ssep[tab_] := StringJoin @ Table["_"<>ToString[tab[[i]]],{i,Length[tab]}]
tab2ssepComma[tab_] := StringJoin[ToString[tab[[1]]], Table[","<>ToString[tab[[i]]],{i,Min[2,Length[tab]], Length[tab]}] ]
str2n[str_] := FromDigits[#, 2] & @ Table[StringTake[str, {i}] // ToExpression, {i, StringLength[str]}]

pid := "_pid"<>ToString[$ProcessID]<>"_kid"<>ToString[$KernelID]
execPrefix = "~/bin/";

Print["module TileBasedOptim loaded."];

(*------------------------- constants -------------------------*)
phi = N[GoldenRatio,16];
oneoverphi = 1/phi;
oneoverphi2 = oneoverphi^2;

mxRot0 =   {{1,0}, {0,1}};
mxRot90 =  {{0, -1}, {1, 0}};
mxRot180 = {{-1,0}, {0,-1}};
mxRot270 = {{0, 1}, {-1, 0}};


(* debugging flags *)
showTileType = 1;
showTilefcode = 2;
showSFC = 4;
showValue = 8;
showArrows = 16;
showFrame = 32;
showOrdinalNumber = 64;
showSamplingPt = 128;

(*------------------------- end of constants -------------------------*)

testDyadicPartitioningNDFull[set_,showErrFlag_:True] := Module[{sz,powers,tests,i,tab,dim},
	sz=Length[set];
	dim = Length[First@set];
	powers = Table[2^i,{i,0,Log[2,sz]}];
	tests = Select[Tuples[powers, dim], (Times @@ #) == sz^(dim-1) &];
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	If[showErrFlag, If[And @@ tab == False, Print["testHierarchicalStratified3DStrong: ",Select[{tests,tab}//T,Last[#]==False&]//mf, " -> ", Select[{tests,tab}//T,Last[#]==False&]//Length] ] ];
	And @@ tab
] (* testDyadicPartitioningNDFull *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
tstStarBinary2D[]
*)

tstStarBinary2D[nlevels_:10, dstarFlag_:True] :=
    Module[ {},
        dtab = Table[
			npts = 4^ilevel;
			nstrats = 2^ilevel;
			Print["Processing tstStarBinary3D level ",ilevel, " npts = ",npts, " nstrats = ",nstrats];
			codes = Partition[#,ilevel]& /@ (IntegerDigits[#,2,2*ilevel]& /@ Range[0,npts-1]);
			pts = Table[
					{codex,codey} = codes[[i]];
					ix = FromDigits[#,2]& @ codex;
					iy = FromDigits[#,2]& @ codey;
					ixfrac = FromDigits[#,2]& @ Reverse[codey];
					iyfrac = FromDigits[#,2]& @ Reverse[codex];
					{ix / 2^ilevel + ixfrac / npts, iy / 2^ilevel + iyfrac / npts}//N
				,{i,npts}];
			ipts = Round[ npts pts ];
			(*Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]];*)
			{npts,If[dstarFlag,getStarDiscrepancy[pts],getGeneralizedL2discrepancy[pts]]}
        ,{ilevel,nlevels}];
        If[dstarFlag,
			showStarDisrepancyND[2,dtab,"tstStarBinary2D"];
        ,(*ELSE*)
       		showGeneralizedL2discrepancyND[2,dtab,"GeneralizedL2discrepancy2D"];
        ];
    ] (* tstStarBinary2D *)

tstStarBinary3D[nlevels_:2] :=
    Module[ {},
        dtab = Table[
			npts = 8^ilevel;
			nstrats = 2^ilevel;
			Print["Processing tstStarBinary3D level ",ilevel, " npts = ",npts, " nstrats = ",nstrats];
			pts = Flatten[#,2]& @ Table[
					{codex,codey,codez} = {IntegerDigits[ix,2,ilevel],IntegerDigits[iy,2,ilevel],IntegerDigits[iz,2,ilevel]};
					(*codex2 = Reverse @ Flatten[  T[{codey,codez}]];
					codey2 = Reverse @ Flatten[  T[{codez,codex}]];
					codez2 = Reverse @ Flatten[  T[{codex,codey}]];*)

					codex2 = Reverse @ Flatten[  T[{BitXor@@{codey,codez},1-BitXor@@{codey,codez} }]];
					codey2 = Reverse @ Flatten[  T[{BitXor@@{codez,codex},1-BitXor@@{codez,codex} }]];
					codez2 = Reverse @ Flatten[  T[{BitXor@@{codex,codey},1-BitXor@@{codex,codey} }]];

					ixfrac = FromDigits[#,2]& @ codex2;
					iyfrac = FromDigits[#,2]& @ codey2;
					izfrac = FromDigits[#,2]& @ codez2;
					If[ilevel == 2, Print[{ix,iy,iz} -> mf[{codex,codey,codez}] -> mf[{T[{codey,codez}],T[{codez,codex}],T[{codex,codey}]}] -> mf[{codex2,codey2,codez2}]]];
					{ix / 2^ilevel + ixfrac / npts, iy / 2^ilevel + iyfrac / npts, iz / 2^ilevel + izfrac / npts}//N
				,{ix,0,nstrats-1},{iy,0,nstrats-1},{iz,0,nstrats-1}];
			{npts,getStarDiscrepancy[pts]}
        ,{ilevel,nlevels}];
        showStarDisrepancyND[3,dtab,"tstStarBinary3D",{2,12}];
        
        Graphics3D[{Point /@ pts}]//Print;
        Graphics[{AbsolutePointSize[10],Point /@ (Drop[#,-1]& /@ pts)}, Framed->True,ImageSize->{1024,1024}/2]//Print;
    ] (* tstStarBinary2D *)

(*-------------------------- calculate & show discrepancy --------------------------*)
getGeneralizedL2discrepancy[pts_, dbg_:False] :=
    Module[ {execString,nDims = Length[First@pts],prog,returnCode, discrepancy},
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
    	prog = "getGeneralizedL2Discrepancy";
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  prog<>" -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat -d "<>ToString[nDims]<>" > /dev/null";
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = Import["tmp/res"<>pid<>".dat"][[1,2]];
        Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];
        discrepancy
   ] (* getGeneralizedL2discrepancy *)

getStarDiscrepancy[pts_, dbg_:False] :=
    Module[ {execString,nDims = Length[First@pts],prog,returnCode, discrepancy},
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
    	prog = "getStarDiscrepancy";
    	
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  If[nDims == 2,
        	"StarDiscrepancy_fromfile_2dd -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat  > /dev/null"
        ,(*ELSE*)
        	"getStarDiscrepancy -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat -d "<>ToString[nDims]<>" > /dev/null"
        ];
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = If[nDims == 2, Import["tmp/res"<>pid<>".dat"][[2,2]], Import["tmp/res"<>pid<>".dat"][[1,1]] ];
        (*Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];*)
        discrepancy
   ] (* getStarDiscrepancy *)

getCloseestN2D[n_] := Round[Sqrt[n]]^2
getCloseestNND[nDims_:2, n_] := Round[n^(1/nDims)]^nDims

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeWNStarDiscrepancy[]
makeStratStarDiscrepancy[]
makeSobolDiscrepancy[]
*)
makeSobolDiscrepancy[nlevels_:14, nDims_:2,dbg_:True] :=
    Module[ {},
        dtab = {};
        Do[
			npts = 2^ilevel;
        	Print["Processing makeSobolDiscrepancy level ",ilevel, " npts = ",npts, " nDims = ",nDims];
			execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 0 -o tmp/pts"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	pts = Import["tmp/pts"<>pid<>".dat"];
			ipts = Round[ npts pts ];
			If[dbg,Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
			AppendTo[dtab, {npts,getStarDiscrepancy[pts]} ];
	        (*Export["data_StarDiscrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; *)
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]

makeWNStarDiscrepancy[nlevels_:12, ntrials_:64, nDims_:3] :=
    Module[ {},
        dtab = {};
        Do[
			npts = 2^ilevel;
        	Print["Processing makeWNStarDiscrepancy level ",ilevel, " npts = ",npts, " nDims = ",nDims];
			trials = Parallelize @ Table[
				pts = Table[Table[RandomReal[],{nDims}],{i,npts}];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
	        (*Export["data_StarDiscrepancy/"<>ToString[nDims]<>"D/WN.dat", dtab]; *)
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]

makeStratStarDiscrepancy[nlevels_:12, ntrials_:64, nDims_:3] :=
    Module[ {},
        dtab = {};
        Do[
			npts = getCloseestNND[nDims, 2^ilevel];
			nstrats = npts^(1/nDims);
        	Print["Processing makeStratStarDiscrepancy level ",ilevel, " npts = ",npts, " nDims = ",nDims];
			trials = (*Parallelize @*) Table[
				Switch[nDims
				,2, pts = N @ Table[{Mod[i,nstrats], Quotient[i,nstrats]}/nstrats + {RandomReal[],RandomReal[]}/nstrats,{i,0,npts-1}];
				,3, pts = Flatten[#,2]& @ Table[{ix+RandomReal[],iy+RandomReal[],iz+RandomReal[]}/nstrats,{ix,0,nstrats-1},{iy,0,nstrats-1},{iz,0,nstrats-1}];
				];
				If[dbg,Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
	        Export["data_StarDiscrepancy/"<>ToString[nDims]<>"D/Strat.dat", dtab]; push
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]

showStarDisrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"FiboSFC",powfromto_:{2,16},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancy2DFiboSFC},
		{powfrom,powto}=powfromto;
		fontSz = 20;
        dirDiscrepancy = "data_StarDiscrepancy/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);
        discrepancy2DFiboSFC = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"FiboSFC.dat"], thisDiscrepancy];
        plotLabel = " StarDisrepancy "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol, discrepancy2DFiboSFC},
        	{discrepancyWN, discrepancySobol, discrepancyStrat, discrepancy2DFiboSFC}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol", thisDiscrepancyLabel},
       		{"WN","Sobol", "Jitter", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
        ];
        range = {Min @@ ((Last /@ #) &@discrepancySobol), Max @@((Last /@ #) &@discrepancyWN) };
        p = ListLogLogPlot[ alldata
            ,PlotLegends -> Placed[#,{.4,.1}]& @  {Style[#,fontSz]& /@ legends }
            ,PlotStyle -> colors
            ,Joined->True
            ,FrameTicks->{{Automatic,None},{Table[base^pow,{pow,powfrom,powto,1}],Table[base^pow,{pow,powfrom,powto,1}]}}
            ,FrameStyle->Directive[Black,20]
            ,RotateLabel -> True
            ,PlotMarkers->{{\[FilledCircle],5} }
            ,Frame->True
            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "StarDisrepancy", fontSz] }
            ,ImageSize -> {1024,1024}
            ,PlotRange->{{base^powfrom,base^powto},range}
            ,GridLines->{Table[base^pow,{pow,powfrom,powto,1}],None}
            ,GridLinesStyle->Directive[Darker@Gray, Dashed]
            ,AspectRatio->1
            ,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
            ,PlotLabel -> Style[ plotLabel, Bold, 24]
        ];
        p//Print;
        (*p*)
    ] (* showStarDisrepancyND *)

 showGeneralizedL2discrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"FiboSFC",powfromto_:{2,16},col_:Red] := 
	Module[{(*dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancy2DFiboSFC*)},
		{powfrom,powto}=powfromto;
    	fontSz = 20;
		kPlusMinus = .5;
		base = 2;
    	{powfrom,powto,powstep} = {2,10,2};
    	selPows = 2^Range[powfrom,powto,powstep];
        dirDiscrepancy = "data_discrepancyL2/"<>ToString[nDims]<>"D/";
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"WN.dat"]);
			discrepancyWN = Select[#,MemberQ[selPows,#[[1]]]&]& @ Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"Sobol.dat"]);
			discrepancySobol = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];				
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"Strat.dat"]);
			discrepancyStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];				

        (*discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);*)
        
        discrepancy2DFiboSFC = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"FiboSFC.dat"], thisDiscrepancy];
        plotLabel = " GeneralizedL2discrepancyND "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol, discrepancy2DFiboSFC},
        	{discrepancyWN, discrepancySobol, discrepancyStrat, discrepancy2DFiboSFC}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol", thisDiscrepancyLabel},
       		{"WN","Sobol", "Jitter", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
        ];
        range = {Min @@ ((Last /@ #) &@discrepancySobol), Max @@((Last /@ #) &@discrepancySobol) };
        p = ListLogLogPlot[ alldata
            ,PlotLegends -> Placed[#,{.4,.1}]& @  {Style[#,fontSz]& /@ legends }
            ,PlotStyle -> colors
            ,Joined->True
            ,FrameTicks->{{Automatic,None},{Table[base^pow,{pow,powfrom,powto,1}],Table[base^pow,{pow,powfrom,powto,1}]}}
            ,FrameStyle->Directive[Black,20]
            ,RotateLabel -> True
            ,PlotMarkers->{{\[FilledCircle],5} }
            ,Frame->True
            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "GeneralizedL2discrepancyND", fontSz] }
            ,ImageSize -> {1024,1024}
            ,PlotRange->{{base^powfrom,base^powto},range}
            ,GridLines->{Table[base^pow,{pow,powfrom,powto,1}],None}
            ,GridLinesStyle->Directive[Darker@Gray, Dashed]
            ,AspectRatio->1
            ,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
            ,PlotLabel -> Style[ plotLabel, Bold, 24]
        ];
        p//Print;
        (*p*)
    ] (* showGeneralizedL2discrepancyND *)
 
(*----------------------------- FiboSFC --------------------------------*)
type1 = 1;
type2 = 2;
type3 = 3;
type4 = 4;
type5 = 5;
type6 = 6;
FIBOF[symbols_] := With[ {s = Reverse[symbols]}, Total@Table[Fibonacci[i+1 ] s[[i]], {i, Length[s]}] ]
FIBOFinv[symbols_] := Total@Table[Fibonacci[i+1 ] symbols[[i]], {i, Length[symbols]}]
FIBOFxy[symbols_] := FIBOF/@symbols

FIBOphitab = Table[phi^-i, {i, 32}] // N;
FIBOPhi[s_] := Sum[FIBOphitab[[i]] s[[i]], {i, Length[s]}] 

subdivFiboSFCTiles[tlst_] :=
    Module[ {res={}, tileType,refPt,v1,v2,samplingPt,fcode },
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
    Module[ {sfc={}, tileType,refPt,v1,v2,samplingPt,fcode,delta=7},
    	Do[
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


getFiboSFCTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[5]},tileType,refPt,v1,v2,samplingPt,fcode,cont,sfc,norm1,norm2,k1,k2,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	Do[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/phi^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
			samplingPt = refPt + k1 v1 + k2 v2;
    		If[BitAnd[params,showValue] > 0, AppendTo[gl,{GrayLevel[FIBOPhi[Reverse@fcode]],Polygon@cont}] ];
			AppendTo[gl,Flatten[#,1]& @ {bortedStyle,Line@cont } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FIBOPhi[Reverse@fcode],Bold,14,Red],refPt+(v1+v2)/2,{-1.9,-1}]} ] ];		
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
			If[BitAnd[params,showFrame] > 0, AppendTo[gl, {Arrowheads[1/phi^(6+(Length[fcode]+Mod[Length[fcode],2])/2)],Red,Arrow[{refPt+(norm1+norm2)/10,refPt+norm1/2+norm2/10}],Blue,Arrow[{refPt+(norm1+norm2)/10,refPt+norm2/2+norm1/10}]} ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, Point@samplingPt ] ];
    	,{ind,Length[tlst]}];
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcFiboSFC[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Return[gl]
    ] (* getFiboSFCTilesGL *)


fillSamplingPtsFiboSFCTiles[tlst_] :=
    Module[ {res={}, tileType,refPt,v1,v2,k1,k2,fcode},
    	Do[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			k1 = k2 = FIBOPhi[Reverse@fcode];
			If[Min[v1] < 0, k1 = 1-k1];
			If[Min[v2] < 0, k2 = 1-k2];
   			AppendTo[res,{tileType,refPt,{v1,v2},{k1,k2},fcode}];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* fillSamplingPtsFiboSFCTiles *)

getSamplingPtsFiboSFCTiles[tlst_] :=
    Module[ {tileType,refPt,v1,v2,k1,k2,fcode},
    	Parallelize @ Table[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			(*{k1,k2} = {RandomReal[],RandomReal[]};*)
			refPt + k1 v1 + k2 v2
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsFiboSFCTiles *)

demoFiboSFC[niters_:15] :=
    Module[ {dbg},
    	dbg = False;
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
		Graphics[ getFiboSFCTilesGL[tlst] ]//Print;
		Do[
			tlst = subdivFiboSFCTiles @ tlst;
			tlst = fillSamplingPtsFiboSFCTiles @ tlst;
			flags = If[iter < 10, showTileType+showTilefcode+showSFC+showArrows+showFrame+showOrdinalNumber+showSamplingPt, showSFC];
			Graphics[ getFiboSFCTilesGL[tlst,flags], PlotLabel-> iter ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,niters}];

		Graphics[ getFiboSFCTilesGL[tlst,showValue] ]//Print;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = fillSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
	] (* demoFiboSFC *)

getDiscrepancy2DFiboSFC[niters_:22] :=
    Module[ {npts,pts,dND, tlst},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
        res = Table[
			tlst = subdivFiboSFCTiles @ tlst;
			tlst = fillSamplingPtsFiboSFCTiles @ tlst;
            npts = Length[tlst];
            pts = getSamplingPtsFiboSFCTiles[tlst];
            dND = getGeneralizedL2discrepancy[pts];
			{npts,dND}
		,{iter,niters}];
        Export["data_discrepancyL2/2D/FiboSFC.dat", res]; 
        Print[mf @ res]
    ] (* getDiscrepancy2DFiboSFC *)

(*----------------------------- end of FiboSFC --------------------------------*)




(*----------------------------- base3SFC --------------------------------*)
typeSqURdir = 		1;
typeSqULdir = 		2;
typeHRectURdir = 	3;
typeHRectULdir = 	4;
typeVRectURdir = 	5;
typeVRectULdir = 	6;
typeSqURinv = 		11;
typeSqULinv = 		12;
typeHRectURinv = 	13;
typeHRectULinv = 	14;
typeVRectURinv = 	15;
typeVRectULinv = 	16;

subdivbase3SFCTiles[tlst_] :=
    Module[ {res={}, tileType,refPt,v1,v2,samplingPt,xcode,ycode,fcode },
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,{xcode,ycode},fcode} = tlst[[ind]];
            Switch[tileType
              ,typeSqURdir, 
	                  AppendTo[res,{typeHRectURdir,refPt,			{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeHRectULdir,refPt + 1/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeHRectURdir,refPt + 2/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
              (*,typeSqURdir, 
	                  AppendTo[res,{typeVRectURdir,refPt,			{1/3 v1, v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeVRectULinv,refPt + 1/3 v1,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeVRectURdir,refPt + 2/3 v1,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];*)
              ,typeSqURinv, 
	                  AppendTo[res,{typeHRectURinv,refPt + 2/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeHRectULinv,refPt + 1/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeHRectURinv,refPt,			{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
              ,typeSqULdir, 
	                  AppendTo[res,{typeVRectULdir,refPt + 2/3 v1,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeVRectURinv,refPt + 1/3 v1,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeVRectULdir,refPt,			{1/3 v1, v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
              (*,typeSqULinv, 
	                  AppendTo[res,{typeHRectULinv,refPt + 2/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeHRectURinv,refPt + 1/3 v2,	{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeHRectULinv,refPt,			{v1, 1/3 v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];*)
              ,typeSqULinv, 
	                  AppendTo[res,{typeVRectULinv,refPt,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeVRectURdir,refPt + 1/3 v1,	{1/3 v1, v2},samplingPt, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeVRectULinv,refPt + 2/3 v1,			{1/3 v1, v2},samplingPt, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
              ,typeHRectURdir, 
	                  AppendTo[res,{typeSqURdir,refPt,			{1/3 v1, v2},samplingPt, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULinv,refPt + 1/3 v1,	{1/3 v1, v2},samplingPt, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURdir,refPt + 2/3 v1,	{1/3 v1, v2},samplingPt, {Append[xcode,2],ycode}, Append[fcode,2]} ];
              ,typeHRectULdir, 
	                  AppendTo[res,{typeSqULdir,refPt + 2/3 v1,	{1/3 v1, v2},samplingPt, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURinv,refPt + 1/3 v1,	{1/3 v1, v2},samplingPt, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULdir,refPt,			{1/3 v1, v2},samplingPt, {Append[xcode,2],ycode}, Append[fcode,2]} ];
             ,typeVRectURdir, 
	                  AppendTo[res,{typeSqURdir,refPt,			{v1, 1/3 v2},samplingPt, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULdir,refPt + 1/3 v2,	{v1, 1/3 v2},samplingPt, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURdir,refPt + 2/3 v2,	{v1, 1/3 v2},samplingPt, {Append[xcode,2],ycode}, Append[fcode,2]} ];
             ,typeVRectULinv, 
	                  AppendTo[res,{typeSqULinv,refPt + 2/3 v2,	{v1, 1/3 v2},samplingPt, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURinv,refPt + 1/3 v2,	{v1, 1/3 v2},samplingPt, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULinv,refPt,			{v1, 1/3 v2},samplingPt, {Append[xcode,2],ycode}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivbase3SFCTiles *)

getsfcbase3SFC[tlst_] :=
    Module[ {sfc={}, tileType,refPt,v1,v2,samplingPt,xcode,ycode,fcode,norm1,norm2,delta=5},
    	{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
    	Do[
			{tileType,refPt,{v1,v2},samplingPt,{xcode,ycode},fcode} = tlst[[ind]];
            If[tileType == typeSqURdir || tileType == typeHRectURdir || tileType == typeVRectURdir,
	   			AppendTo[sfc,refPt + (norm1+norm2)/1/3^((Length[fcode]+delta)/2) ];
  				AppendTo[sfc,refPt + v1 + v2 + (-norm1-norm2)/3^((Length[fcode]+delta)/2) ] ;
  			];
			If[tileType == typeSqULdir || tileType == typeHRectULdir || tileType == typeVRectULdir,
	   			AppendTo[sfc,refPt + v1 + (-norm1+norm2)/1/3^((Length[fcode]+delta)/2) ];
  				AppendTo[sfc,refPt + v2 + (norm1-norm2)/3^((Length[fcode]+delta)/2) ];
            ];
            If[tileType == typeSqURinv || tileType == typeHRectURinv || tileType == typeVRectURinv,
  				AppendTo[sfc,refPt + v1 + v2 + (-norm1-norm2)/3^((Length[fcode]+delta)/2) ] ;
	   			AppendTo[sfc,refPt + (norm1+norm2)/1/3^((Length[fcode]+delta)/2) ];
  			];
			If[tileType == typeSqULinv || tileType == typeHRectULinv || tileType == typeVRectULinv,
  				AppendTo[sfc,refPt + v2 + (norm1-norm2)/3^((Length[fcode]+delta)/2) ];
	   			AppendTo[sfc,refPt + v1 + (-norm1+norm2)/1/3^((Length[fcode]+delta)/2) ];
            ];
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcbase3SFC *)


getbase3SFCTilesGL[tlst_,params_:showSFC+showArrows+showTileType] :=
    Module[ {gl={AbsolutePointSize[5]},tileType,refPt,v1,v2,samplingPt,fcode,cont,sfc,norm1,norm2,k1,k2,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	Do[
			{tileType,refPt,{v1,v2},{k1,k2},{xcode,ycode},fcode} = tlst[[ind]];
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
			samplingPt = refPt + k1 v1 + k2 v2;
    		(*If[BitAnd[params,showValue] > 0, AppendTo[gl,{GrayLevel[FIBOPhi[Reverse@fcode]],Polygon@cont}] ];*)
			AppendTo[gl,Flatten[#,1]& @ {Point@(refPt+(norm1+norm2)/20),bortedStyle,Line@cont } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Red],refPt+(v1+v2)/2,{-1.9,-1}]} ] ];		
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, Point@samplingPt ] ];
    	,{ind,Length[tlst]}];
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcbase3SFC[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,(*Arrowheads[1/3^(3+(Length[fcode]+Mod[Length[fcode],2])/2)],*)Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Return[gl]
    ] (* getbase3SFCTilesGL *)

fillSamplingPtsbase3SFCTiles[tlst_] :=
    Module[ {res={}, tileType,refPt,v1,v2,k1,k2,fcode},
    	Do[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			k1 = k2 = base3Phi[Reverse@fcode];
			If[Min[v1] < 0, k1 = 1-k1];
			If[Min[v2] < 0, k2 = 1-k2];
   			AppendTo[res,{tileType,refPt,{v1,v2},{k1,k2},fcode}];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* fillSamplingPtsbase3SFCTiles *)

getSamplingPtsbase3SFCTiles[tlst_] :=
    Module[ {tileType,refPt,v1,v2,k1,k2,fcode},
    	Parallelize @ Table[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			(*{k1,k2} = {RandomReal[],RandomReal[]};*)
			refPt + k1 v1 + k2 v2
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsbase3SFCTiles *)

demobase3SFC[niters_:3, dbg_:False] :=
    Module[ {},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {{},{}} ,{}} };
		Graphics[ getbase3SFCTilesGL[tlst] ]//Print;
		
		Do[
			tlst = subdivbase3SFCTiles @ tlst;
			Graphics[ getbase3SFCTilesGL[tlst], PlotLabel-> iter ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,niters}];
Abort[];

		Graphics[ getbase3SFCTilesGL[tlst,showValue] ]//Print;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivbase3SFCTiles @ tlst;
		tlst = fillSamplingPtsbase3SFCTiles @ tlst;
		Graphics[ getbase3SFCTilesGL[tlst,showSamplingPt] ]//Print;
	] (* demobase3SFC *)

getDiscrepancy2Dbase3SFC[niters_:22] :=
    Module[ {npts,pts,dND, tlst},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
        res = Table[
			tlst = subdivbase3SFCTiles @ tlst;
			tlst = fillSamplingPtsbase3SFCTiles @ tlst;
            npts = Length[tlst];
            pts = getSamplingPtsbase3SFCTiles[tlst];
            dND = getGeneralizedL2discrepancy[pts];
			{npts,dND}
		,{iter,niters}];
        Export["data_discrepancyL2/2D/base3SFC.dat", res]; 
        Print[mf @ res]
    ] (* getDiscrepancy2Dbase3SFC *)

(*----------------------------- end of FiboSFC --------------------------------*)
