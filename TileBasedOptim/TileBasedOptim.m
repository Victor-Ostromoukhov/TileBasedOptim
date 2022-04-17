(* TileBasedOptim.m
   2022-04 vo, based on fibo-hilbert.m (version 2002/12/14)
*)
 
(****************** globals *******************)
SetDirectory[ToFileName[$HomeDirectory,"TileBasedOptim/"]];
SetOptions[Graphics, ImageSize -> {3/2 1024, Automatic}];

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

(*------------------------- constants -------------------------*)
phi = N[GoldenRatio,16];
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
showArrows = 16;
showFrame = 32;
showOrdinalNumber = 64;
showSamplingPt = 128;

(*------------------------- end of constants -------------------------*)
FIBOF[symbols_] := With[ {s = Reverse[symbols]}, Total@Table[Fibonacci[i+1 ] s[[i]], {i, Length[s]}] ]
FIBOFinv[symbols_] := Total@Table[Fibonacci[i+1 ] symbols[[i]], {i, Length[symbols]}]
FIBOFxy[symbols_] := FIBOF/@symbols

FIBOphitab = Table[phi^-i, {i, 32}] // N;
FIBOPhi[s_] := Sum[FIBOphitab[[i]] s[[i]], {i, Length[s]}] 

getGeneralizedL2discrepancy[pts_, dbg_:False] :=
    Module[ {execString,nDims = Length[First@pts],prog,returnCode, discrepancy},
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
    	prog = "getGeneralizedL2Discrepancy_from_file";
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
    	prog = "discrepancy";
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  prog<>" -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat -d "<>ToString[nDims]<>" > /dev/null";
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = Import["tmp/res"<>pid<>".dat"][[1,1]];
        Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];
        discrepancy
   ] (* getStarDiscrepancy *)


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


setSamplingPtsFiboSFCTiles[tlst_] :=
    Module[ {res={}, tileType,refPt,v1,v2,k1,k2,fcode},
    	Do[
			{tileType,refPt,{v1,v2},{k1,k2},fcode} = tlst[[ind]];
			k1 = k2 = FIBOPhi[Reverse@fcode];
			If[Min[v1] < 0, k1 = 1-k1];
			If[Min[v2] < 0, k2 = 1-k2];
   			AppendTo[res,{tileType,refPt,{v1,v2},{k1,k2},fcode}];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* setSamplingPtsFiboSFCTiles *)

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
			tlst = setSamplingPtsFiboSFCTiles @ tlst;
			flags = If[iter < 10, showTileType+showTilefcode+showSFC+showArrows+showFrame+showOrdinalNumber+showSamplingPt, showSFC];
			Graphics[ getFiboSFCTilesGL[tlst,flags], PlotLabel-> iter ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,niters}];

		Graphics[ getFiboSFCTilesGL[tlst,showValue] ]//Print;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
		tlst = subdivFiboSFCTiles @ tlst;
		tlst = setSamplingPtsFiboSFCTiles @ tlst;
		Graphics[ getFiboSFCTilesGL[tlst,showSamplingPt] ]//Print;
	] (* demoFiboSFC *)

getDiscrepancy2DFiboSFC[niters_:22] :=
    Module[ {npts,pts,dND, tlst},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
        res = Table[
			tlst = subdivFiboSFCTiles @ tlst;
			tlst = setSamplingPtsFiboSFCTiles @ tlst;
            npts = Length[tlst];
            pts = getSamplingPtsFiboSFCTiles[tlst];
            dND = getGeneralizedL2discrepancy[pts];
			{npts,dND}
		,{iter,niters}];
        Export["data_discrepancyL2/2D/FiboSFC.dat", res]; 
        Print[mf @ res]
    ] (* getDiscrepancy2DFiboSFC *)

showDisrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"FiboSFC",powfromto_:{2,16},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancyOwen,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors},
		{powfrom,powto}=powfromto;
		base=2;
		fontSz = 20;
        dirDiscrepancy = "data_discrepancyL2/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Drop[#,1]& @ Import[dirDiscrepancy<>"WN.dat"])[[;;,1;;2]];
        discrepancyOwen = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Drop[#,1]& @ Import[dirDiscrepancy<>"OwenPure.dat"])[[;;,1;;2]];
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Drop[#,1]& @ Import[dirDiscrepancy<>"Strat.dat"])[[;;,1;;2]], {discrepancyWN[[1]]} ];
        plotLabel = "CascadedSeq GeneralizedL2Discrepancy "<>ToString[nDims]<>"D"<>" GF"<>ToString[base];  
             
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Drop[#,1]& @ Import[dirDiscrepancy<>"Sobol.dat"])[[;;,1;;2]];
        discrepancy2DFiboSFC = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"FiboSFC.dat"], thisDiscrepancy];
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol, discrepancy2DFiboSFC},
        	{discrepancyWN, discrepancySobol, discrepancyStrat, discrepancy2DFiboSFC}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sob", thisDiscrepancyLabel},
       		{"WN","Sob", "Jit", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
        ];

		(*refND = Last /@ getInterpolatedRefsDiscrepancyNDGFN[nDims,base,{powfrom,powto}];*)
        
        range = {Min @@ ((Last /@ #) &@discrepancyOwen), Max @@((Last /@ #) &@discrepancyWN) };
        p = ListLogLogPlot[ alldata
            ,PlotLegends -> Placed[#,{.4,.1}]& @  {Style[#,fontSz]& /@ legends }
            ,PlotStyle -> colors
            ,Joined->True
            ,FrameTicks->{{Automatic,None},{Table[base^pow,{pow,powfrom,powto,1}],Table[base^pow,{pow,powfrom,powto,1}]}}
            ,FrameStyle->Directive[Black,20]
            ,RotateLabel -> True
            ,PlotMarkers->{{\[FilledCircle],5} }
            ,Frame->True
             ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "discrepancy", fontSz] }
               ,ImageSize -> {1024,1024}
            (*,PlotRange->{{base^powfrom,base^powto},{.0003,.3}}*)
            (*,PlotRange->{{base^powfrom,base^powto},{.1,Min[Last /@ thisDiscrepancy]}}*)
            ,PlotRange->{{base^powfrom,base^powto},range}
            ,GridLines->{Table[base^pow,{pow,powfrom,powto,1}],None}
            ,GridLinesStyle->Directive[Darker@Gray, Dashed]
            ,AspectRatio->1
            ,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
            ,PlotLabel -> Style[ plotLabel, Bold, 24]
        ];
        p//Print;
        p
    ] (* showDisrepancyND *)

testDyadicPartitioningNDFull[set_,showErrFlag_:True] := Module[{sz,powers,tests,i,tab,dim},
	sz=Length[set];
	dim = Length[First@set];
	powers = Table[2^i,{i,0,Log[2,sz]}];
	tests = Select[Tuples[powers, dim], (Times @@ #) == sz^(dim-1) &];
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	If[showErrFlag, If[And @@ tab == False, Print["testHierarchicalStratified3DStrong: ",Select[{tests,tab}//T,Last[#]==False&]//mf, " -> ", Select[{tests,tab}//T,Last[#]==False&]//Length] ] ];
	And @@ tab
] (* testDyadicPartitioningNDFull *)

tstStarBinary[nlevels_:7] :=
    Module[ {},
        dtab = Table[
			npts = 4^ilevel;
			codes = Partition[#,ilevel]& /@ (IntegerDigits[#,2,2*ilevel]& /@ Range[0,npts-1]);
			pts = Table[
					{codex,codey} = codes[[i]];
					ix = FromDigits[#,2]& @ codex;
					iy = FromDigits[#,2]& @ codey;
					ixfrac = FromDigits[#,2]& @ Reverse[codex];
					iyfrac = FromDigits[#,2]& @ Reverse[codey];
					{x,y} = {ix / 2^ilevel + iyfrac / npts, iy / 2^ilevel + ixfrac / npts}//N
				,{i,npts}];
			ipts = Round[ npts pts ];
			Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]];
			{npts,getStarDiscrepancy[pts]}
        ,{ilevel,nlevels}];
        showDisrepancyND[2,dtab,"tstBinary"];
    ]

    
makeSobolDiscrepancy[nlevels_:18] :=
    Module[ {},
        dtab = Table[
			npts = 2^ilevel;
        	Print["Processing makeSobolDiscrepancy level ",ilevel, " npts = ",npts];
			execString = "owen -n "<>ToString[npts]<>" --nd 2 -p 0 -o tmp/pts"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	pts = Import["tmp/pts"<>pid<>".dat"];
			ipts = Round[ npts pts ];
			If[dbg,Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
			{npts,getStarDiscrepancy[pts]}
        ,{ilevel,nlevels}];
        Export["data_StarDiscrepancy/2D/Sobol.dat", dtab]; 
    ]

makeWNStarDiscrepancy[nlevels_:9, ntrials_:64] :=
    Module[ {},
        dtab = Table[
			npts = 4^ilevel;
        	Print["Processing makeWNStarDiscrepancy level ",ilevel, " npts = ",npts];
			trials = Parallelize @ Table[
				pts = Table[{RandomReal[],RandomReal[]},{i,npts}];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			Mean @ trials
        ,{ilevel,nlevels}];
        showDisrepancyND[2,dtab,"tstBinary"];
        Export["data_StarDiscrepancy/2D/WN.dat", dtab]; 
        Print[mf @ res]
    ]

makeStratStarDiscrepancy[nlevels_:9, ntrials_:64] :=
    Module[ {},
        dtab = Table[
			npts = 4^ilevel;
			nstrats = 2^ilevel;
        	Print["Processing makeStratStarDiscrepancy level ",ilevel, " npts = ",npts];
			trials = Parallelize @ Table[
				pts = N @ Table[{Mod[i,nstrats], Quotient[i,nstrats]}/nstrats + {RandomReal[],RandomReal[]}/nstrats,{i,0,npts-1}];
				If[dbg,Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			Mean @ trials
        ,{ilevel,nlevels}];
        showDisrepancyND[2,dtab,"tstBinary"];
        Export["data_StarDiscrepancy/2D/Strat.dat", dtab]; 
        Print[mf @ res]
    ]

 (*
 gitpull
 math
 <<TileBasedOptim/TileBasedOptim.m
 makeWNStarDiscrepancy[]
makeStratStarDiscrepancy[]
makeSobolDiscrepancy[]
 *)