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
showGrayValue = 8;
showArrows = 16;
showFrame = 32;
showOrdinalNumber = 64;
showSamplingPt = 128;
showTilexycodes = 256;

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



testStrat2DGF3[set_] := 
  Module[{sz,tests,i,tab,sz1d},
	sz=Length[set];
	If[NumberQ[Sqrt[sz]], sz1d = Sqrt[sz]; tests = {{sz1d,sz1d}}];
	If[NumberQ[Sqrt[sz/3]], sz1d = Sqrt[sz/3]; tests = {{3 sz1d,sz1d},{sz1d,3 sz1d}}] ;
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	And @@ tab
] (* testStrat2DGF3 *)

testStrat3DGF3[set_] := 
  Module[{sz,tests,i,tab,sz1d},
	sz=Length[set];
	If[NumberQ[sz^(1/3)], sz1d = sz^(1/3); tests = {{sz1d,sz1d,sz1d}}];
	If[NumberQ[(sz/3)^(1/3)], sz1d = (sz/3)^(1/3); tests = {{3 sz1d,sz1d,sz1d},{sz1d,3 sz1d,sz1d},{sz1d,sz1d,3 sz1d}}] ;
	If[NumberQ[(sz/9)^(1/3)], sz1d = (sz/3)^(1/3); tests = {{3 sz1d,3 sz1d,sz1d},{3 sz1d,sz1d,3 sz1d},{sz1d,3 sz1d,3 sz1d}}] ;
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	And @@ tab
] (* testStrat3DGF3 *)

testStrat2DGFN[base_,set_] := 
  Module[{sz,tests,i,tab,sz1d},
	sz=Length[set];
	If[NumberQ[sz^(1/base)], sz1d = sz^(1/base); tests = {{sz1d,sz1d,sz1d}}];
	If[NumberQ[(sz/base)^(1/base)], sz1d = (sz/base)^(1/base); tests = {{base sz1d,sz1d},{sz1d,base sz1d}}] ;
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	And @@ tab
] (* testStrat2DGFN *)

testStrat3DGFN[base_,set_] := 
  Module[{sz,tests,i,tab,sz1d},
	sz=Length[set];
	If[NumberQ[sz^(1/base)], sz1d = sz^(1/base); tests = {{sz1d,sz1d,sz1d}}];
	If[NumberQ[(sz/base)^(1/base)], sz1d = (sz/base)^(1/base); tests = {{base sz1d,sz1d,sz1d},{sz1d,base sz1d,sz1d},{sz1d,sz1d,base sz1d}}] ;
	If[NumberQ[(sz/(base^2))^(1/base)], sz1d = (sz/base)^(1/base); tests = {{base sz1d,base sz1d,sz1d},{base sz1d,sz1d,base sz1d},{sz1d,base sz1d,base sz1d}}] ;
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ set]] == sz,{i,Length[tests]}];
	And @@ tab
] (* testStrat3DGFN *)

testDyadicPartitioningNDFullGFN[base_,ipts_,showErrFlag_:False] := Module[{sz,powers,tests,i,tab,dim},
	sz=Length[ipts];
	dim = Length[First@ipts];
	powers = Table[base^i,{i,0,Log[base,sz]}];
	tests = Select[Tuples[powers, dim], (Times @@ #) == sz^(dim-1) &];
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ ipts]] == sz,{i,Length[tests]}];
	If[showErrFlag, If[And @@ tab == False, Print["testDyadicPartitioningNDFullGFN: ",Select[{tests,tab}//T,Last[#]==False&]//mf, " -> ", Select[{tests,tab}//T,Last[#]==False&]//Length] ] ];
	And @@ tab
] (* testDyadicPartitioningNDFullGFN *)

testDyadicPartitioningNDFullGFNTfactor[base_,ipts_,tfactor_:0, showErrFlag_:False] := Module[{sz,powers,tests,i,tab,dim},
	sz=Length[ipts];
	dim = Length[First@ipts];
	powers = Table[base^i,{i,0,Log[base,sz]}];
	tests = Select[Tuples[powers, dim], (Times @@ #) == sz^(dim-1) base^(tfactor) &];
	tab = Table[Length[Union[Quotient[#, tests[[i]]] & /@ ipts]] == sz/base^(tfactor),{i,Length[tests]}];
	If[showErrFlag, If[And @@ tab == False, Print["testDyadicPartitioningNDFullGFN: ",Select[{tests,tab}//T,Last[#]==False&]//mf, " -> ", Select[{tests,tab}//T,Last[#]==False&]//Length] ] ];
	And @@ tab
] (* testDyadicPartitioningNDFullGFN *)

getTfactor[base_,ipts_, showErrFlag_:False] :=
    Module[ {tfactor},
    	Do[
    		If[testDyadicPartitioningNDFullGFNTfactor[base,ipts,tf,showErrFlag],
	    		tfactor = tf;
	    		Break[]
    		];
    	,{tf,0,10}];
    	Return[tfactor];
    ] (* testDyadicPartitioningNDFullGFN *)

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
    Module[ {res={}, tileType,sind,samplingPt,refPt,v1,v2,xcode,ycode,fcode },
    	Table[
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
            Switch[tileType
              ,typeSqURdir, 
	                  AppendTo[res,{typeHRectURdir,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeHRectULdir,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeHRectURdir,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
              ,typeSqURinv, 
	                  AppendTo[res,{typeHRectURinv,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeHRectULinv,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeHRectURinv,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,2]} ];
              ,typeSqULdir, 
	                  AppendTo[res,{typeVRectULdir,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeVRectURinv,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeVRectULdir,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,2]} ];
              ,typeSqULinv, 
	                  AppendTo[res,{typeVRectULinv,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeVRectURdir,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeVRectULinv,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,2]} ];
              ,typeHRectURdir, 
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,2]} ];
              ,typeHRectURinv, 
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,2]} ];
              ,typeHRectULdir, 
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,2]} ];
              ,typeHRectULinv, 
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt,			{1/3 v1, v2}, {Append[xcode,0],ycode}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt + 1/3 v1,	{1/3 v1, v2}, {Append[xcode,1],ycode}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt + 2/3 v1,	{1/3 v1, v2}, {Append[xcode,2],ycode}, Append[fcode,2]} ];
             ,typeVRectURdir, 
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
             ,typeVRectURinv, 
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,2]} ];
             ,typeVRectULdir, 
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURdir,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULdir,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,2]} ];
             ,typeVRectULinv, 
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt + 2/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,2]}, Append[fcode,0]} ];
	                  AppendTo[res,{typeSqURinv,sind,samplingPt,refPt + 1/3 v2,	{v1, 1/3 v2}, {xcode,Append[ycode,1]}, Append[fcode,1]} ];
	                  AppendTo[res,{typeSqULinv,sind,samplingPt,refPt,			{v1, 1/3 v2}, {xcode,Append[ycode,0]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivbase3SFCTiles *)

getsfcbase3SFC[tlst_] :=
    Module[ {sfc={}, tileType,sind,refPt,v1,v2,samplingPt,xcode,ycode,fcode,norm1,norm2,delta=3},
    	{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
    	Do[
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
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


getbase3SFCTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[5]},tileType,sind,samplingPt,refPt,v1,v2,xcode,ycode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcbase3SFC[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,(*Arrowheads[1/3^(3+(Length[fcode]+Mod[Length[fcode],2])/2)],*)Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			fcodelen = Length[fcode];
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
			(*samplingPt = refPt + k1 v1 + k2 v2;*)
    		If[BitAnd[params,showGrayValue] > 0, AppendTo[gl,{GrayLevel[FromDigits[Reverse@fcode,3]/3^fcodelen],Polygon@cont}] ];
			AppendTo[gl,Flatten[#,1]& @ {(*Point@(refPt+(norm1+norm2)/20),*)bortedStyle,Line@cont } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[fcode,3],Bold,14,Red],refPt+(v1+v2)/2,{-1.9,-1}]} ] ];		
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
			If[BitAnd[params,showTilexycodes] > 0, AppendTo[gl, {Text[Style[tab2snosep@xcode,Bold,14,Red],refPt+(v1+v2)/2,{1,1}], Text[Style[tab2snosep@ycode,Bold,14,Blue],refPt+(v1+v2)/2,{-1,1}]} ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Point@samplingPt,Text[Style[ind,Bold,14,Blue], samplingPt,{-1.2,-1.2}]} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getbase3SFCTilesGL *)

selectbase3SFCTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsbase3SFCTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,sind,samplingPt,refPt,v1,v2,xcode,ycode,fcode,v,indVect,nsubdivs,m},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
     		nsubdivs = Length[xcode] + Length[ycode];
			v = Join[xcode,ycode];
			m = If[Length[xcode] == Length[ycode],
				mxInv
			,(*ELSE*)
				If[tileType == typeHRectURdir || tileType == typeHRectULdir || tileType == typeHRectURinv || tileType == typeHRectULinv, mxInvH, mxInvV]
			];
			indVect = Mod[#,3]& /@ (m.v);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,sind,samplingPt,refPt,{v1,v2},fcode}] ];
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsbase3SFCTiles *)


getSamplingPtsbase3SFCTiles[tlst_] :=
    Module[ {tileType,sind,samplingPt,refPt,v1,v2,xcode,ycode,fcode},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			(*{k1,k2} = {RandomReal[],RandomReal[]};*)
			refPt + k1 v1 + k2 v2
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsbase3SFCTiles *)

demobase3SFC[innsubdivs_:6, dbg_:False] :=
    Module[ {},
    	nsubdivs = innsubdivs;
		tlst = {{typeSqURdir,0,{0,0},{0,0}, {{1,0},{0,1}}, {{},{}} ,{}} };
		(*Graphics[ getbase3SFCTilesGL[tlst] ]//Print;*)
		
		Do[
			tlst = subdivbase3SFCTiles @ tlst;
			flags = If[iter <= 4, showSFC+showArrows+showTilexycodes+showTileType, showSFC];
			Graphics[ getbase3SFCTilesGL[tlst,flags], PlotLabel-> iter, ImageSize -> {1024,1024} ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,nsubdivs}];
		mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_000001.dat"];
		mxInv = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; nsubdivs/2, ;; nsubdivs]], mxTab[[2, ;; nsubdivs/2, ;; nsubdivs]]] ;
		mxInvH = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; nsubdivs/2, ;; nsubdivs + 1]], mxTab[[2, ;; nsubdivs/2 + 1, ;; nsubdivs + 1]]] ;
		mxInvV = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; nsubdivs/2 + 1, ;; nsubdivs + 1]], mxTab[[2, ;; nsubdivs/2, ;; nsubdivs + 1]]] ;
		tlst = fillSamplingPtsbase3SFCTiles[tlst,mxTab,mxInv,mxInvH,mxInvV];
		Graphics[ {getbase3SFCTilesGL[tlst,showSFC+showSamplingPt]}, PlotLabel-> nsubdivs, ImageSize -> {1024,1024}3/2 ]//Print;

		seltlst = selectbase3SFCTiles[tlst, .33333333];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .33333333 ]//Print;
		seltlst = selectbase3SFCTiles[tlst, .38];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .38 ]//Print;
		seltlst = selectbase3SFCTiles[tlst, .5];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .5 ]//Print;
		seltlst = selectbase3SFCTiles[tlst, .6];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .6 ]//Print;
		seltlst = selectbase3SFCTiles[tlst, .7];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .7 ]//Print;
		seltlst = selectbase3SFCTiles[tlst, .8];
		Graphics[ getbase3SFCTilesGL[seltlst,showGrayValue], PlotLabel-> .8 ]//Print;
	] (* demobase3SFC *)

getDiscrepancy2Dbase3SFC[niters_:4] :=
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

makeMatBuilderMatrices0m2net2D[] :=
    Module[ {},
    	If[ !FileExistsQ["MatBuilder_matrices/"], CreateDirectory["MatBuilder_matrices/"] ];
		Do[
			execString = "testCplex -i MatBuilder_profiles/2D_0m2net.txt -o MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>".dat --seed "<>ToString[RandomInteger[2^16] ]<>" > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	Print[execString -> returnCode];
			mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>".dat"];
        	(*Do[
				mxInv = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; niters/2, ;; niters]], mxTab[[2, ;; niters/2, ;; niters]]] ;
				mxInvH = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; niters/2, ;; niters + 1]], mxTab[[2, ;; niters/2 + 1, ;; niters + 1]]] ;
				mxInvV = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; niters/2 + 1, ;; niters + 1]], mxTab[[2, ;; niters/2, ;; niters + 1]]] ;
        	,{ilevel,15}];*)

		,{i,16}];
    ] (* makeMatBuilderMatrices0m2net2D *)

enumerate0m2netSolutions[] :=
    Module[ {},
    	npts = 3^6;
		tab = Union @ (Union /@ Table[
			execString = "matrixSampler -d 2 -b 3 -n "<>ToString[npts]<>" -i MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>".dat --size 10 -o tmp/pts_"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	Print[i -> execString -> returnCode];
        	prevpts = pts;
        	pts = Round[npts Import["tmp/pts_"<>pid<>".dat"] ];
        	testStrat2DGF3[pts]//Print;
        	pts
		,{i,256}]);
		Print[Length[tab] ];
		
		Graphics[{Red, Point /@ prevpts, Blue, Point /@ pts}]
    ] (* makeMatBuilderMatrices *)
    
readMatBuilderMatrix[fname_,nDims_:2] :=
    Module[ {data,mxsz},
    	data = Import[fname];
    	mxsz = Length[data[[1]]];
    	Table[
    			Table[data[[isz+(idim-1)*(mxsz+1)]],{isz,mxsz}]
    		,{idim,nDims}]
    ]
 
(*----------------------------- end of base3SFC --------------------------------*)
