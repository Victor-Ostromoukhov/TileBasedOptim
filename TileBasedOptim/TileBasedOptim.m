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
showLightGrayTile = 512;
showFcodeInvNumber = 64;

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

getL2discrepancy[pts_, dbg_:False] :=
    Module[ {execString,nDims = Length[First@pts],prog,returnCode, discrepancy},
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
        prog =  Switch[nDims
        	,2, "L2Discrepancy_fromfile_2dd"
        	,3, "L2Discrepancy_fromfile_3dd"
        	,4, "L2Discrepancy_fromfile_4dd"
        ];
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  prog<>" -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat > /dev/null";
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = Import["tmp/res"<>pid<>".dat"][[2,2]];
        Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];
        discrepancy
   ] (* getL2discrepancy *)

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
makeSobolStarDiscrepancy[]
makeSobolGeneralizedL2Discrepancy[]



makeWNL2Discrepancy[]
makeStratL2Discrepancy[]

makeOwenL2Discrepancy[]

makeSobolL2Discrepancy[]

makeWNStarDiscrepancy[]
makeStratStarDiscrepancy[]
makeSobolDiscrepancy[]
*)


makeSobolL2Discrepancy[nlevels_:14, nDims_:2,dbg_:False] :=
    Module[ {},
        dtab = {};
        nptsMax = 2^nlevels;
        Do[
			npts = inpts;
			execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 0 -o tmp/pts"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	pts = Import["tmp/pts"<>pid<>".dat"];		
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getL2discrepancy[pts];
        	Print["Processing makeSobolL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; 
        ,{inpts,nptsMax}];
        Print[mf @ dtab]
    ] (* makeSobolL2Discrepancy *)

makeSobolGeneralizedL2Discrepancy[nlevels_:12, nDims_:2,dbg_:False] :=
    Module[ {},
        nptsMax = 2^nlevels;
        dtab = Parallelize @ Table[
			npts = inpts;
			execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 0 -o tmp/pts"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	pts = Import["tmp/pts"<>pid<>".dat"];		
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getGeneralizedL2discrepancy[pts];
        	Print["Processing makeSobolGeneralizedL2Discrepancy " -> {npts,d}];
			{npts,d} 
        ,{inpts,nptsMax}];
	    Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; 
        Print[mf @ dtab]
    ] (* makeSobolL2Discrepancy *)

makeSobolStarDiscrepancy[nlevels_:12, nDims_:2,dbg_:False] :=
    Module[ {},
        nptsMax = 2^nlevels;
        dtab = Parallelize @ Table[
			npts = inpts;
			execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 0 -o tmp/pts"<>pid<>".dat > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	pts = Import["tmp/pts"<>pid<>".dat"];		
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getStarDiscrepancy[pts];
        	Print["Processing makeSobolStarDiscrepancy " -> {npts,d}];
			{npts,d}
        ,{inpts,nptsMax}];
	    Export["data_StarDiscrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; 
        Print[mf @ dtab]
    ] (* makeSobolL2Discrepancy *)

makeOwenL2Discrepancy[nlevels_:14, nDims_:2,dbg_:False] :=
    Module[ {},
        dtab = {};
        nptsMax = 2^nlevels;
        Do[
			npts = inpts;
        	d = Mean @ (Parallelize @ Table[
				execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 1 --max_tree_depth_32_flag 0 -s "<>ToString[RandomInteger[2^16]]<>" -o tmp/pts"<>pid<>".dat > /dev/null";
        		returnCode = Run[execPrefix<>execString];
        		pts = Import["tmp/pts"<>pid<>".dat"];		
        		getL2discrepancy[pts]
        	,{16}]);
        	Print["Processing makeOwenL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Owen.dat", dtab]; 
        ,{inpts,nptsMax}];
        Print[mf @ dtab]
    ] (* makeOwenL2Discrepancy *)


makeWNL2Discrepancy[nlevels_:14, ntrials_:64, nDims_:2] :=
    Module[ {},
        dtab = {};
        Do[
			npts = 2^ilevel;
			trials = Parallelize @ Table[
				pts = Table[Table[RandomReal[],{nDims}],{i,npts}];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
        	Print["Processing makeWNL2Discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/WN.dat", dtab]; 
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]

makeStratL2Discrepancy[nlevels_:14, ntrials_:64, nDims_:2] :=
    Module[ {},
        dtab = {};
        Do[
			npts = getCloseestNND[nDims, 2^ilevel];
			nstrats = npts^(1/nDims);
			trials = Parallelize @ Table[
				Switch[nDims
				,2, pts = N @ Table[{Mod[i,nstrats], Quotient[i,nstrats]}/nstrats + {RandomReal[],RandomReal[]}/nstrats,{i,0,npts-1}];
				,3, pts = Flatten[#,2]& @ Table[{ix+RandomReal[],iy+RandomReal[],iz+RandomReal[]}/nstrats,{ix,0,nstrats-1},{iy,0,nstrats-1},{iz,0,nstrats-1}];
				];
				{npts,getStarDiscrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
        	Print["Processing makeStratL2Discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Strat.dat", dtab]; push
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]


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

showStarDisrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"Base3SFC2D",powfromto_:{2,10},col_:Red] := 
	Module[{(*dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D*)},
		{powfrom,powto}=powfromto;
		fontSz = 20;
        dirDiscrepancy = "data_StarDiscrepancy/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);
        discrepancyBase3SFC2D = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D.dat"], thisDiscrepancy];
        plotLabel = " StarDisrepancy "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol, discrepancyBase3SFC2D},
        	{discrepancyWN, discrepancySobol, discrepancyStrat, discrepancyBase3SFC2D}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol", thisDiscrepancyLabel},
       		{"WN","Sobol", "Jitter", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
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
        Export["StarDiscrepancy_Base3SFC2D.pdf",p];
    ] (* showStarDisrepancyND *)

 showGeneralizedL2discrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"Base3SFC2D",powfromto_:{2,10},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D,kPlusMinus,data},
		{powfrom,powto}=powfromto;
    	fontSz = 20;
		kPlusMinus = .5;
		base = 2;

        dirDiscrepancy = "data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/";
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"WN.dat"]);
			discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);				
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"Strat.dat"]);
			discrepancyStrat = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];				
        discrepancyBase3SFC2D = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D.dat"], thisDiscrepancy];
        plotLabel = " GeneralizedL2discrepancyND "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol, discrepancyBase3SFC2D},
        	{discrepancyWN, discrepancySobol, discrepancyStrat, discrepancyBase3SFC2D}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol", thisDiscrepancyLabel},
       		{"WN","Sobol", "Jitter", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
        ];
        range = {Min @@ ((Last /@ #) &@discrepancyBase3SFC2D), Max @@((Last /@ #) &@discrepancyBase3SFC2D) };
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
        Export["GeneralizedL2Discrepancy_Base3SFC2D.pdf",p];
        (*P*)
    ] (* showGeneralizedL2discrepancyND *)
 
showL2discrepancyND[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"Base3SFC2D",powfromto_:{2,10},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D},
		{powfrom,powto}=powfromto;
		fontSz = 20;
        dirDiscrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);
        discrepancyOwen = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Owen.dat"]);
        discrepancyBase3SFC2D = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D.dat"], thisDiscrepancy];
        plotLabel = " L2-Discrepancy "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyBase3SFC2D},
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyStrat, discrepancyBase3SFC2D}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol","Owen", thisDiscrepancyLabel},
       		{"WN","Sobol","Owen", "Jitter", thisDiscrepancyLabel}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {Green,AbsoluteThickness[3]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} }
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
            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "L2-Discrepancy", fontSz] }
            ,ImageSize -> {1024,1024}
            ,PlotRange->{{base^powfrom,base^powto},range}
            ,GridLines->{Table[base^pow,{pow,powfrom,powto,1}],None}
            ,GridLinesStyle->Directive[Darker@Gray, Dashed]
            ,AspectRatio->1
            ,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
            ,PlotLabel -> Style[ plotLabel, Bold, 24]
        ];
        p//Print;
        Export["L2Discrepancy_Base3SFC2D.pdf",p];
        (*p*)
    ] (* showL2discrepancyND *)

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

readMatBuilderInvMatrices[fname_,nDims_:2,nlevels_:16] :=
    Module[ {data,mxsz,res,curIndex},
    	data = Import[fname];
    	curIndex = 1;
    	res = Table[{},{nlevels}];
    	Table[
     		mxsz = ilevel;
    		If[OddQ[ilevel],
    			AppendTo[res[[ilevel]], data[[curIndex;;curIndex+mxsz-1]] ];
    			curIndex += mxsz;
    			AppendTo[res[[ilevel]], data[[curIndex;;curIndex+mxsz-1]] ];
    			curIndex += mxsz;
    		,(*ELSE*)
    			AppendTo[res[[ilevel]], data[[curIndex;;curIndex+mxsz-1]] ];
    			curIndex += mxsz;
    		];
    		If[dbg,Print[ilevel -> (mf /@ res[[ilevel]]) ] ];
    	,{ilevel,nlevels}];
    	Return[res]
    ] (* readMatBuilderInvMatrices *)

(*----------------------------- Base3SFC2D --------------------------------*)

typeSq = 	1;
typeHRec = 	2;
typeVRec = 	3;


mxRot0 =   {{1,0}, {0,1}};
mxRot90 =  {{0, -1}, {1, 0}};
mxRot180 = {{-1,0}, {0,-1}};
mxRot270 = {{0, 1}, {-1, 0}};

subdivBase3SFC2DTiles[tlst_] :=
    Module[ {res={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,dxy },
    	Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2} = {v1,v2};
            Switch[tileType
              ,typeSq, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{},{2}}, {{},{1}}, {{},{0}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} }
              		];
					AppendTo[res,{typeHRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,					{v1, v2/3}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeVRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 + 1/3 v2,	{mxRot90.v1/3, mxRot90.v2},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeHRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, v2/3}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
              ,typeHRec, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{},{2}}, {{},{1}}, {{},{0}} }
              		];
					AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,				{1/3 v1, v2}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 1/3 v1 + v2,	{mxRot270.v1/3, mxRot270.v2},	{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v1,		{1/3 v1, v2}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
              ,typeVRec, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{},{2}}, {{},{1}}, {{},{0}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} }
              		];
					AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,				{v1, 1/3 v2}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 1/3 v2 + v1,	{mxRot90.v1, mxRot90.v2/3},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeSq,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, 1/3 v2}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3SFC2DTiles *)

demoBase3SFC2D[innsubdivs_:6, dbg_:False] :=
    Module[ {},
    	nsubdivs = innsubdivs;
		tlst = {{typeSq,0,{0,0}, {0,0},{{1,0},{0,1}}, {0,0},{{1,0},{0,1}}, {{},{}} ,{}} };
		Graphics[ getBase3SFC2DTilesGL[tlst,showSFC+showTilexycodes+showTileType+showArrows], PlotLabel-> 0, ImageSize -> {1024,1024} ]//Print;
		Do[
			tlst = subdivBase3SFC2DTiles @ tlst;
			flags = If[iter <= 4, showSFC+showTilexycodes+showTileType+showArrows, showSFC];
			Graphics[ getBase3SFC2DTilesGL[tlst,flags], PlotLabel-> iter, ImageSize -> {1024,1024} ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,nsubdivs}];
	] (* demoBase3SFC2D *)


getsfcBase3SFC2D[tlst_] :=
    Module[ {sfc={}, tileType,sind,refPt,v1,v2,samplingPt,prevrefPt,prevv1,prevv2,xcode,ycode,fcode,norm1,norm2,factor=4,delta=6},
    	Do[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
	    	{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/factor);
	   		AppendTo[sfc,refPt + (norm1+norm2)/1/3^((Length[fcode]+delta)/factor) ];
  			AppendTo[sfc,refPt + v1 + v2 + (-norm1-norm2)/3^((Length[fcode]+delta)/factor) ] ;
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcBase3SFC2D *)

(*getsfcBase3SFC2D[tlst_] :=
    Module[ {sfc={}, tileType,sind,refPt,v1,v2,samplingPt,prevrefPt,prevv1,prevv2,xcode,ycode,fcode},
    	Do[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			Switch[tileType
			,typeSq,
	   			AppendTo[sfc,refPt + (v1+v2)/2 ];
	   		,typeHRec,
	   			AppendTo[sfc,refPt + v1 1/6 + v2 1/2 ];
	   			AppendTo[sfc,refPt + v1 5/6 + v2 1/2 ];
	   		,typeVRec,
	   			AppendTo[sfc,refPt + v1 1/2 + v2 1/6 ];
	   			AppendTo[sfc,refPt + v1 1/2 + v2 5/6 ];
			];
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcBase3SFC2D *)
*)

getBase3SFC2DTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[10]},tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcBase3SFC2D[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,(*Arrowheads[1/3^(3+(Length[fcode]+Mod[Length[fcode],2])/2)],*)Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			fcodelen = Length[fcode];
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
    		If[BitAnd[params,showGrayValue] > 0, AppendTo[gl,{GrayLevel[FromDigits[Reverse@fcode,3]/3^fcodelen],Polygon@cont}] ];
    		If[BitAnd[params,showLightGrayTile] > 0, AppendTo[gl,{LightGray,Polygon@cont}] ];
			AppendTo[gl,Flatten[#,1]& @ {(*Point@(refPt+(norm1+norm2)/20),*)bortedStyle,Line@cont } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[fcode,3],Bold,14,Red],refPt+(v1+v2)/2,{-1.9,-1}]} ] ];
			If[BitAnd[params,showFcodeInvNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Black],refPt+(v1+v2)/2,{1.9,-1}]} ] ];
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
			If[BitAnd[params,showTilexycodes] > 0, AppendTo[gl, {Text[Style[tab2snosep@xcode,Bold,14,Red],refPt+(v1+v2)/2,{1,1}], Text[Style[tab2snosep@ycode,Bold,14,Blue],refPt+(v1+v2)/2,{-1,1}]} ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt,Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Blue], samplingPt,{-1.2,-1.2}]} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3SFC2DTilesGL *)

selectBase3SFC2DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3SFC2DTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,v,indVect,nsubdivs,m},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
     		nsubdivs = Length[xcode] + Length[ycode];
			v = Join[xcode,ycode];
			m = If[Length[xcode] == Length[ycode],
				mxInv
			,(*ELSE*)
				If[Max@(Abs@(First /@ {v1, v2})) > Max@(Abs@(Last /@ {v1, v2})), mxInvH, mxInvV]
			];
			indVect = Mod[#,3]& /@ (m.v);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC2DTiles *)


getSamplingPtsBase3SFC2DTiles[tlst_] :=
    Module[ {tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			(*{k1,k2} = {RandomReal[],RandomReal[]};*)
			samplingPt
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC2DTiles *)


getDiscrepancy2DBase3SFC2D[niters_:4] :=
    Module[ {npts,pts,dND, tlst},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
        res = Table[
			tlst = subdivBase3SFC2DTiles @ tlst;
			tlst = fillSamplingPtsBase3SFC2DTiles @ tlst;
            npts = Length[tlst];
            pts = getSamplingPtsBase3SFC2DTiles[tlst];
            dND = getGeneralizedL2discrepancy[pts];
			{npts,dND}
		,{iter,niters}];
        Export["data_discrepancyL2/2D/Base3SFC2D.dat", res]; 
        Print[mf @ res]
    ] (* getDiscrepancy2DBase3SFC2D *)


makeMatBuilderMatrices0m2net2D[] :=
    Module[ {},
    	If[ !FileExistsQ["MatBuilder_matrices/"], CreateDirectory["MatBuilder_matrices/"] ];
    	nlevels = 16;
		Do[
			execString = "testCplex -i MatBuilder_profiles/2D_0m2net.txt -o MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>".dat --seed "<>ToString[RandomInteger[2^16] ]<>" > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	Print[execString -> returnCode];
			mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>".dat"];
			Print[mf /@ mxTab];
			mxInvTab = {};
        	Do[
				If[ilevel != 0,
					mxInv = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/2, ;; ilevel]], mxTab[[2, ;; ilevel/2, ;; ilevel]]];
					Print[ilevel -> mf[mxInv] ];
					AppendTo[mxInvTab,mxInv];
        		];
				If[ilevel != nlevels, 
					mxInvH = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/2, ;; ilevel + 1]], mxTab[[2, ;; ilevel/2 + 1, ;; ilevel + 1]]] ;
					mxInvV = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/2 + 1, ;; ilevel + 1]], mxTab[[2, ;; ilevel/2, ;; ilevel + 1]]] ;
					Print[ilevel+1 -> mf[mxInvH] -> mf[mxInvV] ] ;
					AppendTo[mxInvTab,mxInvH];
					AppendTo[mxInvTab,mxInvV];
				];
        	,{ilevel,0,nlevels,2}];
			Export["MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>"_inv.dat", Flatten[#,1]& @ mxInvTab ];
		,{i,16}];
    ] (* makeMatBuilderMatrices0m2net2D *)

exportSelectionBase3SFC2D[fname_, seltlst_] :=
Module[{newtlst,tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
	newtlst = Flatten /@ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = seltlst[[ind]];			
			{tileType,sind,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2},N@refPt,N@{v1,v2},{xcode,ycode},fcode}
		,{ind,Length[seltlst]}];
	Export[fname,newtlst];
] (* exportSelectionBase3SFC2D *)

prepOptimDataBase3SFC2D[innlevels_:6, dbg_:True] :=
    Module[ {},
    	setNo = 1;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	nlevels = innlevels;
    	If[ !FileExistsQ["optim_data_2D/"], CreateDirectory["optim_data_2D/"] ];
    	If[ !FileExistsQ["optim_figs_2D/"], CreateDirectory["optim_figs_2D/"] ];
		mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat"];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		tlst = {{typeSq,0,{0,0}, {0,0},{{1,0},{0,1}}, {0,0},{{1,0},{0,1}}, {{},{}} ,{}} };
		Do[
			tlst = subdivBase3SFC2DTiles @ tlst;
			If[EvenQ[ilevel], mxInv = mxInvTab[[ilevel,1]] ];
			If[OddQ[ilevel],{mxInvH, mxInvV} = mxInvTab[[ilevel]] ];
			tlst = fillSamplingPtsBase3SFC2DTiles[tlst,mxTab,mxInv,mxInvH,mxInvV];
			(*Graphics[ {getBase3SFC2DTilesGL[tlst,showFcodeInvNumber+showTilefcode]}, PlotLabel-> nsubdivs, ImageSize -> {1024,1024} ]//Print;*)
			Do[
				seltlst = selectBase3SFC2DTiles[tlst, iOrdinalAbsolute/3^ilevel];
				fname = "optim_data_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
				exportSelectionBase3SFC2D[fname,seltlst];
				If[dbg,
				p = Graphics[ Append[background,#]& @ getBase3SFC2DTilesGL[seltlst,showLightGrayTile+showSamplingPt], PlotLabel-> iOrdinalAbsolute ];
					p//Print;
					Export["optim_figs_2D/2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
				];
			,{iOrdinalAbsolute,3^(ilevel-1)+1,3^ilevel}];
		,{ilevel,nlevels}];
	] (* prepOptimDataBase3SFC2D *)


(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeBase3SFC2DL2Discrepancy[]
*)

makeBase3SFC2DL2Discrepancy[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "optim_output_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getL2discrepancy[pts];
        	Print["Processing makeBase3SFC2DL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2Discrepancy *)


showBase3SFC2DOptimImprovement[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
         setNo = 1;
			npts = iOrdinalAbsolute = 3^6;
			inipts = Import["optim_output_2D/2D_0m2net_set_1_level_"<>ToString[iOrdinalAbsolute]<>".dat"][[;;,3;;4]];
			optimpts = Import["optim_data_2D/2D_0m2net_set_1_level_"<>ToString[iOrdinalAbsolute]<>".dat"][[;;,3;;4]];
			pairs = {inipts,optimpts}//T;
			p = Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],
				{Red,Point[#[[1]]],Blue,Point[#[[2]]],Yellow,Line[{#[[1]],#[[2]]}]}&/@pairs
				}, ImageSize->{1024,1024}, PlotLabel->ToString[npts]<>" pts   Red:MatBuilder Blue:optimized"] ;
			p//Print;
			Export["OptimImprovement.png",p];
    ] (* showBase3SFC2DOptimImprovement *)

makeBase3SFC2DGeneralizedL2Discrepancy[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "optim_output_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getGeneralizedL2discrepancy[pts];
        	Print["Processing makeBase3SFC2DGeneralizedL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab];
		p = showGeneralizedL2discrepancyND[nDims,dtab,"Base3SFC2D",{2,10}] ;
    ] (* makeBase3SFC2DGeneralizedL2Discrepancy *)

makeBase3SFC2DStarDiscrepancy[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "optim_output_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getStarDiscrepancy[pts];
        	Print["Processing makeBase3SFC2DStarDiscrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_StarDiscrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab];
		p = showStarDisrepancyND[nDims,dtab,"Base3SFC2D",{2,10}] ;
    ] (* makeBase3SFC2DGeneralizedL2Discrepancy *)


makeBase3SFC2DGeneralizedL2Discrepancy[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "optim_output_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getGeneralizedL2discrepancy[pts];
        	Print["Processing makeBase3SFC2DGeneralizedL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab];
		p = showGeneralizedL2discrepancyND[nDims,dtab,"Base3SFC2D",{2,10}] ;
    ] (* makeBase3SFC2DGeneralizedL2Discrepancy *)

(*----------------------------- end of Base3SFC2D --------------------------------*)



(*----------------------------- Base3SFC3D --------------------------------*)

typeXCube = 	1;
typeYCube = 	2;
typeZCube = 	3;
typeXflatPara = 	11;
typeYflatPara = 	12;
typeZflatPara = 	13;
typeXPara = 	21;
typeYPara = 	22;
typeZPara = 	23;

mxRotX180 = {{1,0,0}, {0,-1,0}, {0,0,-1}};
mxRotY180 = {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
mxRotZ180 = {{-1,0,0}, {0,-1,0}, {0,0,1}};
mxRotZ180 = {{-1,0,0}, {0,-1,0}, {0,0,1}};

mxRotX90 =  {{-1,0,0}, {0,0,-1}, {0,1,0}};
mxRotX270 =  {{-1,0,0}, {0,0,1}, {0,-1,0}};
mxRotY90 =  {{0,0,-1}, {0,-1,0}, {1,0,0}};
mxRotY270 =  {{0,0,1}, {0,-1,0}, {-1,0,0}};

mxRotZ270 = {{0,1,0}, {-1,0,0}, {0,0,-1}};
mxRotZ90 =  {{0,-1,0}, {1,0,0}, {0,0,-1}};

typeCube = 	1;
typeFlatPara = 	11;
typeLongPara = 	21;

subdivBase3SFC3DTiles[tlst_] :=
    Module[ {res={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2,prevv3} = {v1,v2,v3};
            Switch[tileType
              ,typeCube, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mxRotY90.v1/3, mxRotY90.v2, mxRotY90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},									{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeZflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mxRotZ90.v1/3,mxRotZ90.v2,mxRotZ90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeXPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mxRotZ270.v1/3,mxRotZ270.v2,mxRotZ270.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeYPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,	{mxRotZ90.v1,mxRotZ90.v2/3,mxRotZ90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3SFC3DTiles *)


(*subdivBase3SFC3DTiles[tlst_] :=
    Module[ {res={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2,prevv3} = {v1,v2,v3};
            Switch[tileType
              ,typeXCube, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
					AppendTo[res,{typeXflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeXflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mxRotX180.v1/3, mxRotX180.v2, mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},									{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeZCube, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, 1/3 v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
(*					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mxRotZ180.v1, mxRotZ180.v2, mxRotZ180.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeZflatPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},									{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
*)              ,typeXflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v1/3+v3,{mxRotX180.v1/3,mxRotX180.v2,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeZflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mxRotZ90.v1/3,mxRotZ90.v2,mxRotZ90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];

					(*AppendTo[res,{typeHRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,					{v1, v2/3}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeVRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 + 1/3 v2,	{mxRot90.v1/3, mxRot90.v2},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeHRec,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, v2/3}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];*)

              ,typeYflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					(*AppendTo[res,{typeYPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mxRotX180.v1,mxRotX180.v2/3,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];*)
					AppendTo[res,{typeXPara,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeXPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
					AppendTo[res,{typeXCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mxRotX180.v1/3,mxRotX180.v2,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeYPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
					AppendTo[res,{typeYCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeXCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,	{mxRotY180.v1,mxRotY180.v2/3,mxRotY180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeYCube,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3SFC3DTiles *)
*)

demoBase3SFC3D[innsubdivs_:6, dbg_:False] :=
    Module[ {},
    	nsubdivs = innsubdivs;
		tlst = {{typeCube,0,{0,0,0}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {{},{},{}} ,{}} };
		Graphics3D[ getBase3SFC3DTilesGL[tlst,showSFC+showArrows+showTilexycodes+showTileType], PlotLabel-> iter, ImageSize -> {1024,1024}3/2  ]//Print;
		If[dbg, tlst//mf//Print];
		Do[
			tlst = subdivBase3SFC3DTiles @ tlst;
			flags = showSFC+showArrows+showTilexycodes+showTileType+showArrows;
			flags = If[iter <= 5, showSFC+showArrows+showTileType, showSFC];
			Graphics3D[ getBase3SFC3DTilesGL[tlst,flags], PlotLabel-> iter, ImageSize -> {1024,1024} 3/2 ]//Print;
			If[dbg, tlst//mf//Print];
		,{iter,nsubdivs}];
	] (* demoBase3SFC3D *)

getsfcBase3SFC3D[tlst_] :=
    Module[ {sfc={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,norm1,norm2,norm3,k=6,delta=12},
    	Do[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
	    	{norm1,norm2,norm3}={v1/Norm[v1],v2/Norm[v2],v3/Norm[v3]}/3^((Length[fcode]+Mod[Length[fcode],3])/k );
	   		AppendTo[sfc,refPt + (norm1+norm2+norm3)/3^((Length[fcode]+delta)/k  ) ];
  			AppendTo[sfc,refPt + v1 + v2 + v3 - (norm1+norm2+norm3)/3^((Length[fcode]+delta)/k ) ] ;
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcBase3SFC3D *)


getBase3SFC3DTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[10]},tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcBase3SFC3D[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			fcodelen = Length[fcode];
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2, refPt+v2+v3,refPt+v3+v1+v2,refPt+v3+v1,refPt+v3,refPt,refPt+v3,refPt+v3+v2,refPt+v2,refPt,refPt+v1,refPt+v3+v1,refPt+v3+v1+v2,refPt+v1+v2};
			AppendTo[gl,Flatten[#,1]& @ {bortedStyle,Line@cont } ];
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2+v3)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[fcode,3],Bold,14,Red],refPt+(v1+v2+v3)/2,{-1.9,-1}]} ] ];
			If[BitAnd[params,showFcodeInvNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Black],refPt+(v1+v2+v3)/2,{1.9,-1}]} ] ];
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2+v3)/2,{0,1}]} ] ];
			If[BitAnd[params,showTilexycodes] > 0, AppendTo[gl, {
				Text[Style[tab2snosep@xcode ,Bold,14,Black]	,refPt+(v1+v2+v3)/2,{3,-1}],
				Text[Style[tab2snosep@ycode ,Bold,14,Red]	,refPt+(v1+v2+v3)/2,{0,-1}],Text[Style[tab2snosep@zcode ,Bold,14,Blue]	,refPt+(v1+v2+v3)/2,{-3,-1}] 
			}  ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt,Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Blue], samplingPt,{-1.2,-1.2}]} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3SFC3DTilesGL *)

selectBase3SFC3DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3SFC3DTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,v,indVect,nsubdivs,m},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
     		nsubdivs = Length[xcode] + Length[ycode] + Length[zcode];
			v = Join[xcode,ycode,zcode];
			m = If[Length[xcode] == Length[ycode],
				mxInv
			,(*ELSE*)
				If[Max@(Abs@(First /@ {v1, v2, v3})) > Max@(Abs@(Last /@ {v1, v2})), mxInvH, mxInvV]
			];
			indVect = Mod[#,3]& /@ (m.v);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC3DTiles *)


getSamplingPtsBase3SFC3DTiles[tlst_] :=
    Module[ {tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode},
    	Parallelize @ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			(*{k1,k2} = {RandomReal[],RandomReal[]};*)
			samplingPt
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC3DTiles *)


getDiscrepancy2DBase3SFC3D[niters_:4] :=
    Module[ {npts,pts,dND, tlst},
		tlst = tlst = {{typeCube,0,{0,0,0}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {{},{},{}} ,{}} };
        res = Table[
			tlst = subdivBase3SFC3DTiles @ tlst;
			tlst = fillSamplingPtsBase3SFC3DTiles @ tlst;
            npts = Length[tlst];
            pts = getSamplingPtsBase3SFC3DTiles[tlst];
            dND = getGeneralizedL2discrepancy[pts];
			{npts,dND}
		,{iter,niters}];
        Export["data_discrepancyL2/2D/Base3SFC3D.dat", res]; 
        Print[mf @ res]
    ] (* getDiscrepancy2DBase3SFC3D *)

exportSelectionBase3SFC3D[fname_, seltlst_] :=
Module[{newtlst,tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode},
	newtlst = Flatten /@ Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = seltlst[[ind]];			
			{tileType,sind,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2,prevv3},N@refPt,N@{v1,v2,v3},{xcode,ycode,zcode},fcode}
		,{ind,Length[seltlst]}];
	Export[fname,newtlst];
] (* exportSelectionBase3SFC3D *)

prepOptimDataBase3SFC3D[innlevels_:2, dbg_:False] :=
    Module[ {},
    	setNo = 1;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	nlevels = innlevels;
    	If[ !FileExistsQ["optim_data_3D/"], CreateDirectory["optim_data_3D/"] ];
    	If[dbg && !FileExistsQ["optim_figs_3D/"], CreateDirectory["optim_figs_3D/"] ];
		(*mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat"];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];*)
		tlst = {{typeCube,0,{0,0,0}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {0,0,0},{{1,0,0},{0,1,0},{0,0,1}}, {{},{},{}} ,{}} };
		Do[
			tlst = subdivBase3SFC3DTiles @ tlst;
			If[EvenQ[ilevel], mxInv = mxInvTab[[ilevel,1]] ];
			If[OddQ[ilevel],{mxInvH, mxInvV} = mxInvTab[[ilevel]] ];
			tlst = fillSamplingPtsBase3SFC3DTiles[tlst,mxTab,mxInv,mxInvH,mxInvV];
			Do[tlst[[i,3]] = {Random[],Random[],Random[]},{i,Length[tlst]}];
			If[dbg, Graphics3D[ {getBase3SFC3DTilesGL[tlst,showFcodeInvNumber+showTilefcode]}, PlotLabel-> nsubdivs, ImageSize -> {1024,1024} ]//Print];
			Do[
				seltlst = selectBase3SFC3DTiles[tlst, ii/3^ilevel];
				fname = "optim_data_3D/2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[ii]<>".dat";
				exportSelectionBase3SFC3D[fname,seltlst];
			,{ii,3^(ilevel-1)+1,3^ilevel}];
		,{ilevel,nlevels}];
	] (* prepOptimDataBase3SFC3D *)

makeMatBuilderMatrices0m2net3D[generateMatBuilderMatrix_:False] :=
    Module[ {},
    	If[ !FileExistsQ["MatBuilder_matrices/"], CreateDirectory["MatBuilder_matrices/"] ];
    	nlevels = 16;
		Do[
			If[generateMatBuilderMatrix,
				execString = "testCplex -i MatBuilder_profiles/3D_0m2net.txt -o MatBuilder_matrices/3D_0m2net_"<>i2s[i]<>".dat --seed "<>ToString[RandomInteger[2^16] ]<>" > /dev/null";
    	    	returnCode = Run[execPrefix<>execString];
        		Print[execString -> returnCode];
			];
			mxTab = readMatBuilderMatrix["MatBuilder_matrices/3D_0m2net_"<>i2s[i]<>".dat",3];
			Print[mf /@ mxTab];
			mxInvTab = {};
        	Do[
				If[ilevel != 0,
					mxInv = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/3, ;; ilevel]], mxTab[[2, ;; ilevel/3, ;; ilevel]], mxTab[[3, ;; ilevel/3, ;; ilevel]]];
					Print[ilevel -> mf[mxInv] ];
					AppendTo[mxInvTab,mxInv];
        		];
				If[ilevel != nlevels, 
					mxInvH = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/2, ;; ilevel + 1]], mxTab[[2, ;; ilevel/2 + 1, ;; ilevel + 1]]] ;
					mxInvV = Inverse[#,Modulus->3]& @ Join[mxTab[[1, ;; ilevel/2 + 1, ;; ilevel + 1]], mxTab[[2, ;; ilevel/2, ;; ilevel + 1]]] ;
					Print[ilevel+1 -> mf[mxInvH] -> mf[mxInvV] ] ;
					AppendTo[mxInvTab,mxInvH];
					AppendTo[mxInvTab,mxInvV];
				];
        	,{ilevel,0,nlevels,3}];
			Export["MatBuilder_matrices/3D_0m2net_"<>i2s[i]<>"_inv.dat", Flatten[#,1]& @ mxInvTab ];
		,{i,1}];
    ] (* makeMatBuilderMatrices0m2net3D *)

(*----------------------------- end of Base3SFC3D --------------------------------*)
