(* TileBasedOptim.m
   2022-04 vo, based on fibo-hilbert.m (version 2002/12/14)
   
makeWNL2Discrepancy[]
makeStratL2Discrepancy[]
makeSobolL2Discrepancy[]
makeOwenL2Discrepancy[]

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

absoluteTime = Floor[1000000*AbsoluteTime[]];

first := First
second[x_]:= If[Length[x] > 1, x[[2]], First[x] ] (* like First *)
third[x_]:= If[Length[x] > 2, x[[3]], First[x] ] (* like First *)
fourth[x_]:= If[Length[x] > 3, x[[4]], First[x] ] (* like First *)
last := Last
slightlyDarker[color_] := Darker[color,1/5]

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
showFcodeInvNumber = 1024;
showBasicVectors = 2048;

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

(* <<<<<<<<<<<<<<<<<<<<<< L2Discrepancy

gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeWNL2Discrepancy[]
makeStratL2Discrepancy[]
makeSobolL2Discrepancy[]
makeOwenL2Discrepancy[]

*)

makeSobolL2Discrepancy[nlevels_:14, nDims_:3,dbg_:False] :=
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

makeOwenL2Discrepancy[nlevels_:14, nDims_:3,dbg_:False] :=
    Module[ {},
		If[ Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2]];
        dtab = {};
        nptsMax = 2^nlevels;
        Do[
			npts = inpts;
        	d = Mean @ (Parallelize @ Table[
				execString = "owen -n "<>ToString[npts]<>" --nd "<>ToString[nDims]<>" -p 1 --max_tree_depth_32_flag 0 -s "<>ToString[RandomInteger[2^16]]<>" -o tmp/pts"<>pid<>".dat > /dev/null";
        		returnCode = Run[execPrefix<>execString];
        		pts = Import["tmp/pts"<>pid<>".dat"];		
        		getL2discrepancy[pts]
        	,{64}]);
        	Print["Processing makeOwenL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Owen.dat", dtab]; 
        ,{inpts,nptsMax}];
        Print[mf @ dtab]
    ] (* makeOwenL2Discrepancy *)

makeWNL2Discrepancy[nlevels_:14, ntrials_:64, nDims_:3] :=
    Module[ {},
        dtab = {};
        Do[
			npts = 2^ilevel;
			trials = Parallelize @ Table[
				pts = Table[Table[RandomReal[],{nDims}],{i,npts}];
				{npts,getL2discrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
        	Print["Processing makeWNL2Discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/WN.dat", dtab]; 
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]


makeStratL2Discrepancy[nlevels_:14, ntrials_:64, nDims_:3] :=
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
				{npts,getL2discrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
        	Print["Processing makeStratL2Discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Strat.dat", dtab]; push
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]
(*
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
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D,discrepancyOwen},
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

showL2discrepancyNDMatBuilderOnly[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"Base3SFC2D (ref)",powfromto_:{2,10},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D,discrepancyOwen},
		{powfrom,powto}=powfromto;
		fontSz = 20;
        dirDiscrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);
        discrepancyOwen = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Owen.dat"]);
        discrepancyBase3SFC2DRef = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D.dat"], thisDiscrepancy];
        discrepancyBase3SFC2D = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D_MatBuilderOnly.dat"], thisDiscrepancy];
        plotLabel = " L2-Discrepancy "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyBase3SFC2DRef, discrepancyBase3SFC2D},
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyStrat, discrepancyBase3SFC2DRef, discrepancyBase3SFC2D}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol","Owen", thisDiscrepancyLabel,"Base3SFC2D + MatBuilderOnly"},
       		{"WN","Sobol","Owen", "Jitter", thisDiscrepancyLabel,"Base3SFC2D + MatBuilderOnly"}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {Green,AbsoluteThickness[3]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]}, {Black,AbsoluteThickness[3]} }
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
        Export["L2Discrepancy_Base3SFC2D_MatBuilderOnly.dat.pdf",p];
        (*p*)
    ] (* showL2discrepancyNDMatBuilderOnly *)


showL2discrepancyNDExperiment2Tanguy[nDims_:2,thisDiscrepancy_:{},thisDiscrepancyLabel_:"Base3SFC2D (based on MatBuilder)",powfromto_:{2,10},col_:Red] := 
	Module[{dirDiscrepancy,discrepancyWN,discrepancyStrat,discrepancySobol,plotLabel,legends,alldata,p,powfrom,powto,fontSz,range,colors,base=2,discrepancyBase3SFC2D,discrepancyOwen},
		{powfrom,powto}=powfromto;
		fontSz = 20;
        dirDiscrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
        discrepancyWN = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"WN.dat"]);
        discrepancyStrat = Select[#,1<=#[[1]]<=base^(powto+1)&]& @ If[FileExistsQ[dirDiscrepancy<>"Strat.dat"], (Import[dirDiscrepancy<>"Strat.dat"]), {discrepancyWN[[1]]} ];
        discrepancySobol = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Sobol.dat"]);
        discrepancyOwen = Select[#,base^powfrom<=#[[1]]<=base^(powto+1)&]& @ (Import[dirDiscrepancy<>"Owen.dat"]);
        discrepancyBase3SFC2DRef = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D.dat"], thisDiscrepancy];
        discrepancyBase3SFC2D = If[thisDiscrepancy==={}, Import[dirDiscrepancy<>"Base3SFC2D_Experiment2Tanguy.dat"], thisDiscrepancy];
        plotLabel = " L2-Discrepancy "<>ToString[nDims]<>"D";               
        alldata = If[Length[discrepancyStrat] == 1,
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyBase3SFC2DRef, discrepancyBase3SFC2D},
        	{discrepancyWN, discrepancySobol,discrepancyOwen, discrepancyStrat, discrepancyBase3SFC2DRef, discrepancyBase3SFC2D}
        ];
        legends = If[Length[discrepancyStrat] == 1,
        	{"WN","Sobol","Owen", thisDiscrepancyLabel,"Base3SFC2D + random init"},
       		{"WN","Sobol","Owen", "Jitter", thisDiscrepancyLabel,"Base3SFC2D + random init"}
        ];
        colors = If[Length[discrepancyStrat] == 1,
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]} },
        	{ {Red,Dotted,AbsoluteThickness[10]}, {Gray,AbsoluteThickness[3]}, {Green,AbsoluteThickness[3]}, {Blue,Dotted,AbsoluteThickness[10]}, {col,AbsoluteThickness[3]}, {Black,AbsoluteThickness[3]} }
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
        Export["L2Discrepancy_Base3SFC2D_Experiment2Tanguy.dat.pdf",p];
        (*p*)
    ] (* showL2discrepancyNDExperiment2Tanguy *)



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
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt,Text[Style[sind (*FromDigits[Reverse@fcode,3]*),Bold,10,Blue], samplingPt,{-1.1,-1.1}]} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3SFC2DTilesGL *)

selectBase3SFC2DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3SFC2DTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,v,indVect,nsubdivs,m,matBuilderIndex},
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
			matBuilderIndex = FromDigits[#,3]& @ (Reverse @ indVect);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode}
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


(*getDiscrepancy2DBase3SFC2D[niters_:4] :=
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
*)

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

selectBase3SFC2DTilesMatBuilderOnly[tlst_,intensityInt_] := Select[tlst, second[#] < intensityInt & ]


(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
prepOptimDataBase3SFCMatBuilderOnly2D[6]
*)

prepOptimDataBase3SFCMatBuilderOnly2D[innlevels_:6, dbg_:True] :=
    Module[ {},
    	setNo = 1;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	nlevels = innlevels;
    	If[ !FileExistsQ["optim_data_2D_MatBuilderOnly/"], CreateDirectory["optim_data_2D_MatBuilderOnly/"] ];
    	If[ !FileExistsQ["optim_figs_2D_MatBuilderOnly/"], CreateDirectory["optim_figs_2D_MatBuilderOnly/"] ];
		mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat"];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		tlst = {{typeSq,0,{0,0}, {0,0},{{1,0},{0,1}}, {0,0},{{1,0},{0,1}}, {{},{}} ,{}} };
		Do[
			tlst = subdivBase3SFC2DTiles @ tlst;
			If[EvenQ[ilevel], mxInv = mxInvTab[[ilevel,1]] ];
			If[OddQ[ilevel],{mxInvH, mxInvV} = mxInvTab[[ilevel]] ];
			tlst = fillSamplingPtsBase3SFC2DTiles[tlst,mxTab,mxInv,mxInvH,mxInvV];
			(*Graphics[ {getBase3SFC2DTilesGL[tlst,showFcodeInvNumber+showTilefcode]}, PlotLabel-> nsubdivs, ImageSize -> {1024,1024} ]//Print;*)
			Parallelize @ Do[
				seltlst = selectBase3SFC2DTilesMatBuilderOnly[tlst, iOrdinalAbsolute ];
				fname = "optim_data_2D_MatBuilderOnly/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
				exportSelectionBase3SFC2D[fname,seltlst];
				Print[Length[seltlst] -> " Exporting " -> fname];
				If[dbg,
				p = Graphics[ Append[background,#]& @ getBase3SFC2DTilesGL[seltlst,showLightGrayTile+showSamplingPt], PlotLabel-> iOrdinalAbsolute ];
					p//Print;
					Export["optim_figs_2D_MatBuilderOnly/2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
				];
			,{iOrdinalAbsolute,3^(ilevel-1)+1,3^ilevel}];
		,{ilevel,nlevels}];
	] (* prepOptimDataBase3SFC2D *)

prepOptimDataBase3SFC2DExperiment1Tanguy[innlevels_:6, dbg_:True] :=
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
			Do[{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[i]];
				prevrefPt = refPt;
				{prevv1,prevv2} = {v1,v2};
				tlst[[i]] = {tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode};
			,{i,Length[tlst]}];
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
	] (* prepOptimDataBase3SFC2DExperiment1Tanguy *)

prepOptimDataBase3SFC2DExperiment2Tanguy[innlevels_:6, dbg_:True] :=
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
			Do[{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[i]];
				samplingPt = prevrefPt + RandomReal[] prevv1 + RandomReal[] prevv2;
				tlst[[i]] = {tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode};
				,{i,Length[tlst]}];

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
	] (* prepOptimDataBase3SFC2DExperiment2Tanguy *)


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

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeBase3SFC2DL2DiscrepancyExperimentTanguy[]
*)
makeBase3SFC2DL2DiscrepancyExperimentTanguy[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "~tanguy/optim_output_2D_1K_Exp2/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>"Experiment2.dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getL2discrepancy[pts];
        	Print["Processing makeBase3SFC2DL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D_Experiment2Tanguy.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2DiscrepancyExperimentTanguy *)

makeBase3SFC2DL2DiscrepancyMatBuilderOnly[dbg_:False] :=
    Module[ {},
    	nDims = 2;
        dtab = {};
        nptsMax = 3^6;
        setNo = 1;
        
        Do[
			npts = iOrdinalAbsolute;
			fname = "optim_output_2D_MatBuilderOnly/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
			pts = Import[fname][[;;,3;;4]];
			If[dbg, ipts = Round[ npts pts ];
				Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
        	d = getL2discrepancy[pts];
        	Print["Processing makeBase3SFC2DL2Discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D_MatBuilderOnly.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2DiscrepancyMatBuilderOnly *)



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

rotmx3D[theta_, {x_,y_,z_}] := Module[{ux,uy,uz},{ux,uy,uz}={x,y,z}/Norm[{x,y,z}]; Cos[theta] IdentityMatrix[3] + Sin[theta] {{0, -uz, uy}, {uz, 0, -ux}, {-uy, ux, 0}} + (1-Cos[theta]) {{ux^2, ux uy, ux uz},{ux uy, uy^2, uy uz},{ux uz, uy uz, uz^2}} ]

permut12[{v1_,v2_,v3_}] := {v2,v1,v3}
permut23[{v1_,v2_,v3_}] := {v1,v3,v2}
permut13[{v1_,v2_,v3_}] := {v3,v2,v1}

typeCubeRight = 1;	(* Right-hand Coordinate Systems, XYZ *)
typeCubeLeft = 	2;	(* Left-hand Coordinate Systems, YZZ *)

typeParaXflatRight = 	11;
typeParaXflatLeft = 	12;
typeParaYflatRight = 	13;
typeParaYflatLeft = 	14;
typeParaZflatRight = 	15;
typeParaZflatLeft = 	16;

typeParaXlongRight = 	21;
typeParaXlongLeft = 	22;
typeParaYlongRight = 	23;
typeParaYlongLeft = 	24;
typeParaZlongRight = 	25;
typeParaZlongLeft = 	26;

subdivBase3SFC3DTiles[tlst_] :=
    Module[ {res={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2,prevv3} = {v1,v2,v3};
            Switch[tileType
              ,typeCubeRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaZflatRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaZflatLeft, sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	permut12 @ {mx.v1, mx.v2, mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaZflatRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeCubeLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXflatLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 					{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx =  rotmx3D[Pi, v1];
					AppendTo[res,{typeParaXflatRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v3+v1/3,	permut23 @ {mx.v1/3, mx.v2, mx.v3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXflatLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXflatLeft,
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaYlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 						{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeParaYlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3+v1/3+v2,(permut23 @ {mx.v1/3,mx.v2,mx.v3}),{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaYlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3}, 				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaXlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,permut12 @ {mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 		{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatLeft, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,permut23 @ {mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2, v3/3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3/3+v1+v2,permut12 @ {mx.v1,mx.v2,mx.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v3,		{v1, v2, v3/3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,permut23 @ {mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3SFC3DTiles *)

debugBase3SFC3D[innsubdivs_:1, dbg_:False] :=
    Module[ {},
    	nsubdivs = innsubdivs;
    	
    	flags = showSFC+showArrows+showTileType+showBasicVectors;
		tlst = {{typeCubeRight,0,{0,0,0}, {},{{},{},{}}, 		{0,0,0},	{{1,0,0},{0,1,0},{0,0,1}}, {{},{},{}},{}},
				{typeCubeLeft,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{0,1,0},{1,0,0},{0,0,1}}, {{},{},{}},{}},
				{typeParaXlongRight,0,{0,0,0}, {},{{},{},{}}, 		{0,0,0},	{{1,0,0},{0,1,0}/3,{0,0,1}/3}, {{},{},{}},{}},
				{typeParaZlongRight,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{1,0,0}/3,{0,1,0}/3,{0,0,1}}, {{},{},{}},{}},
				{typeParaYlongLeft,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{0,1,0}/3,{1,0,0},{0,0,1}/3}, {{},{},{}},{}},
				{typeParaXlongLeft,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{0,1,0},{1,0,0}/3,{0,0,1}/3}, {{},{},{}},{}},
				{typeParaXflatLeft,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{0,1,0},{1,0,0},{0,0,1}/3}, {{},{},{}},{}},
				{typeParaYflatRight,0,{0,0,0}, {},{{},{},{}},		{0,0,0},	{{1,0,0},{0,1,0}/3,{0,0,1}}, {{},{},{}},{}},
				{typeParaZflatRight,0,{0,0,0}, {},{{},{},{}}, 		{0,0,0},	{{1,0,0},{0,1,0},{0,0,1}/3}, {{},{},{}},{}},
				{typeParaZflatLeft,0,{0,0,0}, {},{{},{},{}}, 		{0,0,0},	{{0,1,0},{1,0,0},{0,0,1}/3}, {{},{},{}},{}} };

		Table[
				Graphics3D[ getBase3SFC3DTilesGL[{tlst[[idemo]]}, flags], PlotLabel-> 0, ImageSize -> {1024,1024}/2 , PlotRange->{{-.25,1.25},{-.25,1.25},{-.25,1.25}} ]
		,{idemo,Length[tlst]}]//Print;		
		Do[
			flags = If[iter <= 5, showSFC+showArrows+showTileType+showBasicVectors, showSFC];
			Table[
				tlst[[idemo]] = subdivBase3SFC3DTiles @ {tlst[[idemo]]};
				Graphics3D[ getBase3SFC3DTilesGL[tlst[[idemo]],flags], PlotLabel-> iter, ImageSize -> {1024,1024}/2 , PlotRange->{{-.25,1.25},{-.25,1.25},{-.25,1.25}} ]
			,{idemo,Length[tlst]}]//Print;
		,{iter,nsubdivs}];
	] (* debugBase3SFC3D *)

demoBase3SFC3D[innsubdivs_:1, dbg_:False] :=
    Module[ {},
    	nsubdivs = innsubdivs;
    	
    	flags = showSFC+showArrows+showTileType+showOrdinalNumber+showBasicVectors;
		tlst = {{typeCubeRight,0,{0,0,0}, {},{{},{},{}}, 		{0,0,0},	{{1,0,0},{0,1,0},{0,0,1}}, {{},{},{}},{}} };
		Graphics3D[ getBase3SFC3DTilesGL[tlst, flags], PlotLabel-> 0, ImageSize -> {1024,1024} , PlotRange->{{0,1},{0,1},{0,1}} ] //Print;		
		Do[
			flags = If[iter <= 5, showSFC+showArrows+showTileType+showOrdinalNumber, showSFC];
			tlst = subdivBase3SFC3DTiles @ tlst;
			Graphics3D[ getBase3SFC3DTilesGL[tlst, flags], PlotLabel-> iter, ImageSize -> {1024,1024} , PlotRange->{{0,1},{0,1},{0,1}} ] //Print;	
		,{iter,nsubdivs}];
	] (* demoBase3SFC3D *)

getBase3SFC3DTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[10]},tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={GrayLevel[.6],AbsoluteThickness[5]}},
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
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[ind,Bold,14,Red],refPt+(v1+v2+v3)/2,{-1.9,-1}]} ] ];
			If[BitAnd[params,showFcodeInvNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Black],refPt+(v1+v2+v3)/2,{1.9,-1}]} ] ];
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2+v3)/2,{0,1}]} ] ];
			If[BitAnd[params,showTilexycodes] > 0, AppendTo[gl, {
				Text[Style[tab2snosep@xcode ,Bold,14,Black]	,refPt+(v1+v2+v3)/2,{3,-1}],
				Text[Style[tab2snosep@ycode ,Bold,14,Red]	,refPt+(v1+v2+v3)/2,{0,-1}],Text[Style[tab2snosep@zcode ,Bold,14,Blue]	,refPt+(v1+v2+v3)/2,{-3,-1}] 
			}  ] ];
			If[BitAnd[params,showBasicVectors] > 0, AppendTo[gl, {AbsoluteThickness[3],Red,Arrow[{refPt,refPt+v1}],Green,Arrow[{refPt,refPt+v2}],Blue,Arrow[{refPt,refPt+v3}]} ]];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt,Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Blue], samplingPt,{-1.2,-1.2}]} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3SFC3DTilesGL *)

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
              (*,typeCubeYXZ, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeZflatParaX,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 									{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeZflatParaY,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mxRotY90.v1, mxRotY90.v2/3, mxRotY90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeZflatParaX,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
*)

(*
subdivBase3SFC3DTiles[tlst_] :=
    Module[ {res={}, tileType,sind,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2,prevv3} = {v1,v2,v3};
            Switch[tileType
              ,typeCubeRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaZflatRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaZflatRight, sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mx.v1, mx.v2, mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaZflatRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeCubeLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXflatLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 					{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeParaXflatLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v3+v1/3,	{mx.v1/3, mx.v2, mx.v3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXflatLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXflatLeft,
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaYlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 						{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaYlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,{1,1,1}{mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaYlongLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,{mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 		{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2, v3/3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3/3+v1+v2,	{mx.v1,mx.v2,mx.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v3,		{v1, v2, v3/3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,	{mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,sind,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3SFC3DTiles *)
*)


(*----------------------------------------- MSE *)
(*


gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nintegrands = 16 1024;
nDims = 2;

nPointsets = 2;
integrandType = 1;
makeMSEref[1, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];
integrandType = 2;
makeMSEref[1, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];

nPointsets = 64;
integrandType = 1;
makeMSEref[19, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];
integrandType = 2;
makeMSEref[19, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];

nPointsets = 64;
integrandType = 1;
makeMSEref[10, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];
integrandType = 2;
makeMSEref[10, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];

nPointsets = 64;
integrandType = 1;
makeMSEref[11, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];
integrandType = 2;
makeMSEref[11, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];

nPointsets = 32;
integrandType = 1;
makeMSEref[208, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];
integrandType = 2;
makeMSEref[208, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];

*)
makeMSEref[inpointsetTypes_:10, nTrialsMSE_:1024, powParams_:{2,18,1}, inIntegrandType_:1, innDims_:2, nIntegrands_:1024, dbg_:False] :=
    Module[ {},
    	firstDim = 0;
    	(*If[ Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2]];*)
    	nDims = innDims;
    	integrandType = inIntegrandType;
		integrandTypeLabel = Switch[integrandType,  1,"Ellipses", 2,"SoftEllipses", 3,"Rectangles", 4,"SoftRectangles" ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#Analytical	#MSE	#NbPtsets	#Nbintegrands\n";
		fnameLabel = integrandTypeLabel ;
		nPointsets 	= nTrialsMSE ;
        {powfrom,powto,powstep} = powParams;

		dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>fnameLabel<>"/";
        If[ !FileExistsQ[dirMSE], CreateDirectory[dirMSE] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
		pointsetType = inpointsetTypes;
		pointsetLabel = Switch[pointsetType 
			(* sequence, ND *)   ,1,"Sobol" ,2,"Halton" ,3,"Faure" ,4,"Niederreiter" ,5,"SobolPlusPlus",6,"SobolGlobal",7,"OwenGlobal"
			(* pointsets, ND *)  ,10,"WN" ,11,"Strat" ,12,"OwenPlus" ,13,"Rank1Lattice" ,14,"DartThrowing"  ,15,"RegGrid" ,16,"SOT" ,17,"SOTPlusLloyd",18,"OwenPureFirstDim2",19,"OwenPure" 
			(* ExtensibleLattices, ND *),20,"ExtensibleLatticeType0" ,21,"ExtensibleLatticeType1" ,22,"ExtensibleLatticeType2" ,23,"ExtensibleLatticeType3"
			(* pointsets, 2D only *) ,200,"HexGrid"  ,201,"Hammersley" ,202,"LarcherPillichshammer" ,203,"NRooks" ,204,"BNOT" ,205,"CMJ" ,206,"BNLDS" 
									 ,207,"PMJ" ,208,"PMJ02" ,209,"LDBN" ,210,"Penrose" ,211,"Fattal" ,212, "HexGridTorApprox"
	    	(* pointsets, 3D only *) ,300,"BCC", 301,"FCC", 302,"HCP", 303,"WeairePhelan"
    		(* 4D Variants of Sobol,OwenPlus,OwenPure *) ,400,"Sobol2356" ,401,"OwenPlus2356" ,402,"OwenPure2356",403,"Sobol2367" ,404,"OwenPlus2367" ,405,"OwenPure2367"
			(* uniformND *) ,500,"UniformND" ,501,"UniformNDwithoutSobol"
			(* uniformND *) ,600,"zsampler",601,"morton",602,"morton01"
			(* pointsets, SobolShiftedKx *) ,701,"SobolShifted1x",702,"SobolShifted2x",703,"SobolShifted3x",704,"SobolShifted4x"
			(* pointsets, OwenPlusShiftedKx *) ,801,"OwenPlusShifted1x",802,"OwenPlusShifted2x",803,"OwenPlusShifted3x",804,"OwenPlusShifted4x"
			(* pointsets, OwenShiftedKx *) ,900,"OwenMicroShift",999,"OwenMicroShiftGlobal",901,"OwenShifted1x",902,"OwenShifted2x",903,"OwenShifted3x",904,"OwenShifted4x"
	    	,_, "unknown" 
		];
    	If[pointsetLabel == "SOT", {powfrom,powto,powstep} = Switch[nDims,2,{2,17,1},3,{2,16,1},4,{2,17,1}] ];
		Print[pointsetLabel,{powfrom,powto,powstep} -> " makeMSEref from ",2^powfrom," to ",2^powto];
		dataMSE = {};
		Do[	
     		If[pointsetLabel == "SOT" && nDims == 4 && iptsPow == 17, nPointsets = 1 ]; (* Only 1 available *)
   			npts = getRealNPts[nDims, 2^iptsPow, pointsetType];
    		resFname = pointsetLabel<>"_"<>fnameLabel<>".dat";
			mseTab = ( Parallelize @  Table[   				
				ptsfname = "tmp/pts_"<>ToString[iPointSet]<>".dat";
				msefname = "tmp/mse"<>pid<>".dat";
				Switch[pointsetLabel					
				,"UniformND",
					getUniformND[nDims, npts, ptsfname];
				,"OwenPureFirstDim2", execString = "owen --firstDim 2 --nDims "<>ToString[nDims]<>" -p 1 -t 32 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"HexGrid", 
					pts = makeHexGridPts[npts];
		     		Export[ptsfname,pts]
				,"SOT", navailable = Switch[nDims,2,64,3,64,4,16,_,64];
					ptsfname = "_pointsets_SobolPlusPlus/"<>ToString[nDims]<>"D/SOT/SOT_"<>i2s[npts,8]<>"/SOT_"<>i2s[npts,8]<>"_"<>i2s[Mod[iPointSet,navailable,1],6]<>".dat";
				,"Rank1Lattice", 
					dir="_pointsets_SobolPlusPlus/"<>ToString[nDims]<>"D/"<>pointsetLabel<>"/"<>pointsetLabel<>"_"<>i2s[npts,8]<>"/";
					ptsfname = dir<>pointsetLabel<>"_"<>i2s[npts,8]<>"_"<>i2s[iPointSet,6]<>".dat";
				,"PMJ02", 
					ptsfname = "pointsets/"<>ToString[nDims]<>"D/pmj02/pmj02_set"<>i2s[iPointSet-1,5]<>"_"<>i2s[npts,8]<>".dat";
				,"Sobol", execString = "owen --nDims "<>ToString[nDims]<>" -o "<>ptsfname<>" -n "<>ToString[npts]<>" --permut 0  > /dev/null";
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"OwenPure", execString = "owen --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   					res = Run[execPrefix<>execString];
			     	If[dbg, Print[execString -> res] ];
					,"SobolShifted1x", execString = "owen --start "<>ToString[1*npts]<>" --nDims "<>ToString[nDims]<>" -p 0 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"SobolShifted2x", execString = "owen --start "<>ToString[2*npts]<>" --nDims "<>ToString[nDims]<>" -p 0 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"SobolShifted3x", execString = "owen --start "<>ToString[3*npts]<>" --nDims "<>ToString[nDims]<>" -p 0 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"SobolShifted4x", execString = "owen --start "<>ToString[4*npts]<>" --nDims "<>ToString[nDims]<>" -p 0 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
(*>>>>>>>>>>>>>*)   ,"OwenMicroShift", execString = "owen -f "<>ToString[firstDim]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
			     		data = If[2/npts  <= #[[1]] < 1-2/npts && 2/npts  <= #[[2]] < 1-2/npts, Plus[#, {RandomReal[],RandomReal[]}/npts ], #] & /@ Import[ptsfname];
			     		Export[ptsfname, data];
(*>>>>>>>>>>>>>*)   ,"OwenMicroShiftGlobal", execString = "owen -f "<>ToString[firstDim]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
			     		delta = 1/npts Table[RandomReal[], {nDims}];
			     		data = Plus[#, delta] & /@ Import[ptsfname];
			     		Export[ptsfname, data];
					,"OwenShifted1x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[1*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenShifted2x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[2*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenShifted3x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[3*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenShifted4x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[4*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenPlus", execString = "owen --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenPlusShifted1x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[1*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenPlusShifted2x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[2*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenPlusShifted3x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[3*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
					,"OwenPlusShifted4x", execString = "owen -f "<>ToString[firstDim]<>" --start "<>ToString[4*npts]<>" --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"Sobol2356", execString = "owen --nDims 4 -d data/sobol_init_2356.dat -f 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" --permut 0  > /dev/null";
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"OwenPlus2356", execString = "owen --nDims 4 -d data/sobol_init_2356.dat -f 1 -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"OwenPure2356", execString = "owen --nDims 4 -d data/sobol_init_2356.dat -f 1 -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"Sobol2367", execString = "owen --nDims 4 -d data/sobol_init_2367.dat -f 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" --permut 0  > /dev/null";
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"OwenPlus2367", execString = "owen --nDims 4 -d data/sobol_init_2367.dat -f 1 -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"OwenPure2367", execString = "owen --nDims 4 -d data/sobol_init_2367.dat -f 1 -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 0 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"WN", 
					pts = getWN[nDims, npts];
		     		Export[ptsfname,pts];
				,"Strat", (* something goes wrong in Stratified_3dd *)
					pts = getStratND[nDims, npts];
		     		Export[ptsfname,pts]
				,"OwenPlusTree16bits", execString = "owen --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" -t 16 > /dev/null";
   						res = Run[execPrefix<>execString];
			     		If[dbg, Print[execString -> res] ];
				,"ExtensibleLatticeType0", execString = Switch[nDims
					,2,"getExtensibleLattice2D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 0 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,3,"getExtensibleLattice3D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 0 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,4,"getExtensibleLattice4D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 0 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,6,"getExtensibleLattice6D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 0 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					];
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"ExtensibleLatticeType1", execString = Switch[nDims
					,2,"getExtensibleLattice2D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 1 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,3,"getExtensibleLattice3D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 1 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,4,"getExtensibleLattice4D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 1 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,6,"getExtensibleLattice6D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 1 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					];
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"ExtensibleLatticeType2", execString = Switch[nDims
					,2,"getExtensibleLattice2D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 2 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,3,"getExtensibleLattice3D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 2 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,4,"getExtensibleLattice4D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 2 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,6,"getExtensibleLattice6D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 2 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					];
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,"ExtensibleLatticeType3", execString = Switch[nDims
					,2,"getExtensibleLattice2D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 3 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,3,"getExtensibleLattice3D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 3 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,4,"getExtensibleLattice4D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 3 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					,6,"getExtensibleLattice6D -o "<>ptsfname<>" -n "<>ToString[npts]<>" -l 3 --seed "<>ToString[RandomInteger[{0, 2^31}]]
					];
					res = Run[execPrefix<>execString];
		     		If[dbg, Print[execString -> res] ];
				,_, Print["makeMSEref ",pointsetType -> pointsetLabel, " not implemented yet"]; Abort[];
					];
 				execString = "new_integrateND_from_file --nintegrands "<>ToString[nIntegrands]<>" -i "<>ptsfname<>" -o "<>msefname<>" --integrandType "<>ToString[integrandType]<>" --nDims "<>ToString[nDims]<>" > /dev/null";
				res = Run[execPrefix<>execString];
		     	mse = Last @ (Flatten @ Import[msefname]);
		     	If[dbg, Print[execString -> res -> mse] ];
		     	If[!NumberQ[mse], Abort[] ];
		     	(*If[ pointsetLabel != "SOT" && pointsetLabel != "Rank1Lattice"  && pointsetLabel != "PMJ02" && FileExistsQ[ptsfname], DeleteFile[ptsfname] ];*)
		     	If[ FileExistsQ[msefname], DeleteFile[msefname] ];
				(*Run["rm -rf "<>ptsfname<>" "<>msefname ];*)
				If[dbg, Print[{pointsetLabel,integrandTypeLabel}," ",{npts,iPointSet,mse} ] ];
     			mse
     		,{iPointSet,nPointsets}]);

			If[dbg, Print[Sort@mseTab] ];
			
     		mseMean = Mean @ mseTab;
     		mseVariance = Variance @ mseTab;
     		{mseMin,mseMax} = {Min@mseTab, Max@mseTab};
     		AppendTo[dataMSE,{Round[npts],mseMean,mseVariance,mseMin,mseMax,0,0,nPointsets,nIntegrands}];

     		Print["makeMSEref: ",integrandTypeLabel -> pointsetLabel," ",ToString[nDims]<>"D" -> nTrialsMSE ->  Last[ dataMSE[[;;,1;;2]] ] -> dirMSE];
   				Export[dirMSE<>resFname,header,"TEXT"];
 				Export["tmp/tmpdat"<>pid<>".dat",dataMSE];
 				Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFname];
				Print[dirMSE<>resFname, " written."];
		,{iptsPow,powfrom,powto,powstep}];	 (* available from 1K *)
        Run["rm -rf tmp/" ];
   ] (* makeMSEref *)

getCloseestN2D[n_] := Round[Sqrt[n]]^2
getCloseestNND[nDims_:2, n_] := Round[n^(1/nDims)]^nDims

getRealNPts[nDims_:2, npts_:16, pointsetType_:10] :=
    Switch[pointsetType
    ,11, getCloseestNND[nDims, npts]    (* Strat *)
    ,15, getCloseestNND[nDims, npts]    (* RegGrid *)
    ,777,	First @ getOmegaApproxRealNpts[nDims, npts]
    (*,200, getHexGridRealNpts[nDims, npts]*)    (* HexGrid *)
    ,212, getHexGridTorApproxRealNpts[nDims, npts]    (* HexGridTorApproxReal *)
    ,209, getLDBNRealNpts[nDims, npts]    (* LDBN *)
    ,_, npts
    ]
(*		pointsetLabel = Switch[pointsetType 
				(* sequence, ND *)   ,1,"Sobol" ,2,"Halton" ,3,"Faure" ,4,"Niederreiter" ,5,"SobolPlusPlus",6,"SobolGlobal",7,"OwenGlobal"
			(* pointsets, ND *)  ,10,"WN" ,11,"Strat" ,12,"OwenPlus" ,13,"Rank1Lattice" ,14,"DartThrowing"  ,15,"RegGrid" ,16,"SOT" ,17,"SOTPlusLloyd",18,"OwenPureFirstDim2",19,"OwenPure" 
			(* ExtensibleLattices, ND *),20,"ExtensibleLatticeType0" ,21,"ExtensibleLatticeType1" ,22,"ExtensibleLatticeType2" ,23,"ExtensibleLatticeType3"
			(* pointsets, 2D only *) ,200,"HexGrid"  ,201,"Hammersley" ,202,"LarcherPillichshammer" ,203,"NRooks" ,204,"BNOT" ,205,"CMJ" ,206,"BNLDS" 
									 ,207,"PMJ" ,208,"PMJ02" ,209,"LDBN" ,210,"Penrose" ,211,"Fattal" ,212, "HexGridTorApprox"
	    	(* pointsets, 3D only *) ,300,"BCC", 301,"FCC", 302,"HCP", 303,"WeairePhelan"
    		(* 4D Variants of Sobol,OwenPlus,OwenPure *) ,400,"Sobol2356" ,401,"OwenPlus2356" ,402,"OwenPure2356",403,"Sobol2367" ,404,"OwenPlus2367" ,405,"OwenPure2367"
			(* uniformND *) ,500,"UniformND" ,501,"UniformNDwithoutSobol"
			(* uniformND *) ,600,"zsampler",601,"morton",602,"morton01"
			(* pointsets, SobolShiftedKx *) ,701,"SobolShifted1x",702,"SobolShifted2x",703,"SobolShifted3x",704,"SobolShifted4x"
			(* pointsets, OwenPlusShiftedKx *) ,801,"OwenPlusShifted1x",802,"OwenPlusShifted2x",803,"OwenPlusShifted3x",804,"OwenPlusShifted4x"
			(* pointsets, OwenShiftedKx *) ,900,"OwenMicroShift",999,"OwenMicroShiftGlobal",901,"OwenShifted1x",902,"OwenShifted2x",903,"OwenShifted3x",904,"OwenShifted4x"
	    	,_, "unknown" 
		];
*)
getStrat2D[npts_:256] :=
    Block[ {nstrats,xshift,yshift},
    	nstrats = Sqrt[npts];
    	Flatten[#,1]& @ Parallelize @ (Table[
    		{xshift,yshift} = {RandomReal[], RandomReal[]}/nstrats;
   			{(ix-1)/nstrats + xshift,(iy-1)/nstrats + yshift}
    	,{ix,nstrats},{iy,nstrats}]) //N
    ] (* getStrat2D *)

getWN[nDims_:3,npts_:512] := Table[Table[RandomReal[],{nDims}] ,{npts}]
	
getStratND[nDims_:3,npts_:512] :=
    Block[ {nstrats,xshift,yshift,ushift,vshift,sshift,tshift},
    	nstrats = npts^(1/nDims);	(* suppose that npts is already appropriate, passed through getRealNPts[] *)
    	Switch[nDims
    	,1, Flatten[#,1]& @ (Table[
    			xshift = RandomReal[]/nstrats;
   				(ix-1)/nstrats + xshift
    		,{ix,nstrats}]) //N
    	,2, Flatten[#,1]& @ (Table[
    			{xshift,yshift} = {RandomReal[], RandomReal[]}/nstrats;
   				{(ix-1)/nstrats + xshift,(iy-1)/nstrats + yshift}
    		,{ix,nstrats},{iy,nstrats}]) //N
    	,3, Flatten[#,2]& @ (Table[
    			{xshift,yshift,ushift} = {RandomReal[], RandomReal[], RandomReal[]}/nstrats;
   				{(ix-1)/nstrats + xshift,(iy-1)/nstrats + yshift, (iu-1)/nstrats + ushift}
    		,{ix,nstrats},{iy,nstrats},{iu,nstrats}]) //N
    	,4, Flatten[#,3]& @ (Table[
    			{xshift,yshift,ushift,vshift} = Table[RandomReal[],{4}]/nstrats;
   				{(ix-1)/nstrats + xshift,(iy-1)/nstrats + yshift, (iu-1)/nstrats + ushift, (iv-1)/nstrats + vshift}
    		,{ix,nstrats},{iy,nstrats},{iu,nstrats},{iv,nstrats}]) //N
    	,6, Flatten[#,5]& @ (Table[
    			{xshift,yshift,ushift,vshift,sshift,tshift} = Table[RandomReal[],{6}]/nstrats;
   				{(ix-1)/nstrats + xshift,(iy-1)/nstrats + yshift, (iu-1)/nstrats + ushift, (iv-1)/nstrats + vshift, (is-1)/nstrats + sshift, (it-1)/nstrats + tshift}
    		,{ix,nstrats},{iy,nstrats},{iu,nstrats},{iv,nstrats},{is,nstrats},{it,nstrats}]) //N
    	]
    ] (* getStratND *)



showstdRefMSEandDiscrepancy[] := {showstdRefMSE[], showstdRefDiscrepancy[]}

showstdRefMSE[inlbl_:"Gauss2D_pbnot_nIntegrands4096_kmag100_setno0"] :=
    Module[ {fontSz=20(*,powfrom,powto,powstep,kPlusMinus,data,ffitpow10,fitpow10,ffitpow15,fitpow15,ffitpow20,fitpow20,ffitpow30,fitpow30,delta,plotLabel,legends,alldata,fnameLabel*)},
		lbl = inlbl;
		kPlusMinus = .5;
    	{powfrom,powto,powstep} = {2,16,1};

		integrandTypeLabel = "Gauss";
		nDims = 2;
		(*Manipulate[*)
			fnameLabel = integrandTypeLabel ;
	        plotLabel = "Ref MSE "<>ToString[nDims]<>"D   integrandType = "<>integrandTypeLabel;
	        plotLabel = "Ref MSE "<>ToString[nDims]<>"D   integrandType = ";
			dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>fnameLabel<>"/";
			Switch[nDims
			,6,
				data = (Drop[#,1]& @ Import[dirMSE<>"WN_"<>fnameLabel<>".dat"]);
				mseWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				data = (Drop[#,1]& @ Import[dirMSE<>"OwenPure_"<>fnameLabel<>".dat"]);
				mseOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			    alldata = {mseWN,mseOwen01Pure};
		        legends = {"WN","Owen"};
			,_,
				data = (Drop[#,1]& @ Import[dirMSE<>"WN_"<>fnameLabel<>".dat"]);
				mseWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				data = (Drop[#,1]& @ Import[dirMSE<>"Strat_"<>fnameLabel<>".dat"]);
				mseStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				data = (Drop[#,1]& @ Import[dirMSE<>"OwenPure_"<>fnameLabel<>".dat"]);
				mseOwenPure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				data = (Drop[#,1]& @ Import[dirMSE<>"Sobol_"<>fnameLabel<>".dat"]);
				mseSobol = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				data = (Drop[#,1]& @ Import[dirMSE<>"PMJ02_"<>fnameLabel<>".dat"]);
				msePMJ02 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
				
			    alldata = {mseWN,mseStrat,mseOwenPure,msePMJ02} ;
		        legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat", "Owen", "PMJ02"} ] ];
			];
	        
			p = ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[10]},
							{Blue,AbsoluteThickness[10]},
							{Black,AbsoluteThickness[12]},
							{Cyan,AbsoluteThickness[8]},
							{Red,AbsoluteThickness[6]},
							{Darker@Green,AbsoluteDashing[{10, 5}], AbsoluteThickness[6]},
							{Yellow,AbsoluteThickness[6]},
							
							{Blue,Dotted,AbsoluteThickness[2]},

							{Blue,Dotted,AbsoluteThickness[5]},
							{Black,Dotted,AbsoluteThickness[5]},
							{Red,Dotted,AbsoluteThickness[5]},

							{Red,AbsoluteThickness[10]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[2^pow,{pow,powfrom,powto,2}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "MSE", fontSz] }
		           		,ImageSize -> {1024,1024}
		            	,PlotRange->{{2^powfrom,2^powto},{Min[First /@ Flatten[(second /@ mseOwenPure)]], Max[First /@ (second /@ mseOwenPure)] (*Min[ First /@ Flatten[ (second /@ #)& /@ ( alldata)] ]*) }} (*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[2^pow,{pow,powfrom,powto,2}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel<>lbl, Bold, 24] 
		            ];
			(*,Control[{{nDims,2},{2,3,4,6}}]
			,Control[{{integrandTypeLabel,"Gauss"},{"Heaviside", "Gauss" }}]
         ]*)
         p//Print;
         Export["p_"<>lbl<>".png",p]
     ] (* showstdRefMSE *)

(*
rotmx = RandomVariate[CircularRealMatrixDistribution[2], 1]
*)