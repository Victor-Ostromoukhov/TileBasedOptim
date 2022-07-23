(* TileBasedOptim.m
   2022-04 vo, based on fibo-hilbert.m (version 2002/12/14)

makeMSEref[]
showstdRefMSE[]

makeDiscrepancyRef[]
showstdRefDiscrepancy[]

makeOptimL2discrepancy[]  ->  getL2discrepancy[]
showOptimL2discrepancy[]

makeOptimMSESeq[]  ->  getMSE[]
makeOptimMSEPointSets[]
showstdOptimMSE[] :

=========================== prep Optim Data
prepOptimDataSequences[]
prepOptimDataPointsets[]
prepOptimDataPointsetsVarSize[]

=========================== prepSoftEllipses2D[] & prepHeavisideND[]

=========================== MatBuilder
makeMatBuilderMatrices0m2net2D[]
before S22 submission: testCplex and matrixSampler
suppl  S22 submission: MatBuilder and sampler
getMatBuiderPtsND[]

=========================== lois MCQMC July 2022
loismakeMatBuilderMatrices[]
loismakeL2discrepancy[]

*)
 
(****************** globals *******************)
SetDirectory[ToFileName[$HomeDirectory,"TileBasedOptim/"]];
SetOptions[Graphics, ImageSize -> {1024, Automatic}];

epsilon = 10^-10.;
eps = 10^-6.;

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
showMatBuilderIndex = 4096;
showPrevRect = 8192;

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
    	prog = "getGeneralizedL2discrepancy";
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  prog<>" -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat -d "<>ToString[nDims]<>" > /dev/null";
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = Import["tmp/res"<>pid<>".dat"][[1,2]];
        Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];
        discrepancy
   ] (* getGeneralizedL2discrepancy *)

getL2discrepancy[pts_, inptsfname_:"", innDims_:2, dbg_:False] :=
    Module[ {execString,nDims,ptsfname,prog,returnCode, discrepancy},
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
        If[inptsfname == "", (* given : pts *)
         	nDims = Length[First@pts];
        	ptsfname = "tmp/tmp"<>pid<>".dat";
        	Export["tmp/tmp"<>pid<>".dat",N[pts]];      	
         ,(*ELSE  given : inptsfname & innDims *)
        	nDims = innDims;
        	ptsfname = inptsfname;
        ];
        prog =  Switch[nDims
        	,2, "L2discrepancy_fromfile_2dd"
        	,3, "L2discrepancy_fromfile_3dd"
        	,4, "L2discrepancy_fromfile_4dd"
        ];
        execString =  prog<>" -i "<>ptsfname<>" -o tmp/res"<>pid<>".dat > /dev/null";
        returnCode = Run[execPrefix<>execString];
        If[dbg, Print[execString -> returnCode ] ];
        discrepancy = Import["tmp/res"<>pid<>".dat"][[2,2]];
        If[inptsfname == "",  (* given : pts *)
        	Run["rm tmp/tmp"<>pid<>".dat tmp/res"<>pid<>".dat"];
        ,(*ELSE  given : inptsfname & innDims *)
        	Run["rm tmp/res"<>pid<>".dat"];
        ];
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
makeSobolGeneralizedL2discrepancy[]



makeWNL2discrepancy[]
makeStratL2discrepancy[]

makeOwenL2discrepancy[]

makeSobolL2discrepancy[]

makeWNStarDiscrepancy[]
makeStratStarDiscrepancy[]
makeSobolDiscrepancy[]
*)


makeSobolGeneralizedL2discrepancy[nlevels_:12, nDims_:2,dbg_:False] :=
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
        	Print["Processing makeSobolGeneralizedL2discrepancy " -> {npts,d}];
			{npts,d} 
        ,{inpts,nptsMax}];
	    Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; 
        Print[mf @ dtab]
    ] (* makeSobolL2discrepancy *)

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
    ] (* makeSobolL2discrepancy *)


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
        Export["GeneralizedL2discrepancy_Base3SFC2D.pdf",p];
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
        Export["L2discrepancy_Base3SFC2D.pdf",p];
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
        Export["L2discrepancy_Base3SFC2D_MatBuilderOnly.dat.pdf",p];
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
        Export["L2discrepancy_Base3SFC2D_Experiment2Tanguy.dat.pdf",p];
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
    Module[ {res={}, tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,dxy },
    	Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2} = {v1,v2};
            Switch[tileType
              ,typeSq, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{},{2}}, {{},{1}}, {{},{0}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} }
              		];
					AppendTo[res,{typeHRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,					{v1, v2/3}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeVRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 + 1/3 v2,	{mxRot90.v1/3, mxRot90.v2},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeHRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, v2/3}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
              ,typeHRec, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{},{2}}, {{},{1}}, {{},{0}} }
              		];
					AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,				{1/3 v1, v2}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 1/3 v1 + v2,	{mxRot270.v1/3, mxRot270.v2},	{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v1,		{1/3 v1, v2}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
              ,typeVRec, 
              		dxy = Which[
              			v1[[1]] > 0 && v2[[2]] > 0, {{{},{0}}, {{},{1}}, {{},{2}} },
              			v1[[1]] < 0 && v2[[2]] < 0, {{{},{2}}, {{},{1}}, {{},{0}} },
              			v1[[2]] > 0 && v2[[1]] < 0, {{{2},{}}, {{1},{}}, {{0},{}} },
              			v1[[2]] < 0 && v2[[1]] > 0, {{{0},{}}, {{1},{}}, {{2},{}} }
              		];
					AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,				{v1, 1/3 v2}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 1/3 v2 + v1,	{mxRot90.v1, mxRot90.v2/3},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeSq,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, 1/3 v2}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];
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
    Module[ {sfc={}, tileType,matBuilderIndex,refPt,v1,v2,samplingPt,prevrefPt,prevv1,prevv2,xcode,ycode,fcode,norm1,norm2,factor=4,delta=6},
    	Do[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
	    	{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/factor);
	   		AppendTo[sfc,refPt + (norm1+norm2)/1/3^((Length[fcode]+delta)/factor) ];
  			AppendTo[sfc,refPt + v1 + v2 + (-norm1-norm2)/3^((Length[fcode]+delta)/factor) ] ;
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcBase3SFC2D *)

(*getsfcBase3SFC2D[tlst_] :=
    Module[ {sfc={}, tileType,matBuilderIndex,refPt,v1,v2,samplingPt,prevrefPt,prevv1,prevv2,xcode,ycode,fcode},
    	Do[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
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
    Module[ {gl={AbsolutePointSize[10]},tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcBase3SFC2D[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,(*Arrowheads[1/3^(3+(Length[fcode]+Mod[Length[fcode],2])/2)],*)Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
			fcodelen = Length[fcode];
			samplingPt = samplingPt; (*/3^fcodelen;*)
			{norm1,norm2}={v1/Norm[v1],v2/Norm[v2]}/3^((Length[fcode]+Mod[Length[fcode],2])/2);
			cont = {refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt};
    		If[BitAnd[params,showGrayValue] > 0, AppendTo[gl,{GrayLevel[FromDigits[Reverse@fcode,3]/3^fcodelen],Polygon@cont}] ];
    		If[BitAnd[params,showLightGrayTile] > 0, AppendTo[gl,{LightGray,Polygon@cont}] ];
			AppendTo[gl,Flatten[#,1]& @ {(*Point@(refPt+(norm1+norm2)/20),*)bortedStyle,Line@cont } ];
			If[BitAnd[params,showMatBuilderIndex] > 0, AppendTo[gl, {Text[Style[matBuilderIndex,Bold,14,Blue],refPt+(v1+v2)/2 ]} ] ];		
			If[BitAnd[params,showTileType] > 0, AppendTo[gl, {Text[Style[tileType,Bold,14,Blue],refPt+(v1+v2)/2,{1.9,-1}]} ] ];		
			If[BitAnd[params,showOrdinalNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[fcode,3],Bold,14,Red],refPt+(v1+v2)/2,{-1.9,-1}]} ] ];
			If[BitAnd[params,showFcodeInvNumber] > 0, AppendTo[gl, {Text[Style[FromDigits[Reverse@fcode,3],Bold,14,Black],refPt+(v1+v2)/2,{1.9,-1}]} ] ];
			If[BitAnd[params,showTilefcode] > 0, AppendTo[gl, {Text[Style[tab2snosep@fcode,Bold,14,Gray],refPt+(v1+v2)/2,{0,1}]} ] ];
			If[BitAnd[params,showTilexycodes] > 0, AppendTo[gl, {Text[Style[tab2snosep@xcode,Bold,14,Red],refPt+(v1+v2)/2,{1,1}], Text[Style[tab2snosep@ycode,Bold,14,Blue],refPt+(v1+v2)/2,{-1,1}]} ] ];
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt,Text[Style[matBuilderIndex (*FromDigits[Reverse@fcode,3]*),Bold,10,Blue], samplingPt,{-1.1,-1.1}]} ] ];
			If[BitAnd[params,showPrevRect] > 0, AppendTo[gl, {Red, Line@{prevrefPt,prevrefPt+prevv1,prevrefPt+prevv1+prevv2,prevrefPt+prevv2,prevrefPt}} ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3SFC2DTilesGL *)

selectBase3SFC2DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3SFC2DTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,v,indVect,nsubdivs,m},
    	Parallelize @ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
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
			If[dbg, Print[i -> {tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC2DTiles *)


getSamplingPtsBase3SFC2DTiles[tlst_] :=
    Module[ {tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
    	Parallelize @ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[ind]];
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
    	nlevels = 19;
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
		,{i,256}];
    ] (* makeMatBuilderMatrices0m2net2D *)


prepOptimDataBase3SFC2D[innlevels_:6, dbg_:True] :=
    Module[ {},
    	setNo = 1;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	nlevels = innlevels;
    	If[ !FileExistsQ["optim_input_2D/"], CreateDirectory["optim_input_2D/"] ];
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
				fname = "optim_input_2D/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
				exportSelectionBase3SFC2D[fname,seltlst];
				If[dbg,
				p = Graphics[ Append[background,#]& @ getBase3SFC2DTilesGL[seltlst,showLightGrayTile+showSamplingPt], PlotLabel-> iOrdinalAbsolute ];
					p//Print;
					Export["optim_figs_2D/2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
				];
			,{iOrdinalAbsolute,3^(ilevel-1)+1,3^ilevel}];
		,{ilevel,nlevels}];
	] (* prepOptimDataBase3SFC2D *)

exportSelectionBase3SFC2D[fname_, seltlst_] :=
Module[{newtlst,tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
	newtlst = Flatten /@ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = seltlst[[ind]];			
			{tileType,matBuilderIndex,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2},N@refPt,N@{v1,v2},{xcode,ycode},fcode}
		,{ind,Length[seltlst]}];
	Export[fname,newtlst];
] (* exportSelectionBase3SFC2D *)
(*
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
			Do[{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[i]];
				prevrefPt = refPt;
				{prevv1,prevv2} = {v1,v2};
				tlst[[i]] = {tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode};
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
			Do[{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = tlst[[i]];
				samplingPt = prevrefPt + RandomReal[] prevv1 + RandomReal[] prevv2;
				tlst[[i]] = {tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode};
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
*)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeBase3SFC2DL2discrepancy[]
*)

makeBase3SFC2DL2discrepancy[dbg_:False] :=
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
        	Print["Processing makeBase3SFC2DL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2discrepancy *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeBase3SFC2DL2discrepancyExperimentTanguy[]
*)
makeBase3SFC2DL2discrepancyExperimentTanguy[dbg_:False] :=
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
        	Print["Processing makeBase3SFC2DL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D_Experiment2Tanguy.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2discrepancyExperimentTanguy *)

makeBase3SFC2DL2discrepancyMatBuilderOnly[dbg_:False] :=
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
        	Print["Processing makeBase3SFC2DL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D_MatBuilderOnly.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab]
    ] (* makeBase3SFC2DL2discrepancyMatBuilderOnly *)



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

makeBase3SFC2DGeneralizedL2discrepancy[dbg_:False] :=
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
        	Print["Processing makeBase3SFC2DGeneralizedL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab];
		p = showGeneralizedL2discrepancyND[nDims,dtab,"Base3SFC2D",{2,10}] ;
    ] (* makeBase3SFC2DGeneralizedL2discrepancy *)

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
    ] (* makeBase3SFC2DGeneralizedL2discrepancy *)


makeBase3SFC2DGeneralizedL2discrepancy[dbg_:False] :=
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
        	Print["Processing makeBase3SFC2DGeneralizedL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_GeneralizedL2discrepancy/"<>ToString[nDims]<>"D/Base3SFC2D.dat", dtab]; 
        ,{iOrdinalAbsolute,2,nptsMax}];
        Print[mf @ dtab];
		p = showGeneralizedL2discrepancyND[nDims,dtab,"Base3SFC2D",{2,10}] ;
    ] (* makeBase3SFC2DGeneralizedL2discrepancy *)

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
    Module[ {res={}, tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
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
					AppendTo[res,{typeParaZflatRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaZflatLeft, matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	permut12 @ {mx.v1, mx.v2, mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaZflatRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeCubeLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXflatLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 					{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx =  rotmx3D[Pi, v1];
					AppendTo[res,{typeParaXflatRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v3+v1/3,	permut23 @ {mx.v1/3, mx.v2, mx.v3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXflatLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXflatLeft,
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaYlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 						{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeParaYlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3+v1/3+v2,(permut23 @ {mx.v1/3,mx.v2,mx.v3}),{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaYlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3}, 				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaXlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,permut12 @ {mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 		{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatLeft, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,permut23 @ {mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2, v3/3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3/3+v1+v2,permut12 @ {mx.v1,mx.v2,mx.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v3,		{v1, v2, v3/3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,permut13 @ {mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,permut23 @ {mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
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
    Module[ {gl={AbsolutePointSize[10]},tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,cont,sfc,norm1,norm2,fcodelen,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={GrayLevel[.6],AbsoluteThickness[5]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcBase3SFC3D[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
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
    Module[ {sfc={}, tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,norm1,norm2,norm3,k=6,delta=12},
    	Do[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
	    	{norm1,norm2,norm3}={v1/Norm[v1],v2/Norm[v2],v3/Norm[v3]}/3^((Length[fcode]+Mod[Length[fcode],3])/k );
	   		AppendTo[sfc,refPt + (norm1+norm2+norm3)/3^((Length[fcode]+delta)/k  ) ];
  			AppendTo[sfc,refPt + v1 + v2 + v3 - (norm1+norm2+norm3)/3^((Length[fcode]+delta)/k ) ] ;
    	,{ind,Length[tlst]}];
    	Return[sfc]
    ] (* getsfcBase3SFC3D *)


selectBase3SFC3DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3SFC3DTiles[tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,v,indVect,nsubdivs,m},
    	Parallelize @ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
     		nsubdivs = Length[xcode] + Length[ycode] + Length[zcode];
			v = Join[xcode,ycode,zcode];
			m = If[Length[xcode] == Length[ycode],
				mxInv
			,(*ELSE*)
				If[Max@(Abs@(First /@ {v1, v2, v3})) > Max@(Abs@(Last /@ {v1, v2})), mxInvH, mxInvV]
			];
			indVect = Mod[#,3]& /@ (m.v);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3SFC3DTiles *)


getSamplingPtsBase3SFC3DTiles[tlst_] :=
    Module[ {tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode},
    	Parallelize @ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
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
Module[{newtlst,tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode},
	newtlst = Flatten /@ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = seltlst[[ind]];			
			{tileType,matBuilderIndex,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2,prevv3},N@refPt,N@{v1,v2,v3},{xcode,ycode,zcode},fcode}
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
    Module[ {res={}, tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2,prevv3} = {v1,v2,v3};
            Switch[tileType
              ,typeXCube, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
					AppendTo[res,{typeXflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeXflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mxRotX180.v1/3, mxRotX180.v2, mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},									{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeZCube, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
					AppendTo[res,{typeZflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, 1/3 v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
(*					AppendTo[res,{typeZflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mxRotZ180.v1, mxRotZ180.v2, mxRotZ180.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeZflatPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},									{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
*)              ,typeXflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeYPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v1/3+v3,{mxRotX180.v1/3,mxRotX180.v2,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeYPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeZflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeXPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mxRotZ90.v1/3,mxRotZ90.v2,mxRotZ90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];

					(*AppendTo[res,{typeHRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,					{v1, v2/3}, 					{Join[xcode,dxy[[1,1]]],Join[ycode,dxy[[1,2]]]}, Append[fcode,0]} ];
	                AppendTo[res,{typeVRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 + 1/3 v2,	{mxRot90.v1/3, mxRot90.v2},		{Join[xcode,dxy[[2,1]]],Join[ycode,dxy[[2,2]]]}, Append[fcode,1]} ];
	                AppendTo[res,{typeHRec,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt + 2/3 v2,		{v1, v2/3}, 					{Join[xcode,dxy[[3,1]]],Join[ycode,dxy[[3,2]]]}, Append[fcode,2]} ];*)

              ,typeYflatPara, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
					AppendTo[res,{typeXPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 										{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					(*AppendTo[res,{typeYPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mxRotX180.v1,mxRotX180.v2/3,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];*)
					AppendTo[res,{typeXPara,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeXPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
					AppendTo[res,{typeXCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeYCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mxRotX180.v1/3,mxRotX180.v2,mxRotX180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeXCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeYPara, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
					AppendTo[res,{typeYCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},								{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeXCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,	{mxRotY180.v1,mxRotY180.v2/3,mxRotY180.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeYCube,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
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
					AppendTo[res,{typeZflatParaX,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 									{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					AppendTo[res,{typeZflatParaY,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mxRotY90.v1, mxRotY90.v2/3, mxRotY90.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeZflatParaX,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},								{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
*)

(*
subdivBase3SFC3DTiles[tlst_] :=
    Module[ {res={}, tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,prevv3,refPt,v1,v2,v3,xcode,ycode,zcode,fcode,dxyz },
    	Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,{v1,v2,v3},{xcode,ycode,zcode},fcode} = tlst[[ind]];
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
					AppendTo[res,{typeParaZflatRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaZflatRight, matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2+v3/3,	{mx.v1, mx.v2, mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaZflatRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, 1/3 v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeCubeLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{},{2}}, {{},{},{1}}, {{},{},{0}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{},{0}}, {{},{},{1}}, {{},{},{2}} }
              		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXflatLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1/3, v2, v3}, 					{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeParaXflatLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2+v3+v1/3,	{mx.v1/3, mx.v2, mx.v3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXflatLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1 2/3,	{v1/3, v2, v3},				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2/3, v3}, 				{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v2/3+v3,{mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2 2/3,	{v1, v2/3, v3}, 			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXflatLeft,
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaYlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 						{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaYlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,{1,1,1}{mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaYlongLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 				{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYflatRight, 
               		dxyz = Which[
             			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{},{0},{}}, {{},{1}x,{}}, {{},{2},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{},{2},{}}, {{},{1},{}}, {{},{0},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{},{0},{}}, {{},{1},{}}, {{},{2},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,	{v1, v2, v3/3}, 			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1+v3/3+v2,{mx.v1,mx.v2,mx.v3/3},{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeParaXlongRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3 2/3,	{v1, v2, v3/3}, 		{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaZlongRight, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2, v3/3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v3];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v3/3+v1+v2,	{mx.v1,mx.v2,mx.v3/3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeRight,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v3,		{v1, v2, v3/3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaYlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1, v2/3, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v2];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v2/3+v1+v3,	{mx.v1,mx.v2/3,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v2,		{v1, v2/3, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
              ,typeParaXlongLeft, 
              		dxyz = Which[
              			v1[[1]] < 0 && v2[[2]] < 0 && v3[[3]] > 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] < 0 && v2[[2]] > 0 && v3[[3]] < 0, {{{2},{},{}}, {{1},{},{}}, {{0},{},{}} },
              			v1[[1]] > 0 && v2[[2]] < 0 && v3[[3]] < 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} },
              			v1[[1]] > 0 && v2[[2]] > 0 && v3[[3]] > 0, {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} }
             		];
             		dxyz = {{{0},{},{}}, {{1},{},{}}, {{2},{},{}} };
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt,				{v1/3, v2, v3},			{Join[xcode,dxyz[[1,1]]],Join[ycode,dxyz[[1,2]]],Join[zcode,dxyz[[1,3]]]}, Append[fcode,0]} ];
					mx = rotmx3D[Pi, v1];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+v1/3+v2+v3,	{mx.v1/3,mx.v2,mx.v3},	{Join[xcode,dxyz[[2,1]]],Join[ycode,dxyz[[2,2]]],Join[zcode,dxyz[[2,3]]]}, Append[fcode,1]} ];
					AppendTo[res,{typeCubeLeft,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2,prevv3},refPt+2/3 v1,		{v1/3, v2, v3},			{Join[xcode,dxyz[[3,1]]],Join[ycode,dxyz[[3,2]]],Join[zcode,dxyz[[3,3]]]}, Append[fcode,2]} ];
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
nintegrands = 256 1024;
nDims = 2;
Do[
	nPointsets = 16;                                                                                                                                                                                        
	makeMSEref[11, nPointsets, {2,14,1}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[12, nPointsets, {2,14,1/4.}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[19, nPointsets, {2,14,1/4.}, integrandType, nDims, nintegrands];                                                                                                                               
,{integrandType,1,2}]


gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nintegrands = 256 1024;
nDims = 2;
Do[
	nPointsets = 1024;                                                                                                                                                                                        
	makeMSEref[10, nPointsets, {2,14,1}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[11, nPointsets, {2,14,1}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[12, nPointsets, {2,14,1/4.}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[19, nPointsets, {2,14,1/4.}, integrandType, nDims, nintegrands];                                                                                                                               
,{integrandType,1,2}]

gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nintegrands = 256 1024;
nDims = 2;
	integrandType=2;                                                                                                                                                                                   
	nPointsets = 16;     
	makeMSEref[12, nPointsets, {1,16,1/8.}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[19, nPointsets, {1,16,1/8.}, integrandType, nDims, nintegrands];                                                                                                                               

	nPointsets = 1;     
	makeMSEref[1, nPointsets, {1,16,1}, integrandType, nDims, nintegrands];                                                                                                                               
	nPointsets = 256;     
	makeMSEref[10, nPointsets, {1,16,1}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[11, nPointsets, {1,16,1}, integrandType, nDims, nintegrands];                                                                                                                               


gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nintegrands = 256 1024;
nDims = 2;
	integrandType=2;                                                                                                                                                                                   
	nPointsets = 64;     
	integrandType=2;                                                                                                                                                                                   
	makeMSEref[501, nPointsets, {1,8,1/9.}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[12, nPointsets, {1,16,1/8.}, integrandType, nDims, nintegrands];                                                                                                                               
	makeMSEref[11, nPointsets, {1,16,1/8.}, integrandType, nDims, nintegrands];                                                                                                                               

gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nintegrands = 256 1024;
nDims = 2;
integrandType=1;
nPointsets=1;
makeMSEref[208, nPointsets, {2,16,1}, integrandType, nDims, nintegrands];                                                                                                                               


*)

makeMSEref[inpointsetTypes_:10, innPointsets_:1024, powParams_:{2,18,1}, inIntegrandType_:1, innDims_:2, nIntegrands_:1024, consecutiveFlag_:False, dbg_:False] :=
    Module[ {},
    	firstDim = 0;
    	(*If[ Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2]];*)
    	nDims = innDims;
    	integrandType = inIntegrandType;
		integrandTypeLabel = Switch[integrandType,  1,"Heaviside", 2,"SoftEllipses", 3,"Rectangles", 4,"Ellipses", 5,"SoftEllipses_noRot" ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#Analytical	#MSE	#NbPtsets	#Nbintegrands\n";
		fnameLabel = integrandTypeLabel ;
		nPointsets 	= innPointsets ;
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
			(* uniformND *) ,500,"MatBuider" ,501,"MatBuiderMaxDepth"
			(* uniformND *) ,600,"zsampler",601,"morton",602,"morton01"
			(* pointsets, SobolShiftedKx *) ,701,"SobolShifted1x",702,"SobolShifted2x",703,"SobolShifted3x",704,"SobolShifted4x"
			(* pointsets, OwenShiftedKx *) ,900,"OwenMicroShift",999,"OwenMicroShiftGlobal",901,"OwenShifted1x",902,"OwenShifted2x",903,"OwenShifted3x",904,"OwenShifted4x"
	    	,_, "unknown" 
		];
    	If[pointsetLabel == "SOT", {powfrom,powto,powstep} = Switch[nDims,2,{2,17,1},3,{2,16,1},4,{2,17,1}] ];
    	base = If[pointsetType==500 || pointsetType==501, 3, 2];
		{nPtsfrom,nPtsto} = {base^powfrom, base^powto};
		Print[pointsetLabel,{powfrom,powto,powstep} -> " makeMSEref from ",nPtsfrom," to ",nPtsto];
		dataMSE = {};
		{iCounterfrom,iCounterto,iCounterstep} = If[consecutiveFlag, {nPtsfrom,nPtsto,1}, {powfrom, powto, powstep}];
		Do[	
     		If[pointsetLabel == "SOT" && nDims == 4 && iptsPow == 17, nPointsets = 1 ]; (* Only 1 available *)
     		nptsTarget = If[consecutiveFlag, iCounter, base^iCounter];
   			npts = Round[nptsTarget];
   			If[!consecutiveFlag, npts = getRealNPts[nDims, npts, pointsetType] ];
    		resFname = If[consecutiveFlag, pointsetLabel<>"_"<>fnameLabel<>"_consecutive.dat", pointsetLabel<>"_"<>fnameLabel<>".dat"];
    		
			mseTab = ( Parallelize @  Table[   				
				ptsfname = "tmp/pts_"<>ToString[iPointSet]<>pid<>".dat";
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
				,"OwenPlus", execString = "owen --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
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
				,"MatBuider", 
					npts = Round[base^iCounter];
					infname = "MatBuilder_matrices/2D_0m2net_"<>i2s[RandomInteger[{1,16}]]<>".dat";
					owenFlag = False;
					depth = Ceiling[Log[base, npts]];
					(*Print["Processing MatBuider @ ", iCounter -> npts -> {infname,owenFlag,depth}];*)
					pts = getMatBuiderPtsND[npts, infname, owenFlag,depth, 2, base]; (* 3^19=1162261467 *)
		     		Export[ptsfname,pts];
				,"MatBuiderMaxDepth", 
					npts = Round[base^iCounter];
					infname = "MatBuilder_matrices/2D_0m2net_"<>i2s[RandomInteger[{1,16}]]<>".dat";
					owenFlag = True;
					depth = 19;
					(*Print["Processing MatBuiderMaxDepth @ ", iCounter -> npts -> {infname,owenFlag,depth}];*)
					pts = getMatBuiderPtsND[npts, infname, owenFlag,depth, 2, base]; (* 3^19=1162261467 *)
		     		Export[ptsfname,pts];
				,"WN", 
					pts = getWN[nDims, npts];
		     		Export[ptsfname,pts];
				,"Strat", (* something goes wrong in Stratified_3dd *)
					pts = getStratND[nDims, npts];
		     		Export[ptsfname,pts]
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
     		mseVariance = If[nPointsets == 1, 0 , Variance @ mseTab];
     		{mseMin,mseMax} = {Min@mseTab, Max@mseTab};
     		AppendTo[dataMSE,{Round[npts],mseMean,mseVariance,mseMin,mseMax,0,0,nPointsets,nIntegrands}];

			If[!consecutiveFlag || (consecutiveFlag &&  Mod[iCounter,100] == 99 ) || iCounter == iCounterto,
     			Print["makeMSEref: ",integrandTypeLabel -> pointsetLabel," ",ToString[nDims]<>"D" -> nPointsets ->  Last[ dataMSE[[;;,1;;2]] ] -> dirMSE];
   				Export[dirMSE<>resFname,header,"TEXT"];
 				Export["tmp/tmpdat"<>pid<>".dat",dataMSE];
 				Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFname];
				Print[dirMSE<>resFname, " written."];
			];
		,{iCounter, iCounterfrom,iCounterto,iCounterstep}]; 
        Run["rm -rf tmp/" ];
   ] (* makeMSEref *)

getCloseestN2D[n_] := Round[Sqrt[n]]^2
getCloseestNND[nDims_:2, n_] := Round[n^(1/nDims)]^nDims

getRealNPts[nDims_:2, npts_:16, pointsetType_:10] :=
    Switch[pointsetType
    ,11, getCloseestNND[nDims, npts]    (* Strat *)
    ,15, getCloseestNND[nDims, npts]    (* RegGrid *)
    ,777,	First @ getOmegaApproxRealNpts[nDims, npts]
    ,212, getHexGridTorApproxRealNpts[nDims, npts]    (* HexGridTorApproxReal *)
    ,209, getLDBNRealNpts[nDims, npts]    (* LDBN *)
    ,_, npts
    ]

getWN[nDims_:3,npts_:512] := Table[Table[RandomReal[],{nDims}] ,{npts}]


getMatBuiderPtsND[npts_:3^2,infname_:"MatBuilder_matrices/2D_0m2net_000001.dat",owenFlag_:True,depth_:19,nDims_:2,base_:3,inseed_:0] := (* 3^19=1162261467 *)
    Block[ {execString,outfname,res,seed},
    	seed = If[inseed == 0, RandomInteger[2^16], inseed];
		If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
    	outfname = "tmp/pts"<>pid<>".dat";
		execString = "sampler --nDims "<>ToString[nDims]<>" -o "<>outfname<>" -n "<>ToString[npts]<>" -i "<>infname
			<>If[owenFlag," --owen  --depth "<>ToString[depth], " " ]
			<>" -p "<>ToString[base]<>" --matrixSize "<>ToString[19]<>" --seed "<>ToString[seed]<>" > /dev/null";
		Print[execString];
		res = Run[execPrefix<>execString];
		If[res != 0, Print[execString -> res] ];
		Import[outfname]
    ]

		(*execString = "matrixSampler --nDims "<>ToString[nDims]<>" -o "<>outfname<>" -n "<>ToString[npts]<>" -i "<>infname
			<>If[owenFlag," -p  --depth "<>ToString[depth], " " ]
			<>" -b "<>ToString[base]<>" --size "<>ToString[19]<>" --seed "<>ToString[RandomInteger[2^16] ]<>" > /dev/null";*)

tstMatBuider[dbg_:True] :=
    Module[ {},
    	base = 3;
    	frame={Blue,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}] };
    	Do[
			npts = Round[base^iCounter];
			vticks = Table[Line[{{0,(i-1)/base^Floor[(iCounter+1)/2]},{1,(i-1)/base^Floor[(iCounter+1)/2]}}],{i,base^Floor[(iCounter+1)/2]}];
			hticks = Table[Line[{{(i-1)/base^Floor[(iCounter)/2],0},{(i-1)/base^Floor[(iCounter)/2],1}}],{i,base^Floor[(iCounter)/2]}];
			infname = "MatBuilder_matrices/2D_0m2net_"<>i2s[RandomInteger[{1,16}]]<>".dat";
			owenFlag = True;
			depth = 19;
			pts = getMatBuiderPtsND[npts, infname, owenFlag,depth, 2, base]; (* 3^19=1162261467 *)
			ipts = Floor[npts pts];
			Print["Processing MatBuiderMaxDepth @ ", iCounter -> npts -> {infname,owenFlag,depth} -> getTfactor[base,ipts]];
			If[dbg, Print[Graphics[{frame,Cyan,vticks,hticks,Black,AbsolutePointSize[5],Point/@pts},PlotLabel->getTfactor[base,ipts],ImageSize->3/2{1024,1024}] ] ];
    	,{iCounter,1,10,1}]	
]
    
getStratND[nDims_:3,npts_:512] :=
    Block[ {nstrats,xshift,yshift,ushift,vshift,sshift,tshift,nstratsAsked,res},
    	nstratsAsked = npts^(1/nDims);	(* suppose that npts is already appropriate, passed through getRealNPts[] *)
    	nstrats = If[IntegerQ[nstratsAsked], nstratsAsked, Ceiling[nstratsAsked] ];
    	res = Switch[nDims
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
    	];
    	If[nstratsAsked == nstrats, res, RandomSample[res][[;;npts]] ]
    ] (* getStratND *)

showstdRefMSE[] :=
    Module[ {powfrom,powto,powstep,kPlusMinus,data,plotLabel,legends,alldata,fnameLabel,dirMSE},
    	consecutiveFlag = False;
		fontSz = 14;
		kPlusMinus = .5;
    	{powfrom,powto,powstep} = {2,16,1};

		nDims = 2;
		(*integrandTypeLabel = "Heaviside";*)
		
		Manipulate[
			fnameLabel = integrandTypeLabel ;
	        plotLabel = "Ref MSE "<>ToString[nDims]<>"D   integrandType = "<>integrandTypeLabel;
			dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>fnameLabel<>"/";

			data = (Drop[#,1]& @ Import[dirMSE<>"WN_"<>fnameLabel<>".dat"]);
			mseWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirMSE<>"Strat_"<>fnameLabel<>".dat"]);
			mseStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirMSE<>"OwenPure_"<>fnameLabel<>".dat"]);
			mseOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			mseOwen01PureRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirMSE<>"OwenPlus_"<>fnameLabel<>".dat"]);
			mseOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			mseOwenPlusRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];

			(*data = (Drop[#,1]& @ Import[dirMSE<>"Sobol_"<>fnameLabel<>".dat"]);
			mseSobol01 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)
			(*data = (Drop[#,1]& @ Import[dirMSE<>"PMJ02_"<>fnameLabel<>".dat"]);
			msePMJ02 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)

		    alldata = {mseWN, mseStrat, mseOwen01Pure,  mseOwenPlus} ;
	        legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat", "Owen", "OwenPlus32" } ] ];
	        
			ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[2]},
							{Blue,AbsoluteThickness[2]},
							{Black,AbsoluteThickness[2]},
							{Red,AbsoluteThickness[2]},
							{Cyan,AbsoluteThickness[2]},
							{Darker@Green,AbsoluteThickness[2]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[2^pow,{pow,powfrom,powto,2}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "MSE", fontSz] }
		           		,ImageSize -> {1024,1024}
		            	(*,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ mseOwenPlusRaw), Min @@ (second /@ mseOwenPlusRaw) }} *)(*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[2^pow,{pow,powfrom,powto,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel, Bold, 24] 
		            ]			
			(*,Control[{{consecutiveFlag,False},{True,False}}]*)
			,Control[{{integrandTypeLabel,"Heaviside"},{"SoftEllipses", "Heaviside"(*, "Ellipses", "Rectangles", "SoftEllipses_noRot"*) }}]
         ]
     ] (* showstdRefMSE *)



(*----------------------------------------- L2discrepancy *)
(*


gitpull
math
<<TileBasedOptim/TileBasedOptim.m
nDims = 2;
	DiscrepancyType = 1;                                                                                                                                                                                     
	nPointsets = 256;   
	makeDiscrepancyRef[11, nPointsets, {2,16,1/2.}, DiscrepancyType, nDims];                                                                                                                               

	makeDiscrepancyRef[10, nPointsets, {2,16,1/4.}, DiscrepancyType, nDims];                                                                                                                               
	nPointsets = 16;   
	makeDiscrepancyRef[12, nPointsets, {2,16,1/4.}, DiscrepancyType, nDims];                                                                                                                               
	makeDiscrepancyRef[19, nPointsets, {2,16,1/4.}, DiscrepancyType, nDims];                                                                                                                               

*)
makeDiscrepancyRef[inpointsetTypes_:10, innPointsets_:1024, powParams_:{2,18,1}, inDiscrepancyType_:1, innDims_:2, consecutiveFlag_:False, dbg_:False] :=
    Module[ {},
    	firstDim = 0;
    	(*If[ Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2]];*)
    	nDims = innDims;
    	DiscrepancyType = inDiscrepancyType;
		DiscrepancyTypeLabel = Switch[DiscrepancyType,  1,"L2Discrepancy", 2,"StarDiscrepancy", 3,"GeneralizedL2Discrepancy" ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
		nPointsets 	= innPointsets ;
        {powfrom,powto,powstep} = powParams;

		dirDiscrepancy = "data_"<>DiscrepancyTypeLabel<>"/"<>ToString[nDims]<>"D/";
        If[ !FileExistsQ[dirDiscrepancy], CreateDirectory[dirDiscrepancy] ];
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
		{nPtsfrom,nPtsto} = {2^powfrom, 2^powto};
		Print[pointsetLabel,{powfrom,powto,powstep} -> " makeDiscrepancyRef from ",nPtsfrom," to ",nPtsto];
		dataDiscrepancy = {};
		{iCounterfrom,iCounterto,iCounterstep} = If[consecutiveFlag, {nPtsfrom,nPtsto,1}, {powfrom, powto, powstep}];
		Do[	
     		If[pointsetLabel == "SOT" && nDims == 4 && iptsPow == 17, nPointsets = 1 ]; (* Only 1 available *)
     		nptsTarget = If[consecutiveFlag, iCounter, 2^iCounter];
   			(*npts = If[consecutiveFlag, nptsTarget, getRealNPts[nDims, pointsetLabel, pointsetType] ];*)
   			npts = Round[nptsTarget];
   			If[!consecutiveFlag, npts = getRealNPts[nDims, npts, pointsetType] ];

    		resFname = If[consecutiveFlag, pointsetLabel<>"_"<>DiscrepancyTypeLabel<>"_consecutive.dat", pointsetLabel<>"_"<>DiscrepancyTypeLabel<>".dat"];
			DiscrepancyTab = ( Parallelize @  Table[   				
				ptsfname = "tmp/pts_"<>ToString[iPointSet]<>pid<>".dat";
				Discrepancyfname = "tmp/d"<>pid<>".dat";
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
				,"OwenPlus", execString = "owen --nDims "<>ToString[nDims]<>" -p 1 -o "<>ptsfname<>" -n "<>ToString[npts]<>" -s "<>ToString[RandomInteger[2^31]]<>" --max_tree_depth_32_flag 1 > /dev/null";
   					res = Run[execPrefix<>execString];
			     	If[dbg, Print[execString -> res] ];
				,"WN", 
					pts = getWN[nDims, npts];
		     		Export[ptsfname,pts];
				,"Strat", (* something goes wrong in Stratified_3dd *)
					pts = getStratND[nDims, npts];
		     		Export[ptsfname,pts]
				,_, Print["makeDiscrepancyRef ",pointsetType -> pointsetLabel, " not implemented yet"]; Abort[];
					];
 		     	Discrepancy = Switch[DiscrepancyType,  1, getL2discrepancy[pts,ptsfname] ];
		     	If[dbg, Print[execString -> res -> Discrepancy] ];
		     	If[!NumberQ[Discrepancy], Abort[] ];
		     	(*If[ pointsetLabel != "SOT" && pointsetLabel != "Rank1Lattice"  && pointsetLabel != "PMJ02" && FileExistsQ[ptsfname], DeleteFile[ptsfname] ];*)
		     	If[ FileExistsQ[Discrepancyfname], DeleteFile[Discrepancyfname] ];
				(*Run["rm -rf "<>ptsfname<>" "<>Discrepancyfname ];*)
				If[dbg, Print[{pointsetLabel,DiscrepancyTypeLabel}," ",{npts,iPointSet,Discrepancy} ] ];
     			Discrepancy
     		,{iPointSet,nPointsets}]);

			If[dbg, Print[Sort@DiscrepancyTab] ];
			
     		DiscrepancyMean = Mean @ DiscrepancyTab;
     		DiscrepancyVariance = If[nPointsets == 1, 0 , Variance @ DiscrepancyTab];
     		{DiscrepancyMin,DiscrepancyMax} = {Min@DiscrepancyTab, Max@DiscrepancyTab};
     		AppendTo[dataDiscrepancy,{Round[npts],DiscrepancyMean,DiscrepancyVariance,DiscrepancyMin,DiscrepancyMax,0,0,nPointsets,0}];

			If[!consecutiveFlag || (consecutiveFlag &&  Mod[iCounter,100] == 99 ) || iCounter == iCounterto,
     			Print["makeDiscrepancyRef: ",DiscrepancyTypeLabel -> pointsetLabel," ",ToString[nDims]<>"D" -> nPointsets ->  Last[ dataDiscrepancy[[;;,1;;2]] ] -> dirDiscrepancy];
   				Export[dirDiscrepancy<>resFname,header,"TEXT"];
 				Export["tmp/tmpdat"<>pid<>".dat",dataDiscrepancy];
 				Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirDiscrepancy<>resFname];
				Print[dirDiscrepancy<>resFname, " written."];
			];
		,{iCounter, iCounterfrom,iCounterto,iCounterstep}]; 
        (*Run["rm -rf tmp/" ];*)
   ] (* makeDiscrepancyRef *)

showstdRefDiscrepancy[] :=
    Module[ {powfrom,powto,powstep,kPlusMinus,data,plotLabel,legends,alldata,dirDiscrepancy},
    	consecutiveFlag = False;
		fontSz = 14;
		kPlusMinus = .5;
    	{powfrom,powto,powstep} = {2,16,1};

    	DiscrepancyType = 1;
		DiscrepancyTypeLabel = Switch[DiscrepancyType,  1,"L2Discrepancy", 2,"StarDiscrepancy", 3,"GeneralizedL2Discrepancy" ];

		nDims = 2;
		
		(*Manipulate[*)
	        plotLabel = "Ref "<>DiscrepancyTypeLabel<>" "<>ToString[nDims]<>"D  "<>DiscrepancyTypeLabel;
			dirDiscrepancy = "data_"<>DiscrepancyTypeLabel<>"/"<>ToString[nDims]<>"D/";

			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"WN_"<>DiscrepancyTypeLabel<>".dat"]);
			DiscrepancyWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"Strat_"<>DiscrepancyTypeLabel<>".dat"]);
			DiscrepancyStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"OwenPure_"<>DiscrepancyTypeLabel<>".dat"]);
			DiscrepancyOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			DiscrepancyOwen01PureRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirDiscrepancy<>"OwenPlus_"<>DiscrepancyTypeLabel<>".dat"]);
			DiscrepancyOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			DiscrepancyOwenPlusRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];

			(*data = (Drop[#,1]& @ Import[dirDiscrepancy<>"Sobol_"<>fnameLabel<>".dat"]);
			DiscrepancySobol01 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)
			(*data = (Drop[#,1]& @ Import[dirDiscrepancy<>"PMJ02_"<>fnameLabel<>".dat"]);
			DiscrepancyPMJ02 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)

			(*limit = Table[npts = 2^i; {npts, k npts^-3 (Log[npts])^},{i,Length[data]}];*)

		    alldata = {DiscrepancyWN, DiscrepancyStrat, DiscrepancyOwen01Pure,  DiscrepancyOwenPlus} ;
	        legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat", "Owen", "OwenPlus32" } ] ];
	        
			ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[2]},
							{Blue,AbsoluteThickness[2]},
							{Black,AbsoluteThickness[2]},
							{Red,AbsoluteThickness[2]},
							{Cyan,AbsoluteThickness[2]},
							{Darker@Green,AbsoluteThickness[2]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[2^pow,{pow,powfrom,powto,2}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ DiscrepancyTypeLabel, fontSz] }
		           		,ImageSize -> {1024,1024}
		            	,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ DiscrepancyOwenPlusRaw), Min @@ (second /@ DiscrepancyOwenPlusRaw) }} (*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[2^pow,{pow,powfrom,powto,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel, Bold, 24] 
		            ]			
			(*,Control[{{consecutiveFlag,False},{True,False}}]*)
			(*,Control[{{integrandTypeLabel,"Heaviside"},{"SoftEllipses", "Heaviside", "Ellipses", "Rectangles", "SoftEllipses_noRot" }}]
         ]*)
     ] (* showstdRefDiscrepancy *)

(* <<<<<<<<<<<<<<<<<<<<<< consecutive L2discrepancy

gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeWNL2discrepancy[]
makeStratL2discrepancy[]
makeSobolL2discrepancy[]
makeOwenL2discrepancy[]

*)

makeSobolL2discrepancy[nlevels_:14, nDims_:3,dbg_:False] :=
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
        	Print["Processing makeSobolL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Sobol.dat", dtab]; 
        ,{inpts,nptsMax}];
        Print[mf @ dtab]
    ] (* makeSobolL2discrepancy *)

makeOwenL2discrepancy[nlevels_:14, nDims_:3,dbg_:False] :=
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
        	Print["Processing makeOwenL2discrepancy " -> {npts,d}];
			AppendTo[dtab, {npts,d} ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Owen.dat", dtab]; 
        ,{inpts,nptsMax}];
        Print[mf @ dtab]
    ] (* makeOwenL2discrepancy *)

makeWNL2discrepancy[nlevels_:14, ntrials_:64, nDims_:3] :=
    Module[ {},
        dtab = {};
        Do[
			npts = 2^ilevel;
			trials = Parallelize @ Table[
				pts = Table[Table[RandomReal[],{nDims}],{i,npts}];
				{npts,getL2discrepancy[pts]}
			,{itrail,ntrials}];
			AppendTo[dtab, Mean @ trials ];
        	Print["Processing makeWNL2discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/WN.dat", dtab]; 
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]


makeStratL2discrepancy[nlevels_:14, ntrials_:64, nDims_:3] :=
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
        	Print["Processing makeStratL2discrepancy level ",ilevel -> mf[dtab] ];
	        Export["data_L2discrepancy/"<>ToString[nDims]<>"D/Strat.dat", dtab]; push
        ,{ilevel,nlevels}];
        Print[mf @ dtab]
    ]
(*
*)


(*------------------------------------------ prepOptimDataBase3Simple2D ---------------------------------------------------*)

typeSimpleHSq = 	1;
typeSimpleVSq = 	2;
typeSimpleH1Rect = 	3;
typeSimpleH2Rect = 	4;
typeSimpleV1Rect = 	5;
typeSimpleV2Rect = 	6;


subdivBase3Simple2DTiles[tlst_] :=
    Module[ {res={}, tileType,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2, matBuilderIndex  },
    	Table[
			{tileType,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2}, matBuilderIndex } = tlst[[ind]];
			prevrefPt = refPt; {prevv1,prevv2} = {v1,v2};
            Switch[tileType
              ,typeSimpleHSq, 
					AppendTo[res,{typeSimpleH1Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt,	{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleH2Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2/3,	{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleH1Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2 2/3,	{v1, v2/3} ,matBuilderIndex } ];
              ,typeSimpleVSq, 
					AppendTo[res,{typeSimpleV1Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt,	{v1/3, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleV2Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1/3,	{v1/3, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleV1Rect,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 2/3,	{v1/3, v2} ,matBuilderIndex } ];
              ,typeSimpleH1Rect, 
 					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt,			{1/3 v1, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1/3,	{1/3 v1, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 2/3,{1/3 v1, v2} ,matBuilderIndex } ];
              ,typeSimpleH2Rect, 
 					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt,			{1/3 v1, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1/3,	{1/3 v1, v2} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v1 2/3,{1/3 v1, v2} ,matBuilderIndex } ];
              ,typeSimpleV1Rect, 
 					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt,			{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2/3,	{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2 2/3,{v1, v2/3} ,matBuilderIndex } ];
              ,typeSimpleV2Rect, 
 					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt,			{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleVSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2/3,	{v1, v2/3} ,matBuilderIndex } ];
					AppendTo[res,{typeSimpleHSq,samplingPt,prevrefPt,{prevv1,prevv2},refPt + v2 2/3,{v1, v2/3} ,matBuilderIndex } ];
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivBase3Simple2DTiles *)

getBase3Simple2DTilesGL[tlst_,params_:showSFC] :=
    Module[ {gl={AbsolutePointSize[10]},tileType,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,cont,sfc,norm1,norm2,fcodelen,matBuilderIndex,
    		bortedStyle={Cyan,AbsoluteThickness[1]}, sfcStyle={Orange,AbsoluteThickness[3]}},
    	If[BitAnd[params,showSFC] > 0, sfc = getsfcBase3SFC2D[tlst]; 
    		AppendTo[gl,Flatten[#,1]& @ {sfcStyle,Line@sfc}];
    		If[BitAnd[params,showArrows] > 0, AppendTo[gl,Flatten[#,1]& @ {sfcStyle,(*Arrowheads[1/3^(3+(Length[fcode]+Mod[Length[fcode],2])/2)],*)Arrow/@(Partition[#,2]&@sfc)}] ] ];    	
    	Do[
			{tileType,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2}, matBuilderIndex } = tlst[[ind]];
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
			If[BitAnd[params,showSamplingPt] > 0, AppendTo[gl, {Black,Point@samplingPt } ] ];
    	,{ind,Length[tlst]}];
    	Return[gl]
    ] (* getBase3Simple2DTilesGL *)

selectBase3Simple2DTiles[tlst_,intensity_:.8] := Select[tlst, FromDigits[Reverse@Last[#],3]/3^Length[Last[#]] < intensity & ]

fillSamplingPtsBase3Simple2DTiles[ilevel_,tlst_, mxTab_,mxInv_,mxInvH_,mxInvV_] :=
     Module[ {tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode,v,indVect,nsubdivs,m,matBuilderIndex},
    	Parallelize @ Table[
			{tileType,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2}, matBuilderIndex } = tlst[[ind]];
     		nsubdivs = ilevel;
			v = Join[xcode,ycode];
			m = If[Length[xcode] == Length[ycode],
				mxInv
			,(*ELSE*)
				If[Max@(Abs@(First /@ {v1, v2})) > Max@(Abs@(Last /@ {v1, v2})), mxInvH, mxInvV]
			];
			indVect = Mod[#,3]& /@ (m.v);
			matBuilderIndex = FromDigits[#,3]& @ (Reverse @ indVect);
			samplingPt = (FromDigits[#,3]& /@ (Mod[#,3]& /@ {mxTab[[1,;;nsubdivs,;;nsubdivs]].indVect, mxTab[[2,;;nsubdivs,;;nsubdivs]].indVect}) ) / 3^nsubdivs;
			If[dbg, Print[i -> {tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},fcode}] ];
			{tileType,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},matBuilderIndex}
    	,{ind,Length[tlst]}]
    ] (* getSamplingPtsBase3Simple2DTiles *)

exportSelectionBase3Simple2D[fname_, seltlst_] :=
Module[{newtlst,tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
	newtlst = Flatten /@ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = seltlst[[ind]];			
			{tileType,matBuilderIndex,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2},N@refPt,N@{v1,v2},{xcode,ycode},fcode}
		,{ind,Length[seltlst]}];
	Export[fname,newtlst];
] (* exportSelectionBase3Simple2D *)



(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m

Parallelize @ Do[prepOptimDataSequences[10, i, True, False], {i, 65,256}]

*)
prepOptimDataSequences[innoctaves_:5, insetNo_: 1, prevFlag_: True, dbg_:True] :=
    Module[ {},
        (*If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];*)
        
    	owenFlag = True;
    	depth = 19;
    	nDims = 2;
    	base = 3;
    	seed = RandomInteger[2^16];
    	
    	setNo = insetNo;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	frame={AbsoluteThickness[3],Green,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}] };
    	noctaves = innoctaves;
    	tilesDir = If[prevFlag, "Tiles_Seq_PrevLevel/", "Tiles_Seq_CurLevel/" ];
   		tilesDirFigs = If[prevFlag, "Tiles_Seq_PrevLevel_Figs/", "Tiles_Seq_CurLevel_Figs/" ];
    	If[ !FileExistsQ[tilesDir] && !dbg, CreateDirectory[tilesDir] ];
    	If[ !FileExistsQ[tilesDirFigs] && dbg, CreateDirectory[tilesDirFigs] ];
   	
    	mxfname = "MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat";
		mxTab = readMatBuilderMatrix[mxfname];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		
		pts0 = getMatBuiderPtsND[1, mxfname, owenFlag, depth, nDims, base, seed ][[1]];
		tlst = {{0,pts0, {0,0},{{1,0},{0,1}} } };
		
		pts = getMatBuiderPtsND[base^noctaves, mxfname, owenFlag, depth, nDims, base, seed ];
		Do[
			{iOrdinalAbsoluteFrom,iOrdinalAbsoluteTo} = {base^(ioctave-1)+1,base^ioctave};
			If[dbg,
				ngrid = base^Floor[ioctave/2]; ngridstep = 1/ngrid;
				octaveGrid={Cyan,Table[Line[{{i ngridstep,0},{i ngridstep,1}}],{i,0,ngrid}],Table[Line[{{0,i ngridstep},{1,i ngridstep}}],{i,0,ngrid}]};
			];
			Do[
				gl = {};
				npts = base^ioctave;
				prevoctave = ioctave-1;
				{x,y} = npts pts[[iOrdinalAbsolute]];
				If[EvenQ[prevoctave],
					{ix,iy} = Quotient[{x,y}, base^((ioctave+1)/2)];
					{k1,k2} = If[EvenQ[ix+iy], {base^((ioctave-1)/2),base^((ioctave+1)/2)}, {base^((ioctave+1)/2),base^((ioctave-1)/2)}];
					{v1,v2} = {{1/k1,0},{0,1/k2}} ;
					{refx,refy} = Quotient[{x,y}, {k2,k1}] / {k1,k2};
					
					{prevv1,prevv2} = {{1,0},{0,1}} / base^((ioctave-1)/2);
					{prevrefx,prevrefy} = {ix,iy} / base^((ioctave-1)/2);
					If[dbg,
						curRect = {Blue,Dashed,AbsoluteThickness[4],Line[{{prevrefx,prevrefy}, {prevrefx,prevrefy}+prevv1,{prevrefx,prevrefy}+prevv1+prevv2,{prevrefx,prevrefy}+prevv2, {prevrefx,prevrefy}} ]};
						AppendTo[gl,{LightYellow,Rectangle[{refx,refy}, {refx,refy}+v1+v2],Red,Line[{{refx,refy}, {refx,refy}+v1,{refx,refy}+v1+v2,{refx,refy}+v2, {refx,refy}} ] } ];
						AppendTo[gl,curRect];
					];
				,(*ELSE*)
					{dx,dy} = {base^(ioctave/2),base^(ioctave/2)};
					{refx,refy} = Quotient[{x,y}, {dx,dy}] / {dx,dy};
					{v1,v2} = {{1,0},{0,1}} / {dx,dy};
					
					{ix,iy} = Quotient[{x,y}, base^((ioctave+2)/2)];
					{k1,k2} = If[EvenQ[ix+iy], {base^((ioctave)/2),base^((ioctave+2)/2)}, {base^((ioctave+2)/2),base^((ioctave)/2)}];
					{prevv1,prevv2} = 3 {{1/k1,0},{0,1/k2}};
					{prevrefx,prevrefy} =3  Quotient[{x,y}, {k2,k1}] / {k1,k2};
					If[dbg,
						curRect = {Blue,Dashed,AbsoluteThickness[4],Line[{{prevrefx,prevrefy}, {prevrefx,prevrefy}+prevv1,{prevrefx,prevrefy}+prevv1+prevv2,{prevrefx,prevrefy}+prevv2, {prevrefx,prevrefy}} ]};
						AppendTo[gl,{LightYellow,Rectangle[{refx,refy}, {refx,refy}+v1+v2],Red,Line[{{refx,refy}, {refx,refy}+v1,{refx,refy}+v1+v2,{refx,refy}+v2, {refx,refy}} ] } ];
						AppendTo[gl,curRect];
					];
				];
				If[prevFlag, 
					AppendTo[tlst,{iOrdinalAbsolute-1,{x,y}/npts,N@{prevrefx,prevrefy},N@{prevv1,prevv2}}];
				,(*ELSE*)
					AppendTo[tlst,{iOrdinalAbsolute-1,{x,y}/npts,N@{refx,refy},N@{v1,v2}}];
				];
				If[dbg,
					p = Graphics[ {frame,octaveGrid,gl,AbsolutePointSize[5],Point/@pts[[;;iOrdinalAbsolute]],
						Table[Text[Style[i-1,14],pts[[i]],{-1,-1}],{i,iOrdinalAbsolute}]}, PlotLabel-> {ioctave,iOrdinalAbsolute} ];
					p//Print;
					Export[tilesDirFigs<>"2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
				];
			,{iOrdinalAbsolute,iOrdinalAbsoluteFrom,iOrdinalAbsoluteTo}];
		,{ioctave,noctaves}];
		fname = tilesDir<>"2D_0m2net_set_"<>i2s[setNo]<>"_uptoOctave_"<>ToString[noctaves]<>"_seed_"<>ToString[seed]<>".dat";
		If[!dbg, Export[fname,Flatten/@(tlst)]; Print["Writing ",fname," done."] ];
	] (* prepOptimDataSequences *)


(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m

Do[prepOptimDataPointsets[10, i, True, False], {i, 65,256}]

*)
prepOptimDataPointsets[innoctaves_:5, insetNo_: 1, prevFlag_: True, dbg_:True] :=
    Module[ {},
        (*If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];*)
        
    	owenFlag = True;
    	depth = 19;
    	nDims = 2;
    	base = 3;
    	seed = RandomInteger[2^16];
    	
    	setNo = insetNo;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	frame={AbsoluteThickness[3],Green,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}] };
    	noctaves = innoctaves;
    	tilesDir = If[prevFlag, "Tiles_Pointsets_PrevLevel/SetNo_"<>i2s[setNo]<>"_seed_"<>ToString[seed]<>"/", "Tiles_Pointsets_CurLevel/SetNo_"<>i2s[setNo]<>"_seed_"<>ToString[seed]<>"/" ];
   		tilesDirFigs = If[prevFlag, "Tiles_Pointsets_PrevLevel_Figs/SetNo_"<>i2s[setNo]<>"_seed_"<>ToString[seed]<>"/", "Tiles_Pointsets_CurLevel_Figs/SetNo_"<>i2s[setNo]<>"_seed_"<>ToString[seed]<>"/" ];
    	If[ !FileExistsQ[tilesDir] && !dbg, CreateDirectory[tilesDir] ];
    	If[ !FileExistsQ[tilesDirFigs] && dbg, CreateDirectory[tilesDirFigs] ];
    	
    	mxfname = "MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat";
		mxTab = readMatBuilderMatrix[mxfname];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		
		targetList = makeOctavesBaseN[{1,noctaves,1},3];
		
		pts = getMatBuiderPtsND[base^noctaves, mxfname, owenFlag, depth, nDims, base, seed ];
		{iOrdinalAbsoluteFrom,iOrdinalAbsoluteTo} = {1,base^noctaves};
		Do[
			iOrdinalAbsolute = targetList[[ilst]];
			tlst = Parallelize @ Table[
				ioctave = 1 + Floor @ Log[base,iOrdinalAbsolute];
				gl = {};
				npts = base^ioctave;
				prevoctave = ioctave-1;
				{x,y} = npts pts[[iWithinSet]];
				If[EvenQ[prevoctave],
					{ix,iy} = Quotient[{x,y}, base^((ioctave+1)/2)];
					{k1,k2} = If[EvenQ[ix+iy], {base^((ioctave-1)/2),base^((ioctave+1)/2)}, {base^((ioctave+1)/2),base^((ioctave-1)/2)}];
					{v1,v2} = {{1/k1,0},{0,1/k2}} ;
					{refx,refy} = Quotient[{x,y}, {k2,k1}] / {k1,k2};
					
					{prevv1,prevv2} = {{1,0},{0,1}} / base^((ioctave-1)/2);
					{prevrefx,prevrefy} = {ix,iy} / base^((ioctave-1)/2);
					If[dbg && iWithinSet == iOrdinalAbsolute,
						curRect = {Blue,Dashed,AbsoluteThickness[4],Line[{{prevrefx,prevrefy}, {prevrefx,prevrefy}+prevv1,{prevrefx,prevrefy}+prevv1+prevv2,{prevrefx,prevrefy}+prevv2, {prevrefx,prevrefy}} ]};
						AppendTo[gl,{LightYellow,Rectangle[{refx,refy}, {refx,refy}+v1+v2],Red,Line[{{refx,refy}, {refx,refy}+v1,{refx,refy}+v1+v2,{refx,refy}+v2, {refx,refy}} ] } ];
						AppendTo[gl,curRect];
					];
				,(*ELSE*)
					{dx,dy} = {base^(ioctave/2),base^(ioctave/2)};
					{refx,refy} = Quotient[{x,y}, {dx,dy}] / {dx,dy};
					{v1,v2} = {{1,0},{0,1}} / {dx,dy};
					
					{ix,iy} = Quotient[{x,y}, base^((ioctave+2)/2)];
					{k1,k2} = If[EvenQ[ix+iy], {base^((ioctave)/2),base^((ioctave+2)/2)}, {base^((ioctave+2)/2),base^((ioctave)/2)}];
					{prevv1,prevv2} = 3 {{1/k1,0},{0,1/k2}};
					{prevrefx,prevrefy} =3  Quotient[{x,y}, {k2,k1}] / {k1,k2};
					If[dbg && iWithinSet == iOrdinalAbsolute,
						curRect = {Blue,Dashed,AbsoluteThickness[4],Line[{{prevrefx,prevrefy}, {prevrefx,prevrefy}+prevv1,{prevrefx,prevrefy}+prevv1+prevv2,{prevrefx,prevrefy}+prevv2, {prevrefx,prevrefy}} ]};
						AppendTo[gl,{LightYellow,Rectangle[{refx,refy}, {refx,refy}+v1+v2],Red,Line[{{refx,refy}, {refx,refy}+v1,{refx,refy}+v1+v2,{refx,refy}+v2, {refx,refy}} ] } ];
						AppendTo[gl,curRect];
					];
				];
				If[prevFlag, 
					{iWithinSet-1,{x,y}/npts,N@{prevrefx,prevrefy},N@{prevv1,prevv2}}
				,(*ELSE*)
					{iWithinSet-1,{x,y}/npts,N@{refx,refy},N@{v1,v2}}
				]
			,{iWithinSet,1,iOrdinalAbsolute}];
			If[dbg,
				ngrid = base^Floor[ioctave/2]; ngridstep = 1/ngrid;
				octaveGrid={Cyan,Table[Line[{{i ngridstep,0},{i ngridstep,1}}],{i,0,ngrid}],Table[Line[{{0,i ngridstep},{1,i ngridstep}}],{i,0,ngrid}]};
				p = Graphics[ {frame,octaveGrid,gl,AbsolutePointSize[5],Point/@pts[[;;iOrdinalAbsolute]],
					Table[Text[Style[i-1,14],pts[[i]],{-1,-1}],{i,iOrdinalAbsolute}]}, PlotLabel-> {{ioctave,Floor[ioctave/2]}->iOrdinalAbsolute} ];
				p//Print;
				Export[tilesDirFigs<>"2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
			];
			fname = tilesDir<>"2D_0m2net_set_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".dat";
			If[!dbg, Export[fname,Flatten/@(tlst)]; Print["Writing ",fname," done."] ];
		,{ilst,Length[targetList]}];
	] (* prepOptimDataPointsets *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m

Do[prepOptimDataPointsetsVarSize[mag, 5, i, False], {mag, .7, 2, .01}, {i, 10}]
Do[prepOptimDataPointsetsVarSize[mag, 6, i, False], {mag, .7, 2, .01}, {i, 10}]


*)
prepOptimDataPointsetsVarSize[inmag_:1.0, innoctaves_:5, insetNo_: 1, dbg_:True] :=
    Module[ {},
        (*If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];*)
        mag = inmag;
    	owenFlag = True;
    	depth = 19;
    	nDims = 2;
    	base = 3;
    	seed = RandomInteger[2^16];
    	
    	setNo = insetNo;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	frame={AbsoluteThickness[3],Green,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}] };
    	noctaves = innoctaves;
    	tilesDir = "Tiles_Pointsets_VarSize/SetNo_"<>i2s[setNo]<>"_mag_"<>r2s[mag, 3,2]<>"/" ;
   		tilesDirFigs = "Tiles_Pointsets_VarSize_Figs/SetNo_"<>i2s[setNo]<>"_mag_"<>r2s[mag, 3,2]<>"/" ;
    	If[ !FileExistsQ[tilesDir] && !dbg, CreateDirectory[tilesDir] ];
    	If[ !FileExistsQ[tilesDirFigs] && dbg, CreateDirectory[tilesDirFigs] ];
    	
    	mxfname = "MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat";
		mxTab = readMatBuilderMatrix[mxfname];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		
		targetList = {(*3^(noctaves-2),*) 3^noctaves};
		
		pts = getMatBuiderPtsND[base^noctaves, mxfname, owenFlag, depth, nDims, base, seed ];
		{iOrdinalAbsoluteFrom,iOrdinalAbsoluteTo} = {1,base^noctaves};
		Do[
			iOrdinalAbsolute = targetList[[ilst]];
			tlst = (*Parallelize @*) Table[
				ioctave = Floor @ Log[base,iOrdinalAbsolute];
				gl = {};
				npts = base^ioctave;
				{x,y} = npts pts[[iWithinSet]];
				{ix,iy} = Quotient[{x,y}, base^((ioctave)/2)];
				k = base^((ioctave)/2);
				{v1,v2} = {{1/k,0},{0,1/k}} ;
				{refx,refy} = Quotient[{x,y}, {k,k}] / {k,k};
				center = {refx,refy} + (v1+v2)/2;
				{bottomLeftx,bottomLefty} = center - mag (v1+v2)/2;
				{topRightx,topRighty} = center + mag (v1+v2)/2;
					
				If[OddQ[ioctave],
					{ix,iy} = Quotient[{x,y}, base^((ioctave+1)/2)];
					{k1,k2} = If[EvenQ[ix+iy], {base^((ioctave-1)/2),base^((ioctave+1)/2)}, {base^((ioctave+1)/2),base^((ioctave-1)/2)}];
					{v1,v2} = {{1/k1,0},{0,1/k2}} ;
					{refx,refy} = Quotient[{x,y}, {k2,k1}] / {k1,k2};
					center = {refx,refy} + (v1+v2)/2;
					{bottomLeftx,bottomLefty} = center - mag Sqrt[Det[{v1, v2}]]/2;
					{topRightx,topRighty} = center + mag Sqrt[Det[{v1, v2}]]/2;
				];					
				bottomLeftx = Max[0,bottomLeftx];
				bottomLefty = Max[0,bottomLefty];
				topRightx = Min[1,topRightx];
				topRighty = Min[1,topRighty];
				If[dbg && iWithinSet == iOrdinalAbsolute,
					(*AppendTo[gl,{LightYellow,Rectangle[{refx,refy}, {refx,refy}+v1+v2],Red,Line[{{refx,refy}, {refx,refy}+v1,{refx,refy}+v1+v2,{refx,refy}+v2, {refx,refy}} ] } ];*)
					AppendTo[gl,{LightYellow,Rectangle[{bottomLeftx,bottomLefty}, {topRightx,topRighty}] ,Red,Point@ center} ];
				];
				{iWithinSet-1,{x,y}/npts,N@{bottomLeftx,bottomLefty},N@{{topRightx-bottomLeftx,0},{0,topRighty-bottomLefty}}}
			,{iWithinSet,1,iOrdinalAbsolute}];
			If[dbg,
				ngrid = base^Floor[ioctave/2]; ngridstep = 1/ngrid;
				octaveGrid={Cyan,Table[Line[{{i ngridstep,0},{i ngridstep,1}}],{i,0,ngrid}],Table[Line[{{0,i ngridstep},{1,i ngridstep}}],{i,0,ngrid}]};
				p = Graphics[ {frame,octaveGrid,gl,AbsolutePointSize[5],Point/@pts[[;;iOrdinalAbsolute]],
					Table[Text[Style[i-1,14],pts[[i]],{-1,-1}],{i,iOrdinalAbsolute}]}, PlotLabel-> {{ioctave,Floor[ioctave/2]}->iOrdinalAbsolute} ];
				p//Print;
				Export[tilesDirFigs<>"2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
			];
			fname = tilesDir<>"2D_0m2net_set_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".dat";
			If[!dbg, Export[fname,Flatten/@(tlst)]; Print["Writing ",fname," done."] ];
		,{ilst,Length[targetList]}];
	] (* prepOptimDataPointsetsVarSize *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
	makeOptimMSEVarSize[3^6, {1,1}]
	showOptimMSEVarSize[3^6]
*)
makeOptimMSEVarSize[innpts_:3^5, setFromTo_:{1,10}, suffix_:"VarSize", innDims_:2, dbg_:False] :=
    Module[ {},
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
    	nDims = innDims;
        dtab = {};
        setNo = 1;
    	integrandType = 2;
    	npts=innpts;
    	{setFrom,setTo} = setFromTo;
    	optimType = optimTypeMSEOptimisationSoftEllipses;
		integrandTypeLabel = Switch[integrandType,  1,"Heaviside", 2,"SoftEllipses", 3,"Rectangles", 4,"Ellipses", 5,"SoftEllipses_noRot" ];
		optimTypeL2OptimisationLabel = Switch[optimType
			,optimTypeL2Optimisation,"L2Optimisation"
			,optimTypeMSEOptimisationHeaviside,"MSEOptimisationHeaviside"
			,optimTypeMSEOptimisationSoftEllipses,"MSEOptimisationSoftEllipses"];
        
		dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/";
        If[ !FileExistsQ[dirMSE], CreateDirectory[dirMSE] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
		datamse = {};
   	    resFname = optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_"<>suffix<>"_level_"<>i2s[npts]<>".dat";
   	    
   	    	
		Do[
        	mseTab = Parallelize @ Table[
	        	setNo = isetNo;
	    		dataDir = "Output/Tiles_Pointsets_VarSize/SetNo_"<>i2s[setNo]<>"_mag_"<>r2s[mag, 3,2]<>"/" ;
	       		fname = dataDir<>"2D_0m2net_set_"<>i2s[setNo]<>"_level_"<>i2s[npts]<>".dat";
	       		If[FileExistsQ[fname],
	       			data = Import[fname];
	       			If[Length[data] >= npts,
						pts = data[[;;npts, 2;;3]];
						If[dbg, ipts = Round[ npts pts ];Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
		    	    	mse = getMSE[pts,"",nDims,integrandType];
		    	    	Print[{npts,mag,setNo} -> mse];
						{mag,mse}
	       			,(*ELSE*)
	       				Nothing
	       			]
				,(*ELSE*)
					Print[fname, " does not exist"];
					Nothing
	       		]
			,{isetNo,setFrom,setTo}];
	 		mseMean = Mean @ mseTab;
	 		mseVariance = If[Length[mseTab] <= 1, 0 , Variance @ (Last /@ mseTab)];
	 		{mseMin,mseMax} = {Min@(Last /@ mseTab), Max@(Last /@ mseTab)};
		    Print[iOrdinalAbsolute, " ", resFname  -> mf[{{mseMean,mseVariance},{mseMin,mseMax}}] -> Length[mseTab] ];			
	 		AppendTo[datamse,Flatten @ {mseMean,mseVariance,mseMin,mseMax,0,0,Length[mseTab],0}];	
			Export[dirMSE<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",datamse];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFname];
			Print[dirMSE<>resFname, " written."];
			mseTab
		,{mag, 0.7, 2, 0.01}];
   ] (* makeOptimMSEVarSize *)

 showOptimMSEVarSize[npts_:3^5] :=
    Module[ {},
    	kPlusMinus = 1;
    	suffix = "VarSize";
    	integrandTypeLabel="SoftEllipses";
    	optimTypeL2OptimisationLabel="MSEOptimisationSoftEllipses";
		dirMSE = "data_MSE/2D/"<>integrandTypeLabel<>"/";
		resFname = optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_"<>suffix<>"_level_"<>i2s[npts]<>".dat";
		data = Import[dirMSE<>resFname];
		msedata = Drop[#,1]& @ Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
		msedata//mf//Print;
		ListPlot[msedata, ImageSize -> {2 1024, Automatic}]
		(* 243 pts :
			1.34	(4.3\[PlusMinus]1.1)*10^-10 
			1.45	(4.2\[PlusMinus]1.3)*10^-10
		*)
] (* showOptimMSEVarSize *)
    
selectBase3SFC2DTilesMatBuilderOnly[tlst_,intensityInt_] := Select[tlst, second[#] < intensityInt & ]

prepOptimDataBase3PointSets2DFromMatBuilder[innlevels_:6, dbg_:False] :=
    Module[ {},
    	setNo = 1;
		background = {LightYellow, Polygon[{{0,0},{0,1},{1,1},{1,0},{0,0}}]};
    	nlevels = innlevels;
    	If[ !FileExistsQ["Tiles_PointSets/"], CreateDirectory["Tiles_PointSets/"] ];
    	If[ !FileExistsQ["Tiles_PointSets_Figs/"], CreateDirectory["Tiles_PointSets_Figs/"] ];
		mxTab = readMatBuilderMatrix["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>".dat"];
		mxInvTab = readMatBuilderInvMatrices["MatBuilder_matrices/2D_0m2net_"<>i2s[setNo]<>"_inv.dat"];
		tlst = {{typeSq,0,{0,0}, {0,0},{{1,0},{0,1}}, {0,0},{{1,0},{0,1}}, {{},{}} ,{}} };
		counts = {1,3,4,6,9,13,19,27,39,56,81,117,168,243,350,505,729,1051,1516,2187,3154,4549,6561};
		Do[
			tlst = subdivBase3SFC2DTiles @ tlst;
			If[EvenQ[ilevel], mxInv = mxInvTab[[ilevel,1]] ];
			If[OddQ[ilevel],{mxInvH, mxInvV} = mxInvTab[[ilevel]] ];
			tlst = fillSamplingPtsBase3SFC2DTiles[tlst,mxTab,mxInv,mxInvH,mxInvV];
			Parallelize @ Do[
				seltlst = selectBase3SFC2DTilesMatBuilderOnly[tlst, iOrdinalAbsolute ];
				fname = "Tiles_PointSets/2D_0m2net_set_"<>ToString[setNo]<>"_level_"<>ToString[iOrdinalAbsolute]<>".dat";
				exportSelectionBase3SFC2D[fname,seltlst];
				Print[Length[seltlst] -> " Exporting " -> fname];
				If[dbg,
					p = Graphics[ Append[background,#]& @ getBase3SFC2DTilesGL[seltlst,showLightGrayTile+showMatBuilderIndex+showPrevRect+showSamplingPt], PlotLabel-> iOrdinalAbsolute ];
					(*p//Print;*)
					Export["Tiles_PointSets_Figs/2D_0m2net_"<>i2s[setNo]<>"_level_"<>i2s[iOrdinalAbsolute]<>".png", p];
				];
			,{iOrdinalAbsolute,3^(ilevel-1)+1,3^ilevel}];
		,{ilevel,nlevels}];
	] (* prepOptimDataBase3SFCMatBuilderOnly2D *)


exportSelectionBase3SFC2D[fname_, seltlst_] :=
Module[{newtlst,tileType,matBuilderIndex,samplingPt,prevrefPt,prevv1,prevv2,refPt,v1,v2,xcode,ycode,fcode},
	newtlst = Sort @ (Flatten /@ Table[
			{tileType,matBuilderIndex,samplingPt,prevrefPt,{prevv1,prevv2},refPt,{v1,v2},{xcode,ycode},fcode} = seltlst[[ind]];			
			{matBuilderIndex,N@samplingPt,N@prevrefPt,N@{prevv1,prevv2} }
		,{ind,Length[seltlst]}] );
	Export[fname,newtlst];
] (* exportSelectionBase3SFC2D *)


getMSE[pts_, inptsfname_:"", innDims_:2, inIntegrandType_:2, dbg_:False] :=
    Module[ {(*execString,nDims,ptsfname,integrandType,nIntegrands,res,mse,msefname*)},
		msefname = "tmp/mse"<>pid<>".dat";
    	If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
        If[inptsfname == "", (* given : pts *)
         	nDims = Length[First@pts];
        	ptsfname = "tmp/tmp"<>pid<>".dat";
        	Export["tmp/tmp"<>pid<>".dat",N[pts]];      	
         ,(*ELSE  given : inptsfname & innDims *)
        	nDims = innDims;
        	ptsfname = inptsfname;
        ];
    	integrandType = inIntegrandType;
    	nIntegrands = Switch[inIntegrandType
    		,1,512 1024	(* "Heaviside" *)
   			,2,256 1024	/2	(* "SoftEllipses" *)
   		];
 		execString = "new_integrateND_from_file --nintegrands "<>ToString[nIntegrands]<>" -i "<>ptsfname<>" -o "<>msefname<>" --integrandType "<>ToString[integrandType]<>" --nDims "<>ToString[nDims]<>" > /dev/null";
		res = Run[execPrefix<>execString];
     	mse = Last @ (Flatten @ Import[msefname]);
     	If[dbg, Print[execString -> res -> mse] ];
     	If[!NumberQ[mse], Abort[] ];
        mse
   ] (* getMSE *)


doubleCheck[] :=
    Module[ {},
    	nDims = 2;
    	integrandType = 2;
        npts = 81;
        pts = Import["src/New_Optimize_MSE_2DTiles/Data/Input/Tiles_Seq_PrevLevel/2D_0m2net_set_000001_uptoOctave_8_seed_36827.dat"][[;;npts,2;;3]];
        mse = getMSE[pts,"",nDims,integrandType];
        Print[npts -> mse];
    ]

optimTypeL2Optimisation = 1;
optimTypeMSEOptimisationSoftEllipses = 2;
optimTypeMSEOptimisationHeaviside = 3;
(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
	makeOptimMSESeq[]
	makeOptimMSEPointSets[]

*)

makeOptimMSESeq[optimType_:optimTypeMSEOptimisationSoftEllipses, inIntegrandType_:1, setFromTo_:{1,64}, octaves_:{1,7,1}, suffix_:"Seq_PrevLevel", innDims_:2, dbg_:False] :=
    Module[ {},
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
    	nDims = innDims;
        dtab = {};
        setNo = 1;
    	integrandType = inIntegrandType;
		integrandTypeLabel = Switch[integrandType,  1,"Heaviside", 2,"SoftEllipses", 3,"Rectangles", 4,"Ellipses", 5,"SoftEllipses_noRot" ];
		optimTypeL2OptimisationLabel = Switch[optimType
			,optimTypeL2Optimisation,"L2Optimisation"
			,optimTypeMSEOptimisationHeaviside,"MSEOptimisationHeaviside"
			,optimTypeMSEOptimisationSoftEllipses,"MSEOptimisationSoftEllipses"];
        
		dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/";
        If[ !FileExistsQ[dirMSE], CreateDirectory[dirMSE] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
		datamse = {};

		{octaveFrom,octaveTo,octaveStep} = octaves;
		counters = makeOctavesBaseN[{octaveFrom,octaveTo,octaveStep}];
   	    resFname = optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_"<>suffix<>".dat";
   	    
		dataDir = "Output/Tiles_"<>suffix<>"/";		
		files = FileNames["*.dat",{dataDir}];
		If[setFromTo != {1,1},
   	    	{setFrom,setTo} = setFromTo;
   	    ,(*ELSE*)
  	    	{setFrom,setTo} = {1,Length[files]};
		];

		bestCandidates = {};
		
        mseTabAll = Table[
			npts = counters[[iOrdinalAbsolute]];
	        mseTab = Parallelize @ Table[
	        	setNo = isetNo;
	       		fname = Switch[optimType
	       			,optimTypeMSEOptimisationSoftEllipses,
	       			files[[setNo]]
	       		];
	       		(*Print["Pricessing ",npts," pts "->fname ];*)
	       		If[FileExistsQ[fname],
	       			data = Import[fname];
	       			If[Length[data] >= npts,
						pts = data[[;;npts, 2;;3]];
						If[dbg, ipts = Round[ npts pts ];Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
		    	    	mse = getMSE[pts,"",nDims,integrandType];
						{npts,mse}
	       			,(*ELSE*)
	       				Nothing
	       			]
				,(*ELSE*)
					Print[fname, " does not exist"];
					Nothing
	       		]
        	,{isetNo,setFrom,setTo}];
	 		mseMean = Mean @ mseTab;
	 		mseVariance = If[Length[mseTab] <= 1, 0 , Variance @ (Last /@ mseTab)];
	 		{mseMin,mseMax} = {Min@(Last /@ mseTab), Max@(Last /@ mseTab)};
		    Print[iOrdinalAbsolute, " ", resFname  -> mf[{{mseMean,mseVariance},{mseMin,mseMax}}] -> Length[mseTab] ];			
	 		AppendTo[datamse,Flatten @ {mseMean,mseVariance,mseMin,mseMax,0,0,Length[mseTab],0}];	
			Export[dirMSE<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",datamse];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFname];
			Print[dirMSE<>resFname, " written."];
			mseTab
        ,{iOrdinalAbsolute,Length[counters]}];
   ] (* makeOptimMSESeq *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
	makeOptimMSEPointSets[]

*)

makeOptimMSEPointSets[optimType_:optimTypeMSEOptimisationSoftEllipses, inIntegrandType_:1, setFromTo_:{1,19}, octaves_:{1,7,1}, suffix_:"Pointsets_PrevLevel", innDims_:2, dbg_:False] :=
    Module[ {},
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
    	nDims = innDims;
        dtab = {};
        setNo = 1;
    	integrandType = inIntegrandType;
		integrandTypeLabel = Switch[integrandType,  1,"Heaviside", 2,"SoftEllipses", 3,"Rectangles", 4,"Ellipses", 5,"SoftEllipses_noRot" ];
		optimTypeL2OptimisationLabel = Switch[optimType
			,optimTypeL2Optimisation,"L2Optimisation"
			,optimTypeMSEOptimisationHeaviside,"MSEOptimisationHeaviside"
			,optimTypeMSEOptimisationSoftEllipses,"MSEOptimisationSoftEllipses"];
        
		dirMSE = "data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/";
        If[ !FileExistsQ[dirMSE], CreateDirectory[dirMSE] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
		datamse = {};

		{octaveFrom,octaveTo,octaveStep} = octaves;
		counters = makeOctavesBaseN[{octaveFrom,octaveTo,octaveStep}];
   	    resFname = optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_"<>suffix<>".dat";
   	    
		dataDir = "Output/Tiles_"<>suffix<>"/";		
		dirs = FileNames["SetNo*",{dataDir}];
		If[setFromTo != {1,1},
   	    	{setFrom,setTo} = setFromTo;
   	    ,(*ELSE*)
  	    	{setFrom,setTo} = {1,Length[dirs]};
		];
		
        mseTabAll = Table[
			npts = counters[[iOrdinalAbsolute]];
	        mseTab = (*Parallelize @*) Table[
	        	setNo = isetNo;
				dir = dirs[[setNo]];
				fnames = FileNames["*"<>i2s[npts]<>".dat",{dir}];
	       		If[ fnames =!= {},
	       			data = Import[fnames[[1]] ];
					{fnames[[1]] -> Length[data]}//Print;
	       			If[Length[data] >= npts,
						pts = data[[;;npts, 2;;3]];
						If[dbg, ipts = Round[ npts pts ];Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
		    	    	mse = getMSE[pts,"",nDims,integrandType];
		       			Print["Processing ",npts," pts "->fnames[[1]] -> mse];
						{npts,mse}
	       			,(*ELSE*)
	       				Nothing
	       			]
				,(*ELSE*)
					Nothing
	       		]
        	,{isetNo,setFrom,setTo}];
	 		mseMean = Mean @ mseTab;
	 		mseVariance = If[Length[mseTab] <= 1, 0 , Variance @ (Last /@ mseTab)];
	 		{mseMin,mseMax} = {Min@(Last /@ mseTab), Max@(Last /@ mseTab)};
		    Print[iOrdinalAbsolute, " ", resFname  -> mf[{{mseMean,mseVariance},{mseMin,mseMax}}] -> Length[mseTab] ];			
	 		AppendTo[datamse,Flatten @ {mseMean,mseVariance,mseMin,mseMax,0,0,Length[mseTab],0}];	
			Export[dirMSE<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",datamse];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFname];
			Print[dirMSE<>resFname, " written."];
			mseTab
        ,{iOrdinalAbsolute,Length[counters]}];
   ] (* makeOptimMSE *)

        (*Do[
	        thebest = Intersection @@ bestCandidates[[ibest;;]];
	        If[Length[thebest] > 0, Break[] ];
        ,{ibest,Length[bestCandidates]}];
		datamse = {};
        Do[
        	mseTab = mseTabAll[[iOrdinalAbsolute , thebest]];
	 		mseMean = Mean @ mseTab;
	 		mseVariance = If[Length[mseTab] <= 1, 0 , Variance @ (Last /@ mseTab)];
	 		{mseMin,mseMax} = {Min@(Last /@ mseTab), Max@(Last /@ mseTab)};
		    Print[iOrdinalAbsolute, " ", resFnameBest  -> mf[{{mseMean,mseVariance},{mseMin,mseMax}}] -> Length[mseTab] ];
	 		AppendTo[datamse,Flatten @ {mseMean,mseVariance,mseMin,mseMax,0,0,Length[mseTab],0}];	
			Export[dirMSE<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",datamse];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirMSE<>resFnameBest];
			Print[dirMSE<>resFnameBest, " written."];
        ,{iOrdinalAbsolute,Length[counters]}];*)
         

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
makeOptimL2discrepancy[optimTypeL2Optimisation]
*)
makeOptimL2discrepancy[optimType_:optimTypeL2Optimisation, setFromTo_:{1,5}, innDims_:2, dbg_:False] :=
    Module[ {},
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
    	nDims = innDims;
        dtab = {};
        ordinalAbsoluteMax = 3^6;
		optimTypeL2OptimisationLabel = Switch[optimType,  1,"L2Optimisation",  2,"MSEOptimisationSoftEllipses",  3,"MSEOptimisationHeaviside" ];
   		dirL2discrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
        If[ !FileExistsQ[dirL2discrepancy], CreateDirectory[dirL2discrepancy] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
   	    resFname = optimTypeL2OptimisationLabel<>".dat";	    
   	    {setFrom,setTo} = setFromTo;
		dataDiscrepancy = {};
		dirDiscrepancy = "data_L2Discrepancy/"<>ToString[nDims]<>"D/";
        Do[
			npts = iOrdinalAbsolute;
	        DiscrepancyTab = (Parallelize @ Table[
	       		fname = Switch[optimType
	       			,optimTypeL2Optimisation,
	       			"src/Optimize_L2Discrepancy_2DTiles_Noise_Cancelling/Repetitions/Repetition_"<>ToString[setNo]<>"/Output/level_"<>ToString[iOrdinalAbsolute]<>".dat"
	       		];
	       		If[FileExistsQ[fname],
					pts = Import[fname][[;;,2;;3]];
					If[dbg, ipts = Round[ npts pts ];
						Print[Graphics[{{Cyan,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},AbsolutePointSize[10],Point/@pts}, ImageSize->{1024,1024}/2, PlotLabel->{ilevel,npts,testDyadicPartitioningNDFull@ipts}]]];
		        	L2discrepancy = getL2discrepancy[pts,"",nDims]; 
		        	Print[iOrdinalAbsolute, " ", resFname  -> L2discrepancy];
					{npts,L2discrepancy}
				,(*ELSE*)
					Nothing
	       		]
        	,{setNo,setFrom,setTo}]);
	 		DiscrepancyMean = Mean @ DiscrepancyTab;
	 		DiscrepancyVariance = If[Length[DiscrepancyTab] <= 1, 0 , Variance @ (Last /@ DiscrepancyTab)];
	 		{DiscrepancyMin,DiscrepancyMax} = {Min@(Last /@ DiscrepancyTab), Max@(Last /@ DiscrepancyTab)};
	 		AppendTo[dataDiscrepancy,Flatten @ {DiscrepancyMean,DiscrepancyVariance,DiscrepancyMin,DiscrepancyMax,0,0,setTo-setFrom+1,0}];	
			Export[dirDiscrepancy<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",dataDiscrepancy];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirDiscrepancy<>resFname];
			Print[dirDiscrepancy<>resFname, " written."];
	    ,{iOrdinalAbsolute,2,ordinalAbsoluteMax}];
		Run["rm tmp/tmpdat"<>pid<>".dat"];
   ] (* makeOptimL2discrepancy *)

showOptimL2discrepancy[optimType_:optimTypeL2Optimisation, insetNo_:1, innDims_:2, dbg_:False] :=
    Module[ {},
		fontSz = 14;
		kPlusMinus = 1;
    	{powfrom,powto,powstep} = {2,16,1};

    	nDims = innDims;
        dtab = {};
        nptsMax = 653; (*3^6;*)
        setNo = insetNo;
		optimTypeL2OptimisationLabel = Switch[optimType,  1,"L2Optimisation",  2,"MSEOptimisationSoftEllipses",  3,"MSEOptimisationHeaviside" ];
        
		dirL2discrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
   	    resFname = dirL2discrepancy<>optimTypeL2OptimisationLabel<>".dat";
        If[ !FileExistsQ[resFname], Print[resFname," does not exist."]; Abort[]; ];
        data = (Drop[#,1]& @ Import[resFname]);
   		L2discrepancyOptim = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
		(*Print[L2discrepancyOptim//mf];*)
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"WN_L2Discrepancy.dat"]);
			L2discrepancyWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"Strat_L2Discrepancy.dat"]);
			L2discrepancyStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"OwenPure_L2Discrepancy.dat"]);
			L2discrepancyOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			L2discrepancyOwen01PureRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"OwenPlus_L2Discrepancy.dat"]);
			L2discrepancyOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			L2discrepancyOwenPlusRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];

		    alldata = {L2discrepancyWN, L2discrepancyStrat, L2discrepancyOwen01Pure,  L2discrepancyOwenPlus, L2discrepancyOptim} ;
	        legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat", "Owen", "OwenPlus32", optimTypeL2OptimisationLabel} ] ];

		Manipulate[
	        plotLabel = "Optim vs. Ref L2discrepancy "<>ToString[nDims]<>"D"; 
	        
			ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[2]},
							{Blue,AbsoluteThickness[2]},
							{Black,AbsoluteThickness[2]},
							{Red,AbsoluteThickness[2]},
							{Darker@Green,AbsoluteThickness[2]},
							{Cyan,AbsoluteThickness[2]},
							{Gray,AbsoluteThickness[2]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[2^pow,{pow,powfrom,powto,2}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "L2discrepancy", fontSz] }
		           		,ImageSize -> {1024,1024}
		            	(*,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ L2discrepancyOwenPlusRaw), Min @@ (second /@ L2discrepancyOwenPlusRaw) }} *)(*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[2^pow,{pow,powfrom,powto,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel, Bold, 24] 
		            ]			
			,Control[{{optimTypeL2OptimisationLabel,"L2Optimisation"},{"L2Optimisation", "MSEOptimisationSoftEllipses", "MSEOptimisationHeaviside"} } ]
         ]

   ] (* showOptimL2discrepancy *)
optimTypeL2Optimisation = 1;
optimTypeMSEOptimisationSoftEllipses = 2;
optimTypeMSEOptimisationHeaviside = 3;

 
 showstdOptimMSE[octaves_:{1,8,1}] :=
    Module[ {powfrom,powto,powstep,kPlusMinus,data,plotLabel,legends(*,alldata*)},
    	consecutiveFlag = False;
		fontSz = 14;
		kPlusMinus = 1/2.;
    	{powfrom,powto,powstep} = octaves; (* powers of 3 *)

		nDims = 2;
		(*integrandTypeLabel = "Heaviside";*)
		
		optimTypeL2OptimisationLabel = "MSEOptimisationSoftEllipses";
		(*integrandTypeLabel = "Heaviside";
		integrandTypeLabel = "SoftEllipses";*)
		Manipulate[
		
	        plotLabel = "Optim vs. Ref MSE "<>ToString[nDims]<>"D   integrandType = "<>integrandTypeLabel;

			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"WN_"<>integrandTypeLabel<>".dat"]),   #[[1]] <= 3^(powto+1) &];
			mseWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"Strat_"<>integrandTypeLabel<>".dat"]), #[[1]] <= 3^(powto+1) &];
			mseStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data],8}];

			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"PMJ02_"<>integrandTypeLabel<>".dat"]), #[[1]] <= 3^(powto+1) &];
			msePMJ02 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data]}];

			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"OwenPlus_"<>integrandTypeLabel<>".dat"]),  #[[1]] <= 3^(powto+1) &];
			mseOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data],8}];
			mseOwenPlusRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];

			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"MatBuiderMaxDepth_"<>integrandTypeLabel<>".dat"]), #[[1]] <= 3^(powto+1) &];
			mseMatBuiderMaxDepth = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data],9}];


			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Seq_PrevLevel.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOptimSeq = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			(*data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Seq_CurLevel.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOptim2 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)
			(*data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Seq_PrevLevel_step12octave.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOptimSeqstep12 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)
			(*data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Seq_CurLevel_step12octave.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOptim4 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];*)
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Pointsets_PrevLevel.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			(*mseOptimPointsets = Table[{data[[i,1]], data[[i,2]] },{i,Length[data]}];*)
			mseOptimPointsets = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];

			 
		    alldata = {mseWN, mseStrat, mseOwenPlus, mseMatBuiderMaxDepth, msePMJ02, mseOptimSeq, mseOptimPointsets} ;
	        (*legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat",  "OwenMaxDepth", "MatBuiderMaxDepth", "MB++ Seq", "MB++ Pointsets"} ] ];*)
	        legends = {"WN", "Strat",  "OwenMaxDepth", "MatBuiderMaxDepth", "PMJ02", "MB++ Seq", "MB++ Pointsets"} ;
	        
	        
			p = ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[2]},
							{Blue,AbsoluteThickness[2]},
							{Black,AbsoluteThickness[2]},
							{Orange,AbsoluteThickness[4]},
							{Darker@Green,AbsoluteThickness[4]},
							
							{Red,Dashed,AbsoluteThickness[3], Dashed},
							{Blue,Dashed,AbsoluteThickness[3]},
							{Darker@Green,AbsoluteThickness[2], Dashed},
							{Cyan,AbsoluteThickness[3], Dashed},
							{Gray,Dashed,AbsoluteThickness[3], Dashed}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[3^pow,{pow,1,10,1}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "MSE", fontSz] }
		           		,ImageSize -> 3/2 {1024,1024}
		            	(*,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ mseOwenPlusRaw), Min @@ (second /@ mseOwenPlusRaw) }} *)(*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,PlotRange->{{3^powfrom,3^powto}, Automatic } (*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[3^pow,{pow,1,8,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel, Bold, 24] 
		            ]
		            
			(*,Control[{{optimTypeL2OptimisationLabel,"MSEOptimisationSoftEllipses"},{"L2Optimisation","MSEOptimisationHeaviside","MSEOptimisationSoftEllipses"
				(*,"MSEOptimisationHardRectangles","MSEOptimisationSoftRectangles","MSEOptimisationHardEllipses"*)} } ]*)
			,Control[{{integrandTypeLabel,"SoftEllipses"},{"Heaviside", "SoftEllipses"	(*, "Ellipses", "Rectangles", "SoftEllipses_noRot" *)}}]
         ]
			(*Export["p_MSE.pdf",p];
			p//Print;*)
			

     ] (* showstdOptimMSE *)
     
     
(*(*showstdOptimMSE[] :=
    Module[ {powfrom,powto,powstep,kPlusMinus,data,plotLabel,legends,alldata},
    	consecutiveFlag = False;
		fontSz = 14;
		kPlusMinus = 1;
    	{powfrom,powto,powstep} = {2,8,1}; (* powers of 3 *)

		nDims = 2;
		(*integrandTypeLabel = "Heaviside";*)
		
		Manipulate[
	        plotLabel = "Optim vs. Ref MSE "<>ToString[nDims]<>"D   integrandType = "<>integrandTypeLabel;

			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"WN_"<>integrandTypeLabel<>".dat"]),  3^powfrom  <= #[[1]] <= 3^powto &];
			mseWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"Strat_"<>integrandTypeLabel<>".dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"OwenPure_"<>integrandTypeLabel<>".dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			mseOwen01PureRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"OwenPlus_"<>integrandTypeLabel<>".dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			mseOwenPlusRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"MatBuider_"<>integrandTypeLabel<>".dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseMatBuider = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>"MatBuiderMaxDepth_"<>integrandTypeLabel<>".dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseMatBuiderMaxDepth = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];


			data = Select[(Drop[#,1]& @ Import["data_MSE/"<>ToString[nDims]<>"D/"<>integrandTypeLabel<>"/"<>optimTypeL2OptimisationLabel<>"_"<>integrandTypeLabel<>"_Seq_CurLevel.dat"]), 3^powfrom  <= #[[1]] <= 3^powto &];
			mseOptim1 = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			 
		    alldata = {mseWN, mseStrat, mseOwenPlus, mseMatBuiderMaxDepth, mseOptim1} ;
	        legends = Join[ StringJoin[#, (" dims "<>Switch[nDims,2,"01",3,"012",4,"0123"])] & /@ Join[{"WN", "Strat",  "OwenPlus32", "MatBuiderMaxDepth"} ] ];
	        
	        
			ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz]& /@ legends}
						,PlotStyle -> {
							{Green,AbsoluteThickness[2]},
							{Blue,AbsoluteThickness[2]},
							{Black,AbsoluteThickness[2]},
							{Orange,AbsoluteThickness[2]},
							{Darker@Green,AbsoluteThickness[2]},
							{Cyan,AbsoluteThickness[2]},
							{Red,Dashed,AbsoluteThickness[2]},
							{Blue,Dashed,AbsoluteThickness[2]},
							{Gray,Dashed,AbsoluteThickness[2]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[3^pow,{pow,1,10,1}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "MSE", fontSz] }
		           		,ImageSize -> {1024,1024}
		            	(*,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ mseOwenPlusRaw), Min @@ (second /@ mseOwenPlusRaw) }} *)(*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,PlotRange->{{3^powfrom,3^powto}, Automatic } (*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[3^pow,{pow,1,8,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	,PlotLabel -> Style[ plotLabel, Bold, 24] 
		            ]			
			,Control[{{optimTypeL2OptimisationLabel,"MSEOptimisationSoftEllipses"},{"L2Optimisation","MSEOptimisationHeaviside","MSEOptimisationSoftEllipses"
				(*,"MSEOptimisationHardRectangles","MSEOptimisationSoftRectangles","MSEOptimisationHardEllipses"*)} } ]
			,Control[{{integrandTypeLabel,"SoftEllipses"},{"Heaviside", "SoftEllipses"	(*, "Ellipses", "Rectangles", "SoftEllipses_noRot" *)}}]
         ]
     ] (* showstdOptimMSE *)
     *)*)
(*======================================= prepSoftEllipses2D[] & prepHeavisideND[] =======================================*)
getHeavisideND[pt_,{mu_,normVector_}] := If[(pt-mu).normVector > 0, 1, 0]
getMultivariateND[pt_,{mu_,mxCInv_},mulFactor_:1] := Quiet[1./mulFactor Exp[-.5 (pt-mu).mxCInv.(pt-mu)]]

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
prepSoftEllipses2D[0]
prepSoftEllipses2D[1]

*)
prepSoftEllipses2D[setNo_:1, innbatches_:9] :=
    Module[ {},
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
        nDims = 2;
     	maxtime = 10;
        dir = "integrands/";
        If[ !FileExistsQ[dir], CreateDirectory[dir] ];
		{precision,maxRecursion} = {20,10000};

		nbatches = innbatches; 
		nmus1D = 256;
		(* Nintegrands = nbatches*nmus1D*nmus1D = 256K == 524288 *)

		integrandTypeLabel = "SoftEllipses";
		suffix = integrandTypeLabel<>ToString[nDims]<>"D"<>"_setNo"<>ToString[setNo];
        resfname = dir<>suffix<>".dat";
		hppfname = dir<>suffix<>".hpp";		
        Print["prepSoftEllipses2D" -> hppfname];
		cpptype = "t_GaussianStruct2D" ;
		varName = "tab_SoftEllipses2D" ;
        res = {};
        Do[
          	partial = RandomSample @ (Flatten[#,1]& @ (Parallelize @  Table[
		       	While[True,
						rotmx = RandomVariate[CircularRealMatrixDistribution[nDims], 1][[1]];
						sigma = {RandomReal[{.1, 1}],RandomReal[{.1, 1}]} / 2.;
						mu = {(ixmu-1+(RandomReal[]))/nmus1D,(iymu-1+(RandomReal[]))/nmus1D};
						sigmamx = sigma IdentityMatrix[nDims];
		    			mxC =  rotmx.sigmamx;
	        			mxCInv = Inverse[mxC];
	        			mxCInv = T[mxCInv] . mxCInv;
						Off[NIntegrate::slwcon];
						Off[NIntegrate::eincr];
						Off[NIntegrate::precw];
						Off[NIntegrate::maxp];
						Off[NIntegrate::inumr];
						Off[General::stop];
						Off[NIntegrate`SymbolicPiecewiseSubdivision::maxpwc];
						integral = (NIntegrate[getMultivariateND[Table[x[i],{i,nDims}],{mu,mxCInv}], ## , MaxRecursion->maxRecursion, 
								PrecisionGoal->precision, WorkingPrecision->precision, AccuracyGoal->precision] & @@ Table[{x[i],0,1},{i,nDims}]) ;
					Print[suffix -> mf[{{ibatch,nbatches},{ixmu,nmus1D},{iymu,nmus1D}}] -> mf[{mu,sigma}] -> integral];
					If[integral < eps, Print["Bad trial " -> suffix -> mf[{{ibatch,nbatches},{ixsigma,iysigma},{ixmu,iymu}}] -> mf[{mu,sigma}]  -> integral] ];
					If[integral > eps, Break[] ];
	        	];
				Flatten@{integral,mu,mxCInv}
	        ,{ixmu,nmus1D},{iymu,nmus1D}]) );
	        res = Join[res,partial];
			finalLength = Length[res];
			Put[(CForm /@ #) & /@ SetPrecision[res,precision], resfname]; (* e^-10 rather than 2.5*^-10 *) 
			Print["output into ",hppfname];		
	        Run["echo ' "<>cpptype<>" "<>varName<>"["<>ToString[finalLength]<>"] = ' > "<>hppfname ];
	        Run["cat "<> resfname<>" >> "<>hppfname];
	        Run["echo ';' >> "<>hppfname];       
	    ,{ibatch,nbatches}];
        DeleteFile[resfname];
    ] (* prepSoftEllipses2D *)

(*prepSoftEllipses2D[setNo_:1] :=
    Module[ {},
        nDims = 2;
     	maxtime = 10;
        dir = "integrands/";
        If[ !FileExistsQ[dir], CreateDirectory[dir] ];
		{precision,maxRecursion} = {20,10000};

		nsigmas1D = 4; 
		nmus1D = 128;
		nIntegrands = nsigmas1D^2 nmus1D^2;
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];

		integrandTypeLabel = "SoftEllipses";
		suffix = integrandTypeLabel<>"_setNo"<>ToString[setNo];
		hppsuffix = integrandTypeLabel<>ToString[nDims]<>"D"<>"_setNo"<>ToString[setNo];
		cppsuffix = "t_GaussianStruct2D" ;
		varName = "tab_SoftEllipses2D" ;
		nsigmas = nsigmas1D^2;
        res = Flatten[#,1]& @ (Table[
		    {ixsigma,iysigma} = 1+{Quotient[(isigma-1),nsigmas1D],Mod[(isigma-1),nsigmas1D]};	(* 1+ to skip too small sigmas *)
          	partial = Flatten[#,1]& @ (Parallelize @  Table[
		       	While[True,
						rotmx = RandomVariate[CircularRealMatrixDistribution[nDims], 1][[1]];
						sigma = {(ixsigma+(RandomReal[]))/nsigmas1D,(iysigma+(RandomReal[]))/nsigmas1D} / 2.;
						mu = {(ixmu-1+(RandomReal[]))/nmus1D,(iymu-1+(RandomReal[]))/nmus1D};
						sigmamx = sigma IdentityMatrix[nDims];
		    			mxC =  rotmx.sigmamx;
	        			mxCInv = Inverse[mxC];
	        			mxCInv = T[mxCInv] . mxCInv;
						Off[NIntegrate::slwcon];
						Off[NIntegrate::eincr];
						Off[NIntegrate::precw];
						Off[NIntegrate::maxp];
						Off[NIntegrate::inumr];
						Off[General::stop];
						Off[NIntegrate`SymbolicPiecewiseSubdivision::maxpwc];
						integral = (NIntegrate[getMultivariateND[Table[x[i],{i,nDims}],{mu,mxCInv}], ## , MaxRecursion->maxRecursion, 
								PrecisionGoal->precision, WorkingPrecision->precision, AccuracyGoal->precision] & @@ Table[{x[i],0,1},{i,nDims}]) ;
					Print[suffix -> mf[{{isigma,nsigmas},{ixmu,nmus1D},{iymu,nmus1D}}] -> mf[{mu,sigma}] -> integral];
					If[integral < eps, Print["Bad trial " -> suffix -> mf[{{isigma,nsigmas},{ixsigma,iysigma},{ixmu,iymu}}] -> mf[{mu,sigma}]  -> integral] ];
					If[integral > eps, Break[] ];
	        	];
				Flatten@{integral,mu,mxCInv}
	        ,{ixmu,nmus1D},{iymu,nmus1D}]);
	        partial
        ,{isigma,nsigmas}]);
        res = RandomSample @ res;
		finalLength = Length[res];
        resfname = dir<>suffix<>".dat";
        Print["prepIntegrands2D" -> finalLength -> hppfname];
		alldata = (Flatten/@res);
		hppfname = dir<>hppsuffix<>".cpp";		
		Put[(CForm /@ #) & /@ SetPrecision[alldata,precision], resfname]; (* e^-10 rather than 2.5*^-10 *) 
		Print["output into ",hppfname];		
        Run["echo ' "<>cppsuffix<>" "<>varName<>"["<>ToString[finalLength]<>"] = ' >> "<>hppfname ];
        Run["cat "<> resfname<>" >> "<>hppfname];
        Run["echo ';' >> "<>hppfname];       
        DeleteFile[resfname];
    ] (* prepSoftEllipses2D *)
*)
(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
prepHeavisideND[91, 9]
prepHeavisideND[92, 9]

*)

prepHeavisideND[setNo_:1, innbatches_:9] :=
    Module[ {},
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
        nDims = 2;
     	maxtime = 10;
        dir = "integrands/";
        If[ !FileExistsQ[dir], CreateDirectory[dir] ];
		{precision,maxRecursion} = Switch[nDims
			,1,{20,10000}
			,2,{18,10000}
			,3,{18,10000}
			,4,{18,10000}
			,5,{14,10000}
			,6,{11,10000}
			,_,{10,10000}
		];

		nbatches = innbatches; 
		nmus1D = 256;
		(* Nintegrands = nbatches*nmus1D*nmus1D = 256K == 524288 *)
		batchsz = nmus1D^2;
		nIntegrands = nbatches * batchsz;

		integrandTypeLabel = "Heaviside";
		suffix = integrandTypeLabel<>"_setNo"<>ToString[setNo];
		hppsuffix = integrandTypeLabel<>ToString[nDims]<>"D"<>"_setNo"<>ToString[setNo];
		cpptype = "t_Heaviside"<>ToString[nDims]<>"D" ;
		varName = "tab_Heaviside"<>ToString[nDims]<>"D" ;
        resfname = dir<>suffix<>".dat";
		hppfname = dir<>suffix<>".hpp";		

        res = {};
        Do[
          	partial = RandomSample @ (Flatten[#,1]& @ (Parallelize @  Table[
		       	While[True,
		    		muDiscotinuity = {(ixmu-1+(RandomReal[]))/nmus1D,(iymu-1+(RandomReal[]))/nmus1D};
					normVector = If[nDims == 2,
						i = (iymu-1)*nmus1D + (ixmu-1);
		    			j = (ibatch-1)*batchsz + (i-1);
		    			alpha = 2*PI*(j+RandomReal[])/nIntegrands;
		    			{Cos[alpha], Sin[alpha]}
					,(*ELSE*)
						getUniformDirsND[nDims]
					];
					Off[NIntegrate::slwcon];
					Off[NIntegrate::eincr];
					Off[NIntegrate::precw];
					Off[NIntegrate::maxp];
					Off[NIntegrate::inumr];
					Off[General::stop];
					Off[NIntegrate`SymbolicPiecewiseSubdivision::maxpwc];
					integral = If[nDims == 3,
						TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ##, Method->"LocalAdaptive", PrecisionGoal->precision ] & @@ Table[{x[i],0,1},{i,nDims}], maxtime,0]
					,(*ELSE*)
						If[nDims < 8,
							TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ## , MaxRecursion->maxRecursion, 
								PrecisionGoal->precision, WorkingPrecision->precision, AccuracyGoal->precision] & @@ Table[{x[i],0,1},{i,nDims} ], maxtime,0]
						,(*ELSE*)
							TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ##] & @@ Table[{x[i],0,1},{i,nDims} ], maxtime,0]
						]
					];
					Print[suffix -> mf[{{ibatch,nbatches},{ixmu,nmus1D},{iymu,nmus1D}}] -> mf[{muDiscotinuity,normVector}] -> integral];
					If[integral < eps, Print["Bad trial " -> suffix -> mf[{{ibatch,nbatches},{ixsigma,iysigma},{ixmu,iymu}}] -> mf[{muDiscotinuity,normVector}]  -> integral] ];
					If[integral > eps, Break[] ];
	        	];
				Flatten@{integral,muDiscotinuity,normVector}
	        ,{ixmu,nmus1D},{iymu,nmus1D}]) );
	        res = Join[res,partial];
			finalLength = Length[res];
			Put[(CForm /@ #) & /@ SetPrecision[res,precision], resfname]; (* e^-10 rather than 2.5*^-10 *) 
			Print["output into ",hppfname];		
	        Run["echo ' "<>cpptype<>" "<>varName<>"["<>ToString[finalLength]<>"] = ' > "<>hppfname ];
	        Run["cat "<> resfname<>" >> "<>hppfname];
	        Run["echo ';' >> "<>hppfname];       
	    ,{ibatch,nbatches}];
        DeleteFile[resfname];
    ] (* prepHeavisideND *)
(*
prepHeavisideND[innDims_:2, setNo_:1] :=
    Module[ {nIntegrands,nDims,suffix,maxtime,dir,precision,maxRecursion,batchsz,nbatches,res1024,res,trial,finalLength,resfname,alldata,hppfname,integral,muDiscotinuity,normVector,alpha,j,
    	integrandTypeLabel,hppsuffix,cppsuffix,varName},
    	nIntegrands = 9 256 256;
        If[$ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
        nDims = innDims;

		integrandTypeLabel = "Heaviside";
		suffix = integrandTypeLabel<>"_setNo"<>ToString[setNo];
		hppsuffix = integrandTypeLabel<>ToString[nDims]<>"D"<>"_setNo"<>ToString[setNo];
		cppsuffix = "t_Heaviside"<>ToString[nDims]<>"D" ;
		varName = "tab_Heaviside"<>ToString[nDims]<>"D" ;

     	maxtime = If[nDims <= 10, 10, 3600];
        dir = "integrands/";
        If[ !FileExistsQ[dir], CreateDirectory[dir] ];
		{precision,maxRecursion} = Switch[nDims
			,1,{20,10000}
			,2,{18,10000}
			,3,{18,10000}
			,4,{18,10000}
			,5,{14,10000}
			,6,{11,10000}
			,_,{10,10000}
		];
		batchsz = 1024;
		nbatches = nIntegrands/batchsz;
        res = {};
        Do[
          	res1024 = Parallelize @ (Table[
          		trial = 0;
		       	While[True,
		       		trial++;
		    		(*muDiscotinuity = Table[.5,{nDims}] + .1 getUniformDirsND[nDims];*)
		    		(*muDiscotinuity = Table[.5 - (RandomReal[]-.5)/2,{nDims}];*)
		    		muDiscotinuity = Table[RandomReal[],{nDims}];
					normVector = If[nDims == 2,
		    			j = (ibatch-1)*batchsz + (i-1);
		    			alpha = 2*PI*(j+RandomReal[])/nIntegrands;
		    			{Cos[alpha], Sin[alpha]}
					,(*ELSE*)
						getUniformDirsND[nDims]
					];
					Off[NIntegrate::slwcon];
					Off[NIntegrate::eincr];
					Off[NIntegrate::precw];
					Off[NIntegrate::maxp];
					Off[General::stop];
					Off[NIntegrate`SymbolicPiecewiseSubdivision::maxpwc];
					(*integral = (*TimeConstrained[#,maxtime,0] & @*) *)
					integral = If[nDims == 3,
						TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ##, Method->"LocalAdaptive", PrecisionGoal->precision ] & @@ Table[{x[i],0,1},{i,nDims}], maxtime,0]
					,(*ELSE*)
						If[nDims < 8,
							TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ## , MaxRecursion->maxRecursion, 
								PrecisionGoal->precision, WorkingPrecision->precision, AccuracyGoal->precision] & @@ Table[{x[i],0,1},{i,nDims} ], maxtime,0]
						,(*ELSE*)
							TimeConstrained[ NIntegrate[getHeavisideND[Table[x[i],{i,nDims}],{muDiscotinuity,normVector}], ##] & @@ Table[{x[i],0,1},{i,nDims} ], maxtime,0]
						]
					];
					(*Print["trial = ",trial -> {i,j} -> {muDiscotinuity,normVector} -> integral];*)
					If[integral > eps, Break[] ];
	        	];
				Print[suffix -> ibatch,"/",nbatches -> i,"/",batchsz -> integral];
		        {integral,muDiscotinuity,normVector}
	        ,{i,batchsz}]);
	        res = Join[res, res1024];
        ,{ibatch,nbatches}];
		finalLength = Length[res];
        resfname = dir<>suffix<>".dat";
        Print["prepHeavisideND" -> finalLength -> resfname];
 		alldata = RandomSample @ (Flatten/@res);
		hppfname = dir<>hppsuffix<>".hpp";		
		Put[(CForm /@ #) & /@ SetPrecision[alldata,precision], resfname]; (* e^-10 rather than 2.5*^-10 *) 
		Print["output into ",hppfname];		
        Run["echo ' "<>cppsuffix<>" "<>varName<>"["<>ToString[finalLength]<>"] = ' >> "<>hppfname ];
        Run["cat "<> resfname<>" >> "<>hppfname];
        Run["echo ';' >> "<>hppfname];       
        DeleteFile[resfname];
    ] (* prepHeavisideND *)
*)
getUniformDirsND[nDims_:6]:= 
Module[{v0 = Table[1.,{nDims}]/Sqrt[nDims]},
	RandomVariate[CircularRealMatrixDistribution[nDims]].v0
] (* getUniformDirsND *)

    
(*
gitpull
math
<<uniformity/uniformity.m
prepSoftEllipsesND[2, 16 1024, 1024, 10]

prepSoftEllipsesND[2, 16 1024, 1024, 0]
prepSoftEllipsesND[2, 16 1024, 1024, 1]

prepSoftEllipsesND[2, 16 1024, 1024, 2]
prepSoftEllipsesND[2, 16 1024, 1024, 3]
prepSoftEllipsesND[2, 16 1024, 1024, 4]
prepSoftEllipsesND[2, 16 1024, 1024, 5]

prepSoftEllipsesND[innDims_:2, innIntegrands_:16 1024, inbatchsz_:1024, setno_:0, dbg_:False] :=
    Module[ {},
    	smallValue = 1/1000.;
        If[ $ProcessorCount != 10 && Length[Kernels[]] < $ProcessorCount*2, LaunchKernels[$ProcessorCount*2] ];
    	nIntegrands = innIntegrands;
        nDims = innDims;
		suffix = "SoftEllipses2D_nIntegrands"<>ToString[nIntegrands]<>"_setno"<>ToString[setno];
		cppsuffix = "t_GaussianStruct2D";
     	maxtime = If[nDims <= 10, 10, 3600];
        dir = "integrands/";
        If[ !FileExistsQ[dir], CreateDirectory[dir] ];
		{precision,maxRecursion} = {20,10000};
		batchsz = Min[1024,inbatchsz];
		nbatches = nIntegrands/batchsz;
        res = {};
        fIntegrand := getMultivariateND;
        Do[
	        k1k2Tab = smallValue + 1/2. RandomSample @ getStratND[nDims, batchsz];
          	partial = Parallelize @ (Table[
		       	While[True,
					ktab = k1k2Tab[[inbatch]];
					{mu,mxCInv} = getRandMuAndMxCInvNDEllipses[nDims,ktab];
					Off[NIntegrate::slwcon];
					Off[NIntegrate::eincr];
					Off[NIntegrate::precw];
					Off[NIntegrate::maxp];
					Off[General::stop];.
					Off[NIntegrate`SymbolicPiecewiseSubdivision::maxpwc];
					integral = (*TimeConstrained[#,maxtime,0] & @*) If[nDims == 3,
						precision = 17;
						(NIntegrate[fIntegrand[Table[x[i],{i,nDims}],{mu,mxCInv}], ##, Method->"LocalAdaptive", PrecisionGoal->precision ] & @@ Table[{x[i],0,1},{i,nDims}])
					,(*ELSE*)
						If[nDims < 8,
							(NIntegrate[fIntegrand[Table[x[i],{i,nDims}],{mu,mxCInv}], ## , MaxRecursion->maxRecursion, 
								PrecisionGoal->precision, WorkingPrecision->precision, AccuracyGoal->precision] & @@ Table[{x[i],0,1},{i,nDims}])
						,(*ELSE*)
							(NIntegrate[fIntegrand[Table[x[i],{i,nDims}],{mu,mxCInv}], ##] & @@ Table[{x[i],0,1},{i,nDims}])
						]
					];
					If[integral < eps, Print["Bad trial " -> suffix -> ibatch,"/",nbatches -> inbatch,"/",batchsz -> mf[{mu,mxCInv,ktab}] -> integral] ];
					If[integral > eps, Break[] ];
	        	];
				Print[suffix -> ibatch,"/",nbatches -> inbatch,"/",batchsz -> {mu,mxCInv} -> integral];
		        {integral,mu,mxCInv}
	        ,{inbatch,batchsz}]);
	        res = Join[res, partial];
			finalLength = Length[res];
	        resfname = dir<>suffix<>".dat";
	        Print["prepSoftEllipsesND" -> finalLength -> resfname];
			alldata = (Flatten/@res);
			hppfname = dir<>suffix<>".cpp";		
			Put[(CForm /@ #) & /@ SetPrecision[alldata,precision], resfname]; (* e^-10 rather than 2.5*^-10 *)  
			Print["output into ",resfname," and ",hppfname];		
	        Run["echo ' #include \"../Integration/Integration.h\" ' > "<>hppfname ];
	        Run["echo ' "<>cppsuffix<>" tab_SoftEllipses2D["<>ToString[finalLength]<>"] = ' >> "<>hppfname ];
	        (*Run["echo ' "<>cppsuffix<>" tab_"<>suffix<>"[16384] = ' >> "<>hppfname ];*)
	        Run["cat "<> resfname<>" >> "<>hppfname];
	        Run["echo ';' >> "<>hppfname];       
        ,{ibatch,nbatches}]
    ] (* prepSoftEllipsesND *)
*)


makeOctavesBaseN[powParams_:{0,10,1},base_:3] :=
    Module[ {},
        {powfrom,powto,powstep} = powParams;
        tab = Union @ Table[
        		n = Round[base^ipow];
        		(*Print[ipow -> n];*)
        		n
        	,{ipow,powfrom,powto,powstep}];
        Print[tab];
        Export["tab.dat", {tab}];
        tab
    ]
    
    
(*=========================== lois MCQMC July 2022 *)
(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
loismakeMatBuilderMatrices["net_t0"]
loismakeMatBuilderMatrices["net_t1"]
loismakeMatBuilderMatrices["net_t2"]
loismakeMatBuilderMatrices["net_t3"]
loismakeMatBuilderMatrices["net_t4"]
*)

loismakeMatBuilderMatrices[basename_:"net_t0",ntrials_:256] :=
    Module[ {},
    	If[ !FileExistsQ["lois/MatBuilder_matrices/"], CreateDirectory["lois/MatBuilder_matrices/"] ];
    	nlevels = 19;
		fname = "lois/profiles/tvalue/"<>basename<>".txt";
		Parallelize @ Do[
			execString = "matbuilder -i "<>fname<>" -o lois/MatBuilder_matrices/"<>basename<>"_"<>i2s[itrial]<>".dat --seed "<>ToString[RandomInteger[2^16] ]<>" > /dev/null";
        	returnCode = Run[execPrefix<>execString];
        	Print[execString -> returnCode];
(*			mxTab = readMatBuilderMatrix["MatBuilder_matrices/"<>basename<>"_"<>i2s[i]<>".dat"];
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
			Export["MatBuilder_matrices/2D_0m2net_"<>i2s[i]<>"_inv.dat", Flatten[#,1]& @ mxInvTab ];*)
		,{itrial,1,ntrials}];
    ] (* loismakeMatBuilderMatrices *)


(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
loismakeL2discrepancy["net_t0"]
loismakeL2discrepancy["net_t1"]
loismakeL2discrepancy["net_t2"]
loismakeL2discrepancy["net_t3"]
loismakeL2discrepancy["net_t4"]

gitpull
math
<<TileBasedOptim/TileBasedOptim.m
loismakeL2discrepancy["net_t3",{1,10,1},{1,10}]
*)

loismakeL2discrepancy[basename_:"net_t0", octaves_:{1,10,1}, setFromTo_:{1,256}, innDims_:2, dbg_:False] :=
    Module[ {},
       	header = "#Nbpts	#Mean	#Var	#Min	#Max	#VOID	#VOID	#NbPtsets	#VOID\n";
    	nDims = innDims;
    	base = 3;
		owenFlag = True;
        dtab = {};
        ordinalAbsoluteMax = 3^10;
   	    resFname = basename<>".dat";	    
   	    {setFrom,setTo} = setFromTo;
   	    {powfrom,powto,powstep} = octaves; (* powers of 3 *)
		dataDiscrepancy = {};
		dirDiscrepancy = "data_L2Discrepancy/"<>ToString[nDims]<>"D/";
        If[ !FileExistsQ[dirDiscrepancy], CreateDirectory[dirDiscrepancy] ];
        If[ !FileExistsQ["tmp/"], CreateDirectory["tmp/"] ];
        Do[
			npts = base^pow;
	        DiscrepancyTab = (Parallelize @ Table[
	       			mxfname = "lois/MatBuilder_matrices/"<>basename<>"_"<>i2s[setNo]<>".dat";
	       			If[!FileExistsQ[mxfname], 
	       			Print[mxfname," does not exist."]; 
	       			Nothing
	       		,(*ELSE*)
	       			depth = 19;
	       			seed = RandomInteger[2^16];
	       			pts = getMatBuiderPtsND[base^pow, mxfname, owenFlag, depth, nDims, base, seed ];
		        	L2discrepancy = getL2discrepancy[pts,"",nDims]; 
		        	Print[iOrdinalAbsolute, " ", resFname  -> L2discrepancy];
					{npts,L2discrepancy}
	       		]
        	,{setNo,setFrom,setTo}]);
	 		DiscrepancyMean = Mean @ DiscrepancyTab;
	 		DiscrepancyVariance = If[Length[DiscrepancyTab] <= 1, 0 , Variance @ (Last /@ DiscrepancyTab)];
	 		{DiscrepancyMin,DiscrepancyMax} = {Min@(Last /@ DiscrepancyTab), Max@(Last /@ DiscrepancyTab)};
	 		AppendTo[dataDiscrepancy,Flatten @ {DiscrepancyMean,DiscrepancyVariance,DiscrepancyMin,DiscrepancyMax,0,0,setTo-setFrom+1,0}];	
			Export[dirDiscrepancy<>resFname,header,"TEXT"];
			Export["tmp/tmpdat"<>pid<>".dat",dataDiscrepancy];
			Run["cat tmp/tmpdat"<>pid<>".dat >> "<>dirDiscrepancy<>resFname];
			Print[dirDiscrepancy<>resFname, " written."];
	    ,{pow,powfrom,powto,powstep}];
		Run["rm tmp/tmpdat"<>pid<>".dat"];
   ] (* loismakeL2discrepancy *)

(*
gitpull
math
<<TileBasedOptim/TileBasedOptim.m
loisshowL2discrepancy[]
loisshowL2discrepancy["net_t0"]
loisshowL2discrepancy["net_t1"]
loisshowL2discrepancy["net_t2"]
loisshowL2discrepancy["net_t3"]
loisshowL2discrepancy["net_t4"]
*)
loisshowL2discrepancy[lbl_:"", innDims_:2, dbg_:False] :=
    Module[ {},
		fontSz = 36;
		kPlusMinus = .7;
    	{powfrom,powto,powstep} = {2,16,1};

    	nDims = innDims;
		dirL2discrepancy = "data_L2discrepancy/"<>ToString[nDims]<>"D/";
        dtab = {};
        nptsMax = 653; (*3^6;*)
 		(*Print[L2discrepancyOptim//mf];*)
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"WN_L2Discrepancy.dat"]);
			L2discrepancyWN = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"Strat_L2Discrepancy.dat"]);
			L2discrepancyStrat = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"OwenPure_L2Discrepancy.dat"]);
			L2discrepancyOwen01Pure = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data],4}];
			L2discrepancyOwen01PureRaw = Table[{data[[i,1]],  data[[i,2]]},{i,Length[data]}];
			data = (Drop[#,1]& @ Import[dirL2discrepancy<>"OwenPlus_L2Discrepancy.dat"]);
			L2discrepancyOwenPlus = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,1,Length[data],4}];

			If[lbl =!= "",
	   	    	fname = dirL2discrepancy<>lbl<>".dat";
	        	If[ !FileExistsQ[fname], Print[fname," does not exist."]; Abort[]; ];
		       	data = (Drop[#,1]& @ Import[fname]);
	   			L2discrepancyMatBuilder = Table[{data[[i,1]], Around[ data[[i,2]], kPlusMinus Sqrt@data[[i,3]] ] },{i,Length[data]}];
	        	legends = {"WN", "Strattified", "Sobol", "Owen tree depth=32", "MatBuilder, "<>lbl};
	   		,(*ELSE*)
	   			L2discrepancyMatBuilder = {{0,0}};
	        	legends = {"White Noise", "Strat", "Owen", "OwenPlus32"};
			];

		    alldata = {L2discrepancyWN, L2discrepancyStrat, L2discrepancyOwen01Pure,  L2discrepancyOwenPlus, L2discrepancyMatBuilder} ;

	        plotLabel = "L2discrepancy "<>ToString[nDims]<>"D"; 
	        
			g = ListLogLogPlot[ alldata
						,PlotLegends -> Placed[#,{.3,.2}]& @  {Style[#,fontSz,Bold]& /@ legends}
						,PlotStyle -> {
							{Red,AbsoluteThickness[3]},
							{Blue,AbsoluteThickness[3]},
							{Black,AbsoluteThickness[3]},
							{Darker@Green,AbsoluteThickness[3]},
							
							{Orange,Dashed,AbsoluteThickness[12]},
							{Darker@Green,Dashed,AbsoluteThickness[10]},
							{Black,Dashed,AbsoluteThickness[2]},
							{Red,Dashed,AbsoluteThickness[2]},
							{Gray,Dashed,AbsoluteThickness[2]}
						}
						,Joined->True
		            	,FrameTicks->{{Automatic,None},{Table[2^pow,{pow,powfrom,powto,2}],Table[2^pow,{pow,powfrom,powto,2}]}}
			            ,FrameStyle->Directive[Black,20]
			            ,RotateLabel -> True
			            ,PlotMarkers->{{\[FilledCircle],5} }
			            ,Frame->True
		 	            ,FrameLabel-> {Style[ "Number of Samples", fontSz],Style[ "L2discrepancy", fontSz] }
		           		,ImageSize -> {1024,1024}
		            	(*,PlotRange->{{2^powfrom,2^powto},{Max @@ (second /@ L2discrepancyOwenPlusRaw), Min @@ (second /@ L2discrepancyOwenPlusRaw) }} *)(*{{4,2^powto},Automatic}*)	(* {{2^5,2^12},Automatic} *)
		            	,GridLines->{Table[2^pow,{pow,powfrom,powto,1}],None}
		            	,GridLinesStyle->Directive[Darker@Gray, Dashed]
		            	,AspectRatio->1
		            	,InterpolationOrder -> 1, IntervalMarkers -> "Bands", Sequence[PlotTheme -> "Scientific", PlotRange -> All]
		            	(*,PlotLabel -> Style[ plotLabel, Bold, 24] *)
		            ];
		    g//Print;	
		    Export["g_lois_"<>lbl<>".pdf", g];	

   ] (* loisshowL2discrepancy *)
