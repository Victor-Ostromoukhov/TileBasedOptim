(* fibo-hilbert.m
   V.O. version 2002/12/14
*)
 
(****************** params *******************)
SetOptions[Graphics, ImageSize -> { 2 1024, Automatic}];

mf := MatrixForm;

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
(*------------------------- end of constants -------------------------*)

demoFiboSFC[] :=
    Module[ {},
		tlst = {{type1,{0,0}, {{1,0},{0,1}}, {0,0}, {}} };
		tlst = subdivFiboSFCTiles @ tlst;
		tlst//mf//Print;
		Graphics[ getFiboSFCTilesGL[tlst] ]//Print;
]


subdivFiboSFCTiles[tlst_] :=
    Block[ {res={} },
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,fcode} = tlst[[ind]];
            Switch[tileType
              ,1, {t1,t3} = {z1,z3};
                  AppendTo[res,{type4,refPt,{oneoverphi v1,v2},samplingPt,AppendTo[fcode,0]}];
                  AppendTo[res,{type4,refPt+oneoverphi v1+v2,{ (mxRot270.v1),oneoverphi2 (mxRot270.v2)},samplingPt,AppendTo[fcode,1]}];
              ,2, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z1;
              ,3, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z3;
              ,4, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z3;
              ,5, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,6, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor+1);
            ];
    	,{ind,Length[tlst]}];
    	Return[res]
    ] (* subdivFiboSFCTiles *)

getFiboSFCTilesGL[tlst_,params_] :=
    Block[ {gl={},bortedStyle={Cyan,AbsoluteThickness[1]} },
    	Table[
			{tileType,refPt,{v1,v2},samplingPt,fcode} = tlst[[ind]];
			AppendTo[gl,{bortedStyle,Line[{refPt,refPt+v1,refPt+v1+v2,refPt+v2,refPt}] } ];
    	,{ind,Length[tlst]}];
    	Reurn[gl]
    ] (* getFiboSFCTilesGL *)
    
         (*cont = {z0,z1,z2,z3} = getTileShape[fig];
        If[ showdir,
            col = Magenta;
            Switch[tileType
              ,1, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z1;
              ,2, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z1;
              ,3, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z3;
              ,4, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z3;
              ,5, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,6, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z2;
              ,7, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+1);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,8, {t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor+1);
                  end = z2;
              ,9, {t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,10,{t1,t3} = {z1,z3};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,11,{t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
              ,12,{t1,t3} = {z3,z1};
                  du = (t1-z0)/tau^(scalefactor+2);
                  dv = (t3-z0)/tau^(scalefactor);
                  end = z2;
            ];
            If[ 1 <= tileType <= 4,
                AppendTo[gl,{{col,arr[{z0+du+dv,end-du+dv}]},Point[z0+du+dv],Point[end-du+dv]}];,(*ELSE*)
                AppendTo[gl,{{col,arr[{z0+du+dv,end-du-dv}]},Point[z0+du+dv],Point[end-du-dv]}];
            ];
        ]; (* If[showdir, *)
        If[ showRefPt,
            AppendTo[gl,{refPtCol,Point[refPt]}];
        ];
        If[ showSamplingPt,
            AppendTo[gl,{samplingPtCol,Point[getSamplingPt[fig]]}];
        ];
        If[ showTileShape,
            AppendTo[gl,{tileShapeCol,Thickness[tileShapeTh],Line@@{Append[cont,First[cont]]}}];
            AppendTo[border,{tileShapeCol,Thickness[tileShapeTh],Line@@{Append[cont,First[cont]]}}];
        ];
        If[ showThVal,
            AppendTo[gl,{GrayLevel[thval],Polygon@@{Append[cont,First[cont]]}}];
        ];
        If[ showThValTxt,
            AppendTo[gl,Text[ToString[thval],refPt,{-1,-1}] ]
        ];
        If[ labelOrdinalNumber,
            AppendTo[gl,Text[ToString[i],(z0+z2)/2,{-1,-1}] ]
        ];
        If[ labelTileTypes,
            AppendTo[gl,{Red,Text[ToString[inflationRules[[tileType,1]]],(z0+z2)/2,{1,1}]} ]
        ];
        Return[gl]*)

