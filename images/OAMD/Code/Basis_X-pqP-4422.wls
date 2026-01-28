<<FeynCalc`

(* Definitions *)

(* Charge conjugation matrix *)
(*Cm=I.GA[2].GA[0];*) (*Gamma4_E=Gamma0_M; I*Gamma2_E=Gamma2_M*)

(* respective projectors *)
\[Lambda]p[t_]:=(1+GS[t])/2;
\[Lambda]m[t_]:=(1-GS[t])/2;

(* Omega_r *)
\[Omega]p[t_]:={\[Lambda]p[t].GA5 ,\[Lambda]p[t]};
\[Omega]m[t_]:={\[Lambda]m[t].GA5 ,\[Lambda]p[t]};

C\[Omega]p[t_]:={FCCCT[GA5 ].FCCCT[\[Lambda]p[t] ] ,FCCCT[\[Lambda]p[t] ] };
C\[Omega]m[t_]:={FCCCT[GA5 ].FCCCT[\[Lambda]m[t] ] ,FCCCT[\[Lambda]p[t] ] };

(* Gammametrix *)
Gm[p_,q_,i_]:={1,GS[p],GS[q],GS[p].GS[q]}[[i]];

(* Set S,P,A,V *)
Sp[p_,q_,P_,i_,j_]:=Fpair[Gm[p,q,i].\[Omega]p[P][[1]], Gm[p,q,j].\[Omega]p[P][[2]]];
Sm[p_,q_,P_,i_,j_]:=Fpair[Gm[p,q,i].\[Omega]m[P][[1]], Gm[p,q,j].\[Omega]m[P][[2]]];

CSp[p_,q_,P_,i_,j_]:=Fpair[Calc[C\[Omega]p[P][[1]].FCCCT[Gm[p,q,i] ] ], Calc[C\[Omega]p[P][[2]]].FCCCT[Gm[p,q,j] ] ];
CSm[p_,q_,P_,i_,j_]:=Fpair[Calc[C\[Omega]m[P][[1]].FCCCT[Gm[p,q,i] ] ], Calc[C\[Omega]m[P][[2]]].FCCCT[Gm[p,q,j] ] ];

Pp[p_,q_,P_,i_,j_]:=Fpair[GA5.Gm[p,q,i].\[Omega]p[P][[1]], GA5.Gm[p,q,j].\[Omega]p[P][[2]]];
Pm[p_,q_,P_,i_,j_]:=Fpair[GA5.Gm[p,q,i].\[Omega]m[P][[1]], GA5.Gm[p,q,j].\[Omega]m[P][[2]]];

CPp[p_,q_,P_,i_,j_]:=Fpair[Calc[C\[Omega]p[P][[1]].FCCCT[Gm[p,q,i] ].FCCCT[GA5] ], Calc[C\[Omega]p[P][[2]].FCCCT[Gm[p,q,j] ].FCCCT[GA5] ] ];
CPm[p_,q_,P_,i_,j_]:=Fpair[Calc[C\[Omega]m[P][[1]].FCCCT[Gm[p,q,i] ].FCCCT[GA5] ], Calc[C\[Omega]m[P][[2]].FCCCT[Gm[p,q,j] ].FCCCT[GA5] ] ];

(* Base of 64 sets *)

X[p_,q_,P_] := Module[{array},
    array = ConstantArray[0, {4,4,2,2}]; 
    array[[1, 1, 1, 1]] = Sp[p, q, P, 1, 1]; array[[1, 1, 1, 2]] = Sm[p, q, P, 1, 1];
    array[[1, 2, 1, 1]] = Sp[p, q, P, 1, 2]; array[[1, 2, 1, 2]] = Sm[p, q, P, 1, 2];
    array[[1, 3, 1, 1]] = Sp[p, q, P, 1, 3]; array[[1, 3, 1, 2]] = Sm[p, q, P, 1, 3];
    array[[1, 4, 1, 1]] = Sp[p, q, P, 1, 4]; array[[1, 4, 1, 2]] = Sm[p, q, P, 1, 4];

    array[[1, 1, 2, 1]] = Pp[p, q, P, 1, 1]; array[[1, 1, 2, 2]] = Pm[p, q, P, 1, 1];
    array[[1, 2, 2, 1]] = Pp[p, q, P, 1, 2]; array[[1, 2, 2, 2]] = Pm[p, q, P, 1, 2];
    array[[1, 3, 2, 1]] = Pp[p, q, P, 1, 3]; array[[1, 3, 2, 2]] = Pm[p, q, P, 1, 3];
    array[[1, 4, 2, 1]] = Pp[p, q, P, 1, 4]; array[[1, 4, 2, 2]] = Pm[p, q, P, 1, 4];

    array[[2, 1, 1, 1]] = 1/Sqrt[3]*Sp[p, q, P, 3, 3]+1/Sqrt[3]*Sp[p, q, P, 4, 4]+1/Sqrt[3]*Pp[p, q, P, 2, 2];
    array[[2, 1, 1, 2]] = 1/Sqrt[3]*Sm[p, q, P, 3, 3]+1/Sqrt[3]*Sm[p, q, P, 4, 4]-1/Sqrt[3]*Pm[p, q, P, 2, 2];
    array[[2, 2, 1, 1]] = 1/Sqrt[3]*Sp[p, q, P, 3, 4]-1/Sqrt[3]*Sp[p, q, P, 4, 3]-1/Sqrt[3]*Pp[p, q, P, 2, 1];
    array[[2, 2, 1, 2]] = 1/Sqrt[3]*Sm[p, q, P, 3, 4]-1/Sqrt[3]*Sm[p, q, P, 4, 3]+1/Sqrt[3]*Pm[p, q, P, 2, 1];
    array[[2, 3, 1, 1]] = 1/Sqrt[3]*Sp[p, q, P, 3, 1]-1/Sqrt[3]*Sp[p, q, P, 4, 2]+1/Sqrt[3]*Pp[p, q, P, 2, 4];
    array[[2, 3, 1, 2]] = 1/Sqrt[3]*Sm[p, q, P, 3, 1]-1/Sqrt[3]*Sm[p, q, P, 4, 2]-1/Sqrt[3]*Pm[p, q, P, 2, 4];
    array[[2, 4, 1, 1]] = 1/Sqrt[3]*Sp[p, q, P, 3, 2]+1/Sqrt[3]*Sp[p, q, P, 4, 1]-1/Sqrt[3]*Pp[p, q, P, 2, 3];
    array[[2, 4, 1, 2]] = 1/Sqrt[3]*Sm[p, q, P, 3, 2]+1/Sqrt[3]*Sm[p, q, P, 4, 1]+1/Sqrt[3]*Pm[p, q, P, 2, 3];

    array[[2, 1, 2, 1]] = 1/Sqrt[3]*Pp[p, q, P, 3, 3]+1/Sqrt[3]*Pp[p, q, P, 4, 4]+1/Sqrt[3]*Sp[p, q, P, 2, 2];
    array[[2, 1, 2, 2]] = 1/Sqrt[3]*Pm[p, q, P, 3, 3]+1/Sqrt[3]*Pm[p, q, P, 4, 4]-1/Sqrt[3]*Sm[p, q, P, 2, 2];
    array[[2, 2, 2, 1]] = 1/Sqrt[3]*Pp[p, q, P, 3, 4]-1/Sqrt[3]*Pp[p, q, P, 4, 3]-1/Sqrt[3]*Sp[p, q, P, 2, 1];
    array[[2, 2, 2, 2]] = 1/Sqrt[3]*Pm[p, q, P, 3, 4]-1/Sqrt[3]*Pm[p, q, P, 4, 3]+1/Sqrt[3]*Sm[p, q, P, 2, 1];
    array[[2, 3, 2, 1]] = 1/Sqrt[3]*Pp[p, q, P, 3, 1]-1/Sqrt[3]*Pp[p, q, P, 4, 2]+1/Sqrt[3]*Sp[p, q, P, 2, 4];
    array[[2, 3, 2, 2]] = 1/Sqrt[3]*Pm[p, q, P, 3, 1]-1/Sqrt[3]*Pm[p, q, P, 4, 2]-1/Sqrt[3]*Sm[p, q, P, 2, 4];
    array[[2, 4, 2, 1]] = 1/Sqrt[3]*Pp[p, q, P, 3, 2]+1/Sqrt[3]*Pp[p, q, P, 4, 1]-1/Sqrt[3]*Sp[p, q, P, 2, 3];
    array[[2, 4, 2, 2]] = 1/Sqrt[3]*Pm[p, q, P, 3, 2]+1/Sqrt[3]*Pm[p, q, P, 4, 1]+1/Sqrt[3]*Sm[p, q, P, 2, 3];

    array[[3, 1, 1, 1]] =  2/Sqrt[6]*Sp[p, q, P, 3, 3]-1/Sqrt[6]*Sp[p, q, P, 4, 4]-1/Sqrt[6]*Pp[p, q, P, 2, 2];
    array[[3, 1, 1, 2]] =  2/Sqrt[6]*Sm[p, q, P, 3, 3]-1/Sqrt[6]*Sm[p, q, P, 4, 4]+1/Sqrt[6]*Pm[p, q, P, 2, 2];
    array[[3, 2, 1, 1]] =  1/Sqrt[6]*Sp[p, q, P, 3, 4]-1/Sqrt[6]*Sp[p, q, P, 4, 3]+2/Sqrt[6]*Pp[p, q, P, 2, 1];
    array[[3, 2, 1, 2]] =  1/Sqrt[6]*Sm[p, q, P, 3, 4]-1/Sqrt[6]*Sm[p, q, P, 4, 3]-2/Sqrt[6]*Pm[p, q, P, 2, 1];
    array[[3, 3, 1, 1]] =  2/Sqrt[6]*Sp[p, q, P, 3, 1]+1/Sqrt[6]*Sp[p, q, P, 4, 2]-1/Sqrt[6]*Pp[p, q, P, 2, 4];
    array[[3, 3, 1, 2]] =  2/Sqrt[6]*Sm[p, q, P, 3, 1]+1/Sqrt[6]*Sm[p, q, P, 4, 2]+1/Sqrt[6]*Pm[p, q, P, 2, 4];
    array[[3, 4, 1, 1]] =  2/Sqrt[6]*Sp[p, q, P, 4, 1]-1/Sqrt[6]*Sp[p, q, P, 3, 2]+1/Sqrt[6]*Pp[p, q, P, 2, 3];
    array[[3, 4, 1, 2]] =  2/Sqrt[6]*Sm[p, q, P, 4, 1]-1/Sqrt[6]*Sm[p, q, P, 3, 2]-1/Sqrt[6]*Pm[p, q, P, 2, 3];

    array[[3, 1, 2, 1]] =  2/Sqrt[6]*Pp[p, q, P, 3, 3]-1/Sqrt[6]*Pp[p, q, P, 4, 4]-1/Sqrt[6]*Sp[p, q, P, 2, 2];
    array[[3, 1, 2, 2]] =  2/Sqrt[6]*Pm[p, q, P, 3, 3]-1/Sqrt[6]*Pm[p, q, P, 4, 4]+1/Sqrt[6]*Sm[p, q, P, 2, 2];
    array[[3, 2, 2, 1]] =  1/Sqrt[6]*Pp[p, q, P, 3, 4]-1/Sqrt[6]*Pp[p, q, P, 4, 3]+2/Sqrt[6]*Sp[p, q, P, 2, 1];
    array[[3, 2, 2, 2]] =  1/Sqrt[6]*Pm[p, q, P, 3, 4]-1/Sqrt[6]*Pm[p, q, P, 4, 3]-2/Sqrt[6]*Sm[p, q, P, 2, 1];
    array[[3, 3, 2, 1]] =  2/Sqrt[6]*Pp[p, q, P, 3, 1]+1/Sqrt[6]*Pp[p, q, P, 4, 2]-1/Sqrt[6]*Sp[p, q, P, 2, 4];
    array[[3, 3, 2, 2]] =  2/Sqrt[6]*Pm[p, q, P, 3, 1]+1/Sqrt[6]*Pm[p, q, P, 4, 2]+1/Sqrt[6]*Sm[p, q, P, 2, 4];
    array[[3, 4, 2, 1]] =  2/Sqrt[6]*Pp[p, q, P, 4, 1]-1/Sqrt[6]*Pp[p, q, P, 3, 2]+1/Sqrt[6]*Sp[p, q, P, 2, 3];
    array[[3, 4, 2, 2]] =  2/Sqrt[6]*Pm[p, q, P, 4, 1]-1/Sqrt[6]*Pm[p, q, P, 3, 2]-1/Sqrt[6]*Sm[p, q, P, 2, 3];

    array[[4, 1, 1, 1]] =  1/Sqrt[2]*Sp[p, q, P, 4, 4]-1/Sqrt[2]*Pp[p, q, P, 2, 2];
    array[[4, 1, 1, 2]] =  1/Sqrt[2]*Sm[p, q, P, 4, 4]+1/Sqrt[2]*Pm[p, q, P, 2, 2];
    array[[4, 2, 1, 1]] =  1/Sqrt[2]*Sp[p, q, P, 3, 4]+1/Sqrt[2]*Sp[p, q, P, 4, 3];
    array[[4, 2, 1, 2]] =  1/Sqrt[2]*Sm[p, q, P, 3, 4]+1/Sqrt[2]*Sm[p, q, P, 4, 3];
    array[[4, 3, 1, 1]] = -1/Sqrt[2]*Sp[p, q, P, 4, 2]-1/Sqrt[2]*Pp[p, q, P, 2, 4];
    array[[4, 3, 1, 2]] = -1/Sqrt[2]*Sm[p, q, P, 4, 2]+1/Sqrt[2]*Pm[p, q, P, 2, 4];
    array[[4, 4, 1, 1]] =  1/Sqrt[2]*Sp[p, q, P, 3, 2]+1/Sqrt[2]*Pp[p, q, P, 2, 3];
    array[[4, 4, 1, 2]] =  1/Sqrt[2]*Sm[p, q, P, 3, 2]-1/Sqrt[2]*Pm[p, q, P, 2, 3];

    array[[4, 1, 2, 1]] =  1/Sqrt[2]*Pp[p, q, P, 4, 4]-1/Sqrt[2]*Sp[p, q, P, 2, 2];
    array[[4, 1, 2, 2]] =  1/Sqrt[2]*Pm[p, q, P, 4, 4]+1/Sqrt[2]*Sm[p, q, P, 2, 2];
    array[[4, 2, 2, 1]] =  1/Sqrt[2]*Pp[p, q, P, 3, 4]+1/Sqrt[2]*Pp[p, q, P, 4, 3];
    array[[4, 2, 2, 2]] =  1/Sqrt[2]*Pm[p, q, P, 3, 4]+1/Sqrt[2]*Pm[p, q, P, 4, 3];
    array[[4, 3, 2, 1]] = -1/Sqrt[2]*Pp[p, q, P, 4, 2]-1/Sqrt[2]*Sp[p, q, P, 2, 4];
    array[[4, 3, 2, 2]] = -1/Sqrt[2]*Pm[p, q, P, 4, 2]+1/Sqrt[2]*Sm[p, q, P, 2, 4];
    array[[4, 4, 2, 1]] =  1/Sqrt[2]*Pp[p, q, P, 3, 2]+1/Sqrt[2]*Sp[p, q, P, 2, 3];
    array[[4, 4, 2, 2]] =  1/Sqrt[2]*Pm[p, q, P, 3, 2]-1/Sqrt[2]*Sm[p, q, P, 2, 3];

    array 
];

CX[p_,q_,P_] := Module[{array},
    array = ConstantArray[0, {4,4,2,2}]; 
    array[[1, 1, 1, 1]] = CSp[-p, -q, -P, 1, 1]; array[[1, 1, 1, 2]] = CSm[-p, -q, -P, 1, 1];
    array[[1, 2, 1, 1]] = CSp[-p, -q, -P, 1, 2]; array[[1, 2, 1, 2]] = CSm[-p, -q, -P, 1, 2];
    array[[1, 3, 1, 1]] = CSp[-p, -q, -P, 1, 3]; array[[1, 3, 1, 2]] = CSm[-p, -q, -P, 1, 3];
    array[[1, 4, 1, 1]] = CSp[-p, -q, -P, 1, 4]; array[[1, 4, 1, 2]] = CSm[-p, -q, -P, 1, 4];

    array[[1, 1, 2, 1]] = CPp[-p, -q, -P, 1, 1]; array[[1, 1, 2, 2]] = CPm[-p, -q, -P, 1, 1];
    array[[1, 2, 2, 1]] = CPp[-p, -q, -P, 1, 2]; array[[1, 2, 2, 2]] = CPm[-p, -q, -P, 1, 2];
    array[[1, 3, 2, 1]] = CPp[-p, -q, -P, 1, 3]; array[[1, 3, 2, 2]] = CPm[-p, -q, -P, 1, 3];
    array[[1, 4, 2, 1]] = CPp[-p, -q, -P, 1, 4]; array[[1, 4, 2, 2]] = CPm[-p, -q, -P, 1, 4];

    array[[2, 1, 1, 1]] = 1/Sqrt[3]*CSp[-p, -q, -P, 3, 3]+1/Sqrt[3]*CSp[-p, -q, -P, 4, 4]+1/Sqrt[3]*CPp[-p, -q, -P, 2, 2];
    array[[2, 1, 1, 2]] = 1/Sqrt[3]*CSm[-p, -q, -P, 3, 3]+1/Sqrt[3]*CSm[-p, -q, -P, 4, 4]-1/Sqrt[3]*CPm[-p, -q, -P, 2, 2];
    array[[2, 2, 1, 1]] = 1/Sqrt[3]*CSp[-p, -q, -P, 3, 4]-1/Sqrt[3]*CSp[-p, -q, -P, 4, 3]-1/Sqrt[3]*CPp[-p, -q, -P, 2, 1];
    array[[2, 2, 1, 2]] = 1/Sqrt[3]*CSm[-p, -q, -P, 3, 4]-1/Sqrt[3]*CSm[-p, -q, -P, 4, 3]+1/Sqrt[3]*CPm[-p, -q, -P, 2, 1];
    array[[2, 3, 1, 1]] = 1/Sqrt[3]*CSp[-p, -q, -P, 3, 1]-1/Sqrt[3]*CSp[-p, -q, -P, 4, 2]+1/Sqrt[3]*CPp[-p, -q, -P, 2, 4];
    array[[2, 3, 1, 2]] = 1/Sqrt[3]*CSm[-p, -q, -P, 3, 1]-1/Sqrt[3]*CSm[-p, -q, -P, 4, 2]-1/Sqrt[3]*CPm[-p, -q, -P, 2, 4];
    array[[2, 4, 1, 1]] = 1/Sqrt[3]*CSp[-p, -q, -P, 3, 2]+1/Sqrt[3]*CSp[-p, -q, -P, 4, 1]-1/Sqrt[3]*CPp[-p, -q, -P, 2, 3];
    array[[2, 4, 1, 2]] = 1/Sqrt[3]*CSm[-p, -q, -P, 3, 2]+1/Sqrt[3]*CSm[-p, -q, -P, 4, 1]+1/Sqrt[3]*CPm[-p, -q, -P, 2, 3];

    array[[2, 1, 2, 1]] = 1/Sqrt[3]*CPp[-p, -q, -P, 3, 3]+1/Sqrt[3]*CPp[-p, -q, -P, 4, 4]+1/Sqrt[3]*CSp[-p, -q, -P, 2, 2];
    array[[2, 1, 2, 2]] = 1/Sqrt[3]*CPm[-p, -q, -P, 3, 3]+1/Sqrt[3]*CPm[-p, -q, -P, 4, 4]-1/Sqrt[3]*CSm[-p, -q, -P, 2, 2];
    array[[2, 2, 2, 1]] = 1/Sqrt[3]*CPp[-p, -q, -P, 3, 4]-1/Sqrt[3]*CPp[-p, -q, -P, 4, 3]-1/Sqrt[3]*CSp[-p, -q, -P, 2, 1];
    array[[2, 2, 2, 2]] = 1/Sqrt[3]*CPm[-p, -q, -P, 3, 4]-1/Sqrt[3]*CPm[-p, -q, -P, 4, 3]+1/Sqrt[3]*CSm[-p, -q, -P, 2, 1];
    array[[2, 3, 2, 1]] = 1/Sqrt[3]*CPp[-p, -q, -P, 3, 1]-1/Sqrt[3]*CPp[-p, -q, -P, 4, 2]+1/Sqrt[3]*CSp[-p, -q, -P, 2, 4];
    array[[2, 3, 2, 2]] = 1/Sqrt[3]*CPm[-p, -q, -P, 3, 1]-1/Sqrt[3]*CPm[-p, -q, -P, 4, 2]-1/Sqrt[3]*CSm[-p, -q, -P, 2, 4];
    array[[2, 4, 2, 1]] = 1/Sqrt[3]*CPp[-p, -q, -P, 3, 2]+1/Sqrt[3]*CPp[-p, -q, -P, 4, 1]-1/Sqrt[3]*CSp[-p, -q, -P, 2, 3];
    array[[2, 4, 2, 2]] = 1/Sqrt[3]*CPm[-p, -q, -P, 3, 2]+1/Sqrt[3]*CPm[-p, -q, -P, 4, 1]+1/Sqrt[3]*CSm[-p, -q, -P, 2, 3];

    array[[3, 1, 1, 1]] =  2/Sqrt[6]*CSp[-p, -q, -P, 3, 3]-1/Sqrt[6]*CSp[-p, -q, -P, 4, 4]-1/Sqrt[6]*CPp[-p, -q, -P, 2, 2];
    array[[3, 1, 1, 2]] =  2/Sqrt[6]*CSm[-p, -q, -P, 3, 3]-1/Sqrt[6]*CSm[-p, -q, -P, 4, 4]+1/Sqrt[6]*CPm[-p, -q, -P, 2, 2];
    array[[3, 2, 1, 1]] =  1/Sqrt[6]*CSp[-p, -q, -P, 3, 4]-1/Sqrt[6]*CSp[-p, -q, -P, 4, 3]+2/Sqrt[6]*CPp[-p, -q, -P, 2, 1];
    array[[3, 2, 1, 2]] =  1/Sqrt[6]*CSm[-p, -q, -P, 3, 4]-1/Sqrt[6]*CSm[-p, -q, -P, 4, 3]-2/Sqrt[6]*CPm[-p, -q, -P, 2, 1];
    array[[3, 3, 1, 1]] =  2/Sqrt[6]*CSp[-p, -q, -P, 3, 1]+1/Sqrt[6]*CSp[-p, -q, -P, 4, 2]-1/Sqrt[6]*CPp[-p, -q, -P, 2, 4];
    array[[3, 3, 1, 2]] =  2/Sqrt[6]*CSm[-p, -q, -P, 3, 1]+1/Sqrt[6]*CSm[-p, -q, -P, 4, 2]+1/Sqrt[6]*CPm[-p, -q, -P, 2, 4];
    array[[3, 4, 1, 1]] =  2/Sqrt[6]*CSp[-p, -q, -P, 4, 1]-1/Sqrt[6]*CSp[-p, -q, -P, 3, 2]+1/Sqrt[6]*CPp[-p, -q, -P, 2, 3];
    array[[3, 4, 1, 2]] =  2/Sqrt[6]*CSm[-p, -q, -P, 4, 1]-1/Sqrt[6]*CSm[-p, -q, -P, 3, 2]-1/Sqrt[6]*CPm[-p, -q, -P, 2, 3];

    array[[3, 1, 2, 1]] =  2/Sqrt[6]*CPp[-p, -q, -P, 3, 3]-1/Sqrt[6]*CPp[-p, -q, -P, 4, 4]-1/Sqrt[6]*CSp[-p, -q, -P, 2, 2];
    array[[3, 1, 2, 2]] =  2/Sqrt[6]*CPm[-p, -q, -P, 3, 3]-1/Sqrt[6]*CPm[-p, -q, -P, 4, 4]+1/Sqrt[6]*CSm[-p, -q, -P, 2, 2];
    array[[3, 2, 2, 1]] =  1/Sqrt[6]*CPp[-p, -q, -P, 3, 4]-1/Sqrt[6]*CPp[-p, -q, -P, 4, 3]+2/Sqrt[6]*CSp[-p, -q, -P, 2, 1];
    array[[3, 2, 2, 2]] =  1/Sqrt[6]*CPm[-p, -q, -P, 3, 4]-1/Sqrt[6]*CPm[-p, -q, -P, 4, 3]-2/Sqrt[6]*CSm[-p, -q, -P, 2, 1];
    array[[3, 3, 2, 1]] =  2/Sqrt[6]*CPp[-p, -q, -P, 3, 1]+1/Sqrt[6]*CPp[-p, -q, -P, 4, 2]-1/Sqrt[6]*CSp[-p, -q, -P, 2, 4];
    array[[3, 3, 2, 2]] =  2/Sqrt[6]*CPm[-p, -q, -P, 3, 1]+1/Sqrt[6]*CPm[-p, -q, -P, 4, 2]+1/Sqrt[6]*CSm[-p, -q, -P, 2, 4];
    array[[3, 4, 2, 1]] =  2/Sqrt[6]*CPp[-p, -q, -P, 4, 1]-1/Sqrt[6]*CPp[-p, -q, -P, 3, 2]+1/Sqrt[6]*CSp[-p, -q, -P, 2, 3];
    array[[3, 4, 2, 2]] =  2/Sqrt[6]*CPm[-p, -q, -P, 4, 1]-1/Sqrt[6]*CPm[-p, -q, -P, 3, 2]-1/Sqrt[6]*CSm[-p, -q, -P, 2, 3];

    array[[4, 1, 1, 1]] =  1/Sqrt[2]*CSp[-p, -q, -P, 4, 4]-1/Sqrt[2]*CPp[-p, -q, -P, 2, 2];
    array[[4, 1, 1, 2]] =  1/Sqrt[2]*CSm[-p, -q, -P, 4, 4]+1/Sqrt[2]*CPm[-p, -q, -P, 2, 2];
    array[[4, 2, 1, 1]] =  1/Sqrt[2]*CSp[-p, -q, -P, 3, 4]+1/Sqrt[2]*CSp[-p, -q, -P, 4, 3];
    array[[4, 2, 1, 2]] =  1/Sqrt[2]*CSm[-p, -q, -P, 3, 4]+1/Sqrt[2]*CSm[-p, -q, -P, 4, 3];
    array[[4, 3, 1, 1]] = -1/Sqrt[2]*CSp[-p, -q, -P, 4, 2]-1/Sqrt[2]*CPp[-p, -q, -P, 2, 4];
    array[[4, 3, 1, 2]] = -1/Sqrt[2]*CSm[-p, -q, -P, 4, 2]+1/Sqrt[2]*CPm[-p, -q, -P, 2, 4];
    array[[4, 4, 1, 1]] =  1/Sqrt[2]*CSp[-p, -q, -P, 3, 2]+1/Sqrt[2]*CPp[-p, -q, -P, 2, 3];
    array[[4, 4, 1, 2]] =  1/Sqrt[2]*CSm[-p, -q, -P, 3, 2]-1/Sqrt[2]*CPm[-p, -q, -P, 2, 3];

    array[[4, 1, 2, 1]] =  1/Sqrt[2]*CPp[-p, -q, -P, 4, 4]-1/Sqrt[2]*CSp[-p, -q, -P, 2, 2];
    array[[4, 1, 2, 2]] =  1/Sqrt[2]*CPm[-p, -q, -P, 4, 4]+1/Sqrt[2]*CSm[-p, -q, -P, 2, 2];
    array[[4, 2, 2, 1]] =  1/Sqrt[2]*CPp[-p, -q, -P, 3, 4]+1/Sqrt[2]*CPp[-p, -q, -P, 4, 3];
    array[[4, 2, 2, 2]] =  1/Sqrt[2]*CPm[-p, -q, -P, 3, 4]+1/Sqrt[2]*CPm[-p, -q, -P, 4, 3];
    array[[4, 3, 2, 1]] = -1/Sqrt[2]*CPp[-p, -q, -P, 4, 2]-1/Sqrt[2]*CSp[-p, -q, -P, 2, 4];
    array[[4, 3, 2, 2]] = -1/Sqrt[2]*CPm[-p, -q, -P, 4, 2]+1/Sqrt[2]*CSm[-p, -q, -P, 2, 4];
    array[[4, 4, 2, 1]] =  1/Sqrt[2]*CPp[-p, -q, -P, 3, 2]+1/Sqrt[2]*CSp[-p, -q, -P, 2, 3];
    array[[4, 4, 2, 2]] =  1/Sqrt[2]*CPm[-p, -q, -P, 3, 2]-1/Sqrt[2]*CSm[-p, -q, -P, 2, 3];

    array 
];

Basis=X[p,q,P];
CBasis=CX[p,q,P];

(* Flatten:{1111, 1112, 1121, 1122, 1211, 1212,
1221, 1222, 1311, 1312, 1321, 1322,
1411, 1412, 1421, 1422, 2111, 2112,
2121, 2122, 2211, 2212, 2221, 2222,
2311, 2312, 2321, 2322, 2411, 2412,
2421, 2422, 3111, 3112, 3121, 3122,
3211, 3212, 3221, 3222, 3311, 3312,
3321, 3322, 3411, 3412, 3421, 3422,
4111, 4112, 4121, 4122, 4211, 4212,
4221, 4222, 4311, 4312, 4321, 4322,
4411, 4412, 4421, 4422} *)

Basis64=Flatten[Basis];
CBasis64=Flatten[CBasis];

(* Rules *)

rules={
    Pair[Momentum[x_], Momentum[y_] ] :> If[x === y, 1, 0]
  };
  
mergeFpair =
   Times[pre___, Fpair[a_, b_], Fpair[c_, d_], post___] :> 
   Times[pre, Fpair[Dot[a, c], Dot[b, d] ], post];

traceFpair =
  Fpair[a_, b_] :> TR@DotSimplify[a]*TR@DotSimplify[b];

(* Check for orthogonality relations *)
(* CloseKernels[];
LaunchKernels[64];
DistributeDefinitions[Basis,CBasis,rules];
ParallelEvaluate[<< FeynCalc`];
Table[Print[""], {5}];
Result = 0.25*ParallelMap[(Expand[#]//.mergeFpair//.traceFpair//.rules) &, Outer[Times, CBasis64, Basis64] ];
Print[Result];
Paths=FileNameJoin[{"Results","Orthogonality.wls"}];
If[$OperatingSystem === "Windows", 
    Paths = StringReplace[Paths, "\\" -> "\\\\"]
];
Export[Paths, Result];
nonZeroPositions = Position[Result, _?(# != 0 &)];
nonZeroValues = Extract[Result, nonZeroPositions];
Print[nonZeroPositions,"\n", nonZeroValues]; *)
