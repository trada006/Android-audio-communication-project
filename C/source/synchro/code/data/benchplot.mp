prologues := 2;
filenametemplate "%j-%c.mps";

input graph;
input sarith;

%input TEX; 
%TEXPRE("%&latex" & char(10) & 
%"\documentclass[12pt]{article}\begin{document}"); 
%TEXPOST("\end{document}"); 

verbatimtex
%&latex
\documentclass[8pt]{extarticle}
  \begin{document}
  etex

init_numbers(btex $-$ etex,btex $1$ etex, btex ${\times}10$ etex, btex ${}^-$ etex, btex ${}^2$ etex);


picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;

picture dd; dd :=
dashpattern(on 4bp off 2bp on 1bp off 2bp);

picture diamond;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm rotated 45;
diamond := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture circle; 
draw fullcircle scaled 0.19cm;
circle := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture cross; 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
cross := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture dot; 
fill fullcircle scaled 0.15cm;
dot := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture times; times := cross rotated 45;


%%%%%%%%%%%%%
%%%%%Legend picture

picture legend;
sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

ls:=0;
label.rt(btex Nonlinearity then Massey etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex Pilot and data (no optimiser) etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex Pilot and data etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=3;
%label.rt(btex Direct etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;

legend := currentpicture; currentpicture := ORIGINAL_PICTURE;

plotheight = 6cm;
plotwidth = 7.65cm;

legendx = plotwidth - 5cm;
legendy = plotheight - 1.2cm;


%%%%%%%%%%%%
%%%%Begin plot 1, with 5 pilots

beginfig(1);
draw begingraph(plotwidth, plotheight);
 setcoords(log,linear);
 glabel.bot(btex L etex, OUT) shifted (-1.8cm,1.0cm);
  glabel.lft(btex $\mu$s etex, OUT) rotated 90 shifted (1.05cm,0.0cm);

  %gdraw "Directbench" dashed evenly; 
  gdraw "Recursivebench"; gdraw "Recursivebench" plot(dot);
  gdraw "RecursiveTabulatedMasseybench"; gdraw "RecursiveTabulatedMasseybench" plot(circle);
  gdraw "RecursiveNoRefinementbench"; gdraw "RecursiveNoRefinementbench" plot(times);

  setrange(whatever,0, whatever,130);
  %otick.lft(btex $60$ etex, 60);
  %otick.lft(btex $50$ etex, 50);
  %otick.lft(btex $40$ etex, 40);
  %otick.lft(btex $30$ etex, 30);
  %otick.lft(btex $20$ etex, 20);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0$ etex, 0);
  %for u=auto.y: otick.lft(format("%g",u), u); endfor
  %for u=auto.x: otick.bot(format("%g",u), u); endfor
  %otick.bot(btex $20$ etex, 20);
  otick.bot(btex $50$ etex, 50);
  otick.bot(btex $200$ etex, 200);
  otick.bot(btex $500$ etex, 500);
  otick.bot(btex $1500$ etex, 1500);
  otick.bot(btex $3500$ etex, 3500);
  autogrid(,otick.lft);
  
draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight+0.25bp)--(plotwidth+0.25bp,plotheight+0.25bp)--(plotwidth+0.25bp,-10cm)--cycle

endfig;




end;