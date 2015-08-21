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
diamond := image( fill unitsquare shifted(-0.5,-0.5) scaled 0.10cm rotated 45; );

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture circle; 
circle := image( draw fullcircle scaled 0.14cm; );

picture cross;
cross := image( 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
);

picture dot; 
dot := image( fill fullcircle scaled 0.11cm; );

picture times; times := cross rotated 45;

%draws and oval of width w and height h
def oval(expr w, h) = (w/2,0)...(0,h/2)...(-w/2,0)...(0,-h/2)...cycle  enddef;

picture nsmall; 
nsmall := image( label(btex $N=75$ etex, (0,0)) rotated -23; );
picture nlarge; 
nlarge := image( label(btex $N=2555$ etex, (0,0)) rotated -23; );

%%%%%%%%%%%%%
%%%%%Legend picture

picture legendtwo;

sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

%ls:=1;
%label.rt(btex $L=40$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;

ls:=0;
label.rt(btex All symbols known etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) withcolor(0.8,0.8,0.8);
draw diamond shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $\hat{\tau}$ (new) etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) withcolor(0.8,0.8,0.8);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=4;
label.rt(btex Pilots only etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;

ls:=2;
label.rt(btex Square then Massey etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) withcolor(0.8,0.8,0.8);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=3;
label.rt(btex Absolute value then Massey etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) withcolor(0.8,0.8,0.8);
draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=5;
label.rt(btex Modified CRB etex, (1cm,sh+ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);

legendtwo := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 10cm;
plotwidth = 17cm;
legendx = 0.1cm;
legendy = 0.1cm;
clipboundary = -3cm;

vardef ymin = "5e-7" enddef;
vardef ymax = "2" enddef;
vardef xmin = "-21" enddef;
vardef xmax = "21" enddef;

%only draws plot markers if in plot range.  For some reason metapost doesn't
%do this automatically
def plotdraw(expr name, d) =
  gdata(name, s,
      if ymin Sleq s2 :
      if s2 Sleq ymax :
	if s1 Sleq xmax :
	  if xmin Sleq s1 :
	    glabel(d,s1,s2);
	  fi; fi; fi; fi;
    )
enddef;

def toplabelledelipse(expr labname, xscale, yscale, x, y) =
  draw fullcircle xscaled xscale yscaled yscale shifted (x,y);
  label.top(labname, (x,y+yscale/2));
enddef;

def botlabelledelipse(expr labname, xscale, yscale, x, y) =
  draw fullcircle xscaled xscale yscaled yscale shifted (x,y);
  label.bot(labname, (x,y-yscale/2));
enddef;


%%%%%%%%%%%%
%%%%Begin plot 1, with 5 pilots

beginfig(1);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0,1.1cm);
  glabel.lft(btex MSE etex, OUT) rotated 90 shifted (1.3cm,-1.6cm);

  
%some light grey lines first
  gdraw "recursiveP15L75" withcolor(0.8,0.8,0.8);
  gdraw "slnP15L75" withcolor(0.8,0.8,0.8);
  gdraw "avnP15L75" withcolor(0.8,0.8,0.8);
  gdraw "recursiveP511L2555" withcolor(0.8,0.8,0.8);
  gdraw "slnP511L2555" withcolor(0.8,0.8,0.8);
  gdraw "avnP511L2555" withcolor(0.8,0.8,0.8);

  %draw pilot only curves
  gdraw "pilotonlyP15L75"  dashed evenly;
  gdraw "pilotonlyP511L2555"  dashed evenly;

  %draw all symbols known curves
  gdraw "allpilotsP15L75" withcolor(0.8,0.8,0.8);
  gdraw "allpilotsP511L2555" withcolor(0.8,0.8,0.8);

  %crbs
  gdraw "crbbeta3L75"; %glabel.lft(nsmall, 25);
  gdraw "crbbeta3L2555"; %glabel.lft(nlarge, 25);

  %set ranges and ticks
  setrange(xmin,ymin,xmax,ymax);
  %otick.lft(btex $10^3$ etex, 1000);
  %otick.lft(btex $100$ etex, 100);
  otick.lft(btex $10$ etex, 10);
  otick.lft(btex $1$ etex, 1);
  otick.lft(btex $0.1$ etex, 0.1);
  otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-3}$ etex, 0.001);
  otick.lft(btex $10^{-4}$ etex, "1e-4");
  otick.lft(btex $10^{-5}$ etex, "1e-5");
  otick.lft(btex $10^{-6}$ etex, "1e-6");
  %otick.lft(btex $10^{-7}$ etex, "1e-7");
  %for u=auto.y: otick.lft(format("%g",u), u); endfor
  %for u=auto.x: otick.bot(format("%g",u), u); endfor
  autogrid(otick.bot,);

  %draw dots %need to manual clipping, using gdata for this
  plotdraw("recursiveP15L75", dot); 
  plotdraw("slnP15L75", circle);
  plotdraw("avnP15L75", cross);
  %gdraw "recursiveP255L1275" plot(dot);
  %gdraw "slnP255L1275" plot(circle);
  %gdraw "avnP255L1275" plot(cross);
  plotdraw("recursiveP511L2555",dot);
  plotdraw("slnP511L2555",circle);
  %gdata("recursiveP511L2555", s, if ymin Sleq s2  : glabel(dot,s1,s2); fi; ); %gdraw "recursiveP511L2555" plot(dot);
  %gdata("slnP511L2555", s, if ymin Sleq s2  : glabel(circle,s1,s2); fi; ); %gdraw "slnP511L2555" plot(circle);
  plotdraw("avnP511L2555", cross);
  plotdraw("allpilotsP15L75",diamond);
  plotdraw("allpilotsP511L2555",diamond);

draw legendtwo shifted(legendx,legendy);

endgraph;

toplabelledelipse(btex $N=75$ etex, 0.7cm,2cm,13.35cm,4.1cm);
botlabelledelipse(btex $N=2555$ etex, 0.7cm,1.8cm,10.95cm,2.6cm);
  
%draw fullcircle xscaled 0.75cm yscaled 2cm shifted (13.3cm,4.2cm);
%label.top(btex $N=75$ etex, (13.3cm,5.2cm));

%draw fullcircle xscaled 0.75cm yscaled 2cm shifted (11cm,2.5cm);
%label.bot(btex $N=2555$ etex, (11cm,1.5cm));

clip currentpicture to (clipboundary,clipboundary)--(clipboundary,plotheight+0.25bp)--(plotwidth+0.25bp,plotheight+0.25bp)--(plotwidth+0.25bp,clipboundary)--cycle;

endfig;



% % %%%%%%%%%%%%
% % %%%% Small rolloff plot

beginfig(2);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0,1.1cm);
  glabel.lft(btex MSE etex, OUT) rotated 90 shifted (1.3cm,-1.6cm);

  
%some light grey lines first
  gdraw "recursiveP15L75smallrolloff" withcolor(0.8,0.8,0.8);
  gdraw "slnP15L75smallrolloff" withcolor(0.8,0.8,0.8);
  gdraw "avnP15L75smallrolloff" withcolor(0.8,0.8,0.8);
  gdraw "recursiveP511L2555smallrolloff" withcolor(0.8,0.8,0.8);
  gdraw "slnP511L2555smallrolloff" withcolor(0.8,0.8,0.8);
  gdraw "avnP511L2555smallrolloff" withcolor(0.8,0.8,0.8);

  %draw pilot only curves
  gdraw "pilotonlyP15L75smallrolloff"  dashed evenly;
  gdraw "pilotonlyP511L2555smallrolloff"  dashed evenly;

  %draw all symbols known curves
  gdraw "allpilotsP15L75smallrolloff" withcolor(0.8,0.8,0.8); 
  gdraw "allpilotsP511L2555smallrolloff" withcolor(0.8,0.8,0.8); 

  %crbs
  gdraw "crbbeta30L75"; glabel.lft(nsmall, 25);
  gdraw "crbbeta30L2555"; glabel.lft(nlarge, 25);

  %set ranges and ticks
  setrange(xmin,ymin,xmax,ymax);
  %otick.lft(btex $10^3$ etex, 1000);
  %otick.lft(btex $100$ etex, 100);
  otick.lft(btex $10$ etex, 10);
  otick.lft(btex $1$ etex, 1);
  otick.lft(btex $0.1$ etex, 0.1);
  otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-3}$ etex, 0.001);
  otick.lft(btex $10^{-4}$ etex, "1e-4");
  otick.lft(btex $10^{-5}$ etex, "1e-5");
  otick.lft(btex $10^{-6}$ etex, "1e-6");
  %otick.lft(btex $10^{-7}$ etex, "1e-7");
  %for u=auto.y: otick.lft(format("%g",u), u); endfor
  %for u=auto.x: otick.bot(format("%g",u), u); endfor
  autogrid(otick.bot,);

  %draw dots %need to manual clipping, using gdata for this
  plotdraw("recursiveP15L75smallrolloff", dot); 
  plotdraw("slnP15L75smallrolloff", circle);
  plotdraw("avnP15L75smallrolloff", cross);
  plotdraw("recursiveP511L2555smallrolloff",dot);
  plotdraw("slnP511L2555smallrolloff",circle);
  plotdraw("avnP511L2555smallrolloff", cross);
  plotdraw("allpilotsP15L75smallrolloff", diamond);
  plotdraw("allpilotsP511L2555smallrolloff", diamond);
  
draw legendtwo shifted(legendx,legendy);

endgraph;

toplabelledelipse(btex $N=75$ etex, 0.7cm,3cm,13.38cm,5cm);
botlabelledelipse(btex $N=2555$ etex, 0.7cm,2.65cm,12.55cm,2.5cm);

clip currentpicture to (clipboundary,clipboundary)--(clipboundary,plotheight+0.25bp)--(plotwidth+0.25bp,plotheight+0.25bp)--(plotwidth+0.25bp,clipboundary)--cycle;

endfig;



end;