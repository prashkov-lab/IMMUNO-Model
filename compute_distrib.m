function compute_distrib


   nsamp=1000;

mu=log(2)/2.5; tend=6;  x0=mu;

vmax = [2.04095779756017, 1.59692849121587, 1.59692849121587, 1.79741691222017,...
    2.04095779756017, 2.52100523745854, 2.19599418149267, 2.39755841682793, ...
    2.39755841682793, 2.29732154535348, 2.29732154535348, 2.39755841682793, ...
    2.29732154535348, 2.52100523745854];


  xmax_m=zeros(14,2);

  EC50low_m=zeros(14,5);



# initial guess and lower/upper bounds

  X0 = 0.1;

# monotone
UB=20;

  data=dlmread ("2026_myDF_new0agg.dat",' ',1,0);
  b_max=data(1:4:end,3);
  b_M_m=data(3:4:end,1);
  b_M_sd=data(3:4:end,2);
  indexs=data(1:4:end,1)
  p_b_m=data(2:4:end,1);
  p_b_sd=data(2:4:end,2);

  for j= 1:2
     # generate random values randn ~ N(0,1) -> p_b_m+p_b_sd*Z
    b_M=randn(nsamp,1)*b_M_sd(j)+b_M_m(j);

    p_b = randn(nsamp,1)*p_b_sd(j)+p_b_m(j);

    p_b=p_b(b_M>1e-99);
    b_M=b_M(b_M>1e-99);

    b_M=b_M(p_b>1e-99);
    p_b=p_b(p_b>1e-99);

    nsamp0=length(b_M);

      EC50low=zeros(1,nsamp0);

 if (j==1)
   xr=(0.1: 1e3: 1e5);
 else
   xr=(0.001: 0.001: UB);
 endif

        for jj=1:nsamp0
       ff=-h0(xr,j,jj);
##
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
  EC50low(jj) = min(xr(ff>=(vmax(indexs(j))+1)/2));
endif
endfor


    EC50low=EC50low(EC50low>1e-6); lenvec=length(EC50low);
    if (~isempty(EC50low))
      EC50low_m(indexs(j),1)=mean(EC50low);
      EC50low_m(indexs(j),2)=std(EC50low);
        [nn,xx]=  hist(EC50low,35);
  out=[xx',nn'];
  filename=["EC50_" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)
    endif
endfor


 LB= 0.0001;
  UB =15 ;


  % first scenario
clearvars data g_max b_M g_M p_b p_g indexs
  % first scenario

  data=dlmread("2026_myDF_new1agg.dat",' ',1,0);
    indexs=data(1:7:end,1);
  b_max=data(1:7:end,3);

  g_max_m=data(4:7:end,1);
    g_max_sd= data(4:7:end,2);
  b_M_m=data(5:7:end,1);
     b_M_sd=data(5:7:end,2);
  g_M_m=data(6:7:end,1);
    g_M_sd =data(6:7:end,2);

  p_b_m=data(2:7:end,1);
    p_b_sd=data(2:7:end,2);

  p_g_m=data(3:7:end,1);
  p_g_sd=data(3:7:end,2);

  for j=1
    # generate random values randn ~ N(0,1) -> p_b_m+p_b_sd*Z
    g_max = randn(nsamp,1)*g_max_sd(j) + g_max_m(j);
    b_M=randn(nsamp,1)*b_M_sd(j)+b_M_m(j);
    g_M = randn(nsamp,1)*g_M_sd(j)+g_M_m(j);
    p_b = randn(nsamp,1)*p_b_sd(j)+p_b_m(j);
    p_g = ones(nsamp,1)*p_g_m(j);

    b_M=b_M(g_max>1e-99);
    g_M=g_M(g_max>1e-99);
    p_b=p_b(g_max>1e-99);
    p_g=p_g(g_max>1e-99);
    g_max=g_max(g_max>1e-99);

    g_max=g_max(b_M>1e-99);
    g_M=g_M(b_M>1e-99);
    p_b=p_b(b_M>1e-99);
    p_g=p_g(b_M>1e-99);
    b_M=b_M(b_M>1e-99);

    g_max=g_max(g_M>1e-99);
    b_M=b_M(g_M>1e-99);
    p_b=p_b(g_M>1e-99);
    p_g=p_g(g_M>1e-99);
    g_M=g_M(g_M>1e-99);

    g_max=g_max(p_b>1e-99);
    b_M=b_M(p_b>1e-99);
    p_g=p_g(p_b>1e-99);
    g_M=g_M(p_b>1e-99);
    p_b=p_b(p_b>1e-99);

    nsamp0=length(g_max);

      xmax =zeros(1,nsamp0);
      EC50low=zeros(1,nsamp0);

    for jj=1:nsamp0
    [xmax(jj), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h1(x,j,jj), [], [], LB, UB);


    xr=(0: 1e-4: xmax(jj));
       ff=-h1(xr,j,jj);

if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50low(jj) = min(xr(ff>=(vmax(indexs(j))+1)/2));
  endif

endfor
    xmax_m(indexs(j),1)=mean(xmax);

  [nn,xx]=  hist(xmax,35);
  out=[xx',nn'];
  filename=["xmax" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,'delimiter', " ",'precision',4)


 xmax_m(indexs(j),2)=std(xmax);

    EC50low=EC50low(EC50low>1e-6); lenvec=length(EC50low);
    if (~isempty(EC50low))
      EC50low_m(indexs(j),1)=mean(EC50low);
      EC50low_m(indexs(j),2)=std(EC50low);

  [nn,xx]=  hist(EC50low,35);
  out=[xx',nn'];
  filename=["EC50_" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)
    endif


endfor

clearvars g_max b_M g_M p_b p_g b_max indexs

  data=dlmread("2026_myDF_new3agg.dat",' ',1,0);
    indexs=data(1:6:end,1);
  b_max=data(1:6:end,3);
  g_max_m=data(3:6:end,1);
    g_max_sd= data(3:6:end,2);

  b_M_m=data(4:6:end,1);
    b_M_sd=data(4:6:end,2);

  g_M_m=data(5:6:end,1);
    g_M_sd =data(5:6:end,2);

  p_b_m=data(2:6:end,1);
  p_b_sd=data(2:6:end,2);


 for j=2
    # generate random values randn ~ N(0,1) -> p_b_m+p_b_sd*Z
    g_max = randn(nsamp,1)*g_max_sd(j) + g_max_m(j);
    b_M=randn(nsamp,1)*b_M_sd(j)+b_M_m(j);
    g_M = randn(nsamp,1)*g_M_sd(j)+g_M_m(j);
    p_b = randn(nsamp,1)*p_b_sd(j)+p_b_m(j);

    b_M=b_M(g_max>1e-99);
    g_M=g_M(g_max>1e-99);
    p_b=p_b(g_max>1e-99);
    g_max=g_max(g_max>1e-99);

    g_max=g_max(b_M>1e-99);
    g_M=g_M(b_M>1e-99);
    p_b=p_b(b_M>1e-99);
    b_M=b_M(b_M>1e-99);

    g_max=g_max(g_M>1e-99);
    b_M=b_M(g_M>1e-99);
    p_b=p_b(g_M>1e-99);
    g_M=g_M(g_M>1e-99);

    g_max=g_max(p_b>1e-99);
    b_M=b_M(p_b>1e-99);
    g_M=g_M(p_b>1e-99);
    p_b=p_b(p_b>1e-99);

    nsamp0=length(g_max);
    p_g=3*ones(nsamp0,1);

          xmax =zeros(1,nsamp0);
      EC50low=zeros(1,nsamp0);

    for jj=1:nsamp0
    [xmax(jj), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h1(x,j,jj), [], [], LB, UB);


    xr=(0: 1e-4: xmax(jj));
       ff=-h1(xr,j,jj);

if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50low(jj) = min(xr(ff>=(vmax(indexs(j))+1)/2));
  endif


endfor
    xmax_m(indexs(j),1)=mean(xmax);

  [nn,xx]=  hist(xmax,35);
  out=[xx',nn'];
  filename=["xmax" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)

 xmax_m(indexs(j),2)=std(xmax);

    EC50low=EC50low(EC50low>1e-6); lenvec=length(EC50low);
    if (~isempty(EC50low))
      EC50low_m(indexs(j),1)=mean(EC50low);
      EC50low_m(indexs(j),2)=std(EC50low);

  [nn,xx]=  hist(EC50low,35);
  out=[xx',nn'];
  filename=["EC50_" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)
    endif


endfor


## model 2


data=0;
  UB =15;
    data=dlmread ("2026_myDF_new5agg.dat",' ',1,0);
  b_max=data(1:5:end,3);

  g_max_m=data(3:5:end,1);
  g_max_sd=data(3:5:end,2);

  indexs=data(1:5:end,1);

  p_b_m=data(2:5:end,1);
  p_b_sd=data(2:5:end,2);

  b_M_m=data(4:5:end,1);
  b_M_sd=data(4:5:end,2);


  data=dlmread ("2026_myDF_new4agg.dat",' ',1,0);
  b_max=data(1:6:end,3);
  g_max_m=data(4:6:end,1);
  b_M_m=data(5:6:end,1);
  indexs=data(1:6:end,1);
  p_b_m=data(2:6:end,1);
  p_g_m=data(3:6:end,1);
  g_max_sd=data(4:6:end,2);
  b_M_sd=data(5:6:end,2);

  p_b_sd=data(2:6:end,2);
  p_g_sd=data(3:6:end,2);

  for j= (4:5)
    # generate random values randn ~ N(0,1) -> p_b_m+p_b_sd*Z
    g_max = randn(nsamp,1)*g_max_sd(j) + g_max_m(j);
    b_M=randn(nsamp,1)*b_M_sd(j)+b_M_m(j);

    p_b = randn(nsamp,1)*p_b_sd(j)+p_b_m(j);
    p_g =randn(nsamp,1)*p_g_sd(j)+p_g_m(j);


    b_M=b_M(g_max>1e-99);
    p_b=p_b(g_max>1e-99);
    p_g=p_g(g_max>1e-99);
    g_max=g_max(g_max>1e-99);

    g_max=g_max(b_M>1e-99);
    p_b=p_b(b_M>1e-99);
    p_g=p_g(b_M>1e-99);
    b_M=b_M(b_M>1e-99);

    g_max=g_max(p_b>1e-99);
    b_M=b_M(p_b>1e-99);
    p_g=p_g(p_b>1e-99);
    p_b=p_b(p_b>1e-99);

    g_max=g_max(p_g>1e-99);
    b_M=b_M(p_g>1e-99);
    p_b=p_b(p_g>1e-99);
    p_g=p_g(p_g>1e-99);

    nsamp0=length(g_max);
      xmax =zeros(1,nsamp0);
      EC50low=zeros(1,nsamp0);

        for jj=1:nsamp0
    [xmax(jj), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h2(x,j,jj), [], [], LB, UB);

         xr=(0: 1e-4: xmax(jj));
       ff=-h2(xr,j,jj);

if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
  EC50low(jj) = min(xr(ff>=(vmax(indexs(j))+1)/2));
endif

endfor


    xmax_m(indexs(j),1)=mean(xmax);

  [nn,xx]=  hist(xmax,35);
  out=[xx',nn'];
  filename=["xmax" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)

 xmax_m(indexs(j),2)=std(xmax);

    EC50low=EC50low(EC50low>1e-6); lenvec=length(EC50low);
    if (~isempty(EC50low))
      EC50low_m(indexs(j),1)=mean(EC50low);
      EC50low_m(indexs(j),2)=std(EC50low);

  [nn,xx]=  hist(EC50low,35);
  out=[xx',nn'];
  filename=["EC50_" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)
    endif
endfor

clearvars b_M p_b p_g g_max
UB=15;
data=0;
   data=dlmread ("2026_myDF_new7agg.dat",' ',1,0);
   indexs=data(1:5:end,1);
   b_max=data(1:5:end,3);

   p_b_m=data(2:5:end,1);
   p_b_sd=data(2:5:end,2);

  g_max_m=data(3:5:end,1);
  g_max_sd=data(3:5:end,2);

  b_M_m=data(4:5:end,1);
  b_M_sd=data(4:5:end,2);

##
  for j= (6:13)
     # generate random values randn ~ N(0,1) -> p_b_m+p_b_sd*Z
    g_max = randn(nsamp,1)*g_max_sd(j) + g_max_m(j);
    b_M=randn(nsamp,1)*b_M_sd(j)+b_M_m(j);
    p_b = randn(nsamp,1)*p_b_sd(j)+p_b_m(j);

    b_M=b_M(g_max>1e-99);
    p_b=p_b(g_max>1e-99);
    g_max=g_max(g_max>1e-99);

    g_max=g_max(b_M>1e-99);
    p_b=p_b(b_M>1e-99);
    b_M=b_M(b_M>1e-99);

    g_max=g_max(p_b>1e-99);
    b_M=b_M(p_b>1e-99);
    p_b=p_b(p_b>1e-99);

    nsamp0=length(g_max);
    p_g= 4*ones(nsamp0,1);
      xmax =zeros(1,nsamp0);
      EC50low=zeros(1,nsamp0);

xr=(0: 1e-4: UB);

        for jj=1:nsamp0
    [xmax(jj), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h2(x,j,jj), [], [], LB, UB);


       ff=-h2(xr,j,jj);

if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
  EC50low(jj) = min(xr(ff>=(vmax(indexs(j))+1)/2));
endif
endfor

    xmax_m(indexs(j),1)=mean(xmax);

  [nn,xx]=  hist(xmax,35);
  out=[xx',nn'];
  filename=["xmax" num2str(indexs(j)) ".dat"];
  dlmwrite(filename,out,"delimiter", " ",'precision',4)

 xmax_m(indexs(j),2)=std(xmax);

    EC50low=EC50low(EC50low>1e-6); lenvec=length(EC50low);
    if (~isempty(EC50low))
      EC50low_m(indexs(j),1)=mean(EC50low);
      EC50low_m(indexs(j),2)=std(EC50low);
        [nn,xx]=  hist(EC50low,35);
        out=[xx',nn'];
        filename=["EC50_" num2str(indexs(j)) ".dat"];
        dlmwrite(filename,out,"delimiter", " ",'precision',4)
    endif
endfor



  out=[EC50low_m,xmax_m];
dlmwrite('Estimates.dat',out,"-append",'precision','%6.4f')

#
function r = h1(x,j0,j)

  r= -b_max(j0).*x.^p_b(jj)./(b_M(j)+x.^p_b(j) ).*exp(-tend.*(g_max(j).*x.^p_g(j)./(g_M(j) + x.^p_g(j) )+mu)) ...
  + (1+b_max(j0).*x.^p_b(j)./(b_M(j)+x.^p_b(j) )).*exp(-tend.*g_max(j).*x.^p_g(j)./(g_M(j) + x.^p_g(j) ));
    r = -r;
endfunction


 #
function r = h2(x,j0,j)

  r= -b_max(j0).*x.^p_b(j)./(b_M(j)+x.^p_b(j) ).*exp(-tend.*(g_max(j).*x.^p_g(j)+mu)) ...
  + (1+b_max(j0).*x.^p_b(j)./(b_M(j)+x.^p_b(j) )).*exp(-tend.*g_max(j).*x.^p_g(j) );
    r = -r;
endfunction

# monotone
function r = h0(x,j0,j)

  r= -b_max(j0).*x.^p_b(j)./(b_M(j)+x.^p_b(j) ).*exp(-tend.*mu) ...
  + (1+b_max(j0).*x.^p_b(j)./(b_M(j)+x.^p_b(j) ));
    r = -r;
endfunction

endfunction
