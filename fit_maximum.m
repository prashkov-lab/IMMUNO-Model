function out=fit_maximum

mu=log(2)/2.5; tend=6;  x0=mu;

vmax = [2.04095779756017, 1.59692849121587, 1.59692849121587, 1.79741691222017,...
    2.04095779756017, 2.52100523745854, 2.19599418149267, 2.39755841682793, ...
    2.39755841682793, 2.29732154535348, 2.29732154535348, 2.39755841682793, ...
    2.29732154535348, 2.52100523745854];
  xmax=zeros(14,2);
  EC50=zeros(14,2);

## initial guess and lower/upper bounds
   X0 = 0.1;
  LB= 0.0001;
  UB =5 ;


  # first scenario

  data=dlmread("2026_myDF_new1agg.dat",' ',1,0);
  b_max=data(1:7:end,3);
  g_max=data(4:7:end,1);
  b_M=data(5:7:end,1);
  g_M=data(6:7:end,1);
  indexs=data(1:7:end,1);
  n=data(2:7:end,1);
  m=data(3:7:end,1);

  for j=1
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h1_new(x,j), [], [], LB, UB);

     xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h1_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end

  end


  data=dlmread("2026_myDF_new3agg.dat",' ',1,0);
  b_max=data(1:6:end,3);
  g_max=data(3:6:end,1);
  b_M=data(4:6:end,1);
  g_M=data(5:6:end,1);
  indexs=data(1:6:end,1);
  n=data(2:6:end,1);
    len=length(b_max);
  m=3*ones(len,1);


  for j= 2
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h1_new(x,j), [], [], LB, UB);

     xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h1_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end

  end

  data=dlmread("2026_myDF_new3agg2.dat",' ',1,0);
  b_max=data(1:6:end,3);
  g_max=data(3:6:end,1);
  b_M=data(4:6:end,1);
  g_M=data(5:6:end,1);
  indexs=data(1:6:end,1);
  n=data(2:6:end,1);
    len=length(b_max);
  m=4*ones(len,1);

  for j=6
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h1_new(x,j), [], [], LB, UB);

     xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h1_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end

  end

## model 2

  UB =5;

    data=dlmread ("2026_myDF_new4agg.dat",' ',1,0);
  b_max=data(1:6:end,3);
  g_max=data(4:6:end,1);
  b_M=data(5:6:end,1);
  indexs=data(1:6:end,1);
  n=data(2:6:end,1);
  m=data(3:6:end,1);

  for j= (4:5)
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h2_new(x,j), [], [], LB, UB);

         xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h2_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end

end

  data=dlmread ("2026_myDF_new5agg.dat",' ',1,0);
  b_max=data(1:5:end,3);
  g_max=data(3:5:end,1);
  b_M=data(4:5:end,1);
  indexs=data(1:5:end,1);
  n=data(2:5:end,1);
  len=length(b_max);
  m=2*ones(len,1);

  for j= 3
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h2_new(x,j), [], [], LB, UB);

         xr=(0: 1e-6: xmax(indexs(j),1));
       ff=-h2_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end


  end

    data=dlmread ("2026_myDF_new7agg.dat",' ',1,0);
  b_max=data(1:5:end,3);
  g_max=data(3:5:end,1);
  b_M=data(4:5:end,1);
  indexs=data(1:5:end,1);
  n=data(2:5:end,1);
  len=length(b_max);
  m=4*ones(len,1);

  for j= (7:13)
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h2_new(x,j), [], [], LB, UB);

         xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h2_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end
  end

## for monotone
UB=50;
  data=dlmread ("2026_myDF_new0agg.dat",' ',1,0);
  b_max=data(1:4:end,3)
  b_M=data(3:4:end,1);
  indexs=data(1:4:end,1);
  n=data(2:4:end,1);

  for j= (1:2)
    [xmax(indexs(j),1), objfun, INFO, ITER, NF, LAMBDA] = sqp (X0, @(x)h0_new(x,j), [], [], LB, UB);

         xr=(0: 1e-4: xmax(indexs(j),1));
       ff=-h0_new(xr,j);
if (~ isempty(xr(ff>=(vmax(indexs(j))+1)/2)))
    EC50(indexs(j),1) = min(xr(ff>=(vmax(indexs(j))+1)/2));
end
  end


## for saturating
function r = h1_new(x,j)

  r= -b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) ).*exp(-tend.*(g_max(j).*x.^m(j)./(g_M(j) + x.^m(j) )+mu)) ...
  + (1+b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) )).*exp(-tend.*g_max(j).*x.^m(j)./(g_M(j) + x.^m(j) ));
    r = -r;
end


## power law for gamma
function r = h2_new(x,j)

  r= -b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) ).*exp(-tend.*(g_max(j).*x.^m(j)+mu)) ...
  + (1+b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) )).*exp(-tend.*g_max(j).*x.^m(j) );
    r = -r;
end

## monotone
function r = h0_new(x,j)

  r= -b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) ).*exp(-tend.*mu) ...
  + (1+b_max(j).*x.^n(j)./(b_M(j)+x.^n(j) ));
    r = -r;
end


out=[xmax(:,1),EC50];
dlmwrite('Estimates2026.dat',out,'precision','%6.4f')
end
