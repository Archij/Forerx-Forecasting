function ExtMapBuffer1 = JMA(Close, Len, phase)

Bars = length(Close);
ExtMapBuffer1 = zeros(Bars,1);

%variable definitions
jj=0;
ii=0;
series=0;
%----
vv=0;
v1=0;
v2=0;
v3=0;
v4=0;
s8=0;
s10=0;
s18=0;
s20=0;
v5=0;
v6=0;
s28=0;
s30=0;
s38=0;
s40=0;
s48=0;
s50=0;
s58=0;
s60=0;
s68=0;
s70=0;
f8=0;
f10=0;
f18=0;
f20=0;
f28=0;
f30=0;
f38=0;
f40=0;
f48=0;
f50=0;
f58=0;
f60=0;
f68=0;
f70=0;
f78=0;
f80=0;
f88=0;
f90=0;
f98=0;
fA0=0;
fA8=0;
fB0=0;
fB8=0;
fC0=0;
fC8=0;
fD0=0;
f0=0;
fD8=0;
fE0=0;
fE8=0;
fF0=0;
fF8=0;
value2=0;
JMA=0;
prevtime=0;
%----
list = zeros(128,1);
ring1 = zeros(128,1);
ring2 = zeros(11,1);
buffer = zeros(62,1);
%----
s28=63;
s30=64;
for(ii=1:s28)
    list(ii+1)=-1000000;
end
for(ii=s30:127)

    list(ii+1)=1000000;
end
f0=1;

%----
for( shift=1:Bars)
      series=Close(shift);
      if (fF0 < 61)
         fF0= fF0 + 1;
         buffer(fF0+1)=series;
      end %{ main cycle } 
      if (fF0 > 30)
         if (Len < 1.0000000002)
            f80=0.0000000001; %{1.0e-10} 
         else
            f80=(Len - 1)/2.0;
         end
         if (phase < -100)
            f10=0.5;
         else
            if (phase > 100)
               f10=2.5;
            else
               f10=phase/100 + 1.5;
            end
         end
         v1=log(sqrt(f80));
         v2=v1;
         if (v1/log(2.0) + 2.0 < 0.0)
            v3=0;
         else
            v3=v2/log(2.0) + 2.0;
         end
         f98=v3;
%----
         if (0.5<=f98 - 2.0)
            f88=f98 - 2.0;
         else
            f88=0.5;
         end
         f78=sqrt(f80) * f98;
         f90=f78/(f78 + 1.0);
         f80=f80 * 0.9;
         f50=f80/(f80 + 2.0);
%----
         if (f0~=0)
            f0=0;
            v5=0;
            for(ii=1:29)
               if (buffer(ii+2)~=buffer(ii+1))
                  v5=1.0;
               end
            end
            fD8=v5*30.0;
            if (fD8==0)
               f38=series;
            else
               f38=buffer(2);
            end
            f18=f38;
            if (fD8 > 29)
               fD8=29;
            end
         else
            fD8=0;
         end
%----
         for( ii=fD8:-1:0)
            %{ another bigcycle...} 
            value2=31-ii;
            if (ii==0)
               f8=series;
            else
               f8=buffer(value2+1);
            end
            f28=f8 - f18;
            f48=f8 - f38;
            if (abs(f28) > abs(f48))
               v2=abs(f28);
            else
               v2=abs(f48);
            end
            fA0=v2;
            vv=fA0 + 0.0000000001; %{1.0e-10;} 
%----
            if (s48<=1)
               s48=127;
            else
               s48=s48 - 1;
            end
            if (s50<=1)
               s50=10;
            else
               s50=s50 - 1;
            end
            if (s70 < 128)
               s70=s70 + 1;
            end
            s8=s8 + vv - ring2(s50+1);
            ring2(s50+1)=vv;
            if (s70 > 10)
               s20=s8/10;
            else
               s20=s8/s70;
            end
%----
            if (s70 > 127)
               s10=ring1(s48+1);
               ring1(s48+1)=s20;
               s68=64;
               s58=s68;
               while (s68 > 1)
                  if (list(s58+1) < s10)
                     s68=s68*0.5;
                     s58=s58 + s68;
                  else
                     if (list(s58+1)<=s10)
                        s68=1;
                     else
                        s68=s68*0.5;
                        s58=s58 - s68;
                     end
                  end
               end
             else
               ring1(s48+1)=s20;
               if (s28 + s30 > 127)
                  s30=s30 - 1;
                  s58=s30;
               else
                  s28=s28 + 1;
                  s58=s28;
               end
               if (s28 > 96)
                  s38=96;
               else
                  s38=s28;
               end
               if (s30 < 32)
                  s40=32;
               else
                  s40=s30;
               end
            end
%----
            s68=64;
            s60=s68;
            while(s68 > 1)
               if (list(s60+1)>=s20)
                  if (list(s60)<=s20)
                     s68=1;
                  else
                     s68=s68 *0.5;
                     s60=s60 - s68;
                  end
               else
                  s68=s68*0.5;
                  s60=s60 + s68;
               end
               if ((s60==127) && (s20 > list(128)))
                  s60=128;
               end
            end

            if (s70 > 127)
               if (s58>=s60)
                  if ((s38 + 1 > s60) && (s40 - 1 < s60))
                     s18=s18 + s20;
                  else
                     if ((s40 > s60) && (s40 - 1 < s58))
                        s18=s18 + list(s40);
                     end
                  end
               else
                  if (s40>=s60)
                     if ((s38 + 1 < s60) && (s38 + 1 > s58))
                        s18=s18 + list(s38 + 2);
                     end
                  else
                     if (s38 + 2 > s60)
                        s18=s18 + s20;
                     else
                        if ((s38 + 1 < s60) && (s38 + 1 > s58))
                           s18=s18 + list(s38 + 2);
                        end
                     end
                  end
               end

               if (s58 > s60)
                  if ((s40 - 1 < s58) && (s38 + 1 > s58))
                     s18=s18 - list(s58+1);
                  else
                     if ((s38 < s58) && (s38 + 1 > s60))
                        s18=s18 - list(s38+1);
                     end
                  end
               else
                  if ((s38 + 1 > s58) && (s40 - 1 < s58))
                     s18=s18 - list(s58+1);
                  else
                     if ((s40 > s58) && (s40 < s60))
                        s18=s18 - list(s40+1);
                     end
                  end
               end
            end
            if (s58<=s60)
               if (s58>=s60)
                  list(s60+1)=s20;
               else
                  for(jj=s58 + 1:s60 - 1)
                     list(jj)=list(jj+1);
                  end
                  list(s60)=s20;
               end
            else
               for(jj=s58 - 1:-1:s60)
                  list(jj + 2)=list(jj+1);
               end
               list(s60+1)=s20;
            end
            if (s70<=127)
               s18=0;
               for(jj=s40:s38)
                  s18=s18 + list(jj+1);
               end
            end
            f60=s18/(s38 - s40 + 1);
            if (fF8 + 1 > 31)
               fF8=31;
            else
               fF8=fF8 + 1;
            end
%----
            if (fF8<=30)
               if (f28 > 0)
                  f18=f8;
               else
                  f18=f8 - f28 * f90;
               end
               if (f48 < 0)
                  f38=f8;
               else
                  f38=f8 - f48 * f90;
               end
               fB8=series;
%{EasyLanguage does not have "continue" statement} 
               if (fF8~=30)
                  continue;
               end
               if (fF8==30)
                  fC0=series;
                  if (ceil(f78)>=1)
                     v4=ceil(f78);
                  else
                     v4=1;
                  end
                  fE8=ceil(v4);
                  if (floor(f78)>=1)
                     v2=floor(f78);
                  else
                     v2=1;
                  end
                  fE0=ceil(v2);
                  if (fE8==fE0)
                     f68=1;
                  else
                     v4=fE8 - fE0;
                     f68=(f78 - fE0)/v4;
                  end
                  if (fE0<=29)
                     v5=fE0;
                  else
                     v5=29;
                  end
                  if (fE8<=29)
                     v6=fE8;
                  else
                     v6=29;
                  end
                  fA8=(series - buffer(fF0 - v5+1)) * (1 - f68)/fE0 + (series - buffer(fF0 - v6+1)) * f68/fE8;
               end
            else
               if (f98>=power(fA0/f60, f88))
                  v1=power(fA0/f60, f88);
               else
                  v1=f98;
               end
               if (v1 < 1)
                  v2=1;
               else
                  if (f98>=power(fA0/f60, f88))
                     v3=power(fA0/f60, f88);
                  else
                     v3=f98;
                  end
                  v2=v3;
               end
               f58=v2;
               f70=power(f90, sqrt(f58));
               if (f28 > 0)
                  f18=f8;
               else
                  f18=f8 - f28 * f70;
               end
               if (f48 < 0)
                  f38=f8;
               else
                  f38=f8 - f48 * f70;
               end
            end
         end
         if (fF8 > 30)
            f30=power(f50, f58);
            fC0=(1 - f30) * series + f30 * fC0;
            fC8=(series - fC0) * (1 - f50) + f50 * fC8;
            fD0=f10 * fC8 + fC0;
            f20=-f30 * 2;
            f40=f30 * f30;
            fB0=f20 + f40 + 1;
            fA8=(fD0 - fB8) * fB0 + f40 * fA8;
            fB8=fB8 + fA8;
         end
         JMA= fB8;
      end
      if (fF0<=30)
            JMA=0;
      end
%Print ("JMA is " + JMA + " shift is " + shift); 
      ExtMapBuffer1(shift)=JMA;
%----
end

%+------------------------------------------------------------------+

end