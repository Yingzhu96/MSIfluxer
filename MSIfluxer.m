%MSIdata is the raw data
%library is the metabolite isotopologue
%%
%step1
r1=size(MSIdata,1)%the row number of MSIdata
c1=size(MSIdata,2)%the cloumn number of MSIdata
r2=size(library,1)%the row number of library
c2=size(library,2)%the cloumn number of library
MSIfluxer = cell(r2,c1/2+4)
for i=1:r2 %the row number of library
MSIfluxer(i,1)=library(i,1)%the first cloumn is metabolite name
MSIfluxer(i,2)=library(i,2)%the second cloumn is adduct form
MSIfluxer(i,3)=library(i,3)%the third cloumn is isotopologue
MSIfluxer(i,4)=library(i,4)%the fourth cloumn is m/z
end
for j=1:2:c1 %the cloumn number of MSIdata
MSIfluxer(1,(j+1)/2+4)=MSIdata(1,j)
for i=2:r2
MSIfluxer(i,(j+1)/2+4)={0}
end
end
%%
%step2
for a=1:2:c1 %MSIdata m/z column
    for i=2:r2 %the row number of library
        for j=3:r1 %the row number of MSIdata
        theor=library{i,4}%%the fourth cloumn is m/z
        mea=MSIdata{j,a}
        if abs(theor-mea)/theor<5*10^(-6) %m/z error within 5ppm
            MSIfluxer(i,(a+1)/2+4)=MSIdata(j,a+1); 
        end
        end
    end
end
clear a r1 r2 c1 c2 i j theor mea
clc