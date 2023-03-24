%This script downloads netCDF HFR totals from hfrnet, constrains the data
%to a region of interest and safes the data as .mat
%Uchenna Chizaram Nwankwo
clear
close all
clc
%%
hmDir=pwd;
% cd Half

disp('This program downloads HFR totals from hfrnet and saves the result .mat')

disp('Input spatial limits of area of interest')
disp('Input latitude northern limit')
n = true; 
Lt = 90;
while n
    latN=input('Input northern limit of the area of interest (-90 to 90). \n');
    if latN<=Lt && latN>=-Lt
       n = false;
    end
end

disp('Input latitude southern limit')
n = true; 
while n
    latS=input('Input southern limit of the area of interest (-90 to 90). \n');
    if latS<=Lt && latS>=-Lt && latS<latN
       n = false;
    end
end

disp('Input longitude eastern limit')
n = true; 
Ln = 180;
while n
    latE=input('Input eastern limit of the area of interest (-180 to 180). \n');
    if latE<=Ln && latE>=-Ln
       n = false;
    end
end

disp('Input latitude western limit')
n = true; 
while n
    latW=input('Input western limit of the area of interest (-180 to 180). \n');
    if latW<=Ln && latW>=-Ln && latW<latE
       n = false;
    end
end

% latS = 29.02; latN = 29.93;     % Half area of Gulf coast
% lonE = -86.51;  lonW = -87.84;
%%
disp('Enter begin and end dates')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Enter begin date:') 

%Inputting date of interest
n = true;
S = date; yr = str2double(S(end-3:end)); 
while n
    nyrb=input('Input integer value for year of interest. \n');
    if nyrb<=yr && nyrb>0
       n = false;
    end
end

n = true;
while n
    nmonb=input('Input integer value for month of interest. \n');
    if nmonb<=12 && nmonb>=1
       n = false;
    end
end 

n = true;
E = eomday(nyrb,nmonb);
while n
    ndayb=input('Input integer value for end day of interest. \n');
    if ndayb<=E && ndayb>=1
       n = false;
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Enter end date:')

n = true;
while n
    nyre=input('Input integer value for year of interest. \n');
    if nyre<=yr && nyre>0
       n = false;
    end
end 

n = true;
while n
    nmone=input('Input integer value for month of interest. \n');
    if nmone<=12 && nmone>=1
       n = false;
    end
end

n = true;
E = eomday(nyre,nmone);
while n
    ndaye=input('Input integer value for end day of interest. \n');
    if ndaye<=E && ndaye>=1
       n = false;
    end
end
%%
if (datenum([nyrb nmonb ndayb]) > datenum([nyre nmone ndaye]))
    error('Input Error: End date is earlier than begin date')
end
dates=datenum([nyrb nmonb ndayb]) :1: datenum([nyre nmone ndaye]);
dVec=datevec(dates);%Then to datevector
cmon=num2str(dVec(:,2),'%2.2d');
cday=num2str(dVec(:,3),'%2.2d');
chr=num2str((0:1:23)','%2.2d');

U_comp = nan(23,16,length(dVec)*24);%Three dimensional grid 
V_comp = nan(23,16,length(dVec)*24);
T_comp = nan(1,length(dVec)*24);
Udop_comp = nan(23,16,length(dVec)*24); 
Vdop_comp = nan(23,16,length(dVec)*24);
   
%%
reg = {'PRVI', 'USEGC', 'USHI', 'USWC'};
reg_len = length(reg);
n = true;
while n
    dy=input('type\n 1 for Puerto Rico and Virgin Islands\n2 for U.S. East and Gulf Coast\n3 for Hawaii\n 4 for U.S. West Coast.\n');
    if dy<=reg_len && dy>=1
       n = false;
    end
end 
%%
yx = 1;
for i = 1:length(dVec)%Looping through months
    for q = 1:24%looping through hours
        ftpobj = ftp('ftp-oceans.ncei.noaa.gov');%ftp('ftp.nodc.noaa.gov')
        cd(ftpobj,['/pub/data.nodc/ndbc/hfradar/rtv/' num2str(dVec(i,1)) '/' num2str(dVec(i,1)) cmon(i,:) '/' upper(reg{dy}) '/']);
   
        mget(ftpobj, [num2str(dVec(i,1)) cmon(i,:) cday(i,:) chr(q,:) '00_hfr_' lower(reg{dy}) '_6km_rtv_uwls_NDBC.nc' ]);

        file = dir('*.nc');

        filename = file.name;

        T_comp(1,yx)=datenum([dVec(i,1) str2double(cmon(i,:)) str2double(cday(i,:)) str2double(chr(q,:)) 0 0]);
        
        lat=ncread(filename,'lat');

        lon=ncread(filename,'lon');

        % find lat/lon indices for spatial domain of interest
        lats=find((lat > latS) & (lat < latN)); 
        lons=find((lon > lonW) & (lon < lonE));

        % NOTE: To check correct domain identification execute :
        % lat(lats) and/or lon(lons)

        [ilat,~]=size(lats);
        [ilon,~]=size(lons);

        nlon=lons(ilon)-lons(1)+1;% Determine number of lat and lon indices needed to cover domain of interest
        nlat=lats(ilat)-lats(1)+1;

        U_comp(:,:,yx)=flip(ncread(filename,'u',[ lons(1) lats(1) 1 ], [ nlon nlat 1 ] )); 

        V_comp(:,:,yx)=flip(ncread(filename,'v',[ lons(1) lats(1) 1 ], [ nlon nlat 1 ] ));

        Udop_comp(:,:,yx)=flip(ncread(filename,'DOPx',[ lons(1) lats(1) 1 ], [ nlon nlat 1 ] ));

        Vdop_comp(:,:,yx)=flip(ncread(filename,'DOPy',[ lons(1) lats(1) 1 ], [ nlon nlat 1 ] ));

        lat = lat(lats);
        
        lon = lon(lons);
        
        close(ftpobj)
        delete *.nc
        yx = yx+1;
    end
end
hwdir = pwd;
save ([hwdir,'rtv.mat'],'U_comp','V_comp','Udop_comp','Vdop_comp','T_comp','lat','lon')
