%PlotXppaut
%Plot bifurcation diagrams in Matlab that have been saved by XPPAUT
%
%How to use:
%In AUTO (XPPAUT) click-->File--->Write Pts
%and save your bifurcation diagram in a .dat file
%In Matlab command window type: plotxppaut 
%Select .dat file saved by AUTO
%A bifurcation diagram will be plotted
%Cutomize the plot using Matlab figure editor
%Or edit program parameters to suit your needs--very useful for those who
%plot lots of bifurcation diagrams
%
%Version 1.0 Tested Under Matlab 6.1 (R12.1)

%---------------------------------------
%Version 1.0
%Time-stamp: <2005-05-03 17:06:11 msi800> 
%-------------------------------------------
%Mohammad S. Imtiaz
%Room 403
%The Neuroscience Group
%Discipline of Human Physiology
%School of Biomedical Sciences
%Faculty of Health
%The University of Newcastle
%Callaghan, NSW 2308
%Australia
%Tel: +61 02 49217857  (Office)
%Tel: +61 02 49215626  (Lab)
%Fax: +61 02 49217406
%Email: Mohammad.Imtiaz@newcastle.edu.au
%-------------------------------------------
%-------------------------------------------

% clear;clc

%-------------------------------------------------------------------
%       USER DEFINED PARAMETERS -- CHANGE TO CUSTOMIZE OUTPUT PLOT
%-------------------------------------------------------------------
%Colors -- SHOULD BE VALID MATLAB COLORS
C_Ss = 'Black';   %STABLE STEADY STATE
C_Us = 'Black';   %UNSTABLE STEADY STATE
C_Sp = 'Red';     %STABLE PERIODIC ORBIT
C_Up = 'Blue';    %UNSTABLE PERIODIC ORBIT

%LineStyle -- SHOULD BE VALID MATLAB LINESTYLES
Lt_Ss = '-';    %STABLE STEADY STATE
Lt_Us = '--';   %UNSTABLE STEADY STATE
Lt_Sp = '-';    %STABLE PERIODIC ORBIT
Lt_Up = '--';    %UNSTABLE PERIODIC ORBIT

%Linewidths
Lw_Ss = 1.5;   %STABLE STEADY STATE
Lw_Us = 1;     %UNSTABLE STEADY STATE
Lw_Sp = 1;     %STABLE PERIODIC ORBIT
Lw_Up = 1;     %UNSTABLE PERIODIC ORBIT

%Markersize
Ms_Ss = 2;    %STABLE STEADY STATE
Ms_Us = 2;    %UNSTABLE STEADY STATE
Ms_Sp = 3;    %STABLE PERIODIC ORBIT
Ms_Up = 3;    %UNSTABLE PERIODIC ORBIT

%MarkerFaceColor -- SHOULD BE VALID MATLAB COLORS
Mfc_Ss = 'none';   %STABLE STEADY STATE
Mfc_Us = 'none';   %UNSTABLE STEADY STATE
Mfc_Sp = 'Red';    %STABLE PERIODIC ORBIT
Mfc_Up = 'none';   %UNSTABLE PERIODIC ORBIT

%	-----------------------------------------	
% END OF USER DEFINED PARAMETERS
%	-----------------------------------------






[file_in,path] = uigetfile('*.dat','.dat file saved by AUTO (XPPAUT) ');
file_name = [path file_in];

if(file_name)
    st = importdata(file_name);

%     [file_in2,path2] = uigetfile('*.dat','.dat file saved by AUTO (XPPAUT) ');
%     file_name2 = [path2 file_in2];
%     st2 = importdata(file_name2);
  
%   temp=[st; st2];
    temp = st;

%     for i = 1:length(temp)
%         if temp(i, 2) < 1e-3
%             temp(i, 2) = 1e-3;
%         end
%         if temp(i, 3) <1e-3
%             temp(i, 3) = 1e-3;
%         end
%     end

%   threshold_bool = false(size(temp, 1), 1);
%   for i=1:size(temp, 1)
%       if (temp(i, 1) > 0.2) || (temp(i, 1) < 0)
%            threshold_bool(i) = true;
%       end
%   end
%   temp(threshold_bool, :) = [];
  
  temp1=zeros(size(temp,1),8)*NaN;
  
  for n=1:size(temp1,1)
    if(temp(n,4)==1)
      temp1(n,[1 5]) = temp(n,[2 3]);
    end
    if(temp(n,4)==2)
      temp1(n,[2 6]) = temp(n,[2 3]);
    end
    if(temp(n,4)==3)
      temp1(n,[3 7]) = temp(n,[2 3]);
    end
    if(temp(n,4)==4)
      temp1(n,[4 8]) = temp(n,[2 3]);
    end
  end
  
  
% $$$ figure;
%	-----------------------------------------  
  h_SS = plot(temp(:,1),temp1(:,[1 5]));
  set(h_SS,'color',C_Ss)
  set(h_SS,'linestyle',Lt_Ss)
  set(h_SS,'linewidth',Lw_Ss)
  set(h_SS,'Markersize',Ms_Ss)
  set(h_SS,'MarkerFacecolor',Mfc_Ss)
  hold on;

  
%	-----------------------------------------  
  h_US = plot(temp(:,1),temp1(:,[2 6]));
  set(h_US,'color',C_Us)
  set(h_US,'linestyle',Lt_Us)
  set(h_US,'linewidth',Lw_Us)
  set(h_US,'Markersize',Ms_Us)
  set(h_US,'MarkerFacecolor',Mfc_Us)

%	-----------------------------------------  
  

h_SP = plot(temp(:,1),temp1(:,[3 7]));
  set(h_SP,'color',C_Sp)
  set(h_SP,'linestyle',Lt_Sp)
  set(h_SP,'linewidth',Lw_Sp)
  set(h_SP,'Markersize',Ms_Sp)
  set(h_SP,'MarkerFacecolor',Mfc_Sp)
%	-----------------------------------------   
  
  
  h_UP = plot(temp(:,1),temp1(:,[4 8]));
  set(h_UP,'color',C_Up)
  set(h_UP,'linestyle',Lt_Up)
  set(h_UP,'linewidth',Lw_Up)
  set(h_UP,'Markersize',Ms_Up)
  set(h_UP,'MarkerFacecolor',Mfc_Up)
%	-----------------------------------------    
  
  
  
  set(gca, 'FontSize', 26); box on;
%   set(gca, 'YScale', 'log');
%   xlim([0 3000]);
  xlabel('c_{M-IL10}');
  ylabel('CCL2');
%   hold off;
%   title(file_in);
end
% 
% % 
% C_Ss = 'Red';   %STABLE STEADY STATE
% C_Us = 'Red';   %UNSTABLE STEADY STATE
% 
% 
% 
% [file_in,path] = uigetfile('*.dat','.dat file saved by AUTO (XPPAUT) ');
% file_name = [path file_in];
% 
% if(file_name)
%     st = importdata(file_name);
% 
% %     [file_in2,path2] = uigetfile('*.dat','.dat file saved by AUTO (XPPAUT) ');
% %     file_name2 = [path2 file_in2];
% %     st2 = importdata(file_name2);
%   
% %   temp=[st; st2];
%     temp = st;
% 
% %     for i = 1:length(temp)
% %         if temp(i, 2) < 1e-3
% %             temp(i, 2) = 1e-3;
% %         end
% %         if temp(i, 3) <1e-3
% %             temp(i, 3) = 1e-3;
% %         end
% %     end
% 
% %   threshold_bool = false(size(temp, 1), 1);
% %   for i=1:size(temp, 1)
% %       if (temp(i, 1) > 0.2) || (temp(i, 1) < 0)
% %            threshold_bool(i) = true;
% %       end
% %   end
% %   temp(threshold_bool, :) = [];
%   
%   temp1=zeros(size(temp,1),8)*NaN;
%   
%   for n=1:size(temp1,1)
%     if(temp(n,4)==1)
%       temp1(n,[1 5]) = temp(n,[2 3]);
%     end
%     if(temp(n,4)==2)
%       temp1(n,[2 6]) = temp(n,[2 3]);
%     end
%     if(temp(n,4)==3)
%       temp1(n,[3 7]) = temp(n,[2 3]);
%     end
%     if(temp(n,4)==4)
%       temp1(n,[4 8]) = temp(n,[2 3]);
%     end
%   end
%   
%   
% % $$$ figure;
% %	-----------------------------------------  
%   h_SS = plot(temp(:,1),temp1(:,[1 5]));
%   set(h_SS,'color',C_Ss)
%   set(h_SS,'linestyle',Lt_Ss)
%   set(h_SS,'linewidth',Lw_Ss)
%   set(h_SS,'Markersize',Ms_Ss)
%   set(h_SS,'MarkerFacecolor',Mfc_Ss)
%   hold on;
% 
%   
% %	-----------------------------------------  
%   h_US = plot(temp(:,1),temp1(:,[2 6]));
%   set(h_US,'color',C_Us)
%   set(h_US,'linestyle',Lt_Us)
%   set(h_US,'linewidth',Lw_Us)
%   set(h_US,'Markersize',Ms_Us)
%   set(h_US,'MarkerFacecolor',Mfc_Us)
% 
% %	-----------------------------------------  
%   
% 
% h_SP = plot(temp(:,1),temp1(:,[3 7]));
%   set(h_SP,'color',C_Sp)
%   set(h_SP,'linestyle',Lt_Sp)
%   set(h_SP,'linewidth',Lw_Sp)
%   set(h_SP,'Markersize',Ms_Sp)
%   set(h_SP,'MarkerFacecolor',Mfc_Sp)
% %	-----------------------------------------   
%   
%   
%   h_UP = plot(temp(:,1),temp1(:,[4 8]));
%   set(h_UP,'color',C_Up)
%   set(h_UP,'linestyle',Lt_Up)
%   set(h_UP,'linewidth',Lw_Up)
%   set(h_UP,'Markersize',Ms_Up)
%   set(h_UP,'MarkerFacecolor',Mfc_Up)
% %	-----------------------------------------    
%   
%   
%   
%   set(gca, 'FontSize', 26); box on;
% %   set(gca, 'YScale', 'log');
%   xlim([0 2]);
%   xlabel('c_{M-IL10}');
%   ylabel('CCL2');
%   hold off;
% %   title(file_in);
% end