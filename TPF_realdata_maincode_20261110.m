clc;clear;close all;
%% read data
load('data_profile05.mat');
d{1}=data_profile05(:,1);%depth
t1{1}=data_profile05(:,2);% normal temperature profile;
c1{1}=data_profile05(:,3);% normal conductivity profile;
t2{1}=data_profile05(:,4);% spike temperature profile;
c2{1}=data_profile05(:,5);% spike conductivity profile;
%% ²ÎÊı
resolution=0.0004;%sbe911 conductivity resolution 0.0004S/m
N=30;%truncation order
for i=1:1% profile num
    for h=1:numel(N)%truncation order index
        j=1;
        ccop{j}=c2{i};%profile with spike
        pcop{j}=d{i};
        range_all=false(length(pcop),1);
        stend{i}=[];
        while j
            [f_result{j},rr{j}]=FLYNH(pcop{j},ccop{j},N(h));%fit
            R(j)=max(abs(rr{j}));%max residual
            stfr=find(abs(rr{j})==R(j));% index of max residual
            stend{i}=[stend{i};stfr(:)];
            %%%%%%%%%%%%%% replaced with fitted value
            for ri = 1:numel(stfr)
                st = stfr(ri);
                if st == 1
                    % if the index of spike is 1,replaced with its fitted value
                    ccop{j}(st) = f_result{j}(st);
                else
                    % else,replaced with the mean of 12 proceeding fitted value
                    s = max(1, st-12);
                    yu = f_result{j}(s:st-1);
                    ccop{j}(st) = mean(yu);
                end
            end
            %%%%%%%%%%%%%%%%% stopping criteria for fitting 
            if j>=3 && (stend{i}(j)==stend{i}(j-1)) && (stend{i}(j-1))==stend{i}(j-2)
                st_final=j-3;%Final fitting times?
                break;
            end
            if R(j)<=resolution
                break;
            end
            j=j+1;
            ccop{j}=ccop{j-1};
            pcop{j}=pcop{j-1};
            
            %% plot
            if j==2
                figure('color',[1 1 1],'position',[100 100 900 900]);
                plot(c2{i}/10,d{i},'k.');
                set(gca,'ydir','reverse','fontsize',15);
                hold on
                xlabel('Conductivity(S/m)');
                ylabel('Depth(m)');
                title(['Profile' num2str(i,'%02d')]);
            end
            plot(c2{i}(stfr)/10,d{i}(stfr),'ro');
        end
        % % % %  Replace all outliers with the last fitted data
        for ri=1:numel(stend{i})
            ccop{j}(stend{i}(ri))=f_result{j-1}(stend{i}(ri));
        end
        h3yu=plot(ccop{j}(stend{i})/10,d{i}(stend{i}),'bo');
        hyu=plot(f_result{j-1}/10,d{i},'g-','linewidth',1.5);
        numel(unique(stend{i}))
        clear ccop pcop range_all  rr R stfr st_final
    end
    clear idx
end