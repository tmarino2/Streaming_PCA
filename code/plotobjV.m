%% PLOTOBJV(dataname,methods,PCA,k,numpasses,ITER,AVG,PLOTREPORT) plots
% PROGRESS-PER-ITERATION (empirical & population) and PROGRESS-PER-SECOND.
%
%  The inputs are as follows:
%
%  dataname is one of these strings: 'VidTIMIT', 'VidTIMIT2', 'SIM'
%  method is a cell of strings: e.g. {'sgd', 'brand', 'batch', 'truth'}
%
%  PCA is a boolean flag - set to 1 if you want to plot PCA, 0 if CCA
%
%  k a positive integer - denotes the desired RANK
%
%  numpasses is a positive integer - denotes number of passes over the data
%
%  ITERS is an array of positve integers represening which random splits
%  use only (1-1000)
%
%  AVG is either an empty string '' in which all iterations are plotted
%  simultaneously or 'avg' in which case all iterations are avergaged
%
%  PLOTREPORT is a boolean flag - set to 0 if you want to generate a report
%  on the cluster - set to 1 if you also want to plot the report
%
%  The output containing the report is written to ../REPORT as a MAT file
%  in the following format: reportPCA[method=sgd,rank=4,numiter=10].mat
%
%  If PLOTREPORT is set to 1, three plots are generated and written to
%  ../PLOTS as pdf files with names in the following formats:
%
%  iteration[dataname=VidTIMIT,rank=1,numiter=1,...
%  methods=sgd,brand]
%  convergence[dataname=VidTIMIT,rank=1,numiter=1,...
%  methods=sgd,brand].pdf
%
%%

function plotobjV(dataname,methods,k,n,numiter,etas,lambdas,betas,method_names)

if length(etas)==1
    etas=etas*ones(size(methods));
end
if length(lambdas)==1
    lambdas=lambdas*ones(size(methods));
end
if length(betas)==1
    betas=betas*ones(size(methods));
end

if(strcmp(dataname,'MNIST'))
     dnam='mnist';
elseif(strcmp(dataname,'syn_exp_decay'))
     dnam='syn';
else
    dnam='other';
end

%% Setup Figures
col={'k','g','b','c','m','r','y'};
marker={'ks','k^','ks','kd','kh','ko','kx'};
fig11=figure(11); clf; set(fig11,'Position',[2 2 1200 800]);
fig22=figure(22); clf; set(fig22,'Position',[2 2 1200 800]);
fig33=figure(33); clf; set(fig33,'Position',[2 2 1200 800]);

%% File names of figures to be plotted
fnames=cell(3,1);
%hackaries
gap = floor(10*lambdas(1)*2);

fnames{1}=sprintf(['../plots/%s_convergence_%d.pdf'],dnam,k);
fnames{2}=sprintf(['../plots/%s_rank_%d.pdf'],dnam,k);
fnames{3}=sprintf(['../plots/%s_time_%d.pdf'],dnam,k);

LWIDTH=6;
MSIZE=22;

maxrank=0; maxobj=0; maxtime=0; maxruntime = 0; minruntime = 10000000;

method='truth';
reportfile=sprintf(['../page/profile/PCA/%s/%s/%s[',...
    'rank=%d].mat'],method,dataname,method,k);

load(reportfile,'objV');
trueobjV=objV;


%% Plot for each method
for imethod=1:length(methods)
    
    method=methods{imethod};
    eta=etas(imethod);
    lambda=lambdas(imethod);
    beta=betas(imethod);
    
    %% IF reading from or writing to a report file
    if(sum(strcmp(method,{'batch','incremental'})))
        reportfile=sprintf(['../page/reportPCA[data=%s,method=%s,',...
            'rank=%d,numiter=%d].mat'],...
            dataname,method,k,numiter);
    elseif(sum(strcmp(method,{'l2rmsg','l2proj'})))
        reportfile=sprintf(['../page/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,lambda=%f,numiter=%d].mat'],...
            dataname,method,k,eta,lambda,numiter);
    elseif(sum(strcmp(method,{'l1rmsg'})))
        reportfile=sprintf(['../page/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,beta=%f,numiter=%d].mat'],...
            dataname,method,k,eta,beta,numiter);
    elseif(sum(strcmp(method,{'l21rmsg'})))
        reportfile=sprintf(['../page/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,lambda=%f,beta=%f,numiter=%d].mat'],...
            dataname,method,k,eta,lambda,beta,numiter);
    else
        reportfile=sprintf(['../page/reportPCA[data=%s,method=%s,',...
            'rank=%d,eta=%f,numiter=%d].mat'],...
            dataname,method,k,eta,numiter);
    end
    
    %% Fetch data if plotting from a report file
    if(exist(reportfile,'file'))
        load(reportfile,'seq','progress','avgprogress',...
            'rank_of_iters','avgrank','comptime','avgtime');
        maxrank=max(max(avgrank),maxrank);
        maxobj=max(max(avgprogress),maxobj);
        maxtime=max(max(avgtime),maxtime);
    else
        %% Set PAGE path and PAGE prefix
        pagepath=sprintf('../page/profile/PCA/%s/%s/',...
            method,dataname);
        
        
        if(sum(strcmp(method,{'batch','incremental','incfd'})))
            pageprefix=@(method,rank,iter)[pagepath,...
                sprintf('%s[rank=%d,iter=%d].mat',...
                method,rank,iter)];
        elseif(sum(strcmp(method,{'l2rmsg','l2proj'})))
            pageprefix=@(method,rank,iter)[pagepath,...
                sprintf('%s[rank=%d,eta=%f,lambda=%f,iter=%d].mat',...
                method,rank,eta,lambda,iter)];
        elseif(sum(strcmp(method,{'l1rmsg'})))
            pageprefix=@(method,rank,iter)[pagepath,...
                sprintf('%s[rank=%d,eta=%f,beta=%f,iter=%d].mat',...
                method,rank,eta,beta,iter)];
        elseif(sum(strcmp(method,{'l21rmsg'})))
            pageprefix=@(method,rank,iter)[pagepath,...
                sprintf('%s[rank=%d,eta=%f,lambda=%f,beta=%f,iter=%d].mat',...
                method,rank,eta,lambda,beta,iter)];
        else
            pageprefix=@(method,rank,iter)[pagepath,...
                sprintf('%s[rank=%d,eta=%f,iter=%d].mat',...
                method,rank,eta,iter)];
        end
        
        %% Set sequence points based on dataset and datasize
        [seq,L]=equilogseq(n,1);
        
        %% Initialize performance metrics
        progress=zeros(L(2),length(numiter));
        rank_of_iters=zeros(L(2),length(numiter));
        comptime=zeros(L(2),length(numiter));
        
        %% Gather performance metrics
        for iter=1:numiter
            load(pageprefix(method,k,iter),...
                'objV','objVtune','rk','runtime');
            progress(:,iter)=objV(1:L(2));
            rank_of_iters(:,iter)=rk(1:L(2));
            comptime(:,iter)=runtime(1:L(2));
        end
        
        avgprogress=sum(progress,2)./numiter;   
        avgrank=sum(rank_of_iters,2)./numiter;
        avgtime=sum(comptime,2)./numiter;
        save(reportfile,'seq','progress','avgprogress',...
            'rank_of_iters','avgrank','comptime','avgtime');
        maxrank=max(max(avgrank),maxrank);
        maxobj=max(max(avgprogress),maxobj);        
        maxtime=max(max(avgtime),maxtime);
    end
    avgprogress=trueobjV-avgprogress;
    avgprogress(avgprogress<0)=1e-6;
    runtimetotal=cumsum(avgtime);
    minruntime = min( minruntime, runtimetotal( 1 ) );
    maxruntime = max( maxruntime, runtimetotal( end ) );

    fig11=figure(11);
    ignore = semilogx(seq,avgprogress,'Color',...
        col{imethod},'LineWidth',LWIDTH);
    hold on;
    set( get( get( ignore, 'Annotation' ), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
    subseq = mod(round( (1:15) * length(seq) / 15 + imethod),length(seq))+1;
    semilogx(seq(subseq),avgprogress(subseq),marker{imethod},...
        'MarkerFaceColor',col{imethod},'MarkerSize',MSIZE);
    hold on;
    
    fig22=figure(22);
    ignore = semilogx(seq,avgrank,'Color',...
        col{imethod},'LineWidth',LWIDTH);
    hold on;
    set( get( get( ignore, 'Annotation' ), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
    subseq = mod(round( (1:15) * length(seq) / 15 + imethod),length(seq))+1;
    semilogx(seq(subseq),avgrank(subseq),marker{imethod},...
        'MarkerFaceColor',col{imethod},'MarkerSize',MSIZE);
    hold on;
        
    fig33=figure(33);
    ignore = semilogx(runtimetotal,avgprogress,'Color',...
        col{imethod},'LineWidth',LWIDTH);
    hold on;
    set( get( get( ignore, 'Annotation' ), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
    subseq = mod(round( (1:15) * length(seq) / 15 + imethod),length(seq))+1;
    semilogx(runtimetotal(subseq),avgprogress(subseq),marker{imethod},...
        'MarkerFaceColor',col{imethod},'MarkerSize',MSIZE);
    hold on;
    
end


FSIZE1=40; %70;
FSIZE2=50; %70;
FSIZE3=40; %40;

figure(11);
hold on; grid;
xlabel('Iteration','FontSize',FSIZE2);
maxobj=max(maxobj,trueobjV);
axis([0 seq(end) 0 maxobj]);
ylabel('Suboptimality','FontSize',FSIZE2);
set(gca,'FontSize',FSIZE1,'XTick',[10^0 10^1 10^2 10^3 10^4 10^5],...
    'XMinorGrid','off');
set(gca, 'YScale', 'log');
% mylegend=legend(methods,'Location','Southwest');
% set(mylegend,'FontSize',FSIZE3);

figure(22);
hold on; grid; xlabel('Iteration','FontSize',FSIZE2);
axis([0 seq(end) 0 maxrank+1]);
ylabel('Rank of iterates','FontSize',FSIZE2);
set(gca,'FontSize',FSIZE1,'XTick',[10^0 10^1 10^2 10^3 10^4 10^5],...
    'XMinorGrid','off');
% mylegend=legend(methods,'Location','Northwest');
% set(mylegend,'FontSize',FSIZE3);

figure(33);
hold on; grid; 
xlabel('Time','FontSize',FSIZE2);
axis([0 inf 0 inf]);
ylabel('Suboptimality','FontSize',FSIZE2);%,'Interpreter','Latex');
set(gca,'FontSize',FSIZE1,'XTick',[10^0 10^1 10^2 10^3 10^4 10^5],...
    'XMinorGrid','off', 'YScale', 'log');
mylegend=legend(method_names,'Location','Southwest');
set(mylegend,'FontSize',FSIZE3);


%% Create PDFs 
topdf(fig11,fnames{1});
topdf(fig22,fnames{2});
topdf(fig33,fnames{3});
end
