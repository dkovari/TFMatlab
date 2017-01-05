function tracks = find_tracks(P, varargin)
% Input:
%   P: cell array of XY point positions
%       each cell corresponds to a different frame

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'MaxDisp',Inf,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'MinLength',1,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'Memory',1,@(x) isscalar(x)&&isnumeric(x)&&x>=1);
addParameter(p,'GraphicDebug',false,@(x) isscalar(x));
addParameter(p,'ShowProgress',true,@(x) isscalar(x));
addParameter(p,'ExpandingSearch',true,@(x) isscalar(x));


parse(p,varargin{:});

maxdist = p.Results.MaxDisp;
maxmem = p.Results.Memory;

DEBUG = logical(p.Results.GraphicDebug);

USE_LENGTH = p.Results.MinLength>1;

PROG_BAR = p.Results.ShowProgress;

EXP_DISP = p.Results.ExpandingSearch;

nFrames = numel(P);
%check dims
nDims = NaN;
orig_reclimit = get(0,'RecursionLimit');
reclimit = orig_reclimit;
for f=1:nFrames
    if numel(P{f})>1
        if isnan(nDims)
            nDims = size(P{f},2);
        else
            if size(P{f},2)~=nDims
                error('Dimension 2 of points list must be the same for all frames');
            end
        end
        
        %get recusion limit
        reclimit = max(reclimit,size(P{f},2));
    end
end
set(0,'RecursionLimit',2*reclimit);


if DEBUG
    hFig = figure();
    clf;
    gca;
    hold on;
    for f=1:nFrames
        plot(P{f}(:,1),P{f}(:,2),'.','markersize',20);
    end
end

rlink_id = cell(nFrames,1);
rlink_fr = cell(nFrames,1);
rdist = cell(nFrames,1);
partdist = cell(nFrames,1);
partdist_id = cell(nFrames,1);
partdist_calced = cell(nFrames,1);

if DEBUG
    rlink_h = cell(nFrames,1);
end


%initialize arrays
for f=1:nFrames
    % link arrays: make rlind_{f}(n)
    rlink_id{f} = zeros(size(P{f},1),1); %index of best future match
    rlink_fr{f} = zeros(size(P{f},1),1); %frame best match occurs in
    rdist{f} = Inf(size(P{f},1),1); %distance between current location and best future location
    if DEBUG
        rlink_h{f} = cell(size(P{f},1),1);
    end
    
    %partdist arrays: make partdist_{f}{n}
    partdist{f} = cell(size(P{f},1),1); %distance to each particle
    partdist_id{f} = cell(size(P{f},1),1); % index of particle for which distance was calculated
    partdist_calced{f} = cell(size(P{f},1),1); % flag specifying if distance was calculated
    for n=1:size(P{f},1) %make partdist_{f}{n}{f}
        partdist{f}{n} = cell(nFrames,1);
        partdist_id{f}{n} = cell(nFrames,1);
        partdist_calced{f}{n} = false(nFrames,1);
    end
end

%initialize partdist arrays

    %recursive function for finding best particle in future frames
    function find_best_forward_match(f, n, next_f)
        %f = current frame
        %n = current particle index
        %next_f = next frame
        if next_f>nFrames||next_f-f>maxmem
            return;
        end
        
        %calculate dist if needed
        if ~partdist_calced{f}{n}(next_f)
            partdist_calced{f}{n}(next_f) = true;
            %build dist array cause it's empty
            partdist{f}{n}{next_f} = sqrt( sum( ...
                bsxfun(@minus,P{next_f},P{f}(n,:)).^2 ,2) ); %calc dist
            partdist_id{f}{n}{next_f} = 1:size(P{next_f},1);
            
            if EXP_DISP
                ind = find(partdist{f}{n}{next_f} > maxdist*(next_f-f));
            else
                ind = find(partdist{f}{n}{next_f} > maxdist);
            end
            partdist{f}{n}{next_f}(ind) = [];
            partdist_id{f}{n}{next_f}(ind) = [];
        end
  
        if DEBUG
            h_start = plot(P{f}(n,1),P{f}(n,2),'+k','markersize',20);
        end
        
        
        if isempty( partdist{f}{n}{next_f} )
            %no candidates try next frame
           find_best_forward_match(f,n,next_f+1);
        else
            if DEBUG
                h_possib = plot(P{next_f}(partdist_id{f}{n}{next_f},1),...
                    P{next_f}(partdist_id{f}{n}{next_f},2),'o','markersize',25);
            end
            
            %best match in candidates
            [d,l] = min( partdist{f}{n}{next_f} ); %find min in the list of distances
            next_id = partdist_id{f}{n}{next_f}(l);
            
            if rlink_id{next_f}(next_id)==0 %next link hasn't been assigned by others
                rlink_id{next_f}(next_id) = n;
                rlink_fr{next_f}(next_id) = f;
                rdist{next_f}(next_id) = d;
                if DEBUG
                    rlink_h{next_f}{next_id} = plot( [P{f}(n,1);P{next_f}(next_id,1)],[P{f}(n,2);P{next_f}(next_id,2)], '-k');
                end
                
            elseif rdist{next_f}(next_id) > d %this particle is closer than prev assigned
                
                old_f = rlink_fr{next_f}(next_id);
                old_id = rlink_id{next_f}(next_id);
                
                if DEBUG
                    fprintf('{f=%d}{id=%d} link to {Nf=%d}{Nid=%d}\n  replacing {of=%d}{oid=%d}\n',...
                        f,n,next_f,next_id,old_f,old_id);

                    h_new = plot( [P{f}(n,1);P{next_f}(next_id,1)],[P{f}(n,2);P{next_f}(next_id,2)], '-m');
                    try
                        set(rlink_h{next_f}{next_id}, 'color','r');
                    catch
                        warning('could not find line');
                    end

                    pause();
                    try
                        delete(rlink_h{next_f}{next_id});
                    catch
                    end
                    rlink_h{next_f}{next_id} = h_new;
                    set(h_new,'color','k');
                end
                
                rlink_id{next_f}(next_id) = n;
                rlink_fr{next_f}(next_id) = f;
                rdist{next_f}(next_id) = d;


                ind = find(partdist_id{old_f}{old_id}{next_f}==next_id);
                partdist{old_f}{old_id}{next_f}(ind) = [];
                partdist_id{old_f}{old_id}{next_f}(ind) = [];

                if ~isempty(partdist{old_f}{old_id}{next_f})
                    %have the displaced track search again
                    find_best_forward_match(old_f,old_id,next_f);
                else %that was the last particle in the tracks list for next_f, have it search in next_f+1
                    find_best_forward_match(old_f,old_id, next_f+1);
                end

            else %the previous assigment was better than this remove that possibility and search again
                ind = find(partdist_id{f}{n}{next_f}==next_id);
                partdist{f}{n}{next_f}(ind) = [];
                partdist_id{f}{n}{next_f}(ind) = [];
                if ~isempty(partdist_id{f}{n}{next_f})
                    %search again
                    find_best_forward_match(f,n,next_f);
                else
                    %search next frame
                    find_best_forward_match(f,n,next_f+1);
                end  
            end
            if DEBUG
                delete(h_possib);
            end
        end
        if DEBUG
            delete(h_start);
        end
    end

%loop forward over all particle and try to find best match
if PROG_BAR
    hWait = waitbar(0,'Linking Tracks');
end
for f=1:nFrames-1
    for n=1:size(P{f},1)
        find_best_forward_match(f,n,f+1)
    end
    if PROG_BAR
        waitbar(f/(nFrames-1),hWait);
    end
end
if PROG_BAR
    close(hWait);
end


%rearrainge data to make it easier to use

P_list = cell(nFrames,1);
for f=1:nFrames
    P_list{f} = 1:size(P{f},1);
end

if USE_LENGTH
    trk_framecount = [];
end

tracks = {};
    function assemble(trk,id,f)
        tracks{trk}(f,:) = P{f}(id,:);
        
        if USE_LENGTH
            trk_framecount(trk)=trk_framecount(trk)+1;
        end
        
        P_list{f}(P_list{f}==id) = [];
        if rlink_id{f}(id)~=0
            assemble(trk,rlink_id{f}(id),rlink_fr{f}(id));
        end
    end
if PROG_BAR
    hWait = waitbar(0,'Assembling Track List');
end
nT=0;
for f=nFrames:-1:1
    while ~isempty(P_list{f})
        nT = nT+1;
        tracks{nT} = NaN(nFrames,nDims);
        
        if USE_LENGTH
            trk_framecount(nT) = 0;
        end
        
        assemble(nT,P_list{f}(1),f);
    end
    if PROG_BAR
        waitbar((nFrames-f+1)/nFrames,hWait);
    end
end
if PROG_BAR
    close(hWait);
end

disp(numel(tracks))
if USE_LENGTH
    tracks(trk_framecount<p.Results.MinLength) = [];
end

set(0,'RecursionLimit',orig_reclimit);

%end of function
end
