function tracks = find_displacement(P,RefFrame,varargin)


p = inputParser;
p.CaseSensitive = false;
addParameter(p,'MaxDisp',Inf,@(x) isscalar(x)&&isnumeric(x));


parse(p,varargin{:});

maxdist = p.Results.MaxDisp;

%% Useful Variables
nFrames = numel(P);
nDims = size(P{RefFrame},2);
nTracks = size(P{RefFrame},1);
tracks = cell(nTracks,1);
for trk=1:nTracks
    tracks{trk} = NaN(nFrames,nDims);
end


%% initialize arrays
rlink_id = cell(nFrames,1);
rlink_fr = cell(nFrames,1);
rdist = cell(nFrames,1);
partdist = cell(nFrames,1);
partdist_id = cell(nFrames,1);
partdist_calced = cell(nFrames,1);

for f=1:nFrames
    % link arrays: make rlind_{f}(n)
    rlink_id{f} = zeros(size(P{f},1),1); %index of best match
    rlink_fr{f} = zeros(size(P{f},1),1); %frame best match occurs in
    rdist{f} = Inf(size(P{f},1),1); %distance between current location and best future location
    
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

%% Find Best match

%recursive function for finding best particle in future frames
    function find_best_match(f, n, next_f)
        %f = current frame
        %n = current particle index
        %next_f = next frame
        if next_f>nFrames
            return;
        end
        
        %calculate dist if needed
        if ~partdist_calced{f}{n}(next_f)
            partdist_calced{f}{n}(next_f) = true;
            %build dist array cause it's empty
            partdist{f}{n}{next_f} = sqrt( sum( ...
                bsxfun(@minus,P{next_f},P{f}(n,:)).^2 ,2) ); %calc dist
            partdist_id{f}{n}{next_f} = 1:size(P{next_f},1);
            
            ind = find(partdist{f}{n}{next_f} > maxdist);

            partdist{f}{n}{next_f}(ind) = [];
            partdist_id{f}{n}{next_f}(ind) = [];
        end
  
        if isempty( partdist{f}{n}{next_f} )
            %no candidates
            return;
        else
            %best match in candidates
            [d,l] = min( partdist{f}{n}{next_f} ); %find min in the list of distances
            next_id = partdist_id{f}{n}{next_f}(l); %id to try to assign this frame
            
            if rlink_id{next_f}(next_id)==0 %next link hasn't been assigned by others
                rlink_id{next_f}(next_id) = n;
                rlink_fr{next_f}(next_id) = f;
                rdist{next_f}(next_id) = d;
                
            elseif rdist{next_f}(next_id) > d %this particle is closer than prev assigned
                old_f = rlink_fr{next_f}(next_id);
                old_id = rlink_id{next_f}(next_id);
                
                rlink_id{next_f}(next_id) = n;
                rlink_fr{next_f}(next_id) = f;
                rdist{next_f}(next_id) = d;

                ind = find(partdist_id{old_f}{old_id}{next_f}==next_id);
                partdist{old_f}{old_id}{next_f}(ind) = [];
                partdist_id{old_f}{old_id}{next_f}(ind) = [];

                if ~isempty(partdist{old_f}{old_id}{next_f})
                    %have the displaced track search again
                    find_best_match(old_f,old_id,next_f);
                else %that was the last particle in the tracks list
                    return
                end
            else %the previous assigment was better than this remove that possibility and search again
                ind = find(partdist_id{f}{n}{next_f}==next_id);
                partdist{f}{n}{next_f}(ind) = [];
                partdist_id{f}{n}{next_f}(ind) = [];
                if ~isempty(partdist_id{f}{n}{next_f})
                    %search again
                    find_best_match(f,n,next_f);
                else %no more to seach
                    return
                end  
            end

        end

    end

%% Find Best with loop
for f=1:nFrames
    for n=1:nTracks
        find_best_match(RefFrame,n,f);
    end
end

%% Convert to tracks structure
for n=1:nTracks
    for f=1:nFrames
        id = find(rlink_id{f}==n);
        if isempty(id)
            tracks{n}(f,:) = NaN;
        else
            tracks{n}(f,:) = P{f}(id,:);
        end
    end
end

end