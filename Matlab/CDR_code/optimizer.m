
%% Defining ranges of acceptable values
stringerNumber = 1:10;
stringerThickness = .1*exp([1:20]/4.5);
skinThickness = .05*exp([1:20]/6);
frontSparLocation = .1:.01:.35;
backSparLocation = .45:.01:.75;
numRibsValues = 8:30;
sparCapAreas = .02*exp([1:20]/5);

minWeight = inf;
for i = 1:2000;
    %% Making Guesses
    % Number of stringers
    num = length(stringerNumber);
    numTopStringers = stringerNumber(ceil(rand*num));
    numBottomStringers = stringerNumber(ceil(rand*num));
    numNoseTopStringers = stringerNumber(ceil(rand*num));
    numNoseBottomStringers = stringerNumber(ceil(rand*num));

    % Stringer Thicknesses
    num = length(stringerThickness);
    topStringerThick = stringerThickness(ceil(rand*num));
    bottomStringerThick = stringerThickness(ceil(rand*num));
    noseTopStringerThick = stringerThickness(ceil(rand*num));
    noseBottomStringerThick = stringerThickness(ceil(rand*num));

    % Skin Thicknesses
    num = length(skinThickness);
    t_upper = skinThickness(ceil(rand*num));
    t_lower = skinThickness(ceil(rand*num));
    t_upper_front = skinThickness(ceil(rand*num));
    t_lower_front = skinThickness(ceil(rand*num));
    t_frontSpar = skinThickness(ceil(rand*num));
    t_rearSpar = skinThickness(ceil(rand*num));

    % Spar Locations
    frontSpar = frontSparLocation(ceil(rand*length(frontSparLocation)));
    backSpar = backSparLocation(ceil(rand*length(backSparLocation)));

    % Rib Numbers
    numRibs = numRibsValues(ceil(rand*length(numRibsValues)));

    % Spar Cap Areas
    num = length(sparCapAreas);
    sparCapArea1 = sparCapAreas(ceil(rand*num));
    sparCapArea2 = sparCapAreas(ceil(rand*num));
    sparCapArea3 = sparCapAreas(ceil(rand*num));
    sparCapArea4 = sparCapAreas(ceil(rand*num));
    
    % Do not plot
    plotting = 0;

    % Run the wing analysis
    weight = wingAnalysis(numTopStringers, numBottomStringers, numNoseTopStringers, numNoseBottomStringers,...
                                      topStringerThick, bottomStringerThick, noseTopStringerThick, noseBottomStringerThick,...
                                      t_upper, t_lower, t_upper_front, t_lower_front, t_frontSpar, t_rearSpar,...
                                      frontSpar, backSpar,...
                                      numRibs,...
                                      sparCapArea1, sparCapArea2, sparCapArea3, sparCapArea4,...
                                      plotting);
                                  
    if weight < minWeight
        minWeight = weight 

        % Save best values
        b_numTopStringers = numTopStringers;
        b_numBottomStringers = numBottomStringers;
        b_numNoseTopStringers = numNoseTopStringers;
        b_numNoseBottomStringers = numNoseBottomStringers;

        % Stringer Thicknesses
        b_topStringerThick = topStringerThick;
        b_bottomStringerThick = bottomStringerThick;
        b_noseTopStringerThick = noseTopStringerThick;
        b_noseBottomStringerThick = noseBottomStringerThick;

        % Skin Thicknesses
        b_t_upper = t_upper;
        b_t_lower = t_lower;
        b_t_upper_front = t_upper_front;
        b_t_lower_front = t_lower_front;
        b_t_frontSpar = t_frontSpar;
        b_t_rearSpar = t_rearSpar;

        % Spar Locations
        b_frontSpar = frontSpar;
        b_backSpar = backSpar;

        % Rib Numbers
        b_numRibs = numRibs;

        % Spar Cap Areas
        b_sparCapArea1 = sparCapArea1;
        b_sparCapArea2 = sparCapArea2;
        b_sparCapArea3 = sparCapArea3;
        b_sparCapArea4 = sparCapArea4;
    end
end

% Run the wing analysis with best values and plotting
plotting=1;
weight = wingAnalysis(b_numTopStringers, b_numBottomStringers, b_numNoseTopStringers, b_numNoseBottomStringers,...
                                  b_topStringerThick, b_bottomStringerThick, b_noseTopStringerThick, b_noseBottomStringerThick,...
                                  b_t_upper, b_t_lower, b_t_upper_front, b_t_lower_front, b_t_frontSpar, b_t_rearSpar,...
                                  b_frontSpar, b_backSpar,...
                                  b_numRibs,...
                                  b_sparCapArea1, b_sparCapArea2, b_sparCapArea3, b_sparCapArea4,...
                                  plotting);
